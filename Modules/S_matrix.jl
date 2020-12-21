#= Computes a scattering matrix given the type of collision, collision energy,
and magnetic field (+ parameters with default values).
Module also contains the S_output datastructure, which contains scattering
matrices and initial conditions all together
Description last updated 21/12/2020 =#

using UnitfulAtomic, Unitful, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, Channels, matchF, matchK, StateStructures, Solvers

""" Data structure for containing a scattering matrix and the initial conditions
that produced it """
struct S_output
    S :: Array{Complex{Float64},2}
    flag :: String
    lmax :: Int
    ϵ :: Unitful.Energy
    B :: Unitful.BField
end


# generates isOpen, kOpen, lOpen vectors from a list of asymptotic energies
function isklOpen(D∞::Vector{Unitful.Energy}, ϵ::Unitful.Energy, μ::Unitful.Mass, lookup::Union{Vector{asym_αβlml_ket},Vector{scat_αβlml_ket}})
    ksq=(V->uconvert(u"bohr^-2",2*μ*(ϵ-V)/(1u"ħ^2"))).(D∞)
    isOpen=Vector{Bool}([]) # initialise isOpen
    kOpen=Vector{typeof(0e0u"bohr^-1")}([]) # initialise kOpen
    lOpen=Vector{Int}([]) # initialise lOpen
    for i in 1:N
        k²=ksq[i]
        if k² > 0u"bohr^-2" # positive energy ⟺ open channel
            push!(isOpen, true)
            push!(kOpen, uconvert(u"bohr^-1",sqrt(k²)))
            push!(lOpen, lookup[i].l) # comp and channel basis share same |l ml> block diag structure
        else
            push!(isOpen, false)
        end
    end
    return isOpen, kOpen, lOpen
end

""" Calculates S matrix of channel states, saving as an S_output datastructure
    Input: flag∈["3-3","3-4","4-4"], lmax, ϵ~[E], B~[BField]
    Parameters: lhs, mid, rhs, rrhs, spacings between reorthogonalising,, μ~[M]
    Output: S_output object containing the scattering matrix, flag, lmax, ϵ, B"""
function S_matrix(flag::String, lmax::Int, ϵ::Unitful.Energy, B::Unitful.BField;
    lhs::Unitful.Length=3e0u"bohr", mid::Unitful.Length=5e1u"bohr",
    rhs::Unitful.Length=2e2u"bohr", rrhs::Unitful.Length=1e6u"bohr",
    lhs2mid_spacing::Unitful.Length=5e0u"bohr", rhs2mid_spacing::Unitful.Length=1e1u"bohr",
    rhs2rrhs_spacing::Unitful.Length=1e4u"bohr",
    μ::Unitful.Mass=0.5*4.002602u"u")
    # Generate lookup, form channels
    lookup=αβlml_lookup_generator(flag, lmax)
    N=length(lookup) # total number of computational states, incl. |lml>
    P = let P=ch_matrix(lookup,B) # change-of-basis matrix, *from channel to computational basis*
        @assert size(P)==(N,N) "size(P)≠N×N"
        [P 0*I; 0*I P] # changes basis for solution matrices of the form [𝛙; 𝛙']
    end
    # generate 𝐤sq, vector of asymptotic k² values for channels
    H∞=Array{Unitful.Energy,2}(zeros(N,N)u"hartree") # initialise H∞, comp basis asymptotic hamiltonian
    for i=1:N, j=1:N
        H∞[i,j]=αβlml_eval(H_zee,lookup[i],lookup[j],B)+αβlml_eval(H_hfs,lookup[i],lookup[j]) # only H_zee and H_hfs at infinite distance
    end
    D∞ = let P=ch_matrix(lookup,B) # larger-scope P is 2N×2N, for changing wavefunction matrices.
        Vector{Unitful.Energy}(diag(inv(P)*H∞*P)) # change to diagonal (channel) basis
    end
    @assert length(D∞)==N "length(ksq) ≠ length(lookup)" # sanity check
    isOpen, kOpen, lOpen = isklOpen(D∞, ϵ, μ, lookup) # kOpen, lOpen used for K_matrix later
    Nₒ=length(all(isOpen)) # number of open channels
    # construct lhs and rhs initial conditions
    AL = let AL = [fill(0e0u"bohr",N,N); I] # only wavefunction derivatives
        P*AL # change into computational basis
    end
    BR = let BR = let
            BFL = [fill(0.0u"bohr",N,N); I]
            BFR = [Matrix(Diagonal(ones(N))[:,isOpen]u"bohr"); zeros(N,Nₒ)]
            [BFL BFR]
        end
        P*BR # change into computational basis
    end
    # precalculate M_el, M_sd, M_zee, M_Γ coefficient matrices
    M_el = Array{Tuple{Float64,Float64,Float64},2}(undef,N,N)
    M_sd, M_Γ = zeros(N,N), zeros(N,N)
    M_zee = zeros(N,N)u"hartree" # H_zee is entirely precalculated (no radial fn)
    for i=1:N,j=1:N
        M_el[i,j]=αβlml_eval(H_el_coeffs,lookup[i],lookup[j])
        M_sd[i,j]=αβlml_eval(H_sd_coeffs,lookup[i],lookup[j])
        M_Γ[i,j]=αβlml_eval(Γ_GMS_coeffs,lookup[i],lookup[j])
        M_zee[i,j]=αβlml_eval(H_zee,lookup[i],lookup[j],B)
    end
    # initialise locations to reorthogonalise
    lhs2mid_locs = let locs=collect(lhs:lhs2mid_spacing:mid)
        if locs[end]!=mid # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
            push!(locs,mid)
        end
        locs
    end
    rhs2mid_locs = let locs=collect(rhs:rhs2mid_spacing:mid)
        if locs[end]!=mid # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
            push!(locs,mid)
        end
        locs
    end
    rhs2rrhs_locs = let locs=collect(rhs:rhs2rrhs_spacing:rrhs)
        if locs[end]!=rrhs # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
            push!(locs,rrhs)
        end
        locs
    end
    # solve lhs → mid ← rhs
    AR = orth_solver(lookup, AL, ϵ, M_el, M_sd, M_zee, M_Γ, lhs2mid_locs, B, μ)
    BL = orth_solver(lookup, BR, ϵ, M_el, M_sd, M_zee, M_Γ, rhs2mid_locs, B, μ)
    # match to find 𝐅=[𝐆; 𝐆'] at rhs which satisfies both BCs
    F = F_matrix(AL, AR, BL, BR)
    # solve F out to rrhs before matching to bessel functions
    F = orth_solver(lookup, F, ϵ, M_el, M_sd, M_zee, M_Γ, rhs2rrhs_locs, B, μ)
    F = inv(P)*F # change F to channel basis
    F = F[isOpen, :] # delete rows of F corresponding to closed channels
    𝐊 = K_matrix(rrhs, F, kOpen, lOpen)
    @assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
    𝐒 = (I+im*𝐊)*inv(I-im*𝐊)
    S_output(𝐒, flag, lmax, ϵ, B)
end
