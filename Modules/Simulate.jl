#= Computes a scattering matrix given the type of collision, collision energy,
and magnetic field (+ parameters with default values).
Module also contains the S_output datastructure, which contains scattering
matrices and initial conditions all together
Description last updated 21/12/2020 =#

module Simulate
export sim, sim_output

using UnitfulAtomic, Unitful, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, Channels, matchF, matchK, StateStructures, Solvers


""" Produces elastic and ionisation cross sections from scattering matrix,
    kOpen vector and lb n"""
function calc_σ(S, kOpen::Vector{typeof(0e0u"bohr^-1")}, lb::Int)
    @assert mod(length(kOpen), lb)==0 "mod(length(kOpen), lb)≠0"
    nb=div(length(kOpen),lb) # number of blocks
    @assert size(S)==(nb*lb,nb*lb) "Size of S ≠ (no. blocks × length of block)²"
    kᵧ = kOpen[1:lb] # wavenumbers of the different channels
    Tsq = abs2.(I-S) # transmission coefficients, for el cs
    Ssq = abs2.(S) # square of S-matrix, for ion cs
    # initialise cross sections
    σ_el = zeros(lb, lb)u"bohr^2"
    σ_ion = zeros(lb)u"bohr^2"
    prefacs=(x->π/x^2).(kᵧ)
    # fill in elastic
    for i=1:lb, j=1:lb
        σ_sum=0.0
        for kx=0:(nb-1), ky=0:(nb-1)
            σ_sum += Tsq[i+kx*lb, j+ky*lb] # sum over boxes, taking [i,j] coord of each box
        end
        σ_el[i,j]=prefacs[j]*σ_sum
    end
    # fill in inelastic
    for i=1:lb
        σ_sum=0.0
        for k=0:(nb-1) # sum over same channel in diff. boxes
            σ_sum += 1 - sum(Ssq[:, i+k*lb]) # sum down column ↔ all nonunitary outgoing
        end
        σ_ion[i]=prefacs[i]*σ_sum
    end
    return σ_el, σ_ion
end

""" Generates isOpen, kOpen, lOpen vectors from a list of asymptotic energies"""
function isklOpen(D∞::Vector{Unitful.Energy}, ϵ::Unitful.Energy, μ::Unitful.Mass, lookup::Union{Vector{asym_αβlml_ket},Vector{scat_αβlml_ket},Vector{test_ket}})
    ksq=(V->uconvert(u"bohr^-2",2*μ*(ϵ-V)/(1u"ħ^2"))).(D∞)
    isOpen=Vector{Bool}([]) # initialise isOpen
    kOpen=Vector{typeof(0e0u"bohr^-1")}([]) # initialise kOpen
    lOpen=Vector{Int}([]) # initialise lOpen
    for i in 1:length(D∞)
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

""" Calculates σ for a single lookup vector (to be called by sim())"""
function blackbox(lookup::Union{Vector{asym_αβlml_ket},Vector{scat_αβlml_ket}},
    ϵ::Unitful.Energy, B::Unitful.BField,
    lhs::Unitful.Length, mid::Unitful.Length,
    rhs::Unitful.Length,
    lhs2mid_spacing::Unitful.Length, rhs2mid_spacing::Unitful.Length,
    μ::Unitful.Mass)
    ################
    N=length(lookup) # total number of computational states, incl. |lml>
    P, Pinv = P_Pinv(lookup,B) # change-of-basis matrix, *from channel to computational basis*
    # generate 𝐤sq, vector of asymptotic k² values for channels
    H∞=Array{Unitful.Energy,2}(zeros(N,N)u"hartree") # initialise H∞, comp basis asymptotic hamiltonian
    for i=1:N, j=1:N
        H∞[i,j]=αβlml_eval(H_zee,lookup[i],lookup[j],B)+αβlml_eval(H_hfs,lookup[i],lookup[j]) # only H_zee and H_hfs at infinite distance
    end
    D∞ = Vector{Unitful.Energy}(diag(Pinv*H∞*P)) # change to diagonal (channel) basis
    @assert length(D∞)==N "length(ksq) ≠ length(lookup)" # sanity check
    isOpen, kOpen, lOpen = isklOpen(D∞, ϵ, μ, lookup) # kOpen, lOpen used for K_matrix later
    Nₒ=count(isOpen) # number of open channels (not summing over l ml yet)
    Nₒ==0 && return zeros(0,0)u"bohr^2", zeros(0)u"bohr^2", zeros(0,0), zeros(0)u"bohr^-1" # trivial case, no need to look at scattering
    # precalculate M_el, M_sd, M_zee, M_Γ coefficient matrices
    M_el = Array{Vector{Float64}}(undef,N,N)
    M_sd = M_Γ = zeros(N,N)
    M_zee = M_hfs = zeros(N,N)u"hartree" # H_zee and H_hfs are entirely precalculated
    for j=1:N,i=1:N # fill in coefficient arrays
        M_el[i,j]=αβlml_eval(H_el_coeffs,lookup[i],lookup[j])
        M_sd[i,j]=αβlml_eval(H_sd_coeffs,lookup[i],lookup[j])
        M_Γ[i,j]=αβlml_eval(Γ_GMS_coeffs,lookup[i],lookup[j])
        M_zee[i,j]=αβlml_eval(H_zee,lookup[i],lookup[j],B)
        M_hfs[i,j]=αβlml_eval(H_hfs,lookup[i],lookup[j])
    end
    # initialise locations to reorthogonalise
    lhs2mid_locs = let locs=collect(lhs:lhs2mid_spacing:mid)
        if locs[end]!=mid # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
            push!(locs,mid)
        end
        locs
    end
    rhs2mid_locs = let locs=collect(rhs:-rhs2mid_spacing:mid)
        if locs[end]!=mid # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
            push!(locs,mid)
        end
        locs
    end
    # construct lhs and rhs initial conditions
    AL = let AL=[fill(0e0u"bohr",N,N); I]  # all wavefncs vanish, derivs do not
        [P    zeros(N,N)u"bohr";
         zeros(N,N)u"bohr^-1" P]*AL
      end
    @assert size(AL)==(2N,N) "size(AL)≠2N×N" # sanity check
    BR = let BR = let
            BFL = [fill(0.0u"bohr",N,N); I]
            BFR = [Matrix(Diagonal(ones(N))[:,isOpen]u"bohr"); zeros(N,Nₒ)]
            [BFL BFR]
        end
        [P    zeros(N,N)u"bohr";
         zeros(N,N)u"bohr^-1" P]*BR # change into computational basis
    end
    @assert size(BR)==(2N,N+Nₒ) "size(BR)≠2N×(N+Nₒ)" # sanity check
    # solve lhs → mid ← rhs
    AR, AL = orth_solver(lookup, AL, ϵ, M_el, M_sd, M_zee, M_hfs, M_Γ, lhs2mid_locs, μ)
    BL, BR = orth_solver(lookup, BR, ϵ, M_el, M_sd, M_zee, M_hfs, M_Γ, rhs2mid_locs, μ)
    # match to find 𝐅=[𝐆; 𝐆'] at rhs which satisfies both BCs
    F = F_matrix(AL, AR, BL, BR)
    F = [Pinv zeros(N,N)u"bohr";
         zeros(N,N)u"bohr^-1" Pinv]*F # change F to channel basis
    F = F[[isOpen;isOpen], :] # delete rows of F corresponding to closed channels
    𝐊 = K_matrix(rhs, F, kOpen, lOpen) # following Mies (1980)
    @assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
    𝐒 = (I+im*𝐊)*inv(I-im*𝐊) # Scattering matrix
    # calculate cross sections
    lb = let lookupOpen=lookup[isOpen] # lookupOpen is physically meaningless
        findlast(x->x.l==lookupOpen[1].l && x.ml==lookupOpen[1].ml,lookupOpen) # length of a block = number of channels
    end
    σ_el, σ_ion = calc_σ(𝐒, kOpen, lb)
    αβ=unique((x->(x.α,x.β)).(lookup)) # unique atomic configurations
    nαβ=length(αβ) # number of atomic configurations
    Pb = let # change of basis matrix for interpreting the cross sections
        P_open_ch = P[:, isOpen] # change of basis matrix with only open channels
        @assert mod(size(P,1),nαβ)==0 "number of rows in P not divisible by number of unique |αβ>"
        P_open_ch[1:nαβ, 1:lb] # one possibly rectangular block of the change of basis matrix
    end
    return σ_el, σ_ion, Pb, kOpen[1:lb]
end

"""Simulation output struct."""
struct sim_output
    σ_el::Matrix{typeof(0e0u"bohr^2")}
    σ_ion::Vector{typeof(0e0u"bohr^2")}
    P::Matrix{Float64}
    αβ::Vector{Tuple{atom_nos,atom_nos}}
    k::Vector{typeof(0e0u"bohr^-1")}
    coltype::String
    ϵ::Unitful.Energy
    B::Unitful.BField
    lmax::Int
end

""" Runs simulation to give scattering matrices for identical and different lookup vectors.
Input: coltype, lmax, ϵ, B, lhs, mid, rhs, lhs2mid_spacing, rhs2mid_spacing; μ
    Output: S_output containing S_matrices for iden_ and diff_ |αβ⟩, their
    associated CoB matrices and lookup vectors, plus initial conditions"""
function sim(coltype::String, lmax::Int, ϵ::Unitful.Energy, B::Unitful.BField,
    lhs::Unitful.Length, mid::Unitful.Length,
    rhs::Unitful.Length,
    lhs2mid_spacing::Unitful.Length, rhs2mid_spacing::Unitful.Length;
    μ::Unitful.Mass=0.5*4.002602u"u")
    # generate two different lookup vectors
    iden_lookup = αβlml_lookup_generator(coltype, "iden", lmax)
    diff_lookup = αβlml_lookup_generator(coltype, "diff", lmax)
    # generate scattering matrix in each case
    # skip if no symmetric states (3-4 case)
    if length(iden_lookup)==0
        iden_σ_el, iden_σ_ion, iden_P, iden_k = zeros(0,0)u"bohr^2", zeros(0)u"bohr^2", zeros(0,0), zeros(0)u"bohr^-1"
        diff_σ_el, diff_σ_ion, diff_P, diff_k = blackbox(diff_lookup,ϵ,B,lhs,mid,rhs,lhs2mid_spacing,rhs2mid_spacing,μ)
    else
        @assert length(diff_lookup)>0 "length(diff_lookup)!>0" # sanity check
        iden_σ_el, iden_σ_ion, iden_P, iden_k = blackbox(iden_lookup,ϵ,B,lhs,mid,rhs,lhs2mid_spacing,rhs2mid_spacing,μ)
        diff_σ_el, diff_σ_ion, diff_P, diff_k = blackbox(diff_lookup,ϵ,B,lhs,mid,rhs,lhs2mid_spacing,rhs2mid_spacing,μ)
    end
    @assert length(iden_k)+length(diff_k)>0 "No open channels found in iden_ or diff_ lookups" # sanity check
    iden_αβ = unique((x->(x.α,x.β)).(iden_lookup))
    diff_αβ = unique((x->(x.α,x.β)).(diff_lookup))
    σ_el = let
            @assert size(iden_σ_el)[1]==size(iden_σ_el)[2] "iden_σ_el not square" # sanity check
            @assert size(diff_σ_el)[1]==size(diff_σ_el)[2] "diff_σ_el not square" # sanity check
            i = size(iden_σ_el)[1]
            d = size(diff_σ_el)[1]
            [iden_σ_el zeros(i,d)u"bohr^2" # patch together both elastic cs matrices
            zeros(d,i)u"bohr^2" diff_σ_el ]
    end
    σ_ion = vcat(iden_σ_ion, diff_σ_ion) # glue together both ion cs vectors
    P = let # patch together the change-of-basis matrix for interpreting
        if size(iden_P)==(0,0)
            iden_P = zeros(length(iden_αβ), 1) # no open channels, make a zero vector
        elseif size(diff_P)==(0,0)
            diff_P = zeros(length(diff_αβ), 1)  # no open channels, make a zero vector
        end
        iden_m, iden_n = size(iden_P)
        diff_m, diff_n = size(diff_P)
        [iden_P zeros(iden_m,diff_n)
         zeros(diff_m,iden_n) diff_P]
    end
    k = vcat(iden_k, diff_k) # asymptotic wavenumbers of the channels
    αβ=vcat(iden_αβ,diff_αβ) # atomic configurations for reference
    # calculate wavenumbers associated with the channels
    sim_output(σ_el, σ_ion, P, αβ, k, coltype, ϵ, B, lmax)
end

end # module
