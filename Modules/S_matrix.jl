#= Computes a scattering matrix given the type of collision, collision energy,
and magnetic field (+ parameters with default values).
Module also contains the S_output datastructure, which contains scattering
matrices and initial conditions all together
Description last updated 21/12/2020 =#

#module ScatMat

#export S_output, S_matrix

using UnitfulAtomic, Unitful, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, Channels, matchF, matchK, StateStructures, Solvers


# generates isOpen, kOpen, lOpen vectors from a list of asymptotic energies
function isklOpen(D∞::Vector{Unitful.Energy}, ϵ::Unitful.Energy, μ::Unitful.Mass, lookup::Union{Vector{asym_αβlml_ket},Vector{scat_αβlml_ket}})
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

""" Calculates S matrix of channel states, saving as an S_output datastructure
    Input: flag∈["3-3","3-4","4-4"], lmax, ϵ~[E], B~[BField]
    Parameters: lhs, mid, rhs, rrhs, spacings between reorthogonalising,, μ~[M]
    Output: S_output object containing the scattering matrix, flag, lmax, ϵ, B"""
function S_matrix(lookup::Union{Vector{asym_αβlml_ket},Vector{scat_αβlml_ket}},
    ϵ::Unitful.Energy, B::Unitful.BField,
    lhs::Unitful.Length, mid::Unitful.Length,
    rhs::Unitful.Length, rrhs::Unitful.Length,
    lhs2mid_spacing::Unitful.Length, rhs2mid_spacing::Unitful.Length,
    rhs2rrhs_spacing::Unitful.Length,
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
    # precalculate M_el, M_sd, M_zee, M_Γ coefficient matrices
    M_el = Array{Tuple{Float64,Float64,Float64},2}(undef,N,N)
    M_sd, M_Γ = zeros(N,N), zeros(N,N)
    M_zee = zeros(N,N)u"hartree" # H_zee is entirely precalculated (no radial fn)
    for i=1:N,j=1:N # fill in coefficient arrays
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
    rhs2mid_locs = let locs=collect(rhs:-rhs2mid_spacing:mid)
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
    AR, AL = orth_solver(lookup, AL, ϵ, M_el, M_sd, M_zee, M_Γ, lhs2mid_locs, B, μ)
    BL, BR = orth_solver(lookup, BR, ϵ, M_el, M_sd, M_zee, M_Γ, rhs2mid_locs, B, μ)
    # match to find 𝐅=[𝐆; 𝐆'] at rhs which satisfies both BCs
    F = F_matrix(AL, AR, BL, BR)
    # solve F out to rrhs before matching to bessel functions
    F = orth_solver(lookup, F, ϵ, M_el, M_sd, M_zee, M_Γ, rhs2rrhs_locs, B, μ)[1] # [1] bc only need final value
    F = [Pinv zeros(N,N)u"bohr";
         zeros(N,N)u"bohr^-1" Pinv]*F # change F to channel basis
    F = F[[isOpen;isOpen], :] # delete rows of F corresponding to closed channels
    𝐊 = K_matrix(rrhs, F, kOpen, lOpen)
    @assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
    𝐒 = (I+im*𝐊)*inv(I-im*𝐊)
end

#end # module

""" Data structure for containing scattering matrices for a simulation,
 and the initial conditions of that simulation"""
struct S_output
    diff_S :: Matrix{Complex{Float64}}
    iden_S :: Matrix{Complex{Float64}}
    coltype :: String # "3-3" etc
    lmax :: Int
    ϵ :: Unitful.Energy
    B :: Unitful.BField
end


""" Runs simulation to give scattering matrices for identical and different lookup vectors.
    Output: S_output containing S_matrices for iden_ and diff_ |αβ⟩, their
    associated CoB matrices and lookup vectors, plus initial conditions"""
function sim(coltype::String, lmax::Int, ϵ::Unitful.Energy, B::Unitful.BField;
    lhs::Unitful.Length=3e0u"bohr", mid::Unitful.Length=5e0u"bohr",
    rhs::Unitful.Length=2e2u"bohr", rrhs::Unitful.Length=1e4u"bohr",
    lhs2mid_spacing::Unitful.Length=1e0u"bohr", rhs2mid_spacing::Unitful.Length=1e1u"bohr",
    rhs2rrhs_spacing::Unitful.Length=1e2u"bohr",
    μ::Unitful.Mass=0.5*4.002602u"u")
    # generate two different lookup vectors
    iden_lookup = αβlml_lookup_generator(coltype, "iden", lmax)
    diff_lookup = αβlml_lookup_generator(coltype, "diff", lmax)
    # generate scattering matrix in each case
    # skipping calculation if the lookup vectors are empty (lmax=0 can produce this scenario)
    iden_S = length(iden_lookup)>0 ? S_matrix(iden_lookup, ϵ, B, lhs, mid, rhs, rrhs,
    lhs2mid_spacing, rhs2mid_spacing, rhs2rrhs_spacing, μ) : Matrix{Complex{Float64}}(undef,0,0)
    diff_S = length(diff_lookup)>0 ? S_matrix(diff_lookup, ϵ, B, lhs, mid, rhs, rrhs,
    lhs2mid_spacing, rhs2mid_spacing, rhs2rrhs_spacing, μ) : Matrix{Complex{Float64}}(undef,0,0)
    S_output(diff_S, iden_S, coltype, lmax, ϵ, B)
end


#################################Testing########################################
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


coltype="4-4"; lmax=2; ϵ=1e-12u"hartree"; B=0u"T";
lhs=3e0u"bohr"; mid=5e0u"bohr"; rhs=2e2u"bohr"; rrhs=1e4u"bohr";
lhs2mid_spacing=1e0u"bohr"; rhs2mid_spacing=1e1u"bohr"; rhs2rrhs_spacing=1e2u"bohr";
μ=0.5*4.002602u"u";

iden_lookup = αβlml_lookup_generator(coltype, "iden", lmax)
diff_lookup = αβlml_lookup_generator(coltype, "diff", lmax)
lookup=diff_lookup # playing around w/ lookup vec of |α≠β⟩ states

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
# precalculate M_el, M_sd, M_zee, M_Γ coefficient matrices
M_el = Array{Tuple{Float64,Float64,Float64},2}(undef,N,N)
M_sd, M_Γ = zeros(N,N), zeros(N,N)
M_zee = zeros(N,N)u"hartree" # H_zee is entirely precalculated (no radial fn)
for i=1:N,j=1:N # fill in coefficient arrays
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
rhs2mid_locs = let locs=collect(rhs:-rhs2mid_spacing:mid)
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
AR, AL = orth_solver(lookup, AL, ϵ, M_el, M_sd, M_zee, M_Γ, lhs2mid_locs, B, μ)
BL, BR = orth_solver(lookup, BR, ϵ, M_el, M_sd, M_zee, M_Γ, rhs2mid_locs, B, μ)
# match to find 𝐅=[𝐆; 𝐆'] at rhs which satisfies both BCs
F = F_matrix(AL, AR, BL, BR)
# solve F out to rrhs before matching to bessel functions
F = orth_solver(lookup, F, ϵ, M_el, M_sd, M_zee, M_Γ, rhs2rrhs_locs, B, μ)[1] # [1] bc only need final value
F = [Pinv zeros(N,N)u"bohr";
     zeros(N,N)u"bohr^-1" Pinv]*F # change F to channel basis
F = F[[isOpen;isOpen], :] # delete rows of F corresponding to closed channels
𝐊 = K_matrix(rrhs, F, kOpen, lOpen)
@assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
𝐒 = (I+im*𝐊)*inv(I-im*𝐊)
# calculate cross sections
lb = let lookupOpen=lookup[isOpen] # lookupOpen is physically meaningless
    findlast(x->x.l==lookupOpen[1].l && x.ml==lookupOpen[1].ml,lookupOpen) # length of a block = number of channels
end
σ_el, σ_ion = calc_σ(𝐒, kOpen, lb)
