#= Computes a scattering matrix given the type of collision, collision energy,
and magnetic field (+ parameters with default values).
Module also contains the S_output datastructure, which contains scattering
matrices and initial conditions all together
Description last updated 21/12/2020 =#

module Simulate
export sim, sim_output, μcalc

using UnitfulAtomic, Unitful, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, Channels, matchF, matchK, StateStructures, Solvers

const m4He = 4.002602u"u" # mass of m4He
const m3He = 3.0160293u"u" # mass of m3He
function μcalc(coltype::String)
    @assert coltype ∈ ["3-3","3-4","4-4"] "μcalc() recieved an unrecognised coltype"
    coltype=="3-3" && return 0.5*m3He
    coltype=="3-4" && return (m4He*m3He)/(m3He+m4He)
    coltype=="4-4" && return 0.5*m4He
end

# old code, from when I separated iden_ and diff_ lookups
#=
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
=#

# |αβ⟩ ket, sans |l mₗ⟩ quantum numbers
struct αβ_ket
    α :: atom_nos
    β :: atom_nos
    αβ_ket(α,β) = new(α,β)
end
""" function to strip the |l ml> numbers off an αβlml ket """
lml_stripper(ket::Union{asym_αβlml_ket, scat_αβlml_ket}) = αβ_ket(ket.α,ket.β)

"""unique() but using isapprox instead of isequal"""
function unique_approx(A)
    U = []
    for a in A
        if !any(u->isapprox(a, u), U)
            push!(U, a)
        end
    end
    U
end

"""changes vectors of coefficients that differ only by a phase to be equal"""
function ignore_phase!(vec::Vector{Vector{Float64}})
    n = length(vec)
    for i=1:n
        v = vec[i]
        for j=i:n # only iterate forwards to avoid unnecessary double counting (this fn still does unnecessary work but it's acceptable)
            dp = v ⋅ vec[j] # dot product
            if dp≈0 # orthogonal case
                continue # skip this j
            elseif abs(dp)≈1 # identical, up to a phase
                vec[j]=v # overwrite the w vector with v, so that they have the same phase
            else # this should never trigger
                @warn "Two vectors are neither orthogonal nor different up to a phase"
            end
        end
    end
    vec
end

""" Generates index list indicating which rows/columns of S are alike/different,
    by comparing the linear combinations of |αβ⟩ states, sans |lmₗ⟩ numbers.
    Input: Scattering matrix ~ Nₒ×Nₒ; isOpen ~ N; change of basis P matrix ~ N×N; lookup ~ N
    Output: Vector of integers ~ Nₒ, identifying which |αβ⟩ channel each row/col of S corresponds to;
    unique αβ combinations; legend of |αβ⟩ states for understanding the unique combinations"""
function αβ_index(S::Matrix{ComplexF64}, isOpen::Vector{Bool}, P::Matrix{Float64}, lookup::Union{Vector{asym_αβlml_ket}, Vector{scat_αβlml_ket}})
    N=length(isOpen)
    Nₒ=count(isOpen)
    if true # assert statements
        @assert size(S)==(Nₒ,Nₒ) "S !~ Nₒ×Nₒ"
        @assert size(P)==(N,N) "P !~ N×N"
        @assert length(lookup)==N "length(lookup) ≠ length(isOpen)"
    end
    # construction of list of vectors, describing each row/col in |αβ⟩ basis
    αβs = unique(lml_stripper.(lookup)) # unique |αβ⟩ numbers
    nαβ = length(αβs)
    vec = Vector{Vector{Float64}}() # initialise final output
    for j=1:Nₒ # iterate across cols of P
        jth_vec = zeros(nαβ) # initialise
        Pindex=findall(isOpen)[j] # index of corresponding eigenvector in P
        eigenvec=P[:,Pindex] # eigenvector corresponding to this row/col in S
        for i=1:N # iterate down rows of this eigenvector
            eigenvec[i]≈0 && continue # skip zero rows
            this_αβ = lml_stripper(lookup[i])
            this_αβ_index = findfirst(isequal(this_αβ),αβs) # number of this nonzero represented αβ
            @assert !isnothing(this_αβ_index) "Did not find a matching αβ" # sanity check
            jth_vec[this_αβ_index]+=eigenvec[i] # add the presence of this αβ to the jth identifying vector
        end
        push!(vec,jth_vec) # save index vector for this row/col of S
    end
    @assert length(vec)==Nₒ "length(vec) ≠ Nₒ" # sanity check
    # overwrite vectors that only differ by a phase to make them identical
    ignore_phase!(vec)
    # now use that vector of vectors to create an index list
    unq_vecs = unique_approx(vec); nch=length(unq_vecs) # nch = number of different open channels. σ matrices ~ nch × nch
    @assert nch <= Nₒ "Too many different channels detected"
    indices=zeros(Int,Nₒ)
    for j=1:Nₒ # iterate through the rows/columns of S
        indices[j] = findfirst(v->isapprox(vec[j],v), unq_vecs) # isapprox bc floating point errors arise for 3-3
    end
    return indices, unq_vecs, αβs
end

""" Calculates cross sections from S, kOpen, and the indexing vector"""
function calc_σ(S::Matrix{ComplexF64}, kOpen::Vector{typeof(0e0u"bohr^-1")}, index::Vector{Int})
    Tsq = abs2.(I-S) # transmission coefficients, for el cs
    Ssq = abs2.(S) # square of S-matrix, for ion cs
    Nₒ = length(kOpen) # dimension of S matrix
    if true # assert statements
        @assert size(S)==(Nₒ,Nₒ) "S !~ Nₒ×Nₒ, where Nₒ=length(kOpen)"
        @assert length(index)==Nₒ "length(index) ≠ dimension of S"
    end
    nch = length(unique(index)) # number of open channels
    # initialise cross sections
    σ_el = zeros(nch, nch)u"bohr^2"
    σ_ion = zeros(nch)u"bohr^2"
    # set aside k values for each ch
    kᵧ = zeros(nch)u"bohr^-1" # initialise
    for ch in 1:nch # channel number
        kᵧ[ch] = kOpen[findfirst(isequal(ch),index)]
    end
    prefacs=(x->π/x^2).(kᵧ)
    # fill in elastic
    for i=1:nch, j=1:nch # j → i cross section
        σ_sum=0.0
        for row=1:Nₒ, col=1:Nₒ
            index[row]==i || continue # check this row of S matrix is correct channel
            index[col]==j || continue # check this col of S matrix is correct channel
            σ_sum += Tsq[row,col]
        end
        σ_el[i,j]=prefacs[j]*σ_sum
    end
    # fill in inelastic
    for j=1:nch
        σ_sum=0.0
        for col=1:Nₒ
            index[col]==j || continue # check this col of S matrix is correct channel
            σ_sum += 1 - sum(Ssq[:, j]) # sum down column ↔ all nonunitary outgoing
        end
        σ_ion[j]=prefacs[j]*σ_sum
    end
    σ_el, σ_ion
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

"""Simulation output struct."""
struct sim_output
    σ_el::Matrix{typeof(0e0u"bohr^2")}
    σ_ion::Vector{typeof(0e0u"bohr^2")}
    vecs::Vector{Vector{Float64}}
    αβs::Vector{αβ_ket}
    k::Vector{typeof(0e0u"bohr^-1")}
    coltype::String
    ϵ::Unitful.Energy
    B::Unitful.BField
    lmax::Int
end

""" Runs simulation and returns sim_output object"""
function sim(coltype::String, lmax::Integer,
    ϵ::Unitful.Energy, B::Unitful.BField,
    lhs::Unitful.Length, mid::Unitful.Length,
    rhs::Unitful.Length,
    lhs2mid_spacing::Unitful.Length, rhs2mid_spacing::Unitful.Length)
    @assert coltype∈["3-3", "4-4", "3-4"] "coltype not recognised"
    μ=μcalc(coltype)
    lookup=αβlml_lookup_generator(coltype,"all",lmax)
    N=length(lookup) # total number of computational states, incl. |lml>
    P, Pinv = P_Pinv(lookup,B) # change-of-basis matrix, *from channel to computational basis*
    # precalculate M_el, M_sd, M_zee, M_Γ coefficient matrices
    M_el = fill(zeros(3),N,N)
    M_sd, M_Γ = zeros(N,N), zeros(N,N)
    M_zee, M_hfs = zeros(N,N)u"hartree", zeros(N,N)u"hartree" # H_zee and H_hfs are entirely precalculated
    for j=1:N,i=1:N # fill in coefficient arrays
        M_el[i,j]+=αβlml_eval(H_el_coeffs,lookup[i],lookup[j])
        M_sd[i,j]+=αβlml_eval(H_sd_coeffs,lookup[i],lookup[j])
        M_Γ[i,j]+=αβlml_eval(Γ_GMS_coeffs,lookup[i],lookup[j])
        M_zee[i,j]+=αβlml_eval(H_zee,lookup[i],lookup[j],B)
        M_hfs[i,j]+=αβlml_eval(H_hfs,lookup[i],lookup[j])
    end
    # generate 𝐤sq, vector of asymptotic k² values for channels
    H∞ = M_zee .+ M_hfs
    D∞ = Vector{Unitful.Energy}(diag(Pinv*H∞*P)) # change to diagonal (channel) basis
    @assert length(D∞)==N "length(ksq) ≠ length(lookup)" # sanity check
    isOpen, kOpen, lOpen = isklOpen(D∞, ϵ, μ, lookup) # kOpen, lOpen used for K_matrix later
    Nₒ=count(isOpen) # number of open channels (not summing over l ml yet)
    Nₒ==0 && return zeros(0,0)u"bohr^2", zeros(0)u"bohr^2", zeros(0,0), zeros(0)u"bohr^-1" # trivial case, no need to look at scattering
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
    AR, AL = QR_solver(lookup, AL, ϵ, M_el, M_sd, M_zee, M_hfs, M_Γ, lhs2mid_locs, μ)
    BL, BR = QR_solver(lookup, BR, ϵ, M_el, M_sd, M_zee, M_hfs, M_Γ, rhs2mid_locs, μ)
    # match to find 𝐅=[𝐆; 𝐆'] at rhs which satisfies both BCs
    F = F_matrix(AL, AR, BL, BR)
    F = [Pinv zeros(N,N)u"bohr";
         zeros(N,N)u"bohr^-1" Pinv]*F # change F to channel basis
    F = F[[isOpen;isOpen], :] # delete rows of F corresponding to closed channels
    𝐊 = K_matrix(rhs, F, kOpen, lOpen) # following Mies (1980)
    @assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
    𝐒 = (I+im*𝐊)*inv(I-im*𝐊) # Scattering matrix
    # calculate cross sections
    index, unq_vecs, αβs = αβ_index(𝐒, isOpen, P, lookup)
    kSave=[kOpen[findfirst(isequal(j),index)] for j in 1:length(unq_vecs)]
    σ_el, σ_ion = calc_σ(𝐒, kOpen, index)
    return sim_output(σ_el, σ_ion, unq_vecs, αβs, kSave, coltype, ϵ, B, lmax)
end

end # module
