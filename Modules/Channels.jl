#=  Numerical diagonalisation of computational basis to produce channels.
Designed to work with system of separating |α=β⟩ and |α≠β⟩ lookup states.
Description last updated 5/01 =#

module Channels

export P_Pinv

using Unitful, UnitfulAtomic
using LinearAlgebra
# for access to interactions
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, StateStructures

# Old code, from when I separated iden_ and diff_ lookups
#=
""" Input: lookup (all α=β ⊻ all α≠β), B
    Output: P, Pinv --- both ~ n×n where n=length(lookup) --- and P-block --- for interpreting channels."""
function P_Pinv(lookup::Union{Vector{scat_αβlml_ket},Vector{asym_αβlml_ket}},
    B::Unitful.BField)
    # first, sanity check that lookup kets are all α=β or α≠β, alike
    let αβequal=(x->x.α==x.β).(lookup)
        @assert all(αβequal) || all(.!(αβequal)) "Problem with lookup: all α=β !⊻ all α≠β"
    end
    n=length(lookup)
    lb=findlast(x->x.l==lookup[1].l && x.ml==lookup[1].ml,lookup) # length of block
    @assert mod(n,lb)==0 "block length does not divide length of lookup" # sanity check
    nb=div(n,lb) # number of blocks
    # second, check that all of the blocks in the H∞ matrix will be the same
    let onlyαβ=(x->(x.α,x.β)).(lookup) # strip l, ml numbers
        for j=0:(nb-1)
            @assert onlyαβ[1:lb]==onlyαβ[1+j*lb:(j+1)*lb] "lookup[$(1+j*lb):$((j+1)*lb)] does not match block"
        end
    end
    unq=lookup[1:lb] # kets representative of one block
    H∞ = Matrix{Unitful.Energy}(zeros(lb,lb)u"hartree") # initialise
    for i=1:lb, j=1:lb # calculate asymptotic energies
        bra, ket = unq[i], unq[j]
        H∞[i,j] = αβlml_eval(H_zee, bra, ket, B) + αβlml_eval(H_hfs, bra, ket)
    end
    Pb = eigen(austrip.(H∞)).vectors # P-block, a block of eigenvectors
    Pbinv = inv(Pb)
    # now construct P and Pinv, the CoB matrices of size n×n, using the blocks
    P = Matrix{Float64}(zeros(n,n)) # initialise
    Pinv = Matrix{Float64}(zeros(n,n)) # initialise
    for j=0:(nb-1) # fill in blocks
        P[1+j*lb:(j+1)*lb,1+j*lb:(j+1)*lb]=Pb
        Pinv[1+j*lb:(j+1)*lb,1+j*lb:(j+1)*lb]=Pbinv
    end
    P, Pinv
end

# test
function P_Pinv(lookup::Vector{test_ket},B::Unitful.BField)
    n=length(lookup)
    P = Matrix{Float64}(zeros(n,n)) # initialise
    Pinv = Matrix{Float64}(zeros(n,n)) # initialise
    H∞ = Matrix{Unitful.Energy}(zeros(n,n)u"hartree")
    for i=1:n, j=1:n
        bra, ket = lookup[i], lookup[j]
        H∞[i,j] = αβlml_eval(H_hfs,bra,ket) + αβlml_eval(H_zee,bra,ket,B)
    end
    P=eigen(austrip.(H∞)).vectors
    Pinv=inv(P)
    P, Pinv
end
=#
#######################Single lookup code#######################################

""" returns eigenvectors of the l0 and (if applicable) l1 blocks
    Input: lookup vector, B-field
    Output: matrices of l0 eigenvectors and (if applicable) l1 eigenvectors"""
function Peven_Podd(lookup::Union{Array{scat_αβlml_ket,1},Array{asym_αβlml_ket,1}},
                        B::Unitful.BField)
    n=length(lookup) # number of states
    # if lmax=0, only return a single matrix
    if isnothing(findfirst(ϕ->ϕ.l==1,lookup))
        # compute asymptotic hamiltonian for l=0 ⟺ for all of lookup
        H∞=Array{Unitful.Energy,2}(zeros(n,n)u"hartree")
        for i=1:n, j=1:n
            H∞[i,j] =αβlml_eval(H_zee,lookup[i],lookup[j],B)
            H∞[i,j]+=αβlml_eval(H_hfs,lookup[i],lookup[j])
        end
        eig=eigen(austrip.(H∞))
        return eig.vectors
    end
    # lmax ≧ 1 case. First, identify l=1 indices
    l1start=findfirst(ϕ->ϕ.l==1,lookup) # start of l=1 part
    l1end = let l2start=findfirst(ϕ->ϕ.l==2,lookup)
        isnothing(l2start) ? n : l2start-1 # lmax=1 case ⟹ go to end of lookup
    end
    H∞_l0=Array{Unitful.Energy,2}(zeros(l1start-1,l1start-1)u"hartree")
    for i=1:(l1start-1), j=1:(l1start-1)
        H∞_l0[i,j] =αβlml_eval(H_zee,lookup[i],lookup[j],B)
        H∞_l0[i,j]+=αβlml_eval(H_hfs,lookup[i],lookup[j])
    end
    eig_l0=eigen(austrip.(H∞_l0))
    # separate out first third of l=1, corresponding to ml=-l=-1
    n_l1=l1end-l1start+1 # number of l=1 lookup entries
    n_mlm1=Int(n_l1/3) # 1/3 of above = number of mₗ=-1 entries
    H∞_mlm1=Array{Unitful.Energy,2}(zeros(n_mlm1,n_mlm1)u"hartree")
    for i=1:n_mlm1, j=1:n_mlm1
        H∞_mlm1[i,j] =αβlml_eval(H_zee,lookup[i+l1start-1],lookup[j+l1start-1],B)
        H∞_mlm1[i,j]+=αβlml_eval(H_hfs,lookup[i+l1start-1],lookup[j+l1start-1])
    end
    eig_mlm1=eigen(austrip.(H∞_mlm1))
    return eig_l0.vectors, eig_mlm1.vectors
end

"""Return channel-to-computational change of basis matrix, a block-diag matrix
    where each column is an eigenvector (a channel) expressed in the comp basis
    Input: lookup, B
    Output: 𝐏, the change-of-basis matrix, size n×n where n=length(lookup), 𝐏⁻¹"""
function P_Pinv(lookup::Union{Array{scat_αβlml_ket,1},Array{asym_αβlml_ket,1}},
                        B::Unitful.BField)
    Peo=Peven_Podd(lookup,B) # P_even_odd, eigenvecs for the even and odd blocks
    typeof(Peo)==Tuple{Array{Float64,2},Array{Float64,2}} || return Peo, inv(Peo) # lmax=0 case, done immediately
    n=length(lookup)
    lmax=lookup[end].l # lookup vector ordered from l=0:lmax
    𝐏=zeros(Float64,n,n) # initialise
    𝐏inv = zeros(Float64,n,n)
    Peoinv = inv(Peo[1]), inv(Peo[2])
    for l=0:lmax
        𝐌=iseven(l) ? Peo[1] : Peo[2] # change of basis for one ml
        𝐌inv=iseven(l) ? Peoinv[1] : Peoinv[2]
        n_l=size(𝐌,1) # number of states per ml
        start=findfirst(ϕ->ϕ.l==l,lookup) # where the l starts in lookup
        degen=2l+1 # number of different ml's
        for k=0:(degen-1)
            𝐏[(start+k*n_l):(start+(k+1)*n_l-1),(start+k*n_l):(start+(k+1)*n_l-1)]=𝐌
            𝐏inv[(start+k*n_l):(start+(k+1)*n_l-1),(start+k*n_l):(start+(k+1)*n_l-1)]=𝐌inv
        end
    end
    return 𝐏, 𝐏inv
end

end # module
