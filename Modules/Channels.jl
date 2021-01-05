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

end # module
