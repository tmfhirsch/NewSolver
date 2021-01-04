push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")

using Interactions, StateStructures, Channels, Solvers
using Unitful, UnitfulAtomic, LinearAlgebra
using OrdinaryDiffEq

lmax=2

lookup34=αβlml_lookup_generator("3-4",lmax)
lookup33=αβlml_lookup_generator("3-3",lmax)
lookup44=αβlml_lookup_generator("4-4",lmax)

B=0.1u"T"
l44=length(lookup44)
H44=Array{Unitful.Energy,2}(zeros(l44,l44)u"hartree")
for i in 1:l44, j in 1:l44
    bra,ket = lookup44[i],lookup44[j]
    H44[i,j]+=αβlml_eval(H_zee, bra, ket, B)
    H44[i,j]+=αβlml_eval(H_hfs, bra, ket)
end

l33=length(lookup33)
H33=Array{Unitful.Energy,2}(zeros(l33,l33)u"hartree")
for i in 1:l33, j in 1:l33
    bra,ket = lookup33[i],lookup33[j]
    H33[i,j]+=αβlml_eval(H_zee, bra, ket, B)
    H33[i,j]+=αβlml_eval(H_hfs, bra, ket)
end

l34=length(lookup34)
H34=Array{Unitful.Energy,2}(zeros(l34,l34)u"hartree")
for i in 1:l34, j in 1:l34
    bra,ket = lookup34[i],lookup34[j]
    H34[i,j]+=αβlml_eval(H_zee, bra, ket, B)
    H34[i,j]+=αβlml_eval(H_hfs, bra, ket)
end 

find_unq(lookup) = let unq_αβ=unique((x->(x.α,x.β)).(lookup[1:findlast(x->x.l==1,lookup)]))
    @assert all(x->typeof(x)==typeof(lookup[1]), lookup) "lookup entries not all same type"
    [typeof(lookup[1])(x[1],x[2],0,0) for x in unq_αβ]
end
#####################test new diagonalisation routine
function eig(lookup)
    unq=find_unq(lookup)
    nunq=length(unq)
    unqH∞=Array{Unitful.Energy,2}(zeros(nunq,nunq)u"hartree") # init
    for i=1:nunq, j=1:nunq
        bra,ket = unq[i],unq[j]
        unqH∞[i,j] =αβlml_eval(H_zee, bra, ket, B) + αβlml_eval(H_hfs, bra, ket)
    end
    eig=eigen(austrip.(unqH∞))
    eig
end

eig33=eig(lookup33)
eig34=eig(lookup34)
eig44=eig(lookup44)
