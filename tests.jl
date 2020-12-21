push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")

using Interactions, StateStructures, Channels, Solvers
using Unitful, UnitfulAtomic, LinearAlgebra
using OrdinaryDiffEq

lmax=0

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

ch33=ch_matrix(lookup33,B)
ch34=ch_matrix(lookup34,B)
ch44=ch_matrix(lookup44,B)

lhs,rhs=3u"bohr",50u"bohr"
locs=LinRange(lhs,rhs,2)
IC44=[fill(0e0u"bohr",l44,l44); I]

# precalculate coefficient matrices for nondiag interactions
@assert size(IC44,1)==2*l44 "IC does not have 2*length(lookup) rows" # sanity check
M_el_44 = Array{Tuple{Float64,Float64,Float64},2}(undef,l44,l44)
M_sd_44, M_Γ_44 = zeros(l44,l44), zeros(l44,l44)
M_zee_44 = zeros(l44,l44)u"hartree" # H_zee is entirely precalculated
@time for i=1:l44,j=1:l44
    M_el_44[i,j]=αβlml_eval(H_el_coeffs,lookup44[i],lookup44[j])
    M_sd_44[i,j]=αβlml_eval(H_sd_coeffs,lookup44[i],lookup44[j])
    M_Γ_44[i,j]=αβlml_eval(Γ_GMS_coeffs,lookup44[i],lookup44[j])
    M_zee_44[i,j]=αβlml_eval(H_zee,lookup44[i],lookup44[j],B)
end

foo2=orth_solver(lookup44,IC44,1e-12u"hartree",M_el_44,M_sd_44,M_zee_44,M_Γ_44,LinRange(lhs,rhs,2),B=B)
foo3=orth_solver(lookup44,IC44,1e-12u"hartree",M_el_44,M_sd_44,M_zee_44,M_Γ_44,LinRange(lhs,rhs,3),B=B)
foo4=orth_solver(lookup44,IC44,1e-12u"hartree",M_el_44,M_sd_44,M_zee_44,M_Γ_44,LinRange(lhs,rhs,4),B=B)
foo10=orth_solver(lookup44,IC44,1e-12u"hartree",M_el_44,M_sd_44,M_zee_44,M_Γ_44,LinRange(lhs,rhs,10),B=B)
