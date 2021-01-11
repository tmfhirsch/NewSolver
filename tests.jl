push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")

using Interactions, StateStructures, Channels, Solvers
using Unitful, UnitfulAtomic, LinearAlgebra
using OrdinaryDiffEq

##################Testing why solutions are blowing up##########################
coltype="4-4"; lmax=0; ϵ=1e-12u"hartree"; B=0.00u"T"; μ=0.5*4.002602u"u"
lookup=vcat(αβlml_lookup_generator(coltype,"iden",lmax),αβlml_lookup_generator(coltype,"diff",lmax))
n=length(lookup)
locs=[2e2u"bohr",1e4u"bohr"]#
IC=[diagm(ones(n))u"bohr" zeros(n,n)u"bohr"
    zeros(n,n) I]
# taken verbatim from proper code###############################################
M_el = Matrix{Tuple{Float64,Float64,Float64}}(undef,n,n)
M_sd, M_Γ = zeros(n,n), zeros(n,n)
M_zee = zeros(n,n)u"hartree" # H_zee is entirely precalculated (no radial fn)
for i=1:n,j=1:n # fill in coefficient arrays
    M_el[i,j]=αβlml_eval(H_el_coeffs,lookup[i],lookup[j])
    M_sd[i,j]=αβlml_eval(H_sd_coeffs,lookup[i],lookup[j])
    M_Γ[i,j]=αβlml_eval(Γ_GMS_coeffs,lookup[i],lookup[j])
    M_zee[i,j]=αβlml_eval(H_zee,lookup[i],lookup[j],B)
end#############################################################################

sol_Ψ, sol_IC =orth_solver(lookup,IC,ϵ,M_el,M_sd,M_zee,M_Γ,locs,μ)
