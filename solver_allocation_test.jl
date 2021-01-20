#= Sketch for figuring out how to reduce allocations in the solver() function.
The blame seems to rest with the TISE function=#
using Revise

using UnitfulAtomic, Unitful, LinearAlgebra, OrdinaryDiffEq
using StaticArrays
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, Channels, matchF, matchK, StateStructures, Solvers, Simulate

# Callback code by Danny Cocks (edited minorly by Tim Hirsch)
mutable struct RenormCallback
    maxvalsqr::Float64
    transform::Vector{Float64}
    debug_counter::Int
end

function CreateRenormCallback(maxval::Float64, n::Int)
    maxvalsqr=maxval^2
    obj = RenormCallback(maxvalsqr, ones(n), 0)
    condition = (u,t,int) -> any(abs2.(u) .> maxvalsqr)
    DiscreteCallback(condition, obj, save_positions=(false,false))
end

function (obj::RenormCallback)(int)
    maxval = sqrt.(maximum(abs2.(int.u), dims=1))
    int.u ./= maxval
    for i = eachindex(int.sol.u)
        int.sol.u[i] ./= maxval
    end
    @assert size(maxval)==(1,length(obj.transform)) "size(maxval)≠(1,length(obj.transform))" # sanity check
    obj.transform ./= vec(maxval) # can't directly broadcast a vector with a 1xN matrix
    @debug "Renormalised" int.t
    obj.debug_counter+=1
    u_modified!(int,true)
    nothing
end

coltype="4-4"; lmax=0; ϵ=1e-12u"hartree"; B=0.0u"T";
lhs=3e0u"bohr"; mid=5e1u"bohr";
μ=0.5*4.002602u"u";

iden_lookup = αβlml_lookup_generator(coltype, "iden", lmax)
diff_lookup = αβlml_lookup_generator(coltype, "diff", lmax)
#lookup=test_lookup_generator() # playing around w/ lookup vec of |α≠β⟩ states
lookup=diff_lookup

N=length(lookup) # total number of computational states, incl. |lml>
P, Pinv = P_Pinv(lookup,B) # change-of-basis matrix, *from channel to computational basis*
# generate 𝐤sq, vector of asymptotic k² values for channels
H∞=Array{Unitful.Energy,2}(zeros(N,N)u"hartree") # initialise H∞, comp basis asymptotic hamiltonian
for i=1:N, j=1:N
    H∞[i,j]=αβlml_eval(H_zee,lookup[i],lookup[j],B)+αβlml_eval(H_hfs,lookup[i],lookup[j]) # only H_zee and H_hfs at infinite distance
end
D∞ = Vector{Unitful.Energy}(diag(Pinv*H∞*P)) # change to diagonal (channel) basis
@assert length(D∞)==N "length(ksq) ≠ length(lookup)" # sanity check
isOpen, kOpen, lOpen = Simulate.isklOpen(D∞, ϵ, μ, lookup) # kOpen, lOpen used for K_matrix later
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
#const M_el_c = M_el; const M_sd_c = M_sd; const M_Γ_c = M_Γ; const M_zee_c = M_zee;


# construct lhs initial conditions
IC = let AL=[fill(0e0u"bohr",N,N); I]  # all wavefncs vanish, derivs do not
    [P    zeros(N,N)u"bohr";
     zeros(N,N)u"bohr^-1" P]*AL
  end
@assert size(IC)==(2N,N) "size(AL)≠2N×N" # sanity check

# solve lhs → mid ← rhs
# Check length and units of IC
n = length(lookup) # number of channels
@assert size(IC)[1]==2*n "Initial condition has wrong number of channels"
if length(size(IC))==1 # IC a vector
    for i=1:n # check wavefunction entries
        @assert dimension(IC[i])==dimension(1u"m") "IC[$i] not a length"
    end
    for i=(n+1):2*n # check derivative entries
        @assert dimension(IC[i])==dimension(1) "IC[$i] not dimensionless"
    end
elseif length(size(IC))==2 # IC a matrix of initial condition vectors
    for i=1:n, j=1:size(IC)[2] # check wavefunction entries
        @assert dimension(IC[i,j])==dimension(1u"m") "IC[$i,$j] not a length"
    end
    for i=(n+1):2*n, j=1:size(IC)[2] # check derivative entries
        @assert dimension(IC[i,j])==dimension(1) "IC[$i,$j] not dimensionless"
    end
else # IC not a 1D vector or 2D array
    error("IC not 1D or 2D")
end

# strip units from constants
ϵ⁰, μ⁰ = austrip(ϵ), austrip(μ)
ħ⁰ = austrip(1.0u"ħ")
lhs⁰, mid⁰ = austrip(lhs), austrip(mid)
# strip units from IC
IC⁰ = austrip.(complex.(IC))

D_times_u=similar(IC⁰) # preallocate
V = zeros(ComplexF64,n,n)u"hartree" # preallocate
V⁰= zeros(ComplexF64,n,n) # preallocate

# TISE differential equation
function TISE(u,p,x)
    # Construct V(R) matrix
    for i=1:n, j=1:n
        print("l.115"); @time V[i,j] = H_rot(lookup[i],lookup[j], x*1u"bohr", μ) # rotational
        @time V[i,j]+= H_el_radial(M_el[i,j], x*1u"bohr") # electronic
        @time V[i,j]+= M_sd[i,j]*H_sd_radial(x*1u"bohr") # spin-dipole
        @time V[i,j]+= M_zee[i,j] # Zeeman
        @time V[i,j]+= αβlml_eval(H_hfs,lookup[i],lookup[j]) # Hyperfine
        @time V[i,j]+= M_Γ[i,j]*Γ_GMS_radial(x*1u"bohr")
    end
    V⁰.=austrip.(V) # strip units from V
    M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matrix
    D = ([0*I I
          M 0*I])
    mul!(D_times_u,D,u)
    #D*u # ⃗u' = D . ⃗u
end

# solve
prob=ODEProblem(TISE,IC⁰,(lhs⁰,mid⁰))
sol_unitless=solve(prob,Tsit5(),reltol=1e-10,save_start=true,save_end=true,save_everystep=false,dense=false,
callback=CreateRenormCallback(1e5,size(IC,2)))
