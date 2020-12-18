#= Contains orth_solver(), re-orthogonalising numerical integrator of the TISE,
and solver(), the function that does the actual numerical integration 'stints'
Description last updated 18/12/20 =#

module Solvers
export solver, orth_solver

using Unitful, UnitfulAtomic, LinearAlgebra
using OrdinaryDiffEq
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using StateStructures, Interactions

"""Callback function for renormalisation of wavefunction. Code by DC"""
function CreateRenormalisedCallback(maxval=1e5)
    maxvalsqr = maxval^2
    condition = (u,t,int) -> any(abs2.(u) .> maxvalsqr)
    DiscreteCallback(condition, _Renormalise!, save_positions=(true,false))
end
function _Renormalise!(int)
    maxval = sqrt.(maximum(abs2.(int.u), dims=1))
    int.u ./= maxval
    for i = eachindex(int.sol.u)
        int.sol.u[i] ./= maxval
    end
    nothing
end


""" TISE solver
    lookup, the lookup vector; IC, the initial wavefunction matrix (at lhs);
    ϵ, the energy; M_el, M_sd, M_zee, M_Γ, the precalculated coefficient
    matrices for the corresponding non-diagonal interactions;
    lhs and rhs, the start/end points; B and μ, the external magnetic field
    and effective collisional mass.
    Output: Solution to numerically integrating the TISE, from IC at lhs to rhs."""
function solver(lookup::Union{Array{asym_αβlml_ket,1},Array{scat_αβlml_ket,1}},
                IC, ϵ::Unitful.Energy,
                M_el, M_sd, M_zee, M_Γ,
                lhs::Unitful.Length, rhs::Unitful.Length,
                B::Unitful.BField, μ::Unitful.Mass)
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
    lhs⁰, rhs⁰ = austrip(lhs), austrip(rhs)
    # TISE differential equation
    function TISE(u,p,x)
        # Construct V(R) matrix
        V = zeros(ComplexF64,n,n)u"hartree" # initialise
        for i=1:n, j=1:n
            V[i,j] = H_rot(lookup[i],lookup[j], x*1u"bohr", μ) # rotational
            V[i,j]+= H_el_radial(M_el[i,j], x*1u"bohr") # electronic
            V[i,j]+= M_sd[i,j]*H_sd_radial(x*1u"bohr") # spin-dipole
            V[i,j]+= M_zee[i,j] # Zeeman
            V[i,j]+= αβlml_eval(H_hfs,lookup[i],lookup[j]) # Hyperfine
            V[i,j]+= M_Γ[i,j]*Γ_GMS_radial(x*1u"bohr")
        end
        V⁰=austrip.(V) # strip units from V
        M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matix
        D = ([0*I I
              M 0*I])
        D*u # ⃗u' = D . ⃗u
    end
    # strip units from IC
    IC⁰ = austrip.(complex.(IC))
    # solve
    prob=ODEProblem(TISE,IC⁰,(lhs⁰,rhs⁰))
    callback=CreateRenormalisedCallback()
    sol_unitless=solve(prob,Tsit5(),reltol=1e-10,save_start=true,save_end=true,save_everystep=false,dense=false,
    callback=callback)
    # add units back
    units = vcat(fill(1.0u"bohr",n),fill(1.0,n))
    sol = x -> sol_unitless(austrip(x)).*units # TODO see notes 12/10
    return sol
end

"""Re-orthogonalising DE solver. Takes a list of R values and solves from first
to last, re-orthogonalising at every interim value."""
function orth_solver(lookup::Union{Array{asym_αβlml_ket,1},Array{scat_αβlml_ket,1}},
                    IC, ϵ::Unitful.Energy, locs;
                    B::Unitful.BField=0u"T", μ::Unitful.Mass=0.5*4.002602u"u")
    # sanity checks on Rs
    @assert length(locs)>=2 "length(locs) < 2"
    @assert all(x->dimension(x)==dimension(1u"m"),locs) "Not all R∈locs are lengths"
    # precalculate coefficient matrices for nondiag interactions
    n=length(lookup)
    @assert size(IC,1)==2*n "IC does not have 2*length(lookup) rows" # sanity check
    M_el = Array{Tuple{Float64,Float64,Float64},2}(undef,n,n)
    M_sd, M_Γ = zeros(n,n), zeros(n,n)
    M_zee = zeros(n,n)u"hartree" # H_zee is entirely precalculated
    @info "Precalculating matrices"
    @time for i=1:n,j=1:n
        M_el[i,j]=αβlml_eval(H_el_coeffs,lookup[i],lookup[j])
        M_sd[i,j]=αβlml_eval(H_sd_coeffs,lookup[i],lookup[j])
        M_Γ[i,j]=αβlml_eval(Γ_GMS_coeffs,lookup[i],lookup[j])
        M_zee[i,j]=αβlml_eval(H_zee,lookup[i],lookup[j],B)
    end
    # initialise Q, R arrays
    Qs, Rs = [IC], [Matrix(I,n,n)]
    units=vcat(fill(1e0u"bohr",n),fill(1e0,n)) # units, for making Q have units
    @info "Beginning DE solving"
    for k=1:(length(locs)-1)
        lhs, rhs = locs[k], locs[k+1] # left and right bounds
        # solve, last stored Q being the IC
        sol=solver(lookup,Qs[k],ϵ,M_el,M_sd,M_zee,M_Γ,lhs,rhs,B,μ)
        ψ=sol(rhs) # solution evaluated at rhs
        @info "up to here, k=$k"
        ψQR=qr(austrip.(ψ)) # units stripped before QR
        push!(Qs,Matrix(ψQR.Q).*units) # save orthogonalised soln w/ units
        push!(Rs,ψQR.R) # save R
    end
    ψ=Qs[end] # at this point Q[end] is the Q of the final solution
    for R in reverse(Rs)
        ψ=ψ*R # reverse the orthogonalisation process
    end
    return ψ
end

end # module
