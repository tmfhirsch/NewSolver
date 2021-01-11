#= Contains orth_solver(), re-orthogonalising numerical integrator of the TISE,
and solver(), the function that does the actual numerical integration 'stints'
Description last updated 18/12/20 =#

module Solvers
export solver, orth_solver

using Unitful, UnitfulAtomic, LinearAlgebra
using OrdinaryDiffEq
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using StateStructures, Interactions

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
    DiscreteCallback(condition, obj, save_positions=(true,false))
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

# following by Tim Hirsch
""" TISE solver
    lookup, the lookup vector; IC, the initial wavefunction matrix (at lhs);
    ϵ, the energy; M_el, M_sd, M_zee, M_Γ, the precalculated coefficient
    matrices for the corresponding non-diagonal interactions;
    lhs and rhs, the start/end points; B and μ, the external magnetic field
    and effective collisional mass; callback, for storing renormalisations
    Output: Solution to numerically integrating the TISE, from IC at lhs to rhs."""
function solver(lookup::Union{Array{asym_αβlml_ket,1},Array{scat_αβlml_ket,1},Vector{test_ket}},
                IC, ϵ::Unitful.Energy,
                M_el, M_sd, M_zee, M_Γ,
                lhs::Unitful.Length, rhs::Unitful.Length,
                μ::Unitful.Mass,
                callback::DiscreteCallback)
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
    sol_unitless=solve(prob,Tsit5(),reltol=1e-10,save_start=true,save_end=true,save_everystep=false,dense=false,
    callback=callback)
    # add units back
    units = vcat(fill(1.0u"bohr",n),fill(1.0,n))
    sol = x -> sol_unitless(austrip(x)).*units # TODO see notes 12/10
    return sol
end

"""Re-orthogonalising DE solver. Takes a list of R values and solves from first
to last, re-orthogonalising at every interim value."""
function orth_solver(lookup::Union{Array{asym_αβlml_ket,1},Array{scat_αβlml_ket,1},Vector{test_ket}},
                    IC, ϵ::Unitful.Energy,
                    M_el, M_sd, M_zee, M_Γ,
                    locs,
                    μ::Unitful.Mass;
                    maxval=1e5) # maxval defines when to renormalise
    # sanity checks on Rs
    @assert length(locs)>=2 "length(locs) < 2"
    @assert all(x->dimension(x)==dimension(1u"m"),locs) "Not all R ∈ locs are lengths"
    # sanity checks on precalculated matrices
    n=length(lookup); ncols=size(IC,2)
    @assert size(M_el)==size(M_sd)==size(M_zee)==size(M_Γ)==(n,n) "Precalculated matrices not of size n×n"
    # initialise callback
    callback=CreateRenormCallback(maxval,ncols)
    # initialise Q, R arrays
    Qs, Rs = Array{Any}([IC]), Array{Any}([Matrix(I,ncols,ncols)])
    units=vcat(fill(1e0u"bohr",n),fill(1e0,n)) # units, for making Q have units
    for k=1:(length(locs)-1)
        start, finish = locs[k], locs[k+1] # start and finish bounds
        # solve, last stored Q being the IC
        sol=solver(lookup,Qs[k],ϵ,M_el,M_sd,M_zee,M_Γ,start,finish,μ,callback)
        ψ=sol(finish) # solution evaluated at rhs
        @assert !any(abs2.(austrip.(ψ)) .> maxval^2) "renorm didn't work" # debugging nonfunctional renorm
        ψQR=qr(austrip.(ψ)) # units stripped before QR
        push!(Qs,Matrix(ψQR.Q).*units) # save orthogonalised soln w/ units
        # renorm R
        R=ψQR.R
        maxR=maximum(abs.(R), dims=1)
        R ./= maxR # normalise R
        callback.affect!.transform ./= vec(maxR) # save renorm in transform
        push!(Rs,R) # save R
        # debugging
        let Qmax=maximum(sqrt.(abs2.(austrip.(Qs[k]))))
            Rmax=maximum(sqrt.(abs2.(Rs[k])))
            @info "k=$k, maxQ=$Qmax, maxR=$Rmax"
        end
    end
    ψ=Qs[end] # at this point Q[end] is the Q of the final solution
    for R in reverse(Rs)
        ψ=ψ*R # reverse the orthogonalisation process
    end
    IC *= diagm(callback.affect!.transform) # retroactively apply identical renormalisation to IC
    @info "Integrating $(locs[1]) → $(locs[end]), renormalised $(callback.affect!.debug_counter) times"
    @info callback.affect!.transform
    return ψ, IC
end

end # module#
