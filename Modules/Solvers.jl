#= Contains orth_solver(), re-orthogonalising numerical integrator of the TISE,
and solver(), the function that does the actual numerical integration 'stints'
Description last updated 18/12/20 =#

module Solvers
export solver, QR_solver, DC_solver

using Unitful, UnitfulAtomic, LinearAlgebra
using OrdinaryDiffEq
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using StateStructures, Interactions, Potentials

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
    @assert length(obj.transform)==1 || size(maxval)==(1,length(obj.transform)) "length(obj.transform)≠0 ⩓ size(maxval)≠(1,length(obj.transform))" # sanity check
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
                M_el, M_sd, M_zee, M_hfs, M_Γ,
                lhs::Unitful.Length, rhs::Unitful.Length,
                μ::Unitful.Mass,
                callback::DiscreteCallback;
                C=nothing)
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
    IC⁰ = austrip.(complex.(IC))
    # TISE differential equation
    function TISE(du,u,p,x)
        # Construct V(R) matrix
        V = zeros(ComplexF64,n,n)u"hartree" # initialise
        pots = [Singlet(x*1u"bohr"),Triplet(x*1u"bohr"),Quintet(x*1u"bohr")]
        for j=1:n, i=1:n
            V[i,j] = H_rot(lookup[i],lookup[j], x*1u"bohr", μ, C) # rotational
            V[i,j]+= M_el[i,j] ⋅ pots # electronic
            V[i,j]+= M_sd[i,j]*H_sd_radial(x*1u"bohr") # spin-dipole
            V[i,j]+= M_zee[i,j] # Zeeman
            V[i,j]+= M_hfs[i,j] # Hyperfine
            #V[i,j]+= M_Γ[i,j]*Γ_GMS_radial(x*1u"bohr") #TODO turned off for testing
        end
        V⁰=austrip.(V) # strip units from V
        M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matix
        D = ([0*I I
              M 0*I])
        mul!(du,D,u) # pre-allocated du matrix
        #D*u # ⃗u' = D . ⃗u
    end
    # solve
    prob=ODEProblem(TISE,IC⁰,(lhs⁰,rhs⁰))
    sol_unitless=solve(prob,Tsit5(),reltol=1e-10,save_start=true,save_end=true,save_everystep=false,dense=false,
    callback=callback)
    # add units back
    units = vcat(fill(1.0u"bohr",n),fill(1.0,n))
    sol = x -> sol_unitless(austrip(x)).*units
    #println("Integrated $lhs to $rhs, renormalsed $(callback.affect!.debug_counter) times.") # debugging
    return sol
end

"""Re-orthogonalising DE solver. Takes a list of R values and solves from first
to last, re-orthogonalising at every interim value."""
function QR_solver(lookup::Union{Array{asym_αβlml_ket,1},Array{scat_αβlml_ket,1},Vector{test_ket}},
                    IC, ϵ::Unitful.Energy,
                    M_el, M_sd, M_zee, M_hfs, M_Γ,
                    locs,
                    μ::Unitful.Mass;
                    maxval=1e5,
                    C=nothing) # maxval defines when to renormalise
    # sanity checks on Rs
    @assert length(locs)>=2 "length(locs) < 2"
    @assert all(x->dimension(x)==dimension(1u"m"),locs) "Not all R ∈ locs are lengths"
    # sanity checks on precalculated matrices
    n=length(lookup); ncols=size(IC,2)
    @assert size(M_el)==size(M_sd)==size(M_zee)==size(M_Γ)==(n,n) "Precalculated matrices not of size n×n"
    # initialise callback
    callback=CreateRenormCallback(maxval,ncols)
    # initialise Q, R arrays
    Q, Rcum = IC, Matrix(I,ncols,ncols) # Qs and Rs we will use
    units=vcat(fill(1e0u"bohr",n),fill(1e0,n)) # units, for making Q have units
    for k=1:(length(locs)-1)
        start, finish = locs[k], locs[k+1] # start and finish bounds
        # solve, last stored Q being the IC
        sol=solver(lookup,Q,ϵ,M_el,M_sd,M_zee,M_hfs,M_Γ,start,finish,μ,callback,C=C)
        ψ=sol(finish) # solution evaluated at rhs
        @assert !any(abs2.(austrip.(ψ)) .> maxval^2) "renorm didn't work" # debugging
        ψQR=qr(austrip.(ψ)) # units stripped before QR
        Q = Matrix(ψQR.Q).*units # latest Q matrix
        R=ψQR.R
        @debug "Finished integrating $(start) to $(finish), cond(ψ)=$(cond(austrip.(ψ)))"
        Rcum *= inv(R)
        @debug "cond(Rcum)=$(cond(Rcum))"
        #@show (maximum(abs.(R),dims=1)) # debugging
        #@show (maximum(abs.(inv(R)),dims=1)) # debugging
    end
    ψ=Q
    #@show Rcum # debugging
    IC *= Rcum
    # optionally renormalise ψ
    if any(abs2.(austrip.(ψ)) .> maxval^2)
        colmaxes = sqrt.(maximum(abs2.(austrip.(ψ)), dims=1))
        ψ ./= colmaxes
        @assert size(colmaxes)==(1,length(callback.affect!.transform)) "size(maxval)≠(1,length(obj.transform))" # sanity check
        callback.affect!.transform ./= vec(colmaxes) # can't directly broadcast a vector with a 1xN matrix
    end
    IC *= diagm(callback.affect!.transform) # retroactively apply identical renormalisation to IC
    return ψ, IC
end

#################################DC's idea, alternative to QR re-orthogonalisation##########################
""" Generate Hamiltonian at point R given lookup (and B)"""
function ham(lookup::Union{Vector{scat_αβlml_ket},Vector{asym_αβlml_ket}},
    R::Unitful.Length, M_el, M_sd, M_zee, M_hfs, μ::Unitful.Mass)
    n = length(lookup)
    V = zeros(n,n)u"hartree"
    pots = [Singlet(R),Triplet(R),Quintet(R)]
    for j=1:n, i=1:n
        V[i,j] = H_rot(lookup[i],lookup[j], R, μ) # rotational
        V[i,j]+= M_el[i,j] ⋅ pots # electronic
        V[i,j]+= M_sd[i,j]*H_sd_radial(R) # spin-dipole
        V[i,j]+= M_zee[i,j] # Zeeman
        V[i,j]+= M_hfs[i,j] # Hyperfine
    end
    V
end

"""Re-orthogonalising solver using DC's idea (projection onto the subspace of
open channels first, and then closed channels)
    ONLY FOR MID←RHS INTEGRATION! Assumes N+Nₒ columns in the IC"""
function DC_solver(lookup::Union{Array{asym_αβlml_ket,1},Array{scat_αβlml_ket,1},Vector{test_ket}},
                    IC, ϵ::Unitful.Energy,
                    M_el, M_sd, M_zee, M_hfs, M_Γ,
                    locs,
                    μ::Unitful.Mass;
                    maxval=1e5) # maxval defines when to renormalise
    # sanity checks on Rs
    @assert length(locs)>=2 "length(locs) < 2"
    @assert all(x->dimension(x)==dimension(1u"m"),locs) "Not all R ∈ locs are lengths"
    # sanity checks on precalculated matrices
    n=length(lookup); ncols=size(IC,2); nopen=ncols-n; nclosed=n-nopen # assuming RHS BCs of N+Nₒ columns
    @assert size(M_el)==size(M_sd)==size(M_zee)==size(M_Γ)==(n,n) "Precalculated matrices not of size n×n"
    # initialise callback
    callback=CreateRenormCallback(maxval,ncols)
    # initialise Q, R arrays
    ψ, Rcum = IC, Matrix(I,ncols,ncols) # Qs and Rs we will use
    units=vcat(fill(1e0u"bohr",n),fill(1e0,n)) # units, for making sol have units
    for k=1:(length(locs)-1)
        start, finish = locs[k], locs[k+1] # start and finish bounds
        # solve, last stored Q being the IC
        sol=solver(lookup,ψ,ϵ,M_el,M_sd,M_zee,M_hfs,M_Γ,start,finish,μ,callback)
        A=austrip.(sol(finish)) # solution evaluated at end of stint
        @assert !any(abs2.(A) .> maxval^2) "renorm didn't work" # debugging
        opens = let V=ham(lookup,locs[k+1],M_el,M_sd,M_zee,M_hfs,μ) # hamiltonian at this location (stripped, in Eh)
            F = eigen(austrip.(V))
            F.vectors[:,1:nopen]
        end
        Bopen = [opens zeros(n,nopen); # [0 ̃I] expressed in hyperfine basis
                 zeros(n,nopen) opens] # [̃I 0]
        C = Bopen' * A # projection of the solution onto subspace of open channels
        D = A - (Bopen * C) # subtract from the solution the bits in the subspace
        Bclosed = svd(D).U[:,1:nclosed] # solns represented by closed channels
        B = [Bopen Bclosed]
        ψ = B.*units # save solution, ready to go into the next step
        R=A\B
        @debug "Finished integrating $(start) to $(finish), cond(A)=$(cond(A)), cond(B)=$(cond(B)), cond(R)=$(cond(R))"
        Rcum *= R
        @debug "cond(Rcum)=$(cond(Rcum))"
    end
    IC *= Rcum
    # optionally renormalise ψ
    if any(abs2.(austrip.(ψ)) .> maxval^2)
        colmaxes = sqrt.(maximum(abs2.(austrip.(ψ)), dims=1))
        ψ ./= colmaxes
        @assert size(colmaxes)==(1,length(callback.affect!.transform)) "size(maxval)≠(1,length(obj.transform))" # sanity check
        callback.affect!.transform ./= vec(colmaxes) # can't directly broadcast a vector with a 1xN matrix
    end
    IC *= diagm(callback.affect!.transform) # retroactively apply identical renormalisation to IC
    return ψ, IC
end


end # module
