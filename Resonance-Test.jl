using Revise

using UnitfulAtomic, Unitful, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, Channels, matchF, matchK, StateStructures, Solvers, Simulate


#C=-1e-5u"hartree"; # testing numbers
ϵ = 1e-8u"hartree"
B=0e0u"T"; #ϵ=-4.257945202650538e-9u"hartree"; B=0.0005u"T";
lhs=3e0u"bohr"; mid=1.4e1u"bohr"; rhs=1.5e1u"bohr";
lhs2mid_spacing=1e-1u"bohr"; rhs2mid_spacing=1e-1u"bohr";
μ=0.5*4.002602u"u";

function restest(C::Unitful.Energy)
    lookup=test_lookup_generator() # playing around w/ lookup vec of |α≠β⟩ states
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
    isOpen, kOpen, lOpen = Simulate.isklOpen(D∞, ϵ, μ, lookup) # kOpen, lOpen used for K_matrix later
    Nₒ=count(isOpen) # number of open channels (not summing over l ml yet)

    # initialise locations to reorthogonalise
    lhs2mid_locs = let locs=collect(lhs:lhs2mid_spacing:mid)
        if locs[end]!=mid # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
            push!(locs,mid)
        end
        for i=1:(length(locs)-1) # check ascending order
            @assert locs[i+1]>locs[i] "lhs2mid_locs[$i+1] ≦ lhs2mid_locs[$i]"
        end
        locs
    end
    rhs2mid_locs = let locs=collect(rhs:-rhs2mid_spacing:mid)
        if locs[end]!=mid # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
            push!(locs,mid)
        end
        for i=1:(length(locs)-1) # check descending order
            @assert locs[i+1]<locs[i] "rhs2mid_locs[$i+1] ≧ rhs2mid_locs[$i]"
        end
        locs
    end

    # construct lhs and rhs initial conditions
    AL = let AL=[fill((0e0+0e0im)u"bohr",N,N); I]  # all wavefncs vanish, derivs do not
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
    ARsol, ALsol = QR_solver(lookup, AL, ϵ, M_el, M_sd, M_zee, M_hfs, M_Γ, lhs2mid_locs, μ, C=C)
    BLsol, BRsol = QR_solver(lookup, BR, ϵ, M_el, M_sd, M_zee, M_hfs, M_Γ, rhs2mid_locs, μ, C=C)
    # match to find 𝐅=[𝐆; 𝐆'] at rhs which satisfies both BCs
    F = F_matrix(ALsol, ARsol, BLsol, BRsol)
    Fch = [Pinv zeros(N,N)u"bohr";
         zeros(N,N)u"bohr^-1" Pinv]*F # change F to channel basis
    FOpen = Fch[[isOpen;isOpen], :] # delete rows of F corresponding to closed channels
    𝐊 = K_matrix(rhs, FOpen, kOpen, lOpen)
    @assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
    𝐒 = (I+im*𝐊)*inv(I-im*𝐊)
    # Following only useful for testing with realistic kets
    # calculate cross sections
    lb = let lookupOpen=lookup[isOpen] # lookupOpen is physically meaningless
        findlast(x->x.l==lookupOpen[1].l && x.ml==lookupOpen[1].ml,lookupOpen) # length of a block = number of channels
    end
    σ_el, σ_ion = Simulate.calc_σ(𝐒, kOpen, lb)
    return σ_el[1], σ_ion[1]
end

σ_els, σ_ions = Vector{typeof(0e0u"bohr^2")}([]), Vector{typeof(0e0u"bohr^2")}([])
Cs = LinRange(-2e-5,-1e-5,50)u"hartree"
for C in Cs
    σ_el, σ_ion = restest(C)
    push!(σ_els,σ_el); push!(σ_ions,σ_ion);
end

using Plots
plot(Cs./1u"hartree",σ_els./1u"bohr^2",
    title="ϵ=$ϵ, mid=$mid, rhs=$rhs, rhs2mid_spacing=$(rhs2mid_spacing)")
vline!([-1.349522044e-5]) # resonance location as identified in Hons
