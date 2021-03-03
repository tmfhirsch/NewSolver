using Revise

using UnitfulAtomic, Unitful, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, Channels, matchF, matchK, StateStructures, Solvers, Simulate

G = 1e-4u"T" # Gauss
k_desired=1e-4u"bohr^-1"
coltype="3-4"; lmax=1; B=0G #B=12000G; #ϵ=-4.257945202650538e-9u"hartree"; B=0.0005u"T";
lhs=3e0u"bohr"; mid=5e1u"bohr"; rhs=1e4u"bohr";
lhs2mid_spacing=1e9u"bohr"; rhs2mid_spacing=1e1u"bohr";
μ=μcalc(coltype)

@assert coltype∈["3-3", "4-4", "3-4"] "coltype not recognised"
unq_lookup = let lookup=αβlml_lookup_generator(coltype,"all",lmax)
    if typeof(lookup[1])==scat_αβlml_ket
        unique((ket->scat_αβlml_ket(ket.α,ket.β,0,0)).(lookup))
    else # asym case
        unique((ket->asym_αβlml_ket(ket.α,ket.β,0,0)).(lookup))
    end
end
unq_N=length(unq_lookup)
H∞ = Matrix{Unitful.Energy}(zeros(unq_N,unq_N)u"hartree") # intialise
for i=1:unq_N, j=1:unq_N
    bra, ket = unq_lookup[i], unq_lookup[j]
    H∞[i,j]=αβlml_eval(H_zee, bra, ket, B)+αβlml_eval(H_hfs, bra, ket)
end
minV∞=(eigen(austrip.(H∞)).values[1])u"hartree" # lowest energy channel is physically relevant
ϵ = (v->auconvert(v+1u"ħ^2"*k_desired^2/(2*μ)))(minV∞) # energy to make low-energy channel desired wavenumber

#=ϵ = let N=length(lookup)
    H∞ = Matrix{Unitful.Energy}(zeros(N,N)u"hartree") # intialise
    for i=1:N, j=1:N
        bra, ket = lookup[i], lookup[j]
        H∞[i,j]=αβlml_eval(H_zee, bra, ket, B)+αβlml_eval(H_hfs, bra, ket)
    end
    minV∞=(eigen(austrip.(H∞)).values[1])u"hartree" # lowest energy channel is physically relevant
    ϵ=(v->auconvert(v+1u"ħ^2"*k_desired^2/(2*μ)))(minV∞)
    @assert isreal(ϵ) "ϵ is not real"
    real(ϵ)
end=#

lookup = αβlml_lookup_generator(coltype, "all", lmax)

#=
#####################################TEST#######################################
testIC = let IC=vcat([i==3 for i in 1:N],fill(0e0,N))
    [P    zeros(N,N)u"bohr";
     zeros(N,N)u"bohr^-1" P]*(IC.*vcat(fill(1e0u"bohr",N),fill(1e0,N)))
end
testsol=solver(lookup,testIC,ϵ,M_el,M_sd,M_zee,M_hfs,M_Γ,rhs,mid,μ,Solvers.CreateRenormCallback(1e5,size(testIC,2)))(mid)
testk=(v->auconvert(sqrt(2*μ*(ϵ-v))/1u"ħ"))(D∞[1])
println("B=$B, testk=$testk")
@show testsol
testQR=QR_solver(lookup,testIC,ϵ,M_el,M_sd,M_zee,M_hfs,M_Γ,rhs2mid_locs,μ)[1]
#################################################################### New BR code
function closed_column(k,N::Int,j::Int,R::Unitful.Length) # generates column for closed channel in BR
    fac = exp(-imag(k)*R)
    wavfn = [i==j ? 1 : 0e0 for i in 1:N]u"bohr"
    deriv = [i==j ? -1/austrip(imag(k))*1 : 0e0 for i in 1:N]
    vcat(wavfn,deriv)
end
function create_BR(k∞,R::Unitful.Length)
    N=length(k∞)
    isOpen=(k->imag(k)==0u"bohr^-1").(k∞)
    Iₒ=Matrix(Diagonal(ones(N)))[:,isOpen]
    nc = count(.!(isOpen)) # number of closed channels
    Bcc = [zeros(N,nc)u"bohr^-1"; # initialise closed ch. part of BR
           zeros(N,nc)]
    counter=1 #count closed channels
    for i=1:N
        isOpen[i] && continue # skip if open channel
        k=k∞[i]
        Bcc[:,counter]=closed_column(k,N,i,R)
        counter+=1
    end
    @assert counter==nc+1 "Didn't write across all of Bcc" # sanity check
    BL = vcat(hcat(Iₒ.*1u"bohr", Iₒ.*0u"bohr"), # open channels allow ψ and ψ'
              hcat(Iₒ.*0,         Iₒ))
    BR = [BL Bcc]
    @assert size(BR)==(2*N,N+count(isOpen)) "BR not ~ 2N×(N+Nₒ)"
    BR
end
#BR=create_BR((V->uconvert(u"bohr^-1",sqrt(complex(2*μ*(ϵ-V)))/(1u"ħ"))).(D∞), rhs)
=#

N = length(lookup)
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
Nₒ==0 && return zeros(0,0)u"bohr^2", zeros(0)u"bohr^2", zeros(0,0), zeros(0)u"bohr^-1" # trivial case, no need to look at scattering
# initialise locations to reorthogonalise
lhs2mid_locs = let locs=collect(lhs:lhs2mid_spacing:mid)
    if locs[end]!=mid # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
        push!(locs,mid)
    end
    locs
end
rhs2mid_locs = let locs=collect(rhs:-rhs2mid_spacing:mid)
    if locs[end]!=mid # in case the spacing doesn't match up, do an extra, shorter stint to finish at the right location
        push!(locs,mid)
    end
    locs
end
# construct lhs and rhs initial conditions
AL = let AL=[fill(0e0u"bohr",N,N); I]  # all wavefncs vanish, derivs do not
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
AR, AL = QR_solver(lookup, AL, ϵ, M_el, M_sd, M_zee, M_hfs, M_Γ, lhs2mid_locs, μ)
BL, BR = QR_solver(lookup, BR, ϵ, M_el, M_sd, M_zee, M_hfs, M_Γ, rhs2mid_locs, μ)
# match to find 𝐅=[𝐆; 𝐆'] at rhs which satisfies both BCs
F = F_matrix(AL, AR, BL, BR)
F = [Pinv zeros(N,N)u"bohr";
     zeros(N,N)u"bohr^-1" Pinv]*F # change F to channel basis
F = F[[isOpen;isOpen], :] # delete rows of F corresponding to closed channels
𝐊 = K_matrix(rhs, F, kOpen, lOpen) # following Mies (1980)
@assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
𝐒 = (I+im*𝐊)*inv(I-im*𝐊) # Scattering matrix
# calculate cross sections
index, unq_vecs, αβs = Simulate.αβ_index(𝐒, isOpen, P, lookup)
kSave=[kOpen[findfirst(isequal(j),index)] for j in 1:length(unq_vecs)]
σ_el, σ_ion = Simulate.calc_σ(𝐒, kOpen, index)

Simulate.sim_output(σ_el,σ_ion,unq_vecs,αβs,kSave,coltype,ϵ,B,lmax)
