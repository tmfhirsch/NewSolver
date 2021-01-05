#=  Numerical diagonalisation of computational basis to produce channels,
plus related functions
Description last updated 17/12 =#

module Channels

export l0l1_eigenvecs, ch_matrix

using Unitful, UnitfulAtomic
using LinearAlgebra
# for access to interactions
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Interactions, StateStructures

""" returns eigenvectors of the l0 and (if applicable) l1 blocks
    Input: lookup vector, B-field
    Output: matrices of l0 eigenvectors and (if applicable) l1 eigenvectors"""
function l0l1_eigenvecs(lookup::Union{Array{scat_αβlml_ket,1},Array{asym_αβlml_ket,1}},
                        B::Unitful.BField)
    n=length(lookup) # number of states
    # if lmax=0, only return a single matrix
    if isnothing(findfirst(ϕ->ϕ.l==1,lookup))
        # compute asymptotic hamiltonian for l=0 ⟺ for all of lookup
        H∞=Array{Unitful.Energy,2}(zeros(n,n)u"hartree")
        for i=1:n, j=1:n
            H∞[i,j] =αβlml_eval(H_zee,lookup[i],lookup[j],B)
            H∞[i,j]+=αβlml_eval(H_hfs,lookup[i],lookup[j])
        end
        eig=eigen(austrip.(H∞))
        return eig.vectors
    end
    # lmax ≧ 1 case. First, identify l=1 indices
    l1start=findfirst(ϕ->ϕ.l==1,lookup) # start of l=1 part
    l1end = findlast(ϕ->ϕ.l==1,lookup)
    H∞_l0=Array{Unitful.Energy,2}(zeros(l1start-1,l1start-1)u"hartree")
    for i=1:(l1start-1), j=1:(l1start-1)
        H∞_l0[i,j] =αβlml_eval(H_zee,lookup[i],lookup[j],B)
        H∞_l0[i,j]+=αβlml_eval(H_hfs,lookup[i],lookup[j])
    end
    eig_l0=eigen(austrip.(H∞_l0))
    # separate out first third of l=1, corresponding to ml=-l=-1
    n_l1=l1end-l1start+1 # number of l=1 lookup entries
    n_mlm1=Int(n_l1/3) # 1/3 of above = number of mₗ=-1 entries
    H∞_mlm1=Array{Unitful.Energy,2}(zeros(n_mlm1,n_mlm1)u"hartree")
    for i=1:n_mlm1, j=1:n_mlm1
        H∞_mlm1[i,j] =αβlml_eval(H_zee,lookup[i+l1start-1],lookup[j+l1start-1],B)
        H∞_mlm1[i,j]+=αβlml_eval(H_hfs,lookup[i+l1start-1],lookup[j+l1start-1])
    end
    eig_mlm1=eigen(austrip.(H∞_mlm1))
    return eig_l0.vectors, eig_mlm1.vectors
end

"""Return channel-to-computational change of basis matrix, a block-diag matrix
    where each column is an eigenvector (a channel) expressed in the comp basis
    Input: lookup, B
    Output: 𝐏, the change-of-basis matrix, size n×n where n=length(lookup)"""
function ch_matrix(lookup::Union{Array{scat_αβlml_ket,1},Array{asym_αβlml_ket,1}},
                        B::Unitful.BField)
    l0l1=l0l1_eigenvecs(lookup,B)
    typeof(l0l1)==Tuple{Array{Float64,2},Array{Float64,2}} || return l0l1 # lmax=0 case, done immediately
    n=length(lookup)
    lmax=lookup[end].l # lookup vector ordered from l=0:lmax
    𝐏=zeros(Float64,n,n) # initialise
    for l=0:lmax
        𝐌=iseven(l) ? l0l1[1] : l0l1[2] # change of basis for one ml
        n_l=size(𝐌,1) # number of states per ml
        start=findfirst(ϕ->ϕ.l==l,lookup) # where the l starts in lookup
        degen=2l+1 # number of different ml's
        for k=0:(degen-1)
            𝐏[(start+k*n_l):(start+(k+1)*n_l-1),(start+k*n_l):(start+(k+1)*n_l-1)]=𝐌
        end
    end
    return 𝐏
end

end # module
