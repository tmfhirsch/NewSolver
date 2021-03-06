using Revise
using Unitful, UnitfulAtomic
savedir=raw"D:\2021-SummerInternship-Results\15-2-noion-3"

const G = 1e-4u"T"

union!(LOAD_PATH,[raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules"])
using GenerateData

# sketch follow ↓
#no_sims, k_fixed, B_fixed = 11,1e-4u"bohr^-1", 0u"T"
#coltype, lmax, Bmin, Bmax, kmin, kmax = "4-4", 0, 0u"T", 0.1u"T",1e-4u"bohr^-1",1e-2u"bohr^-1"
#gen_diffB_constk_data(savedir,Bmin,Bmax,no_sims,k_fixed,coltype,lmax)
#gen_diffk_constB_data(savedir,kmin,kmax,no_sims,B_fixed,coltype,lmax)
using Revise
using Simulate, StateStructures, HalfIntegers, Plots
""" Scatter plot of elastic cross sections, at different B and constant k"""
function σ_vs_B_plot(el_or_ion::String,dir::String, Bmin::Unitful.BField, Bmax::Unitful.BField,
    k::Union{typeof(0u"bohr^-1"),typeof(0e0u"bohr^-1")}, coltype::String, lmax::Integer)
    @assert el_or_ion ∈ ["el","ion"] "el_or_ion flag ≠ el, ion" # sanity check
    data=load_data(dir, coltype, -(Inf)u"hartree", (Inf)u"hartree", Bmin, Bmax, lmax)
    if length(data)==0
        @warn "No data found matching arguments. Aborting plot."
        return nothing
    end
    σs = zeros(0)u"bohr^2"; Bs = zeros(0)u"T" # initialise
    for d in data
        correctks = findall(x->x≈k, d.k)
        correctks==[] && continue # skip if no channels have correct wavenumber
        for i in correctks # channel i has correct wavenumber
            el_or_ion=="el" ? push!(σs, d.σ_el[i,i]) : push!(σs, d.σ_ion[i])
            push!(Bs, d.B)
        end
    end
    if length(σs)==0
        @warn "Data found, but none with desired k. Aborting plot."
        return nothing
    end
    if el_or_ion=="el"
        plt=scatter(Bs./(1G), σs./(1u"bohr^2"), xlabel="B (G)", ylabel="σₑₗ (a₀²)",
        yscale=:log10, legend=false, title=coltype*" lmax=$lmax")
    else
        plt=scatter(Bs./(1G), σs./(1u"bohr^2"), xlabel="B (G)", ylabel="σᵢₒₙ (a₀²)",
        legend=false, title=coltype*" lmax=$lmax")
    end
    plt
end

function σ_vs_k_plot(el_or_ion::String,dir::String, kmin::Union{typeof(0u"bohr^-1"),typeof(0e0u"bohr^-1")},
    kmax::Union{typeof(0u"bohr^-1"),typeof(0e0u"bohr^-1")}, B::Unitful.BField,
    coltype::String, lmax::Integer)
    @assert el_or_ion ∈ ["el","ion"] "el_or_ion flag ≠ el, ion" # sanity check
    data=load_data(dir, coltype, -(Inf)u"hartree", (Inf)u"hartree", B, B, lmax)
    if length(data)==0
        @warn "No data found matching arguments. Aborting plot."
        return nothing
    end
    σs = zeros(0)u"bohr^2"; ks = zeros(0)u"bohr^-1" # initialise
    for d in data
        correctks = findall(x->kmin<=x<=kmax, d.k) # in desired k bounds
        correctks==[] && continue # skip if no channels have desirable wavenumber
        for i in correctks # channel i has desirable wavenumber
            el_or_ion=="el" ? push!(σs, d.σ_el[i,i]) : push!(σs, d.σ_ion[i])
            push!(ks, d.k[i])
        end
    end
    if length(σs)==0
        @warn "Data found, but none with desired k. Aborting plot."
        return nothing
    end
    if el_or_ion=="el"
        plt=scatter(ks./(1u"bohr^-1"), σs./(1u"bohr^2"), xlabel="k (a₀⁻¹)", ylabel="σₑₗ (a₀²)",
        yscale=:log10,xscale=:log10,legend=false, title=coltype*" lmax=$lmax")
    else
        plt=scatter(ks./(1u"bohr^-1"), σs./(1u"bohr^2"), xlabel="k (a₀⁻¹)", ylabel="σᵢₒₙ (a₀²)",
        xscale=:log10,legend=false, title=coltype*" lmax=$lmax")
    end
    plt
end

#σ_vs_B_plot("el",savedir, Bmin, Bmax, k_fixed, coltype, lmax)
