using Revise
using Unitful, UnitfulAtomic
#savedir=raw"D:\2021-SummerInternship-Results\final\aa-noion"
savedir=raw"C:\Users\hirsc\Documents\Summer-Internship-Results\test-rightmass"
iondir=raw"C:\Users\hirsc\Documents\Summer-Internship-Results\wrongmass\aa-ion"
noiondir=raw"C:\Users\hirsc\Documents\Summer-Internship-Results\wrongmass\aa-noion"

const G = 1e-4u"T"

union!(LOAD_PATH,[raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules"])
using GenerateData, Simulate, StateStructures, HalfIntegers
using Plots
#using PyCall
#pygui(:tk)
#using PyPlot

# sketch follow ↓
#no_sims, k_fixed, B_fixed = 11,1e-4u"bohr^-1", 0u"T"
#coltype, lmax, Bmin, Bmax, kmin, kmax = "4-4", 0, 0u"T", 0.1u"T",1e-4u"bohr^-1",1e-2u"bohr^-1"
#gen_diffB_constk_data(savedir,Bmin,Bmax,no_sims,k_fixed,coltype,lmax)
#gen_diffk_constB_data(savedir,kmin,kmax,no_sims,B_fixed,coltype,lmax)
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
        yscale=:log10, legend=false, title=coltype*" lmax=$lmax",
        markershape=:x, markerstrokewidth=1, markersize=2)
    else
        plt=scatter(Bs./(1G), σs./(1u"bohr^2"), xlabel="B (G)", ylabel="σᵢₒₙ (a₀²)",
        legend=false, title=coltype*" lmax=$lmax",
        markershape=:x, markerstrokewidth=1, markersize=2)
    end
    #@show Bs[findmin(σs)[2]]
    println("Data from $savedir")
    plt
end

function ion_noion_plot(ion_dir::String,noion_dir::String, Bmin::Unitful.BField, Bmax::Unitful.BField,
    k::Union{typeof(0u"bohr^-1"),typeof(0e0u"bohr^-1")}, coltype::String, lmax::Integer)
    ion_data=load_data(ion_dir,coltype, -(Inf)u"hartree", (Inf)u"hartree", Bmin,Bmax, lmax)
    if length(ion_data)==0
        @warn "No ion_data found matching arguments. Aborting plot."
        return nothing
    end
    sort!(ion_data,by=x->x.B)
    noion_data=load_data(noion_dir,coltype, -(Inf)u"hartree", (Inf)u"hartree", Bmin,Bmax, lmax)
    if length(noion_data)==0
        @warn "No noion_data found matching arguments. Aborting plot."
        return nothing
    end
    sort!(noion_data,by=x->x.B)
    ion_σs, ion_Bs = zeros(0)u"bohr^2", zeros(0)u"T"
    for d in ion_data
        correctks = findall(x->x≈k, d.k)
        correctks==[] && continue # skip this data if no correct wavenumbers
        for i in correctks # channel i has correct wavenumber
            push!(ion_σs, d.σ_el[i,i])
            push!(ion_Bs, d.B)
        end
    end
    noion_σs, noion_Bs = zeros(0)u"bohr^2", zeros(0)u"T"
    for d in noion_data
        correctks = findall(x->x≈k, d.k)
        correctks==[] && continue # skip this data if no correct wavenumbers
        for i in correctks # channel i has correct wavenumber
            push!(noion_σs, d.σ_el[i,i])
            push!(noion_Bs, d.B)
        end
    end
    plt=plot(ion_Bs./G, ion_σs./1u"bohr^2", label="With ionisation", grid=false,
    thickness_scaling=1.5, linewidth=2, framestyle=:box,
    xlabel="B (G)", ylabel="σ (a₀²)")
    plot!(plt, noion_Bs./G, noion_σs./1u"bohr^2", label="Without ionisation",
    linewidth=2)
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
        yscale=:log10,xscale=:log10,legend=false, title=coltype*" lmax=$lmax",
        markershape=:x, markerstrokewidth=1, markersize=2)
    else
        plt=scatter(ks./(1u"bohr^-1"), σs./(1u"bohr^2"), xlabel="k (a₀⁻¹)", ylabel="σᵢₒₙ (a₀²)",
        xscale=:log10,legend=false, title=coltype*" lmax=$lmax",
        markershape=:x, markerstrokewidth=1, markersize=2)
    end
    println("Data from $savedir")
    plt
end


using DelimitedFiles
csvdir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Summer Internship\March\2-3\csvs"
""" Saves σ vs B data in two columns in a .csv, for use in python"""
function csv_maker(simdir::String, el_or_ion::String,
    Bmin::Unitful.BField, Bmax::Unitful.BField, k::Union{typeof(0u"bohr^-1"),typeof(0e0u"bohr^-1")},
    coltype::String, lmax::Integer,
    csvdir::String, filename::String)
    @assert el_or_ion ∈ ["el","ion"] "el_or_ion flag ≠ el, ion" # sanity check
    data=load_data(simdir, coltype, -(Inf)u"hartree", (Inf)u"hartree", Bmin, Bmax, lmax)
    if length(data)==0
        @warn "No data found matching arguments. Aborting plot."
        return nothing
    end
    sort!(data,by=x->x.B) # increasing magnetic field order
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
    prevdir=pwd()
    cd(csvdir)
    println("Saving $filename in $csvdir")
    writedlm(filename, ["B (G)" "cs (a_0^2)";
                        Bs./G σs./(1u"bohr^2")], ',',
                        header=true)
    cd(prevdir)
end
