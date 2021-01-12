savedir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Summer Internship\Results\12-1-test2"


"""
Plot diagonal elements of pairwise σ (i.e. elastic cross sections) vs wavenumber
Inputs: kmin~[L]⁻¹, kmax~[L]⁻¹ for plot x axis, B~[BField],lmax;
    Emin~[E]=0Eₕ, Emax~[E]=1Eₕ for finding appropriate data
Output: plot of elastic cross sections vs k"""
function diffk_gam_plot(kmin::typeof(0e0u"bohr^-1"),kmax::typeof(0e0u"bohr^-1"),
    B::Unitful.BField,lmax::Int)
    # load all data with correct B
    datas=load_data("gam",-(Inf)u"hartree",(Inf)u"hartree",B,B,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.ϵ)) # sort by increasing energy
    unq = unqkets(reverse(datas)) # reverse to find uniques starting with high energy data
    pltdata=[] # array to store ([k],[σ}) pairs for the different γ
    pltlabel=label_from_lookup(unq)
    for γ in unq
        datatuple::typeof(([0.0u"bohr^-1"],[0.0u"bohr^2"]))=([],[])
        for d in datas # already sorted datas by energy
            if γ in d.γ_lookup
                dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
                k=k∞(γ,d.ϵ,B) # the asymptotic wavenumber for this particular data
                imag(k)==0.0u"bohr^-1" || continue # don't store if the wavenumber is complex⟺channel closed
                kmin <= real(k) <= kmax || continue # don't store if the wavenumber is out of bounds
                push!(datatuple[1],k) # store wavenumber
                push!(datatuple[2],d.σ[dγ_index,dγ_index]) # store cross section
            end
        end
        push!(pltdata,datatuple)
    end
    println("Minimum values are:")
    println("S, mS = ",(x->"$(x.S), $(x.mS)").(unq))
    println(austrip.([pltdata[i][2][1] for i=1:length(pltdata)]))
    # plot first γ_ket
    plot(austrip.(pltdata[1][1]),austrip.(pltdata[1][2]),xlabel="Wavenumber (a₀⁻¹)", xscale=:log10,
    ylabel="σ (a₀²)",yscale=:log10, minorticks=true, label=pltlabel[1], legend=:outertopright,
    title="B=$B, lmax=$lmax",
    left_margin=5mm,bottom_margin=5mm, top_margin=5mm,
    linewidth=2, grid=false)
    if length(pltdata)>1 # plot rest of the γ_kets
        for i in 2:length(pltdata)
            plot!(austrip.(pltdata[i][1]),austrip.(pltdata[i][2]),label=pltlabel[i])
        end
    end
    hline!([4*pi*austrip((7.54u"nm")^2)],label="S=2 4πa²")
end
