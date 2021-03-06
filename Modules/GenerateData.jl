#= Module with data generation functions, for calling in scripts=#

module GenerateData
export gen_diffE_data, gen_diffB_constk_data, gen_diffk_constB_data, load_data

using Revise
using BSON, Dates, Unitful, UnitfulAtomic, HalfIntegers, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Simulate, StateStructures, Interactions

""" Parameters """
const lhs=3e0u"bohr"; const mid=5e1u"bohr"; const rhs=1e4u"bohr"
const lhs2mid_spacing=1e1u"bohr"; const rhs2mid_spacing=1e1u"bohr"
const μ=0.5*4.002602u"u"
#const rrhs=1e3u"bohr"; const rhs2rrhs_spacing=2e2u"bohr" # test

rhs *= 1.1 # TODO investigating resonance location dependence

const G = 1e-4u"T" # Gauss unit of magnetic flux density

function create_params_str(coltype::String,ϵ::Unitful.Energy,B::Unitful.BField,lmax::Int)
    "_CT"*coltype*
    "_E"*string(ustrip(uconvert(u"hartree",ϵ)))*
    "_B"*string(ustrip(uconvert(u"T",B)))*
    "_lmax"*string(lmax)
end

# saves simulation output in <save_dir>, ϵ in Eh and B in T
function save_output(savedir::String, output::sim_output)
    wd=pwd() # current directory, to move back into at the end
    cd(savedir)
    date_str=string(Dates.today())
    params_str=create_params_str(output.coltype,output.ϵ,output.B,output.lmax)
    save_str=date_str*params_str*".sim"
    bson(save_str, data=output)
    cd(wd)
end


""" generate and save .Smat for different energies (in Eh), logarithmically spaced
    Inputs: log₁₀(austrip(Emin)), log₁₀(austrip(Emax)), n = number of different energies, lmax;
    B=0T magnetic field strength
    Outputs: / (saves files using save_output)"""
function gen_diffE_data(savedir::String,Emin_exp,Emax_exp,n::Integer,coltype::String,lmax::Integer;B=0u"T")
    existingfiles=readdir(savedir)
    ϵ_list=exp10.(LinRange(Emin_exp,Emax_exp,n))u"hartree" # energies
    println("lmax=$lmax, B=$B. Generating S_output for ϵ/Eh= ")
    Threads.@threads for ϵ in ϵ_list
        check_str=create_params_str(coltype,ϵ,B,lmax)
        if any(f->occursin(check_str,f), existingfiles)
            println("ϵ=$(ϵ/1u"hartree")Eh, B=$(B/1e-4u"T")G was already calculated.")
            continue
        end
        # datum was not already calculated, proceeding to simulate.
        println("Simulating CT=$coltype, lmax=$lmax, ϵ=$(ϵ/1u"hartree")Eh, B=$B on thread $(Threads.threadid()).")
        output=sim(coltype,lmax,ϵ,B,lhs,mid,rhs,lhs2mid_spacing,rhs2mid_spacing)
        save_output(savedir, output)
    end
end

"""Constant k, different B data generation.
    Inputs: savedir, Bmin, Bmax, n, k, coltype, lmax"""
function gen_diffB_constk_data(savedir::String, Bmin::Unitful.BField,Bmax::Unitful.BField,n::Integer,
    k::Union{typeof(0e0u"bohr^-1"),typeof(0u"bohr^-1")}, coltype::String, lmax::Int)
    println("Running gen_diffB_constk_data, $coltype collisions with lmax of $lmax.")
    println("Desired k=$k, iterating over $n B-fields from $Bmin to $Bmax.")
    existingfiles=readdir(savedir)
    Bs=LinRange(Bmin,Bmax,n) # different magnetic fields to iterate over
	unq_iden_lookup = let iden_lookup = αβlml_lookup_generator(coltype,"iden",lmax)
		@assert all(x->typeof(x)==scat_αβlml_ket,iden_lookup) || all(x->typeof(x)==asym_αβlml_ket,iden_lookup) "I've wrongly assumed lookup only has one type of ket" # Sanity checl
		if length(iden_lookup)==0
			iden_lookup # return empty vector (3-4 case)
		elseif typeof(iden_lookup[1])==scat_αβlml_ket
			unique((ket->scat_αβlml_ket(ket.α,ket.β,0,0)).(iden_lookup)) # possible unphysical ket, just used to calc H_hfs+H_zee
		else # asym_ket_case
			unique((ket->asym_αβlml_ket(ket.α,ket.β,0,0)).(iden_lookup)) # possible unphysical ket, just used to calc H_hfs+H_zee
		end
	end
	iden_N=length(unq_iden_lookup)
	unq_diff_lookup = let diff_lookup = αβlml_lookup_generator(coltype,"diff",lmax)
		@assert all(x->typeof(x)==scat_αβlml_ket,diff_lookup) || all(x->typeof(x)==asym_αβlml_ket,diff_lookup) "I've wrongly assumed lookup only has one type of ket" # Sanity checl
		if length(diff_lookup)==0
			diff_lookup # return empty vector (3-4 case)
		elseif typeof(diff_lookup[1])==scat_αβlml_ket
			unique((ket->scat_αβlml_ket(ket.α,ket.β,0,0)).(diff_lookup)) # possible unphysical ket, just used to calc H_hfs+H_zee
		else # asym_ket case
			unique((ket->asym_αβlml_ket(ket.α,ket.β,0,0)).(diff_lookup)) # possible unphysical ket, just used to calc H_hfs+H_zee
		end
	end
	diff_N=length(unq_diff_lookup)
	lookup=vcat(unq_iden_lookup,unq_diff_lookup); N=length(lookup)
    Threads.@threads for B in Bs # iterate over different B fields
		println("Calculating B=$B on thread $(Threads.threadid())")
        H∞ = Matrix{Unitful.Energy}(zeros(N,N)u"hartree") # intialise
        for i=1:N, j=1:N
            bra, ket = lookup[i], lookup[j]
            H∞[i,j]=αβlml_eval(H_zee, bra, ket, B)+αβlml_eval(H_hfs, bra, ket)
        end
        minV∞=(eigen(austrip.(H∞)).values[1])u"hartree" # lowest energy channel is physically relevant
		ϵ = (v->auconvert(v+1u"ħ^2"*k^2/(2*μ)))(minV∞) # energy to make low-energy channel desired wavenumber
        @assert dimension(ϵ)==dimension(0u"hartree") "ϵ not an energy" # sanity check
        check_str=create_params_str(coltype,ϵ,B,lmax)
        if any(f->occursin(check_str,f), existingfiles) # check if data already exists
            println("ϵ=$(ϵ/1u"hartree")Eh, B=$(B/1e-4u"T")G was already calculated.")
            continue
        end
		println("Simulating CT=$coltype, lmax=$lmax, ϵ=$(ϵ/1u"hartree")Eh, B=$(B/G) G on thread $(Threads.threadid()).")
        try # try/catch in case of bug with this simulation
            output = sim(coltype, lmax, ϵ, B, lhs, mid, rhs,
            lhs2mid_spacing, rhs2mid_spacing)
            save_output(savedir, output) # save
        catch e
            @warn "Error occured running sim. Skipped to next simulation."
			@show e
        end # try/catch
    end # B
end # function

""" Constant B, different k data generation (logarithmically spaced wavenumbers)
    Inputs: savedir, kmin (exponent), kmax (exponent), n, B, coltype, lmax"""
function gen_diffk_constB_data(savedir::String, kmin::Number, kmax::Number, n::Integer,
    B::Unitful.BField, coltype::String, lmax::Integer)
    existingfiles=readdir(savedir)
    ks=exp10.(LinRange(kmin,kmax,n))u"bohr^-1" # different wavenumbers to iterate over
    unq_iden_lookup = let iden_lookup = αβlml_lookup_generator(coltype,"iden",lmax)
		@assert all(x->typeof(x)==scat_αβlml_ket,iden_lookup) || all(x->typeof(x)==asym_αβlml_ket,iden_lookup) "I've wrongly assumed lookup only has one type of ket" # Sanity checl
		if length(iden_lookup)==0
			iden_lookup # return empty vector (3-4 case)
		elseif typeof(iden_lookup[1])==scat_αβlml_ket
			unique((ket->scat_αβlml_ket(ket.α,ket.β,0,0)).(iden_lookup)) # possible unphysical ket, just used to calc H_hfs+H_zee
		else
			unique((ket->asym_αβlml_ket(ket.α,ket.β,0,0)).(iden_lookup)) # possible unphysical ket, just used to calc H_hfs+H_zee
		end
	end
	iden_N=length(unq_iden_lookup)
	unq_diff_lookup = let diff_lookup = αβlml_lookup_generator(coltype,"diff",lmax)
		@assert all(x->typeof(x)==scat_αβlml_ket,diff_lookup) || all(x->typeof(x)==asym_αβlml_ket,diff_lookup) "I've wrongly assumed lookup only has one type of ket" # Sanity checl
		if length(diff_lookup)==0
			diff_lookup # return empty vector (3-4 case)
		elseif typeof(diff_lookup[1])==scat_αβlml_ket
			unique((ket->scat_αβlml_ket(ket.α,ket.β,0,0)).(diff_lookup)) # possible unphysical ket, just used to calc H_hfs+H_zee
		else
			unique((ket->asym_αβlml_ket(ket.α,ket.β,0,0)).(diff_lookup)) # possible unphysical ket, just used to calc H_hfs+H_zee
		end
	end
	diff_N=length(unq_diff_lookup)
	lookup=vcat(unq_iden_lookup,unq_diff_lookup); N=length(lookup)
    H∞ = Matrix{Unitful.Energy}(zeros(N,N)u"hartree") # intialise
    for j=1:N, i=1:N
        bra, ket = lookup[i], lookup[j]
        H∞[i,j]=αβlml_eval(H_zee, bra, ket, B)+αβlml_eval(H_hfs, bra, ket)
    end
    minV∞=(eigen(austrip.(H∞)).values[1])u"hartree" # lowest energy chanel is physically relevant
    Threads.@threads for k in ks # iterate over desired wavenumbers
		println("Calculating k=$k on thread $(Threads.threadid())")
		ϵ = (v->auconvert(v+1u"ħ^2"*k^2/(2*μ)))(minV∞) # energy to make low-energy channel desired wavenumber
		@assert dimension(ϵ)==dimension(0u"hartree") "ϵ not an energy" # sanity check
        check_str=create_params_str(coltype,ϵ,B,lmax)
        if any(f->occursin(check_str,f), existingfiles) # check if data is already calculated
            println("ϵ=$(ϵ/1u"hartree")Eh, B=$(B/1e-4u"T")G was already calculated.")
            continue
        end
		println("Simulating CT=$coltype, lmax=$lmax, ϵ=$(ϵ/1u"hartree")Eh, B=$(B/G) G on thread $(Threads.threadid()).")
        try # try/catch in case of bug with this iteration.
            output = sim(coltype, lmax, ϵ, B, lhs, mid, rhs,
            lhs2mid_spacing, rhs2mid_spacing)
            save_output(savedir, output) # save
        catch e
            @warn "Error occured running sim(), for ϵ=$ϵ, B=$B. Skipped to next simulation."
			@show e
        end # try/catch
    end # k
end # function

###################################Load data####################################
""" load .sim data with filenames fitting bounds on coltype, ϵ, B, lmax
    Inputs: savedir, coltype, Emin/max ~[E], Bmin/max ~ [Tesla], lmax
    Output: array of outputs found from filename filtering"""
function load_data(dir,coltype::String,Emin::Unitful.Energy,Emax::Unitful.Energy,
    Bmin::Unitful.BField,Bmax::Unitful.BField,lmax::Integer)
    @assert coltype in ["3-3","3-4","4-4"] "invalid coltype"
    # directory change
    wd=pwd() # current working directory
    cd(dir)
    # initialise output
    output=Vector{sim_output}([])
    for f in readdir()
		occursin(".sim",f) || continue
        occursin(coltype,f) || continue
        fE=let
            SIndex=findfirst(isequal('E'),f)+1
            EIndex=findall(isequal('_'),f)[3]-1
            str=f[SIndex:EIndex]
            parse(Float64,str)u"hartree"
        end
        fB=let
            SIndex=findfirst(isequal('B'),f)+1
            EIndex=findall(isequal('_'),f)[4]-1
            str=f[SIndex:EIndex]
            parse(Float64,str)u"T"
        end
        flmax=let
            SIndex=findfirst("lmax",f)[end]+1
            EIndex=findfirst(".sim",f)[1]-1
            str=f[SIndex:EIndex]
            parse(Int,str)
        end
        (Emin<=fE<=Emax)&&(Bmin<=fB<=Bmax)&&(flmax==lmax) || continue
        # within parameters, loading:
        datum=BSON.load(f)[:data]
        push!(output,datum)
    end
    cd(wd)
    # sort data by energy
    sort!(output, by=x->x.ϵ)
    output
end

end # module
