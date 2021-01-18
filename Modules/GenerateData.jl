#= Module with data generation functions, for calling in scripts=#

module GenerateData
export gen_diffE_data, gen_diffB_constk_data, gen_diffk_constB_data, load_data

using Revise
using BSON, Dates, Unitful, UnitfulAtomic, HalfIntegers, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Simulate, StateStructures, Interactions

""" Parameters """
const lhs=3e0u"bohr"; const mid=5e1u"bohr"; const rhs=1e3u"bohr"
const lhs2mid_spacing=1e1u"bohr"; const rhs2mid_spacing=2e9u"bohr" # no mid←rhs orthog
const μ=0.5*4.002602u"u"

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
    iden_lookup=unique(αβlml_lookup_generator(coltype,"iden",lmax))
    iden_N=length(iden_lookup)
    diff_lookup=unique(αβlml_lookup_generator(coltype,"diff",lmax))
    diff_N=length(diff_lookup)
    for (lookup, N) in [(iden_lookup, iden_N), (diff_lookup, diff_N)] # iterating over different groups of channels
		println("Iterating over iden_ or diff_ channels")
		N==0 && continue # skip empty group of channels
        Threads.@threads for B in Bs # iterate over different B fields
			println("Calculating B=$B on thread $(Threads.threadid())")
            H∞ = Matrix{Unitful.Energy}(zeros(N,N)u"hartree") # intialise
            for i=1:N, j=1:N
                bra, ket = lookup[i], lookup[j]
                H∞[i,j]=αβlml_eval(H_zee, bra, ket, B)+αβlml_eval(H_hfs, bra, ket)
            end
            D∞ = unique(eigen(austrip.(H∞)).values)u"hartree" # unique, so as to avoid repeated calculation
			minV∞=minimum(D∞) # lowest energy chanel is physically relevant
			ϵ = (v->auconvert(v+1u"ħ^2"*k^2/(2*μ)))(minV∞) # energy to make low-energy channel desired wavenumber
            @assert dimension(ϵ)==dimension(0u"hartree") "ϵ not an energy" # sanity check
            check_str=create_params_str(coltype,ϵ,B,lmax)
            if any(f->occursin(check_str,f), existingfiles) # check if data already exists
                println("ϵ=$(ϵ/1u"hartree")Eh, B=$(B/1e-4u"T")G was already calculated.")
                continue
            end
			println("Simulating CT=$coltype, lmax=$lmax, ϵ=$(ϵ/1u"hartree")Eh, B=$B on thread $(Threads.threadid()).")
            try # try/catch in case of bug with this simulation
                output = sim(coltype, lmax, ϵ, B, lhs, mid, rhs,
                lhs2mid_spacing, rhs2mid_spacing)
                save_output(savedir, output) # save
            catch e
                @warn "Error occured running sim. Skipped to next simulation."
				@show e
            end # try/catch
        end # B
    end # (lookup, N)
end # function

""" Constant B, different k data generation (logarithmically spaced wavenumbers)
    Inputs: savedir, kmin (exponent), kmax (exponent), n, B, coltype, lmax"""
function gen_diffk_constB_data(savedir::String, kmin::Number, kmax::Number, n::Integer,
    B::Unitful.BField, coltype::String, lmax::Integer)
    existingfiles=readdir(savedir)
    ks=exp10.(LinRange(kmin,kmax,n))u"bohr^-1" # different wavenumbers to iterate over
    iden_lookup=unique(αβlml_lookup_generator(coltype,"iden",lmax))
    iden_N=length(iden_lookup)
    diff_lookup=unique(αβlml_lookup_generator(coltype,"diff",lmax))
    diff_N=length(diff_lookup)
    for (lookup, N) in [(iden_lookup, iden_N), (diff_lookup, diff_N)] # iterating over different groups of channels
		println("Iterating over iden_ or diff_ channels")
		N==0 && continue # skip empty group of channels
        H∞ = Matrix{Unitful.Energy}(zeros(N,N)u"hartree") # intialise
        for i=1:N, j=1:N
            bra, ket = lookup[i], lookup[j]
            H∞[i,j]=αβlml_eval(H_zee, bra, ket, B)+αβlml_eval(H_hfs, bra, ket)
        end
        D∞ = unique(eigen(austrip.(H∞)).values)u"hartree" # unique, so as to avoid repeated calculation
		minV∞=minimum(D∞) # lowest energy chanel is physically relevant
        Threads.@threads for k in ks # iterate over desired wavenumbers
			println("Calculating k=$k on thread $(Threads.threadid())")
			ϵ = (v->auconvert(v+1u"ħ^2"*k^2/(2*μ)))(minV∞) # energy to make low-energy channel desired wavenumber
			@assert dimension(ϵ)==dimension(0u"hartree") "ϵ not an energy" # sanity check
            check_str=create_params_str(coltype,ϵ,B,lmax)
            if any(f->occursin(check_str,f), existingfiles) # check if data is already calculated
                println("ϵ=$(ϵ/1u"hartree")Eh, B=$(B/1e-4u"T")G was already calculated.")
                continue
            end
			println("Simulating CT=$coltype, lmax=$lmax, ϵ=$(ϵ/1u"hartree")Eh, B=$B on thread $(Threads.threadid()).")
            try # try/catch in case of bug with this iteration.
                output = sim(coltype, lmax, ϵ, B, lhs, mid, rhs,
                lhs2mid_spacing, rhs2mid_spacing)
                save_output(savedir, output) # save
            catch e
                @warn "Error occured running sim(), for ϵ=$ϵ, B=$B. Skipped to next simulation."
				@show e
            end # try/catch
        end # k
    end # (lookup, N)
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
