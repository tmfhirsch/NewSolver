#= Module with data generation functions, for calling in scripts=#

module GenerateData
export gen_diffE_data, gen_diffB_constk_data

using Revise
using BSON, Dates, Unitful, UnitfulAtomic, HalfIntegers, LinearAlgebra
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using Simulate, StateStructures, Interactions

""" Parameters """
const lhs=3e0u"bohr"; const mid=5e1u"bohr"; const rhs=2e2u"bohr"; const rrhs=1e3u"bohr"
const lhs2mid_spacing=1e1u"bohr"; const rhs2mid_spacing=1e9u"bohr" # no mid←rhs orthog
const rhs2rrhs_spacing=2e2u"bohr"; const μ=0.5*4.002602u"u"

const G = 1e-4u"T" # Gauss unit of magnetic flux density

# saves simulation output in <save_dir>, ϵ in Eh and B in T
function save_output(savedir::String, output::sim_output)
    wd=pwd() # current directory, to move back into at the end
    cd(savedir)
    date_str=string(Dates.today())
    params_str="_CT"*output.coltype*
    "_E"*string(ustrip(uconvert(u"hartree",output.ϵ)))*
    "_B"*string(ustrip(uconvert(u"T",output.B)))*
    "_lmax"*string(output.lmax)
    save_str=date_str*params_str*".sim"
    bson(save_str, data=output)
    cd(wd)
end


""" generate and save .Smat for different energies (in Eh), logarithmically spaced
    Inputs: log₁₀(austrip(Emin)), log₁₀(austrip(Emax)), n = number of different energies, lmax;
    B=0T magnetic field strength
    Outputs: / (saves files using save_output)"""
function gen_diffE_data(savedir::String,Emin_exp,Emax_exp,n::Integer,coltype::String,lmax::Integer;B=0u"T")
    ϵ_list=exp10.(LinRange(Emin_exp,Emax_exp,n))u"hartree" # energies
    println("lmax=$lmax, B=$B. Generating S_output for ϵ/Eh= ")
    for ϵ in ϵ_list
        println("$(austrip(ϵ)), ")
        output=sim(coltype,lmax,ϵ,B,lhs,mid,rhs,rrhs,lhs2mid_spacing,rhs2mid_spacing,rhs2rrhs_spacing)
        save_output(savedir, output)
    end
end

"""Constant k, different B data generation"""
function gen_diffB_constk_data(savedir::String, Bmin::Unitful.BField,Bmax::Unitful.BField,n::Integer,
    k::Union{typeof(0e0u"bohr^-1"),typeof(0u"bohr^-1")}, coltype, lmax::Int)
    println("Running gen_diffB_constk_data, $coltype collisions with lmax of $lmax.")
    println("Desired k=$k, iterating over $n B-fields from $Bmin to $Bmax.")
    Bs=LinRange(Bmin,Bmax,n) # different magnetic fields to iterate over
    iden_lookup=unique(αβlml_lookup_generator(coltype,"iden",lmax))
    iden_N=length(iden_lookup)
    diff_lookup=unique(αβlml_lookup_generator(coltype,"diff",lmax))
    diff_N=length(diff_lookup)
    for (lookup, N) in [(iden_lookup, iden_N), (diff_lookup, diff_N)] # iterating over different groups of channels
        for B in Bs # iterate over different B fields
            H∞ = Matrix{Unitful.Energy}(zeros(N,N)u"hartree")
            for i=1:N, j=1:N
                bra, ket = lookup[i], lookup[j]
                H∞[i,j]=αβlml_eval(H_zee, bra, ket, B)+αβlml_eval(H_hfs, bra, ket)
            end
            D∞ = unique(trunc.(eigen(austrip.(H∞)).values,digits=8))u"hartree" # unique, so as to avoid repeated calculation
            ϵs = (v->auconvert(v+1u"ħ^2"*k^2/(2*μ))).(D∞)
            @assert all(x->dimension(x)==dimension(0u"hartree"),ϵs) "ϵs not a list of energies" # sanity check
            for ϵ in ϵs # iterating over different channels
                @info ϵ, B
                try # try/catch in case of bug with this iteration.
                    output = sim(coltype, lmax, ϵ, B, lhs, mid, rhs, rrhs,
                    lhs2mid_spacing, rhs2mid_spacing, rhs2rrhs_spacing)
                    save_output(savedir, output) # save
                catch
                    @warn "Error occured running sim(), for ϵ=$ϵ, B=$B. Skipped to next simulation."
                end
            end
        end # B
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
