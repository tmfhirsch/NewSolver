#= Interactions as functions of |Φₐ⟩=|S₁S₂SmSlml⟩ states, radial dist. R and μ=#

module Interactions
export H_rot, H_el_coeffs, H_el_radial, H_sd_coeffs, H_sd_radial, H_zee, H_hfs, Γ_GMS_coeffs, Γ_GMS_radial

using Unitful, UnitfulAtomic

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\NewSolver\Modules")
using StateStructures

# rotational
include("./Interactions/H_rot.jl")

# electronic
using Potentials
using WignerSymbols
include("./Interactions/H_el.jl")

# spin-dipole
using WignerSymbols, Wigner9j
include("./Interactions/H_sd.jl")

# zeeman
include("./Interactions/H_zee.jl")

# hyperfine
using HalfIntegers
include("./Interactions/H_hfs.jl")

# ionisation
include("./Interactions/Ionisation.jl")

end # module
