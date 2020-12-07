#= structs used to describe basis states used for the computation, plus utility
functions, such as lookup generator and symmetrised state evaluator.
Description last edited 7/12/2020=#

module StateStructures
using HalfIntegers

# Hyperfine basis states

# (Sⱼiⱼfₖmₖ), atomic quantum numbers
struct atom_nos
    S :: Int
    i :: HalfInteger
    f :: HalfInteger
    m :: HalfInteger
    atom_nos(S,i,f,m) = abs(S-i)<=f<=(S+i) && abs(m)<=f ? new(S,i,f,m) :
    error("Invalid |S i f m⟩ numbers supplied to atom_nos")
end

# asymmetric state, for evaluating bra-ket expressions
struct asym_αβlml_ket
    α :: atom_nos
    β :: atom_nos
    l :: Int
    ml :: Int
    asym_αβlmₗ_ket(α,β,l,ml) = abs(mₗ)<=l ? new(α,β,l,ml) : error("abs(mₗ)>l")
end

# symmetric state, used to represent scattering states
struct scat_αβlml_ket
    α :: atom_nos
    β :: atom_nos
    l :: Int
    ml :: Int
    asym_αβlmₗ_ket(α,β,l,ml) = α.i==β.i && abs(ml)<=l ? new(α,β,l,ml) : error("i₁≠i₂ || abs(mₗ)>l")
end

# Evaluator of matrix elements for symmetric states, by expanding into asym states
function sym_eval(op, bra::scat_αβlmₗ_ket, ket::scat_αβlmₗ_ket, p...)
    α_,β_ = bra.α,bra.β;        α,β = ket.α,ket.β # atomic numbers
    iₐ_, iᵦ_ = α_.i, β_.i;      iₐ, iᵦ = α.i, β.i # nuclear spins, for Bose/Fermi ±
    l_, ml_ = bra.l, bra.ml;    l,ml = ket.l, ket.ml # angular momenta, for symmetry factor and constructing asym states
    prefac = 1/(2*sqrt(1+(α==β ? 1 : 0))*sqrt(1+(α_==β_ ? 1 : 0))) # normalising prefactor
    t1 = op(asym_αβlml_ket(α_,β_,l_,ml_), asym_αβlmₗ_ket(α,β,l,ml), p...) # first term
    t2 = (-1)^(iₐ_+iᵦ_+l_)*op(asym_αβlml_ket(β_,α_,l_,ml_), asym_αβlmₗ_ket(α,β,l,ml), p...) # second term
    t3 = (-1)^(iₐ+iᵦ+l)*op(asym_αβlml_ket(α_,β_,l_,ml_), asym_αβlmₗ_ket(β,α,l,ml), p...) # third term
    t4 = (-1)^(iₐ_+iᵦ_+l_)*(-1)^(iₐ+iᵦ+l)*op(asym_αβlml_ket(β_,α_,l_,ml_), asym_αβlmₗ_ket(β,α,l,ml), p...) # fourth term
    return prefac*(t1+t2+t3+t4)
end

# Generate all symmetrised scattering states for a particular choice of i

# Generate all asymmetrised scattering states, setting α.i=3, β.i =4

end # module
