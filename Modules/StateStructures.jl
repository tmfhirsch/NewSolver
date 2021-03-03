#= structs used to describe basis states used for the computation, plus utility
functions, such as lookup generator and symmetrised state evaluator.
Description last edited 7/12/2020=#

module StateStructures
using HalfIntegers

export asym_αβlml_ket, scat_αβlml_ket, αβlml_eval, αβlml_lookup_generator, atom_nos
export test_ket, test_lookup_generator
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
    asym_αβlml_ket(α,β,l,ml) = abs(ml)<=l ? new(α,β,l,ml) : error("abs(mₗ)>l")
end

# symmetric state, used to represent scattering states
struct scat_αβlml_ket
    α :: atom_nos
    β :: atom_nos
    l :: Int
    ml :: Int
    scat_αβlml_ket(α,β,l,ml) = α.i==β.i && abs(ml)<=l ? new(α,β,l,ml) : error("i₁≠i₂ || abs(mₗ)>l")
end


# evaluator of matrix elements for symmetric or asymmetric states
function αβlml_eval(op,
    bra::Union{scat_αβlml_ket,asym_αβlml_ket}, ket::Union{scat_αβlml_ket,asym_αβlml_ket},
    p...)
    # sanity check
    @assert typeof(bra)==typeof(ket) "Mismatched asymmetric & symmetric bra/ket"
    # simple case of asymmetric states: passes straight through to the operator
    if typeof(bra)==asym_αβlml_ket
        return op(bra,ket,p...)
    end
    # symmetric state case: expands into asymmetric states
    α_,β_ = bra.α,bra.β;        α,β = ket.α,ket.β # atomic numbers
    iₐ_, iᵦ_ = α_.i, β_.i;      iₐ, iᵦ = α.i, β.i # nuclear spins, for Bose/Fermi ±
    l_, ml_ = bra.l, bra.ml;    l,ml = ket.l, ket.ml # angular momenta, for symmetry factor and constructing asym states
    prefac = 1/(2*sqrt(1+(α==β ? 1 : 0))*sqrt(1+(α_==β_ ? 1 : 0))) # normalising prefactor
    t1 = op(asym_αβlml_ket(α_,β_,l_,ml_), asym_αβlml_ket(α,β,l,ml), p...) # first term
    t2 = (-1)^(iₐ_+iᵦ_+l_) .* op(asym_αβlml_ket(β_,α_,l_,ml_), asym_αβlml_ket(α,β,l,ml), p...) # second term
    t3 = (-1)^(iₐ+iᵦ+l) .* op(asym_αβlml_ket(α_,β_,l_,ml_), asym_αβlml_ket(β,α,l,ml), p...) # third term
    t4 = (-1)^(iₐ_+iᵦ_+l_) .* (-1)^(iₐ+iᵦ+l) .* op(asym_αβlml_ket(β_,α_,l_,ml_), asym_αβlml_ket(β,α,l,ml), p...) # fourth term
    return prefac .* (t1 .+ t2 .+ t3 .+ t4) # broadcast over arrays of linear coupling coefficients
end

# Generate hyperfine scattering states for 3-3, 3-4, and 4-4 collisions
function αβlml_lookup_generator(coltype::String, sameness::String, lmax::Int)
    # sanity checks
    @assert coltype ∈ ["3-3","3-4","4-4"] "coltype ∉ [3-3,3-4,4-4]"
    @assert sameness ∈ ["iden","diff","all"] "sameness ∉ [iden, diff, all]"
    @assert lmax >= 0 "lmax < 0"
    Sₐ,Sᵦ = 1, 1 # He* => atomic spin 1
    if coltype=="3-4" # asymmetric case
        lookup=Array{asym_αβlml_ket,1}([]) # initialise lookup
        iₐ,iᵦ = half(1), 0 # α is ³He, β is ⁴He
        for l=0:lmax
            for ml=-l:l
                for fₐ=abs(Sₐ-iₐ):(Sₐ+iₐ) # iterate atom α
                    for mₐ=(-fₐ):fₐ
                        α=atom_nos(Sₐ,iₐ,fₐ,mₐ)
                        for fᵦ=abs(Sᵦ-iᵦ):(Sᵦ+iᵦ) # iterate atom β
                            for mᵦ=(-fᵦ):fᵦ
                                β=atom_nos(Sᵦ,iᵦ,fᵦ,mᵦ)
                                # only store |α==β⟩ or |α!=β⟩ as desired
                                sameness=="all" || sameness=="diff" || continue # simplifies for 3-4
                                state=asym_αβlml_ket(α,β,l,ml)
                                push!(lookup,state)
                            end # mᵦ
                        end # fᵦ
                    end # mₐ
                end # fₐ
            end # ml
        end # l
    elseif coltype=="4-4" # symmetric ⁴He* case
        lookup=Array{scat_αβlml_ket,1}([])
        iₐ,iᵦ=0,0 # Bosonic ⁴He*
        for l=0:lmax
            for ml=-l:l
                for fₐ=abs(Sₐ-iₐ):(Sₐ+iₐ) # iterate atom α
                    for mₐ=(-fₐ):fₐ
                        α=atom_nos(Sₐ,iₐ,fₐ,mₐ)
                        for fᵦ=abs(Sᵦ-iᵦ):(Sᵦ+iᵦ) # iterate atom β
                            for mᵦ=-fᵦ:fᵦ
                                (fₐ,mₐ)<=(fᵦ,mᵦ) || continue # Double counting check: if fₐ == fᵦ, require mₐ <= mᵦ
                                β=atom_nos(Sᵦ,iᵦ,fᵦ,mᵦ)
                                α==β && mod(iₐ+iᵦ+l,2)==1 && continue # symmetrisation condition
                                # only store |α==β⟩ or |α!=β⟩ as desired
                                sameness=="all" || (sameness=="iden" && α==β) || (sameness=="diff" && α!=β) || continue
                                state=scat_αβlml_ket(α,β,l,ml)
                                push!(lookup,state)
                            end #mᵦ
                        end #fᵦ2
                    end #mₐ
                end #fₐ
            end # ml
        end #l
        # self-check for double counting
        for ket in lookup
            ket.α==ket.β && continue
            ketflip = scat_αβlml_ket(ket.β,ket.α,ket.l,ket.ml)
            @assert ketflip ∉ lookup "Double counted $ket"
        end
    elseif coltype=="3-3" #Symmetric ³He* case
        lookup=Array{scat_αβlml_ket,1}([])
        iₐ,iᵦ=half(1),half(1) # Bosonic ⁴He*
        for l=0:lmax
            for ml=-l:l
                for fₐ=abs(Sₐ-iₐ):(Sₐ+iₐ) # iterate atom α
                    for mₐ=(-fₐ):fₐ
                        α=atom_nos(Sₐ,iₐ,fₐ,mₐ)
                        for fᵦ=abs(Sᵦ-iᵦ):(Sᵦ+iᵦ) # iterate atom β
                            for mᵦ=-fᵦ:fᵦ
                                (fₐ,mₐ)<=(fᵦ,mᵦ) || continue # Double counting check: if fₐ == fᵦ, require mₐ <= mᵦ
                                β=atom_nos(Sᵦ,iᵦ,fᵦ,mᵦ)
                                α==β && mod(iₐ+iᵦ+l,2)==1 && continue # symmetrisation condition
                                # only store |α==β⟩ or |α!=β⟩ as desired
                                sameness=="all" || (sameness=="iden" && α==β) || (sameness=="diff" && α!=β) || continue
                                state=scat_αβlml_ket(α,β,l,ml)
                                push!(lookup,state)
                            end #mᵦ
                        end #fᵦ
                    end #mₐ
                end #fₐ
            end # ml
        end #l
        # self-check for double counting
        for ket in lookup
            ket.α==ket.β && continue
            ketflip = scat_αβlml_ket(ket.β,ket.α,ket.l,ket.ml)
            @assert ketflip ∉ lookup "Double counted $ket"
        end
    end # coltype 'if' statement
    return lookup
end

# Testing
struct test_ket
    x::Bool # identifies as closed or open. TRUE = OPEN, FALSE = CLOSED
    l::Int
    ml::Int
    test_ket(x) = new(x,0,0)
end

function test_lookup_generator()
    lookup=Vector{test_ket}([])
    push!(lookup,test_ket(true))
    push!(lookup,test_ket(false))
    lookup
end

αβlml_eval(op, bra::test_ket, ket::test_ket, p...)=op(bra,ket,p...) # no need for symmetrisation

end # module
