# rotational interaction. Diagonal in symmetric states as well as asymmetric
function H_rot(bra::Union{scat_αβlml_ket,asym_αβlml_ket},
               ket::Union{scat_αβlml_ket,asym_αβlml_ket},
               R::Unitful.Length,μ::Unitful.Mass,
               C) # C is just for testing and can be ignored
    bra == ket || return 0.0u"hartree" # delta function
    l=ket.l
    return uconvert(u"hartree", l*(l+1)*1.0u"ħ^2"/(2*μ*R^2))
end

# testing
function H_rot(bra::test_ket, ket::test_ket, R::Unitful.Length, μ::Unitful.Mass, C::Unitful.Energy)
    bra == ket || return 0u"hartree" # diagonal
    ket.x==false || return 0u"hartree" # open channel, zero potential
    3u"bohr" <= R <= 13u"bohr" || return 0u"hartree" # zero outside well region
    return -1u"hartree" + C # closed channel inside well
end
