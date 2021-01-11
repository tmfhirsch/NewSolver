# rotational interaction. Diagonal in symmetric states as well as asymmetric
function H_rot(bra::Union{scat_αβlml_ket,asym_αβlml_ket},
               ket::Union{scat_αβlml_ket,asym_αβlml_ket},
               R::Unitful.Length,μ::Unitful.Mass)
    bra == ket || return 0.0u"hartree" # delta function
    l=ket.l
    return uconvert(u"hartree", l*(l+1)*1.0u"ħ^2"/(2*μ*R^2))
end

# testing
function H_rot(bra::test_ket, ket::test_ket, R, μ)
    3u"bohr" < R < 50u"bohr" || return 0u"hartree" # well in inner region
    bra != ket && return 2*1e-8u"hartree" # off-diagonals; C = 1e-8u"hartree"
    bra.x == 1 && return -1e-8u"hartree" # negative for x=1
    return 0e0u"hartree" # zero for x=2,3 diagonals
end
