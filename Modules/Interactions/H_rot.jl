# rotational interaction. Diagonal in symmetric states as well as asymmetric
function H_rot(bra::Union{scat_αβlml_ket,asym_αβlml_ket},
               ket::Union{scat_αβlml_ket,asym_αβlml_ket},
               R,μ)
    bra == ket || return 0.0u"hartree" # delta function
    l=ket.l
    return uconvert(u"hartree", l*(l+1)*1.0u"ħ^2"/(2*μ*R^2))
end
