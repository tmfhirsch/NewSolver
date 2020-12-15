const ϵ_hfs=1.519830e-7u"hartree" # from Cocks et al. 2019

# hyperfine energy from Cocks et al. (2019)
function E_hfs(α::atom_nos)
    i,f = α.i,α.f # unpack i, f for quick reference
    @assert i==0 || i==half(1) "i ≠ 0 && i ≠ 1/2" # sanity check on nuclear spin
    i==0 && return 0u"hartree" # ⁴He* has no hyperfine structure
    # now in i==0.5 case
    @assert f==half(1) || f==half(3) "i==1/2 but (f ≠ 1/2 && f ≠ 3/2)" # sanity check
    f == half(1) && return ϵ_hfs # f==0.5, i.e. singlet spin case
    # now in f==1.5 case
    return 0u"hartree" # triplet spin defined as zero energy
end


# hyperfine interaction
function H_hfs(bra::asym_αβlml_ket,ket::asym_αβlml_ket)
    bra == ket || return 0u"hartree"
    return E_hfs(ket.α) + E_hfs(ket.β) # assumed independent of R
end
