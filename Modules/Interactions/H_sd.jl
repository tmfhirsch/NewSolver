
# ξ factor (Cocks et al (2019) equation (7))
const μ_ratio = 1.00115965 # ratio of e⁻ mag mom to Bohr magneton (Cocks 2019)
const α_fs = 0.0072973525693 # fine structure constant (Wikipedia)
const ξ = α_fs^2 * μ_ratio^2 * 1.0u"hartree*bohr^3"

""" Radial factor in H_sd, Vₚ(R)
    Input: R~[L]
    Output: Vₚ(R)~[E]
    Based off Cocks et al (2019) eqn (39)=Beams et al (2006) below (9)"""
function H_sd_radial(R)
    # radial factor
    b=-sqrt(6)*ξ # ħ⁻² not written since it cancels with the TTensor factor
    Vₚ=b/R^3
    return Vₚ
end


# returns coefficient, independent of radial term
function H_sd_coeffs(bra::asym_αβlml_ket,ket::asym_αβlml_ket)
    # apply Delta function (nuclear spins not involved)
    (ket.α.i, ket.β.i)==(bra.α.i, bra.β.i) || return 0
    result=0
    # decouple α in ket
    for ket_mS_α=-ket.α.S:ket.α.S, ket_mi_α=-ket.α.i:ket.α.i
        # skip if the clebschgordan is zero
        clebschgordan(ket.α.S,ket_mS_α, ket.α.i,ket_mi_α, ket.α.f,ket.α.m)==0 && continue
        # decouple β in ket
        for ket_mS_β=-ket.β.S:ket.β.S, ket_mi_β=-ket.β.i:ket.β.i
            # skip if the clebschgordan is zero
            clebschgordan(ket.β.S,ket_mS_β, ket.β.i,ket_mi_β, ket.β.f,ket.β.m)==0 && continue
            # decouple α in bra
            for bra_mS_α=-bra.α.S:bra.α.S, bra_mi_α=-bra.α.i:bra.α.i
                # skip if the clebschgordan is zero
                clebschgordan(bra.α.S,bra_mS_α, bra.α.i,bra_mi_α, bra.α.f,bra.α.m)==0 && continue
                # decouple β in bra
                for bra_mS_β=-bra.β.S:bra.β.S, bra_mi_β=-bra.β.i:bra.β.i
                    # skip if the clebschgordan is zero
                    clebschgordan(bra.β.S,bra_mS_β, bra.β.i,bra_mi_β, bra.β.f,bra.β.m)==0 && continue
                    # apply Delta function (nuclear spins not involved)
                    (ket_mi_α, ket_mi_β)==(bra_mi_α, bra_mi_β) || continue
                    # couple ket S
                    for ket_S=abs(ket.α.S-ket.β.S):(ket.α.S+ket.β.S), ket_mS=-ket_S:ket_S
                        # skip if the clebschgordan is zero
                        clebschgordan(ket.α.S,ket_mS_α, ket.β.S,ket_mS_β, ket_S,ket_mS)==0 && continue
                        # couple bra S
                        for bra_S=abs(bra.α.S-bra.β.S):(bra.α.S+bra.β.S), bra_mS=-bra_S:bra_S
                            # skip if the clebschgordan is zero
                            clebschgordan(bra.α.S,bra_mS_α, bra.β.S,bra_mS_β, bra_S,bra_mS)==0 && continue
                            # have now expanded into the |SmS> basis, where my Honours code can work
                            # coeffs of this sum term
                            coeff=clebschgordan(ket.α.S,ket_mS_α, ket.α.i,ket_mi_α, ket.α.f,ket.α.m)*
                                  clebschgordan(ket.β.S,ket_mS_β, ket.β.i,ket_mi_β, ket.β.f,ket.β.m)*
                                  clebschgordan(ket.α.S,ket_mS_α, ket.β.S,ket_mS_β, ket_S,ket_mS)*
                                  clebschgordan(bra.α.S,bra_mS_α, bra.α.i,bra_mi_α, bra.α.f,bra.α.m)*
                                  clebschgordan(bra.β.S,bra_mS_β, bra.β.i,bra_mi_β, bra.β.f,bra.β.m)*
                                  clebschgordan(bra.α.S,bra_mS_α, bra.β.S,bra_mS_β, bra_S,bra_mS)
                            result += coeff*H_sd_SmS_coeffs(bra.α.S,bra.β.S,bra_S,bra_mS,bra.l,bra.ml,
                                                     ket.α.S,ket.β.S,ket_S,ket_mS,ket.l,ket.ml)
                        end # bra S, mS
                    end # ket S, mS
                end # bra β mS, mi
            end # bra α mS, mi
        end # ket β mS, mi
    end # ket α mS, mi
    result
end


"""Redefined CG, returns 0 for unphysical (jᵢ,mᵢ) combination (instead of error)"""
function clebschgordan_lax(j₁,m₁,j₂,m₂,j₃,m₃=m₁+m₂)
    for (jᵢ,mᵢ) in ((j₁, m₁), (j₂, m₂), (j₃, m₃))
        if !WignerSymbols.ϵ(jᵢ,mᵢ) # unphysical jᵢ,mᵢ entered
            return 0::Int
        end
    end
    return clebschgordan(j₁,m₁,j₂,m₂,j₃,m₃) # if all good, send on to CG calc
end

"""Coefficients for spin-dipole interaction for |Φₐ⟩≡|SmS⟩|lml⟩ kets.
    Input: ⟨Φₐ'|, |Φₐ⟩ numbers
    Output: ⟨Φₐ'|̂H|Φₐ⟩/Vₚ(R) ~ 1"""
function H_sd_SmS_coeffs(Sα_,Sβ_,S_,mS_,l_,ml_,Sα,Sβ,S,mS,l,ml)
    # δ function
    (mS+ml)==(mS_+ml_) || return 0
    term=(-1)^(mS_-mS)
    term*= clebschgordan_lax(S,mS,2,mS_-mS,S_,mS_)
    term*=clebschgordan_lax(l,ml,2,ml_-ml,l_,ml_)
    term*=TTensor(Sα,Sβ,S,Sα_,Sβ_,S_)
    term*=CTensor(l,l_)
    term
end

"""⟨Γ'S'||𝐓²||ΓS⟩ from Beams et al. (2006) eqn (19)"""
function TTensor(S₁,S₂,S,S₁_,S₂_,S_)
    (S₁==S₁_ && S₂==S₂_) || return 0
    𝐓² = sqrt(S₁*(S₁+1)*S₂*(S₂+1))
    𝐓²*= sqrt(5*(2*S₁+1)*(2*S₂+1)*(2*S+1))
    𝐓²*= wigner9j(S₁,S₂,S,1,1,2,S₁,S₂,S_)
    𝐓²
end

"""⟨Γ'S'||𝐂²||ΓS⟩ from Beams et al. (2006) eqn (20)"""
CTensor(l,l_)=sqrt((2*l+1)/(2*l_+1))*clebschgordan_lax(l,0,2,0,l_,0)
