#= Γ_GMS ionisation interaction, diagonal in |SmS⟩ basis and therefore nondiagonal
in hyperfine basis.
Description last updated 18/12/20=#

# coupling coefficient. Multiply by 0.3*exp(-R/1.086) Eh, R in a₀, for value
function Γ_GMS_coeffs(bra::asym_αβlml_ket,ket::asym_αβlml_ket)
    # δ fn on |S₁,S₂,i₁,i₂,l,ml>
    (bra.α.S,bra.β.S, bra.α.i,bra.β.i, bra.l,bra.ml)==(ket.α.S,ket.β.S, ket.α.i,ket.β.i, ket.l,ket.ml) || return 0
    # Passed initial δ fn, time to expand into SmS basis
    result=0 # final result accumulates in this variable
    # expand out ket
    for ket_mS_α=-ket.α.S:ket.α.S, ket_mi_α=-ket.α.i:ket.α.i
        # skip if the clebschgordan is zero
        clebschgordan(ket.α.S,ket_mS_α, ket.α.i,ket_mi_α, ket.α.f,ket.α.m)==0 && continue
        # decouple β
        for ket_mS_β=-ket.β.S:ket.β.S, ket_mi_β=-ket.β.i:ket.β.i
            # skip if the clebschgordan is zero
            clebschgordan(ket.β.S,ket_mS_β, ket.β.i,ket_mi_β, ket.β.f,ket.β.m)==0 && continue
            # couple S_α+S_β
            for ket_S=abs(ket.α.S-ket.β.S):(ket.α.S+ket.β.S), ket_mS=-ket_S:ket_S
                # skip if the clebschgordan is zero
                clebschgordan(ket.α.S,ket_mS_α, ket.β.S,ket_mS_β, ket_S,ket_mS)==0 && continue
                ###################### expand out bra
                # decouple α
                for bra_mS_α=-bra.α.S:bra.α.S, bra_mi_α=-bra.α.i:bra.α.i
                    # skip if the clebschgordan is zero
                    clebschgordan(bra.α.S,bra_mS_α, bra.α.i,bra_mi_α, bra.α.f,bra.α.m)==0 && continue
                    # decouple β
                    for bra_mS_β=-bra.β.S:bra.β.S, bra_mi_β=-bra.β.i:bra.β.i
                        # skip if the clebschgordan is zero
                        clebschgordan(bra.β.S,bra_mS_β, bra.β.i,bra_mi_β, bra.β.f,bra.β.m)==0 && continue
                        # couple S_α+S_β
                        for bra_S=abs(bra.α.S-bra.β.S):(bra.α.S+bra.β.S), bra_mS=-bra_S:bra_S
                            # skip if the clebschgordan is zero
                            clebschgordan(bra.α.S,bra_mS_α, bra.β.S,bra_mS_β, bra_S,bra_mS)==0 && continue
                            # have now expanded into a nontrivial |SmS> bra/ket combination
                            # first, δ fn as H_zee diag in this basis (neglecting nuclear spin => diag in mi_α,mi_β)
                            (bra_S,bra_mS,bra_mi_α,bra_mi_β)==(ket_S,ket_mS,ket_mi_α,ket_mi_β) || continue
                            # coefficients of this sum term
                            coeff=clebschgordan(ket.α.S,ket_mS_α, ket.α.i,ket_mi_α, ket.α.f,ket.α.m)*
                                  clebschgordan(ket.β.S,ket_mS_β, ket.β.i,ket_mi_β, ket.β.f,ket.β.m)*
                                  clebschgordan(ket.α.S,ket_mS_α, ket.β.S,ket_mS_β, ket_S,ket_mS)*
                                  clebschgordan(bra.α.S,bra_mS_α, bra.α.i,bra_mi_α, bra.α.f,bra.α.m)*
                                  clebschgordan(bra.β.S,bra_mS_β, bra.β.i,bra_mi_β, bra.β.f,bra.β.m)*
                                  clebschgordan(bra.α.S,bra_mS_α, bra.β.S,bra_mS_β, bra_S,bra_mS)
                            # imaginary potential added for S=0,1
                            Γ = ket_S∈[0,1] ? 1 : 0
                            result += Γ*coeff
                        end # bra S
                    end # bra β
                end # bra α
            end # ket S
        end # ket β
    end # ket α
    result
end

# radial factor (Cocks et al. 2019)
Γ_GMS_radial(R::Unitful.Length) = -im/2*0.3*exp(-R/1.086u"bohr")*1u"hartree"

# testing
Γ_GMS_coeffs(bra::test_ket, ket::test_ket)=0
