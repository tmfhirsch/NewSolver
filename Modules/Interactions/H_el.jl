""" Coupling function: (bra,ket) -> (c₁,c₂,c₃) coefficients,
    where <bra|H_el|ket> = c₁*¹V(R)+c₃*³V(R)+c₅*⁵V(R)"""
function H_el_coeffs(bra::asym_αβlml_ket,ket::asym_αβlml_ket)
    C¹,C³,C⁵ = 0,0,0 # 'tallying' the different couplings => integers
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
                ######################
                # expand out bra
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
                            # the coefficients are nontrivial. Now we apply the H_el delta function
                            SmS_ket = (ket.α.S,ket.β.S,ket_S,ket_mS, ket.α.i,ket_mi_α,ket.β.i,ket_mi_β, ket.l,ket.ml)
                            SmS_bra = (bra.α.S,bra.β.S,bra_S,bra_mS, bra.α.i,bra_mi_α,bra.β.i,bra_mi_β, bra.l,bra.ml)
                            SmS_ket==SmS_bra || continue # H_el diagonal in SmS basis
                            # nontrivial coefficients, and delta function passed, so we do the calculation
                            coeff=clebschgordan(ket.α.S,ket_mS_α, ket.α.i,ket_mi_α, ket.α.f,ket.α.m)*
                                  clebschgordan(ket.β.S,ket_mS_β, ket.β.i,ket_mi_β, ket.β.f,ket.β.m)*
                                  clebschgordan(ket.α.S,ket_mS_α, ket.β.S,ket_mS_β, ket_S,ket_mS)*
                                  clebschgordan(bra.α.S,bra_mS_α, bra.α.i,bra_mi_α, bra.α.f,bra.α.m)*
                                  clebschgordan(bra.β.S,bra_mS_β, bra.β.i,bra_mi_β, bra.β.f,bra.β.m)*
                                  clebschgordan(bra.α.S,bra_mS_α, bra.β.S,bra_mS_β, bra_S,bra_mS)
                            @assert ket_S ∈ [0,1,2] "ket_S ∉ [0,1,2]" # sanity check
                            if ket_S == 0
                                C¹ += coeff # singlet spin case
                            elseif ket_S == 1
                                C³ += coeff # triplet
                            else
                                C⁵ += coeff # quintetS
                            end # if
                        end # bra S,mS
                    end # bra mS_β, mi_β
                end # bra mS_α, mi_α
            end # ket S,mS
        end # ket mS_β, mi_β
    end # ket mS_α, mi_α
    return (C¹,C³,C⁵)
end

""" Radial function: (c₁,c₃,c₅) -> c₁*¹V(R)+c₃*⁵V(R)+c₃*⁵V(R)"""
function H_el_radial(cs::Tuple{Float64,Float64,Float64}, R::Unitful.Length)
    return cs[1]*Singlet(R) + cs[2]*Triplet(R) + cs[3]*Quintet(R)
end

# test
H_el_coeffs(bra::test_ket, ket::test_ket) = (0,0,0)
