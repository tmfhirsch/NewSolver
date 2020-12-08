""" Coupling function: (bra,ket) -> (c₁,c₂,c₃) coefficients,
    where <bra|H_el|ket> = c₁*¹V(R)+c₂*²V(R)+c₃*³V(R)"""
function H_el_coupling(bra::asym_αβlml_ket,ket::asym_αβlml_ket)
end

""" Radial function: (c₁,c₂,c₃) -> c₁*¹V(R)+c₂*²V(R)+c₃*³V(R)"""
function H_el_radial(cs::Tuple{Float64,Float64,Float64}, R::Unitful.Length)
    return cs[1]*Singlet(R) + cs[2]*Triplet(R) + cs[3]*Quintet(R)
end
