#= Matches solution matrix to bessel functions, producing a reactance matrix 𝐊
following Mies (1980).
Description last updated 21/12/20=#

module matchK
export K_matrix

using SpecialFunctions, Unitful, UnitfulAtomic
""" Solver for 𝐊 matrix (following Mies eqn (3.8)), where 𝐅=𝐉-𝐍𝐊.
    Inputs:
        R~[L] the asymptotic radial distance to match K at;
        𝐅 = [𝐆; 𝐆'] where 𝐆 is matrix of wavefunctions for
    different initial conditions evaluated at R (only open channels);
        𝐤 ~ [L]⁻¹ vector of wavenumber for each channel at rhs;
    **Note: this function assumes it is only being given open channels**
        𝐥 vector of l quantum numbers for each channel.
    **Note: eval, 𝐤, 𝐥 must share the same ordering/number of channels**
    Output:
        𝐊 ~ n×n matrix (n=number of channels considered)"""
function K_matrix(R, 𝐅, 𝐤, 𝐥)
    # match for A, B where G=J.A-N.B
    #construct G, G' matrices to match for A, B with
    n=Int(size(𝐅,1)/2) # n = number of channels. Assumes sol in above form.
    @assert size(𝐅,2) == n "solution matrix not of shape 2n × n"
    G, G⁻ = austrip.(𝐅[1:n,1:n]), 𝐅[n+1:2*n,1:n]
    # solve for A,B
    A, B = zeros(ComplexF64,n,n), zeros(ComplexF64,n,n) # initialise
    for i in 1:n, j in 1:n
        # construct [jᵢ nᵢ; jᵢ' nᵢ'] matrix, here called [bj -bn; bj⁻ -bn⁻]
        # expressions for derivatives (⁻) calculated using Mathematica
        # function form from Mies (A2) which ≡ Cocks et al (59)
        k=𝐤[i]
        l=𝐥[i]
        bj=austrip(sqrt(k)*R)*sphericalbesselj(l,k*R)
        bj⁻=austrip(sqrt(k))*((l+1)*sphericalbesselj(l,k*R)
            -k*R*sphericalbesselj(l+1,k*R))
        bn=austrip(sqrt(k)*R)*sphericalbessely(l,k*R)
        bn⁻=austrip(sqrt(k))*((l+1)*sphericalbessely(l,k*R)
            -k*R*sphericalbessely(l+1,k*R))
        Gᵢⱼ, G⁻ᵢⱼ = G[i,j], G⁻[i,j]
        AB = [bj bn; bj⁻ bn⁻]\[Gᵢⱼ; G⁻ᵢⱼ] #12/10 removed minus sign on 𝐍 # AB≡[Aᵢⱼ; Bᵢⱼ], solve J;-N*AB=G;G⁻
        A[i,j], B[i,j] = AB
    end
    𝐊 = B*inv(A)
    return 𝐊
end

end # module
