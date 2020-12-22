#= F_matrix, for matching AR and BL solution matrices to find a solution matrix
F satisfying inner and outer boundary conditions. Returns F evaluated
at the outer boundary condition.
Note that this code is different to the Honours code in that it does not delete
rows of F corresponding to closed channels (that cannot be done in the computational
basis)
Description last updated 21/12/20=#

module matchF

export F_matrix

using LinearAlgebra, Unitful, UnitfulAtomic

"""Boundary condition matching via QR decomposition
    Inputs: AL, BCs on LHS;
            AR, wavefunction solution to AL at matching location;
            BL, wavefunction solution to BR at matching location;
            BR, BCs on RHS
            ***Assumes same channel ordering in AL,AR,BL,BR***
    Output: F, 2Nₒ×Nₒ matrix of valid (i.e. open channel) wavefunction solutions
    evaluated at the RHS"""
function F_matrix(AL,AR,BL,BR)
    # check on dimensions of inputs
    @assert size(AL)==size(AR) "AL and AR have unlike dimensions"
    @assert size(BL)==size(BR) "BL and BR have unlike dimensions"
    @assert size(AL,1)==size(BL,1) "AL/AR and BL/BR have unlike numbers of rows"
    # numbers of channel, for reference
    N = size(AL,1)÷2 # N channels
    Nₒ = size(BL,2)-N # Nₒ open channels
    # take QR decomposition
    Q = qr(austrip.(permutedims([AR -BL]))).Q
    CD = Q[:,(end-Nₒ+1):end] # 4/09/20 cols of V matching to the zero part of Σ
    # sanity check for linear combinations
    @assert size(CD,1)==2*N+Nₒ "[C; D] doesn't have 2*N+N₀ rows"
    @assert size(CD,2)==Nₒ "[C; D] doesn't have Nₒ columns"
    C = CD[1:N,:] # C is the combination of LHS conditions
    D = CD[(N+1):end,:] # D is the combination of RHS conditions
    # forming F
    F = BR*D # returns satisfying solution at RHS
end

end # module
