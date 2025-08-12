include("AdjMatrix.jl")
include("SpecialMatrices.jl")

function KleeneStar(AdjMat::AdjMatrix.PwrSetRealsAdjMatrix)
    if AdjMat == ZERO_MATRIX
        return IDENTITY_MATRIX
    elseif AdjMat == CONSTANT_MATRIX || AdjMat == IDENTITY_MATRIX
        return AdjMat
    else
        C = AdjMat * AdjMat
        return AdjMat + KleeneStar(C)
    end
end
