
using IntervalSets
include("SemiRingOfLifetimes.jl")

function build_zero_matrix(N::Int = 20)
    cells = fill(SemiRingOfLifetimes.THE_EMPTY_SET, N, N)
    return PwrSetRealsAdjMatrix(cells)
end

function build_const_matrix(N::Int = 20)
    cells = fill(SemiRingOfLifetimes.THE_REALS, N, N)
    return PwrSetRealsAdjMatrix(cells)
end

function build_id_matrix(N::Int = 20)
    cells = [i == j ? SemiRingOfLifetimes.THE_REALS : SemiRingOfLifetimes.THE_EMPTY_SET for i in 1:N, j in 1:N]
    return PwrSetRealsAdjMatrix(cells)
end


const ZERO_MATRIX = build_zero_matrix(20)
const CONSTANT_MATRIX = build_const_matrix(20)
const IDENTITY_MATRIX = build_id_matrix(20)