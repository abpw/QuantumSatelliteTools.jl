using IntervalSets
include("SemiRingOfLifetimes.jl")

# --- helpers: measure/coverage of a (finite union of) intervals ---
# length of a single interval (endpoints' openness doesn't change measure)
_measure(I::AbstractInterval) = rightendpoint(I) - leftendpoint(I)
_measure(::EmptyInterval) = 0.0
# union-of-intervals (IntervalSet stores disjoint pieces)
_measure(S::IntervalSet) = sum(_measure, S.intervals)
# fall back: if anything turns infinite
coverage(s::SetOfReals) = try _measure(s.interval) catch; Inf end

const CAP = 86400.0         # total horizon (e.g., 1 day in seconds)
const THRESH = 0.98            # saturation threshold fraction

struct PwrSetRealsAdjMatrix
    cells::Matrix{SetOfReals}
end

function PwrSetRealsAdjMatrix(N::Integer=20)
    cells = [i==j ? THE_REALS : THE_EMPTY_SET for i in 1:N, j in 1:N]
    PwrSetRealsAdjMatrix(cells)
end

Base.size(M::PwrSetRealsAdjMatrix) = size(M.cells)

function Base.:+(A::PwrSetRealsAdjMatrix, B::PwrSetRealsAdjMatrix)
    @assert size(A) == size(B)
    N, M = size(A)
    C = Matrix{SetOfReals}(undef, N, M)
    @inbounds for i in 1:N, j in 1:M
        s = A.cells[i,j] + B.cells[i,j]
        C[i,j] = coverage(s) >= THRESH*CAP ? THE_REALS : s
    end
    PwrSetRealsAdjMatrix(C)
end


function Base.:*(A::PwrSetRealsAdjMatrix, B::PwrSetRealsAdjMatrix)
    @assert size(A,2) == size(B,1)
    N, K = size(A); _, M = size(B)
    C = Matrix{SetOfReals}(undef, N, M)
    @inbounds for i in 1:N, j in 1:M
        acc = THE_EMPTY_SET
        for k in 1:K
            acc = acc + (A.cells[i,k] * B.cells[k,j])
            if coverage(acc) >= THRESH*CAP
                acc = THE_REALS
                break
            end
        end
        C[i,j] = acc
    end
    PwrSetRealsAdjMatrix(C)
end

# add a closed interval [t1,t2] (or set directly via alt)
function add_interval!(M::PwrSetRealsAdjMatrix, i::Int, j::Int;
                       t1::Real=0, t2::Real=0, alt::Union{Nothing,SetOfReals}=nothing)
    _bounds_check(M, i, j)
    if alt !== nothing
        M.cells[i,j] = alt
    else
        s = M.cells[i,j] + SetOfReals(ClosedInterval(float(t1), float(t2)))
        M.cells[i,j] = coverage(s) >= THRESH*CAP ? THE_REALS : s
    end
    return M
end

# bulk add
function add_intervals!(M::PwrSetRealsAdjMatrix, items)
    for (i, j, t1, t2) in items
        add_interval!(M, i, j; t1=t1, t2=t2)
    end
    return M
end

# pretty print non-empty entries
function Base.show(io::IO, M::PwrSetRealsAdjMatrix)
    N, _ = size(M)
    printed = false
    for i in 1:N, j in 1:N
        if M.cells[i,j] != THE_EMPTY_SET
            printed = true
            println(io, lpad(i,2), " -> ", lpad(j,2), " : ", M.cells[i,j])
        end
    end
    printed || print(io, "<all zero>")
end

# matrix equality (structural)
Base.:(==)(A::PwrSetRealsAdjMatrix, B::PwrSetRealsAdjMatrix) = A.cells == B.cells

# bounds check like your staticmethod
_bounds_check(M::PwrSetRealsAdjMatrix, i::Int, j::Int) =
    (1 <= i <= size(M,1) && 1 <= j <= size(M,2)) ||
    throw(BoundsError(M.cells, (i, j)))