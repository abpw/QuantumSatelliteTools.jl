using IntervalSets

const REALS = IntervalSets.OpenInterval(-Inf, Inf)

struct SetOfReals{I}
    interval::I            # can be an AbstractInterval, IntervalSet, or EmptyInterval
end


Base.:+(a::SetOfReals, b::SetOfReals) = SetOfReals(a.interval + b.interval)

Base.:-(a::SetOfReals, b::SetOfReals) = SetOfReals(a.interval - b.interval)

Base.:*(a::SetOfReals, b::SetOfReals) = SetOfReals(a.interval & b.interval)


Base.:<=(a::SetOfReals, b::SetOfReals) = issubset(a,b)


Base.show(io::IO, s::SetOfReals) = print(io, s.interval)


const THE_REALS = SetOfReals(REALS)
const THE_EMPTY_SET = SetOfReals(EmptyInterval(Int))