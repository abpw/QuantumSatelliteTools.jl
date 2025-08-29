using IntervalSets

REALS = IntervalSets.OpenInterval(-Inf, Inf)

struct SetOfReals{I}
    interval::I            
end


Base.:+(a::SetOfReals, b::SetOfReals) = SetOfReals(a.interval + b.interval)

Base.:-(a::SetOfReals, b::SetOfReals) = SetOfReals(a.interval - b.interval)

Base.:*(a::SetOfReals, b::SetOfReals) = SetOfReals(a.interval & b.interval)


Base.:<=(a::SetOfReals, b::SetOfReals) = issubset(a,b)


Base.show(io::IO, s::SetOfReals) = print(io, s.interval)


THE_REALS = SetOfReals(REALS)
THE_EMPTY_SET = SetOfReals(EmptyInterval(Int))