# PASSED TESTS
using LinearAlgebra
include(raw"C:\\Users\\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterCreating\\DistanceCalc.jl")
include(raw"C:\\Users\\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterCreating\\ObjectiveFunction.jl")

dist_km(city1, city2) = dist(city1, city2) / 1000.0

function shapelyAlg(Cities::Vector{NTuple{2,Float64}}, alpha::Float64 = 0.98)
    n = length(Cities)
    n == 0 && return NTuple{2,Float64}[]
    Phi = Vector{Float64}(undef, n)

    @inbounds for i ∈ 1:n
        ci = Cities[i]
        s = 0.0
        for j in 1:n
            j == i && continue
            s += objective_function(dist_km(ci, Cities[j]), Cities)
        end
        Phi[i] = 0.5 * s
    end
    
    Q = collect(1:n)
    heads = NTuple{2,Float64}[]
    
    while !isempty(Q)
        m = Q[argmax(@view Phi[Q])]
        push!(heads, Cities[m])
        Tm = [i for i in Q if objective_function(dist_km(Cities[i], Cities[m]), Cities) >= alpha]
        Q = setdiff(Q, Tm)
    end
    
    return heads
end