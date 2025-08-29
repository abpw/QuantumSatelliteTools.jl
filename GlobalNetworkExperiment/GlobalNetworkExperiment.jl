module GlobalNetworkExperiment
    include(raw"C:\\Users\\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterCreating\\DistanceCalc.jl")
    include(raw"C:\\Users\\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterCreating\\ObjectiveFunction.jl")
    include(raw"C:\\Users\\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterCreating\\ShapelyAlg.jl")
    include(raw"C:\\Users\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterCreating\\KMeansAlg.jl")
    include(raw"C:\\Users\\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterMerging\\ClusterBridging.jl")

    using NearestNeighbors
    using LinearAlgebra
    using Distances
    using StaticArrays 
end