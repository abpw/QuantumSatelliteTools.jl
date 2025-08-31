include(raw"C:\Users\alexb\ProgrammingProjects\URV2025\QSATJulia\QuantumSatelliteTools.jl\GlobalNetworkExperiment\ClusterCreating\KMeansAlg.jl")
include(raw"C:\Users\alexb\ProgrammingProjects\URV2025\QSATJulia\QuantumSatelliteTools.jl\GlobalNetworkExperiment\ClusterCreating\ShapelyAlg.jl")
include(raw"C:\Users\alexb\ProgrammingProjects\URV2025\QSATJulia\QuantumSatelliteTools.jl\GlobalNetworkExperiment\ClusterMerging\ClusterBridging.jl")

function new_GSes(Cities::Vector{NTuple{2,Float64}}, GS_dist:: Float64 = 450.0)
    # Create clusters using Shapely algorithm
    cluster_heads = shapelyAlg(Cities)
    print(length(cluster_heads))
    # Cluster cities around the heads
    clustered_cities = cluster_cities_fixed_heads(Cities, cluster_heads)
    # Create bridge chains between clusters
    new_nodes = create_bridge_chains(clustered_cities, distance_threshold=GS_dist)
   
    return new_nodes
end