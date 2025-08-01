module QuantumSatelliteTools

    using PyPlot, GeoDatasets
    using SatelliteToolboxTle, SatelliteToolboxPropagators, SatelliteToolboxTransformations, SatelliteToolboxBase
    using SatelliteAnalysis
    using Downloads, CSV, Random, Dates, StaticArrays

    module AstronomyGeometry
        include("AstronomyGeometry.jl")
    end

    module GenerateGroundStations
        include("GenerateGroundStations.jl")
    end
    
    module GenerateTLEs
        include("GenerateTLEs.jl")
    end
    
    module VisualizeMap
        include("VisualizeMap.jl")
    end

end
