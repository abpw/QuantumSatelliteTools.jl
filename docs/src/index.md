```@meta
CurrentModule = QuantumSatelliteTools
```

# QuantumSatelliteTools

**QuantumSatelliteTools.jl** is a Julia package for modeling and simulating quantum communication scenarios involving satellites. It provides geometric, orbital, and (soon) atmospheric tools for analyzing quantum link performance, visibility, and channel quality between satellite and ground nodes.

It makes heavy usage of [SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl), and familiarity with that package is necessary for most usage (i.e. for creating satellite propagator objects).

## Features

- Ground station and satellite geometry utilities (ECEF/ECI conversions, ellipsoid intersections)
- Optical line-of-sight calculations including atmospheric loss modeling
- Satellite visibility, and sun interference tracking

## Getting Started

Add the package using Julia's package manager:
```julia
] add QuantumSatelliteTools
```

**QuantumSatelliteTools.jl** contains several namespaces that each provide their own functionality. See their documentation linked below on use.

## Index

```@index
Pages   = ["AstronomyGeometry.md", "GenerateGroundStations.md", "GenerateSatellites.md", "FreespaceChannels.md", "LossCalculation.md"]
Modules = []
```
