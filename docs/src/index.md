```@meta
CurrentModule = QuantumSatelliteTools
```

# QuantumSatelliteTools

**QuantumSatelliteTools.jl** is a Julia package for modeling and simulating quantum communication scenarios involving satellites. It provides geometric, orbital, and (soon) atmospheric tools for analyzing quantum link performance, visibility, and channel quality between satellite and ground nodes.

Makes heavy usage of [SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl), and familiarity with that package is necessary for most usage (i.e. for creating satellite propagator objects).

## Features

- Ground station and satellite geometry utilities (ECEF/ECI conversions, ellipsoid intersections)
- Optical line-of-sight calculations including (soon) atmospheric loss modeling
- Satellite visibility, and sun interference tracking

## Getting Started

Add the package using Julia's package manager:
```julia
] add QuantumSatelliteTools
```

Import it in your Julia environment:
```julia
using QuantumSatelliteTools
```

## Documentation Structure

- [`QuantumSatelliteTools`](@ref): Main module documentation
- Tutorials and usage examples (coming soon)
- API reference (generated automatically below)

## Index

```@index
```

## API Documentation

```@autodocs
Modules = [QuantumSatelliteTools.AstronomyGeometry, QuantumSatelliteTools.GenerateTLEs, QuantumSatelliteTools.GenerateGroundStations]
Order   = [:type, :function]
```
