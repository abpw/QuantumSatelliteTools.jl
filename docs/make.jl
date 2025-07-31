using QuantumSatelliteTools
using Documenter

DocMeta.setdocmeta!(QuantumSatelliteTools, :DocTestSetup, :(using QuantumSatelliteTools); recursive=true)

makedocs(;
    modules=[QuantumSatelliteTools],
    authors="Albert Williams <albertbpwilliams@gmail.com> and contributors",
    sitename="QuantumSatelliteTools.jl",
    format=Documenter.HTML(;
        canonical="https://abpw.github.io/QuantumSatelliteTools.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/abpw/QuantumSatelliteTools.jl",
    devbranch="main",
)
