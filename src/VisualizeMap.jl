using PyPlot
using GeoDatasets

const 🌎_low = GeoDatasets.landseamask(;resolution='c',grid=5)

function plot_gses(gses::Vector)
    lats = [gs[1]/π*180 for gs in gses]
    lons = [gs[2]/π*180 for gs in gses]
    pcolormesh(🌎_low[1],🌎_low[2],🌎_low[3]')
    scatter(lons, lats, marker="v", color="pink")
    gca().set_aspect(1)
end
