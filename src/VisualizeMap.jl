using GeoDatasets
using PyCall
using PyPlot

const 🌎_low = GeoDatasets.landseamask(; resolution=('c'), grid=5)

function plot_gses(gses::Vector)
    lats = [gs[1] / π * 180 for gs ∈ gses]
    lons = [gs[2] / π * 180 for gs ∈ gses]
    pcolormesh(🌎_low[1], 🌎_low[2], 🌎_low[3]')
    scatter(lons, lats, marker="v", color="pink")
    gca().set_aspect(1)
end


function plot_gs_selections(
    spacing::Int,
    pop_ctr_gses::Vector,
    selected_gs_importances::Dict
)
    # Get CRS and features
    ccrs = pyimport("cartopy.crs")
    cfeat = pyimport("cartopy.feature")

    # Create figure and axis with PlateCarree projection
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # Plot coatlines, borders, and background image
    ax.coastlines(linewidth=0.5)
    ax.add_feature(cfeat.BORDERS, linewidth=0.5)
    ax.stock_img()

    # Plot the population centers
    lats = [gs[1][1] / π * 180 for gs ∈ pop_ctr_gses]
    lons = [gs[1][2] / π * 180 for gs ∈ pop_ctr_gses]
    weights = [gs[2] for gs ∈ pop_ctr_gses]

    max_weight = maximum(weights)
    min_weight = minimum(weights)
    sizes = [(weight - min_weight) / (max_weight - min_weight) * 100 for weight ∈ weights]

    ax.scatter(
        lons,
        lats,
        marker="v",
        color="blue",
        alpha=0.5,
        s=sizes,
        transform=ccrs.PlateCarree(),
        label="Population Centers",
        zorder=5
    )

    # Plot the selected ground stations with importance values
    selected_lats = [gs[1] / π * 180 for gs ∈ keys(selected_gs_importances)]
    selected_lons = [gs[2] / π * 180 for gs ∈ keys(selected_gs_importances)]
    selected_importances = collect(values(selected_gs_importances))

    max_importance = maximum(selected_importances)
    min_importance = minimum(selected_importances)
    selected_sizes = [
        (importance - min_importance) / (max_importance - min_importance) * 100
        for importance ∈ selected_importances
    ]

    ax.scatter(
        selected_lons,
        selected_lats,
        marker="o",
        color="red",
        alpha=0.5,
        s=selected_sizes,
        transform=ccrs.PlateCarree(),
        label="Selected Ground Stations",
        zorder=5
    )

    # Display the map
    plt.title("Ground Stations Selected from $(spacing) km Grid")
    plt.show()
end