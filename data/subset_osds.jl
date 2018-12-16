using JLD
using DataFrames

OSD_file = load(joinpath(Pkg.dir(), "ExoplanetsSysSim", "data", "allosds.jld"))

allosds = OSD_file["allosds"]
periods = OSD_file["periods"]
kepids = OSD_file["kepids"]

OSD_file = 0

cat_ids = load(joinpath(Pkg.dir(), "ExoplanetsSysSim", "data", "q1q17_dr25_gaia_fgk.jld"), "stellar_catalog")[:kepid]

tmp_ind = [findfirst(x -> x == y, kepids) for y in cat_ids]

kepids = kepids[tmp_ind]
allosds = allosds[tmp_ind,:,:]

# allosds = convert(Array{Float32,3}, allosds)

save(joinpath(Pkg.dir(), "ExoplanetsSysSim", "data", "dr25fgk_osds.jld"), "allosds", allosds, "periods", periods, "kepids", kepids)
