using ExoplanetsSysSim
using HDF5, JLD

include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","hsu_etal_2018", "christiansen_func.jl"))

koi_catalog_file_in = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q16_koi_cand.csv")
koi_catalog_file_out = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q16_koi_cand.jld")

koi_catalog, koi_usable = read_koi_catalog(ascii(koi_catalog_file_in))

save(koi_catalog_file_out,"koi_catalog", koi_catalog, "koi_catalog_usable", koi_usable)








 
