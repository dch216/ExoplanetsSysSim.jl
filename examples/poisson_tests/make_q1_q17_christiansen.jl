## ExoplanetsSysSim/examples/hsu_etal_2018/make_q1_q16_christiansen.jl
## (c) 2018 Danley C. Hsu & Eric B. Ford
# Script for creating the Q1-Q16 FGK subset stellar catalog in jld format
#    used in the Hsu et al. 2018 paper

using ExoplanetsSysSim
using JLD

include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","poisson_tests", "christiansen_func.jl"))

stellar_catalog_file_in = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q17_dr25_stellar.csv")
stellar_catalog_file_out = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1q17_dr25_kepler_fgk.jld")

stellar_catalog = setup_star_table_christiansen(ascii(stellar_catalog_file_in))

save(stellar_catalog_file_out,"stellar_catalog", stellar_catalog, "stellar_catalog_usable", ExoplanetsSysSim.StellarTable.usable)
 
