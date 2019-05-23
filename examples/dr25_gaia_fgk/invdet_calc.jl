## ExoplanetsSysSim/examples/hsu_etal_2018/invdet_calc.jl
## (c) 2018 Danley C. Hsu
# Script for producing Q1-Q16 FGK planet candidate occurrence rate estimates
#    using both the inverse detection efficiency and the simple
#    Bayesian methods

using ExoplanetsSysSim
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","poisson_tests", "christiansen_func.jl"))

global sim_param_closure = setup_sim_param_christiansen()
sim_param_closure = set_test_param(sim_param_closure)
add_param_fixed(sim_param_closure,"koi_catalog","q1_q16_koi_cand.csv")

df_star = setup_star_table_christiansen(sim_param_closure)
println("# Finished reading in stellar data")
df_koi,usable_koi = read_koi_catalog(sim_param_closure)
println("# Finished reading in KOI data") 
cat_obs = setup_actual_planet_candidate_catalog(df_star, df_koi, usable_koi, sim_param_closure)

@time inv_det_simp_bayes(cat_obs, sim_param_closure)
