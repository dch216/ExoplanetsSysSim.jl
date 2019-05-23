using ExoplanetsSysSim
using StatsFuns
using JLD
using DataFrames
#import ExoplanetsSysSim.StellarTable.df
#import Compat: UTF8String, ASCIIString

## simulation_parameters
function setup_sim_param_model(args::Vector{String} = Array{String}(0) )   # allow this to take a list of parameter (e.g., from command line)
  sim_param = SimParam()
  #add_param_fixed(sim_param,"max_planets_in_sys",20)
  #add_param_fixed(sim_param,"max_tranets_in_sys",8)
  add_param_fixed(sim_param,"max_planets_in_sys",2)
  add_param_fixed(sim_param,"max_tranets_in_sys",2)
  add_param_fixed(sim_param,"num_targets_sim_pass_one",150969)                      # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
  add_param_fixed(sim_param,"num_kepler_targets",150969)                            # Note this is used for the number of Kepler targets for the observational catalog
  #add_param_fixed(sim_param,"generate_star",generate_star_dumb)
  add_param_fixed(sim_param,"generate_planetary_system", ExoplanetsSysSim.generate_planetary_system_simple)
  #add_param_fixed(sim_param,"generate_planetary_system", ExoplanetsSysSim.generate_planetary_system_uncorrelated_incl)

  # add_param_fixed(sim_param,"generate_kepler_target",ExoplanetsSysSim.generate_kepler_target_simple)
  add_param_fixed(sim_param,"generate_kepler_target",generate_kepler_target_from_table)
  #= 
  add_param_fixed(sim_param,"generate_num_planets",ExoplanetsSysSim.generate_num_planets_poisson)
  add_param_active(sim_param,"log_eta_pl",log(2.0))
  =#
  #=
  add_param_fixed(sim_param,"generate_num_planets",generate_num_planets_poisson_mixture)
  add_param_active(sim_param,"frac_zero_planet",0.3)
  add_param_active(sim_param,"frac_one_planet",0.4)
  add_param_active(sim_param,"frac_multi_planet",0.3)
  add_param_active(sim_param,"eta_pl",4.0)
  =#
  add_param_fixed(sim_param,"generate_num_planets",generate_num_planets_categorical)
  add_param_active(sim_param,"fracs_num_planets",ones(get_int(sim_param,"max_planets_in_sys")+1)/(1+get_int(sim_param,"max_planets_in_sys")))
   
  add_param_fixed(sim_param,"generate_planet_mass_from_radius",ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw)
  add_param_fixed(sim_param,"mr_power_index",2.0)
  add_param_fixed(sim_param,"mr_const",1.0)
  #add_param_fixed(sim_param,"generate_period_and_sizes",ExoplanetsSysSim.generate_period_and_sizes_log_normal)
  #add_param_active(sim_param,"mean_log_planet_radius",log(2.0*ExoplanetsSysSim.earth_radius))
  #add_param_active(sim_param,"sigma_log_planet_radius",log(2.0))
  ##add_param_fixed(sim_param,"sigma_log_planet_radius",log(2.0))
  #add_param_active(sim_param,"mean_log_planet_period",log(10.0))
  #add_param_active(sim_param,"sigma_log_planet_period",log(2.0))
  ##add_param_fixed(sim_param,"sigma_log_planet_period",log(2.0))
  #add_param_fixed(sim_param,"generate_period_and_sizes", ExoplanetsSysSim.generate_period_and_sizes_power_law)
  #add_param_active(sim_param,"power_law_P",0.3)
  #add_param_active(sim_param,"power_law_r",-2.44)
  add_param_fixed(sim_param,"min_period",1.0)
  add_param_fixed(sim_param,"max_period",400.0)
  add_param_fixed(sim_param,"min_radius",0.5*ExoplanetsSysSim.earth_radius)
  add_param_fixed(sim_param,"max_radius",10.0*ExoplanetsSysSim.earth_radius)
  add_param_fixed(sim_param,"generate_e_omega",ExoplanetsSysSim.generate_e_omega_rayleigh)
  add_param_fixed(sim_param,"sigma_hk",0.03)
  #add_param_fixed(sim_param,"sigma_hk_one",0.3)
  #add_param_fixed(sim_param,"sigma_hk_multi",0.03)
  add_param_fixed(sim_param,"sigma_incl",0.0)   # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
  add_param_fixed(sim_param,"calc_target_obs_sky_ave",ExoplanetsSysSim.calc_target_obs_sky_ave)
  add_param_fixed(sim_param,"calc_target_obs_single_obs",ExoplanetsSysSim.calc_target_obs_single_obs)
  add_param_fixed(sim_param,"read_target_obs",ExoplanetsSysSim.simulated_read_kepler_observations)
  add_param_fixed(sim_param,"transit_noise_model",ExoplanetsSysSim.transit_noise_model_fixed_noise)
  #add_param_fixed(sim_param,"transit_noise_model",ExoplanetsSysSim.transit_noise_model_diagonal)
  # add_param_fixed(sim_param,"rng_seed",1234)   # If you want to be able to reproduce simulations
  
  #set_inactive(sim_param,["eta_pl","power_law_P","power_law_r"])
  add_param_fixed(sim_param,"star_table_setup",setup_star_table_christiansen)
  add_param_fixed(sim_param,"stellar_catalog","q1_q17_dr25_stellar.jld")
  p_lim_arr_num = [0.5, 1.25, 2.5, 5., 10., 20., 40., 80., 160., 320.]
  r_lim_arr_num = [1., 1.25, 1.5, 1.75, 2.]
  add_param_fixed(sim_param, "p_lim_arr", p_lim_arr_num)
  add_param_fixed(sim_param, "r_lim_arr", r_lim_arr_num*ExoplanetsSysSim.earth_radius)
  #add_param_fixed(sim_param,"generate_num_planets",generate_num_planets_from_rate_table)
  add_param_fixed(sim_param,"generate_period_and_sizes", generate_period_and_sizes_from_rate_table)
  #=
  rate_tab_init = [0.027 0.149 0.354 1.37 2.56 1.96 2.85 4.77 0.0001 ;
                   0.02 0.063 0.432 0.9 1.72 1.46 2.4 5.42 3.63 ;
                   0.0074 0.032 0.389 0.545 1.32 1.65 1.96 1.69 1.17 ;
                   0.0001 0.017 0.156 0.542 0.777 1.35 0.965 1.22 2.49]*0.01

  add_param_active(sim_param,"obs_par", rate_tab_init)
  =#
  return sim_param
end


## planetary_system
using Distributions
function generate_num_planets_categorical(s::Star, sim_param::SimParam)
  const n::Int64 = get_int(sim_param,"max_tranets_in_sys")
  fracs::Array{Float64,1} = get(sim_param,"fracs_num_planets", ones(n+1)/(n+1))
  fracs = abs.(fracs)
  fracs ./= sum(fracs)
  n = rand(Distributions.Categorical(fracs))-1
end

function generate_num_planets_poisson_mixture(s::Star, sim_param::SimParam)
  frac_zero::Float64 = min(max(get_real(sim_param,"frac_zero_planet"), 0.0),1.0)
  frac_one::Float64  = min(max(get_real(sim_param,"frac_one_planet"), 0.0),1.0)
  frac_multi::Float64  = min(max(get_real(sim_param,"frac_multi_planet"),0.0), 1.0)
  frac_sum = frac_zero+frac_one+frac_multi
  frac_zero /= frac_sum
  frac_one /= frac_sum
  frac_multi /= frac_sum
  u = rand()
  if u <= frac_zero 
     return 0
  elseif u<= frac_zero+frac_one
     return 1
  else
     const lambda::Float64 = get_real(sim_param,"eta_pl")
     const max_planets_in_sys::Int64 = get_int(sim_param,"max_planets_in_sys")
     return ExoplanetsSysSim.generate_num_planets_poisson(lambda,max_planets_in_sys,min_planets=2)
  end
end


function generate_num_planets_from_rate_table(s::Star, sim_param::SimParam)
  const max_planets_in_sys::Int64 = get_int(sim_param,"max_planets_in_sys")
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  lambda = sum_kbn(rate_tab)
  #println("# lambda= ", lambda)
  ExoplanetsSysSim.generate_num_planets_poisson(lambda,max_planets_in_sys)
end

function generate_period_and_sizes_from_rate_table(s::Star, sim_param::SimParam; num_pl::Integer = 1)
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  
  limitP::Array{Float64,1} = get_any(sim_param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(sim_param, "r_lim_arr", Array{Float64,1})

  @assert ((length(limitP)-1) == size(rate_tab, 2))
  @assert ((length(limitRp)-1) == size(rate_tab, 1))
  Plist = Array{Float64}(num_pl)
  Rplist = Array{Float64}(num_pl)
  rate_tab_1d = reshape(rate_tab,length(rate_tab))
  #logmaxcuml = logsumexp(rate_tab_1d)
  #cuml = cumsum_kbn(exp(rate_tab_1d-logmaxcuml))
  maxcuml = sum(rate_tab_1d)
  cuml = cumsum_kbn(rate_tab_1d/maxcuml)

  for n in 1:num_pl
    rollp = rand()
    idx = findfirst(x -> x > rollp, cuml)
    i_idx = (idx-1)%size(rate_tab,1)+1
    j_idx = floor(Int64,(idx-1)//size(rate_tab,1))+1
    ### TODO: Keep uniform log sampling here?
    Rp = exp(rand()*(log(limitRp[i_idx+1]/limitRp[i_idx]))+log(limitRp[i_idx]))
    P = exp(rand()*(log(limitP[j_idx+1]/limitP[j_idx]))+log(limitP[j_idx]))
    #push!(Plist, P)
    #push!(Rplist, Rp)
    if !(0<P<Inf)
       println("# BAD P: ",P," should be in ", limitP[j_idx], " - ", limitP[j_idx+1])
    end
    Plist[n] = P
    Rplist[n] = Rp
  end
  return Plist, Rplist
end


## stellar_table
function setup_star_table_christiansen(sim_param::SimParam; force_reread::Bool = false)
  #global df
  df = ExoplanetsSysSim.StellarTable.df
  if haskey(sim_param,"read_stellar_catalog") && !force_reread
     return df
     #return data
  end
  stellar_catalog_filename = convert(String,joinpath(Pkg.dir("ExoplanetsSysSim"), "data", convert(String,get(sim_param,"stellar_catalog","q1_q17_dr25_stellar.csv")) ) )
  df = setup_star_table_christiansen(stellar_catalog_filename)
  add_param_fixed(sim_param,"read_stellar_catalog",true)
  ExoplanetsSysSim.StellarTable.set_star_table(df)
  return df  
end

function setup_star_table_christiansen(filename::String; force_reread::Bool = false)
  df = ExoplanetsSysSim.StellarTable.df
  if ismatch(r".jld$",filename)
  try 
    data = load(filename)
    df::DataFrame = data["stellar_catalog"]
    ExoplanetsSysSim.StellarTable.set_star_table(df)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in jld format."))
  end
  else
  try 
    df = readtable(filename)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in ascii format."))
  end

  has_mass = ! (isna(df[:mass]) | isna(df[:mass_err1]) | isna(df[:mass_err2]))
  has_radius = ! (isna(df[:radius]) | isna(df[:radius_err1]) | isna(df[:radius_err2]))
  has_dens = ! (isna(df[:dens]) | isna(df[:dens_err1]) | isna(df[:dens_err2]))
  has_rest = ! (isna(df[:rrmscdpp04p5]) | isna(df[:dataspan]) | isna(df[:dutycycle]))
  in_Q1Q12 = []
  for x in df[:st_quarters]
    subx = string(x)
    subx = ("0"^(17-length(subx)))*subx
    indQ = search(subx, '1')
    if ((indQ < 1) | (indQ > 12))
      push!(in_Q1Q12, false)
    else
      push!(in_Q1Q12, true)
    end
  end
  is_FGK = []
  for x in 1:length(df[:teff])
    if ((df[x,:teff] > 4000.0) & (df[x,:teff] < 7000.0) & (df[x,:logg] > 4.0))
      push!(is_FGK, true)
    else
      push!(is_FGK, false)
    end
  end
  is_usable = has_radius & is_FGK & has_mass & has_rest #& in_Q1Q12 # & has_dens
  #if contains(filename,"q1_q12_christiansen.jld")
  if contains(filename,"q1_q12_christiansen")   # TODO: Ask Danely what he's trying to do here.
    is_usable = is_usable & in_Q1Q12
  end
  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  # TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp04p5, :dataspan, :dutycycle ]
  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]
  set_star_table(df)
  end
  return df
end

## summary_statistics

# Compile indices for N-tranet systems
function calc_summary_stats_idx_n_tranets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.cache,"idx_n_tranets")
     return css.cache["idx_n_tranets"]
  end
     max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    
     idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
     for n in 1:max_tranets_in_sys-1
       idx_n_tranets[n] = find(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target )
     end
     idx_n_tranets[max_tranets_in_sys] = find(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target )
     css.cache["idx_n_tranets"] = idx_n_tranets
  return idx_n_tranets 
end

# Count total number of tranets using lists of indices for N-tranet systems
function calc_summary_stats_num_tranets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.stat,"num_tranets")
     return css.stat["num_tranets"]
  elseif haskey(css.cache,"num_tranets")
     css.stat["num_tranets"] = css.cache["num_tranets"]
     return css.stat["num_tranets"]
  end
     idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
     num_tranets = 0
     for n in 1:length(idx_n_tranets)
         num_tranets += n*length(idx_n_tranets[n])
     end
     css.stat["num_tranets"] = num_tranets
  return num_tranets 
end

function calc_summary_stats_num_targets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam ; trueobs_cat::Bool = false)
  if !trueobs_cat
    css.stat["num targets"] = get_int(param,"num_targets_sim_pass_one")
  else
    css.stat["num targets"] = get_int(param,"num_kepler_targets")
  end
end

function calc_summary_stats_num_n_tranet_systems!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  #max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    
  num_n_tranet_systems = map(n->length(idx_n_tranets[n]), 1:length(idx_n_tranets) )
  #for n in 1:length(idx_n_tranets)
  #  num_n_tranet_systems[n] = length(idx_n_tranets[n])
  #end
  css.stat["num n-tranet systems"] = num_n_tranet_systems
end

function calc_summary_stats_duration_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.stat,"duration_ratio_list")
     return css.stat["duration_ratio_list"]
  elseif haskey(css.cache,"duration_ratio_list")
     return css.cache["duration_ratio_list"]
  end
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  @assert length(idx_n_tranets) >= 1 

  # Calculate how many duration ratios there will be & allocate storage
  num_ratios = 0
  for i in 2:length(idx_n_tranets)
     num_ratios += length(idx_n_tranets[i])*(i-1)
  end
  duration_ratio_list = Array{Float64}(num_ratios)
 
  k = 0
  for n in 2:length(idx_n_tranets)         # Loop over number of tranets in system
    for i in idx_n_tranets[n]              # Loop over systems with n tranets
       period_in_sys = Array{Float64}(n)
       duration_in_sys = Array{Float64}(n)
       for j in 1:n                        # Loop over periods within a system
         period_in_sys[j] = cat_obs.target[i].obs[j].period
         duration_in_sys[j] = cat_obs.target[i].obs[j].duration
       end
       perm = sortperm(period_in_sys)
       for j in 1:(n-1)                       # Loop over period ratios within a system
          inv_period_ratio = period_in_sys[perm[j+1]]/period_in_sys[perm[j]]
          if 1<inv_period_ratio<Inf
             k += 1
             duration_ratio_list[k] = duration_in_sys[perm[j]]/duration_in_sys[perm[j+1]] * inv_period_ratio^(1//3)
          end
       end
    end
  end
  resize!(duration_ratio_list,k)
  css.stat["duration_ratio_list"] = duration_ratio_list

  return duration_ratio_list
end

function calc_summary_stats_period_radius_ratios_neighbors_internal!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.stat,"period_ratio_list") && haskey(css.stat,"radius_ratio_list")
     return (css.stat["period_ratio_list"], css.stat["radius_ratio_list"])
  elseif haskey(css.cache,"period_ratio_list") && haskey(css.cache,"radius_ratio_list")
     return (css.cache["period_ratio_list"], css.cache["radius_ratio_list"])
  end
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  @assert length(idx_n_tranets) >= 1 

  # Calculate how many period ratios there will be & allocate storage
  num_ratios = 0
  for i in 2:length(idx_n_tranets)
     num_ratios += length(idx_n_tranets[i])*(i-1)
  end
  period_ratio_list = Array{Float64}(num_ratios)
  radius_ratio_list = Array{Float64}(num_ratios)
 
  k = 0
  for n in 2:length(idx_n_tranets)         # Loop over number of tranets in system
    period_in_sys = Array{Float64}(n)
    #radius_in_sys = Array{Float64}(n)
    depth_in_sys = Array{Float64}(n)
    for i in idx_n_tranets[n]              # Loop over systems with n tranets
       for j in 1:n                        # Loop over periods within a system
         period_in_sys[j] = cat_obs.target[i].obs[j].period
         #radius_in_sys[j] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
         depth_in_sys[j] = cat_obs.target[i].obs[j].depth
       end
       perm = sortperm(period_in_sys)
       for j in 1:(n-1)                       # Loop over period ratios within a system
          k = k+1
          period_ratio_list[k] = period_in_sys[perm[j]]/period_in_sys[perm[j+1]]
          #radius_ratio_list[k] = radius_in_sys[perm[j]]/radius_in_sys[perm[j+1]]
          radius_ratio_list[k] = sqrt(depth_in_sys[perm[j]]/depth_in_sys[perm[j+1]]) 
       end
    end
  end
  css.cache["period_ratio_list"] = period_ratio_list
  css.cache["radius_ratio_list"] = radius_ratio_list

  return (period_ratio_list,radius_ratio_list)
end


function calc_summary_stats_period_radius_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  (period_ratio_list,radius_ratio_list) = calc_summary_stats_period_radius_ratios_neighbors!(css,cat_obs,param)
  css.stat["period_ratio_list"] = period_ratio_list
  css.stat["radius_ratio_list"] = radius_ratio_list
  return (period_ratio_list, radius_ratio_list)
end
function calc_summary_stats_period_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  period_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css,cat_obs,param)[1]
  css.stat["period_ratio_list"] = period_ratio_list
  return period_ratio_list
end
function calc_summary_stats_radius_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  radius_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css,cat_obs,param)[2]
  css.stat["radius_ratio_list"] = radius_ratio_list
  return radius_ratio_list
end

function calc_summary_stats_mean_std_log_period_depth!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  # Allocate arrays to store values for each tranet
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  period_list = zeros(num_tranets)
  depth_list = zeros(num_tranets)
  #weight_list = ones(num_tranets)

  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys") 
  @assert max_tranets_in_sys >= 1
   i = 1   # tranet id
   for targ in cat_obs.target                        # For each target 
     for j in 1:min(length(targ.obs),max_tranets_in_sys)          # For each tranet around that target (but truncated if too many tranets in one system)
         #println("# i= ",i," j= ",j)
         period_list[i] = targ.obs[j].period
         depth_list[i] = targ.obs[j].depth
         #weight_list[i] = 1.0
         i = i+1
      end
   end

  css.cache["P list"] = period_list                                     # We can store whole lists, e.g., if we want to compute K-S distances
  css.cache["depth list"] = depth_list
  #css.cache["weight list"] = weight_list

  idx_good = Bool[ period_list[i]>0.0 && depth_list[i]>0.0 for i in 1:length(period_list) ]
  log_period_list = log10(period_list[idx_good])
  log_depth_list = log10(depth_list[idx_good])
  css.stat["mean log10 P"]  =  mean_log_P = mean(log_period_list)
  css.stat["mean log10 depth"]  =  mean_log_depth = mean(log_depth_list)
  css.stat["std log10 P"]  = std_log_P = stdm(log_period_list,mean_log_P)
  css.stat["std log10 depth"]  = std_log_depth = stdm(log_depth_list,mean_log_depth)

  return (mean_log_P, std_log_P, mean_log_depth, std_log_depth) 
end

function calc_summary_stats_cuml_period_depth!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  # Allocate arrays to store values for each tranet
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  period_list = zeros(num_tranets)
  depth_list = zeros(num_tranets)
  #weight_list = ones(num_tranets)

  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys") 
  @assert max_tranets_in_sys >= 1
  i = 0   # tranet id
   for targ in cat_obs.target                        # For each target 
     for j in 1:min(length(targ.obs),max_tranets_in_sys)          # For each tranet around that target (but truncated if too many tranets in one system)
         i = i+1
         #println("# i= ",i," j= ",j)
         period_list[i] = targ.obs[j].period
         depth_list[i] = targ.obs[j].depth
         #weight_list[i] = 1.0
      end
   end
  resize!(period_list,i)
  resize!(depth_list,i)
  css.stat["P list"] = period_list                                     # We can store whole lists, e.g., if we want to compute K-S distances
  css.stat["depth list"] = depth_list
  #css.cache["weight list"] = weight_list
  #println("# P list = ",period_list)
  return (period_list,depth_list)
end


function calc_summary_stats_obs_binned_rates!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam; trueobs_cat::Bool = false)
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)

  idx_tranets = find(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1}             # Find indices of systems with at least 1 tranet = potentially detectable transiting planet
  css.cache["idx_tranets"] = idx_tranets                                   # We can save lists of indices to summary stats for pass 2, even though we won't use these for computing a distance or probability

  #=
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    # Demo that simulation parameters can specify how to evalute models, too
  @assert max_tranets_in_sys >= 1

  # Count total number of tranets and compile indices for N-tranet systems
  num_tranets = 0
  idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
  for n in 1:max_tranets_in_sys-1
    idx_n_tranets[n] = find(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target )
    num_tranets += n*length(idx_n_tranets[n])
  end
  idx_n_tranets[max_tranets_in_sys] = find(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target )
  css.cache["idx_n_tranets"] = idx_n_tranets

  num_tranets += max_tranets_in_sys*length(idx_n_tranets[max_tranets_in_sys])  # WARNING: this means we need to ignore planets w/ indices > max_tranets_in_sys
  if ( length( find(x::KeplerTargetObs-> length(x.obs) > max_tranets_in_sys, cat_obs.target ) ) > 0)   # Make sure max_tranets_in_sys is at least big enough for observed systems
    warn("Observational data has more transiting planets in one systems than max_tranets_in_sys allows.")
  end
  num_tranets  = convert(Int64,num_tranets)            # TODO OPT: Figure out why isn't this already an Int.  I may be doing something that prevents some optimizations
  css.cache["num_tranets"] = num_tranets                                   

  num_sys_tranets = zeros(max_tranets_in_sys)                           # Since observed data, don't need to calculate probabilities.
  for n in 1:max_tranets_in_sys                                         # Make histogram of N-tranet systems
    num_sys_tranets[n] = length(idx_n_tranets[n])
  end
  css.stat["num_sys_tranets"] = num_sys_tranets
  css.stat["planets detected"] = num_tranets 
  =#

  period_list = zeros(num_tranets)
  radius_list = zeros(num_tranets)
  weight_list = ones(num_tranets)

  n = 1    # tranet id
  for i in 1:length(cat_obs.target)
    for j in 1:num_planets(cat_obs.target[i])
      period_list[n] = cat_obs.target[i].obs[j].period
      radius_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
      n = n+1
    end
  end

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})
  @assert length(limitP)>=2 && length(limitRp)>=2

  np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  np_bin_idx = 1
  for i in 1:(length(limitP)-1)
    P_match = find(x -> ((x > limitP[i]) && (x < limitP[i+1])), period_list)
    for j in 1:(length(limitRp)-1)
      R_match = find(x -> ((x > limitRp[j]) && (x < limitRp[j+1])), radius_list)
      
      bin_match = intersect(P_match, R_match)

      np_bin[np_bin_idx] = sum(weight_list[bin_match])
      np_bin_idx += 1
    end
  end

  #css.stat["planets detected"] = sum(np_bin)
  css.stat["planets table"] = np_bin

  return css
end


function calc_summary_stats_model(cat_obs::KeplerObsCatalog, param::SimParam; trueobs_cat::Bool = false)
  css = CatalogSummaryStatistics()
  calc_summary_stats_num_targets!(css,cat_obs,param,trueobs_cat=trueobs_cat)
  calc_summary_stats_num_tranets!(css,cat_obs,param)
  calc_summary_stats_num_n_tranet_systems!(css,cat_obs,param)
  calc_summary_stats_cuml_period_depth!(css,cat_obs,param)
  #calc_summary_stats_obs_binned_rates!(css,cat_obs,param)
  #calc_summary_stats_mean_std_log_period_depth!(css,cat_obs,param)
  calc_summary_stats_period_ratios_neighbors!(css,cat_obs,param)
  calc_summary_stats_duration_ratios_neighbors!(css,cat_obs,param)
  return css
end

## abc_distance
function calc_distance_model_ks(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  d1 = calc_distance_num_planets(summary1,summary2,sim_param,verbose=verbose)
  d2 = calc_distance_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d3 = calc_distance_num_planets_binned(summary1,summary2,sim_param,verbose=verbose)
  #d4 = calc_distance_mean_std_log_period_depth(summary1,summary2,sim_param,verbose=verbose)
  d5 = calc_distance_ks_period(summary1,summary2,sim_param,verbose=verbose)
  d6 = calc_distance_ks_depth(summary1,summary2,sim_param,verbose=verbose)
  d7 = calc_distance_ks_period_ratios(summary1,summary2,sim_param,verbose=verbose)
  d8 = calc_distance_ks_duration_ratios(summary1,summary2,sim_param,verbose=verbose)
  return vcat(d1, d2, d5, d6, d7, d8)
end

function calc_distance_model_kl(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  d1 = calc_distance_kl_num_planets(summary1,summary2,sim_param,verbose=verbose)
  d2 = calc_distance_kl_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d2 = calc_distance_hellinger_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d3 = calc_distance_num_planets_binned(summary1,summary2,sim_param,verbose=verbose)
  #d4 = calc_distance_mean_std_log_period_depth(summary1,summary2,sim_param,verbose=verbose)
  d5 = calc_distance_kl_period(summary1,summary2,sim_param,verbose=verbose)
  d6 = calc_distance_kl_depth(summary1,summary2,sim_param,verbose=verbose)
  d7 = calc_distance_kl_period_ratios(summary1,summary2,sim_param,verbose=verbose)
  d8 = calc_distance_kl_duration_ratios(summary1,summary2,sim_param,verbose=verbose)
  return vcat(d1, d2, d5, d6, d7, d8)
end
calc_distance_model(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false) = calc_distance_model_kl(summary1,summary2,sim_param,verbose=verbose)

function calc_distance_num_planets(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = haskey(summary1.stat,"num_tranets") ? summary1.stat["num_tranets"] : summary1.stat["expected planets detected"]
    np2 = haskey(summary2.stat,"num_tranets") ? summary2.stat["num_tranets"] : summary2.stat["expected planets detected"]
    #println("np1 = ",np1,", np2 = ",np2)

    dist_np = dist_L1_abs(np1/summary1.stat["num targets"], np2/summary2.stat["num targets"])
    #println("np1 (normalized) = ",np1/summary1.stat["num targets"],", np2 (normalized) = ",np2/summary2.stat["num targets"],", d = ",dist_np)
  return dist_np
end

function calc_distance_num_planets_binned(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = haskey(summary1.stat,"planets table") ? summary1.stat["planets table"] : summary1.stat["expected planets table"]
    np2 = haskey(summary2.stat,"planets table") ? summary2.stat["planets table"] : summary2.stat["expected planets table"]
    #println("np1 = ",np1,", np2 = ",np2)

    dist_np_bin = zeros(length(np1))
    for n in 1:length(np1)
      dist_np_bin[n] = dist_L1_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
      #println("True # [Bin ", n,"] = ",np1[n],", Expected # [Bin ", n,"] = ",np2[n])
    end
    #println("np1 (normalized) = ",np1/summary1.stat["num targets"],", np2 (normalized) = ",np2/summary2.stat["num targets"],", dist = ",dist_np_bin)
  return dist_np_bin
end

function calc_distance_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")    
  d = zeros(max_tranets_in_sys)    
  for n in 1:max_tranets_in_sys   
    d[n] = dist_L1_abs(summary1.stat["num n-tranet systems"][n]/summary1.stat["num targets"], summary2.stat["num n-tranet systems"][n]/summary2.stat["num targets"])
  end
  return d
end

function calc_distance_mean_std_log_period_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    d1 = dist_L1_abs(summary1.stat["mean log10 P"], summary2.stat["mean log10 P"])
    d2 = dist_L1_abs(summary1.stat["mean log10 depth"], summary2.stat["mean log10 depth"])
    d3 = dist_L1_abs(summary1.stat["std log10 P"], summary2.stat["std log10 P"])
    d4 = dist_L1_abs(summary1.stat["std log10 depth"], summary2.stat["std log10 depth"])
    return [d1, d2, d3, d4]
end

# compute supremum of differences between empirical cdfs.
# Borrowed from JuliaStats/HypothesisTests.jl
function ksstats{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})
    n_x, n_y = length(x), length(y)
    sort_idx = sortperm([x; y])
    pdf_diffs = [ones(n_x)/n_x; -ones(n_y)/n_y][sort_idx]
    cdf_diffs = cumsum(pdf_diffs)
    deltap = maximum(cdf_diffs)
    deltan = -minimum(cdf_diffs)
    delta = max(deltap,deltan)
    (n_x, n_y, deltap, deltan, delta)
end

function calc_distance_ks_period(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["P list"] 
  samp2 = summary2.stat["P list"] 
  return ksstats(samp1,samp2)[5]
end

function calc_distance_ks_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["depth list"] 
  samp2 = summary2.stat["depth list"] 
  return ksstats(samp1,samp2)[5]
end

function calc_distance_ks_period_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["period_ratio_list"] 
  samp2 = summary2.stat["period_ratio_list"] 
  return ksstats(samp1,samp2)[5]
end

function calc_distance_ks_duration_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["duration_ratio_list"] 
  samp2 = summary2.stat["duration_ratio_list"] 
  return ksstats(samp1,samp2)[5]
end

# Function for Relative Entropy / K-L divergence
include("kde.jl")

function calc_distance_kl_num_planets(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = (haskey(summary1.stat,"num_tranets") ? summary1.stat["num_tranets"] : summary1.stat["expected planets detected"]) 
    np2 = (haskey(summary2.stat,"num_tranets") ? summary2.stat["num_tranets"] : summary2.stat["expected planets detected"]) 
    ntarg1 = summary1.stat["num targets"]
    ntarg2 = summary2.stat["num targets"]
    mintarg = min(ntarg1,ntarg2)
    alpha1 = 1+np1
    beta1 = 1+ntarg1
    #mu1 = alpha1/beta1
    #sigma1 = mu1 /sqrt(alpha1)
    #dist = alpha+lgamma(alpha)+(1-alpha)*digamma(alpha) -log(beta)
    alpha2 = 1+np2
    beta2 = 1+ntarg2
    #dist -= alpha+lgamma(alpha)+(1-alpha)*digamma(alpha) -log(beta)
    #mu2 = alpha2/beta2
    #sigma2 = mu2 /sqrt(alpha1)
    #dist = 0.5*( (sigma2/sigma1)^2 - 1 + 2*log(sigma1/sigma2) + (mu2-mu1)^2/sigma1^2 )  # Normal approximation, substitute denom for difference in rates
    #= Gamma PDF
    dist = abs( (alpha2-alpha1)*digamma(alpha1) - lgamma(alpha1) + lgamma(alpha2)  ) # Gamma approximation, same beta, why need abs? 
    dist += abs( (alpha1-alpha2)*digamma(alpha2) - lgamma(alpha2) + lgamma(alpha1)  )  # Gamma approximation, same beta, why need abs? 
    print("# a1=",alpha1, " a2=",alpha2, " digamma=",digamma(alpha1), " ", digamma(alpha2)," gamma term=",lgamma(alpha1)-lgamma(alpha2))
    if ntarg1!=ntarg2
       print("# ntarg=",ntarg1," ",ntarg2," dist(before betas)=", dist) 
       dist += alpha2*(log(beta1)-log(beta2)) + alpha1*(beta2-beta1)/beta1  # Gammma approximation, beta terms
       dist += alpha1*(log(beta2)-log(beta1)) + alpha2*(beta1-beta2)/beta2  # Gammma approximation, beta terms
       println(" dist (after betas)=",dist) 
    end
    =#
    ## #= Poisson 
    #dist = abs( alpha2-alpha1+alpha1*log(alpha1/alpha2) )
    #dist += abs( alpha2-alpha1+alpha1*log(alpha1/alpha2) )
    rate1 = alpha1/beta1*mintarg
    rate2 = alpha2/beta2*mintarg
    dist = (rate1-rate2)*log(rate1/rate2) 
    ## =#
    dist /= 2
    return dist
end

function calc_distance_kl_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")   
  #= categorical distribution
  f1sum = sum(summary1.stat["num n-tranet systems"]) # summary1.stat["num targets"]
  f2sum = sum(summary2.stat["num n-tranet systems"]) # summary2.stat["num targets"]
  if !(f1sum>0 && f2sum>0) return 0.0  end 
  d = zeros(max_tranets_in_sys)
  for n in 1:max_tranets_in_sys   
    f1 = summary1.stat["num n-tranet systems"][n]/f1sum
    f2 = summary2.stat["num n-tranet systems"][n]/f2sum
    m = (f1+f2)/2
    if m>zero(m)
    if f1>zero(f1) 
       d[n] += 0.5*f1*log(f1/m)
    end
    if f2>zero(f2) 
       d[n] += 0.5*f2*log(f2/m)
    end
#   else
#       d += Inf
    end
  end
  =#
  # Poisson distributions for each
    ntarg1 = summary1.stat["num targets"]
    ntarg2 = summary2.stat["num targets"]
    #mintarg = min(ntarg1,ntarg2)
  d = zeros(max_tranets_in_sys)
  for n in 1:max_tranets_in_sys   
    f1 = (1+summary1.stat["num n-tranet systems"][n])/(1+ntarg1) # *(1+mintarg)
    f2 = (1+summary2.stat["num n-tranet systems"][n])/(1+ntarg2) # *(1+mintarg)
    d[n] = (f1-f2)*log(f1/f2)
  end
  return d
end
function calc_distance_hellinger_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")    
  f1sum = sum(summary1.stat["num n-tranet systems"]) # summary1.stat["num targets"]
  f2sum = sum(summary2.stat["num n-tranet systems"]) # summary2.stat["num targets"]
  d = 1
  for n in 1:max_tranets_in_sys   
    f1 = summary1.stat["num n-tranet systems"][n]/f1sum
    f2 = summary2.stat["num n-tranet systems"][n]/f2sum
    d -= sqrt(f1*f2)
  end
  #d = sqrt(d)
  return d
end

function calc_distance_kl_period(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  #println("# P list 1: n=",length(summary1.stat["P list"])," min=",minimum(summary1.stat["P list"]), " max=",maximum(summary1.stat["P list"]))
  #println("# P list 2: n=",length(summary2.stat["P list"])," min=",minimum(summary2.stat["P list"]), " max=",maximum(summary2.stat["P list"]))
  samp1 = log(summary1.stat["P list"] )
  samp2 = log(summary2.stat["P list"] )
  calc_kl_distance_ab(samp1,samp2,log(0.5),log(320.) )
end

function calc_distance_kl_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = log(summary1.stat["depth list"] )
  samp2 = log(summary2.stat["depth list"] )
  calc_kl_distance_ab(samp1,samp2,log(0.000025),log(0.025) )
end

function calc_distance_kl_period_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  min_ratios_to_compute_distances = 3
  distance_when_not_enough_ratios = 10.0
  samp1 = summary1.stat["period_ratio_list"] 
  samp2 = summary2.stat["period_ratio_list"] 
  if length(samp1)<min_ratios_to_compute_distances || length(samp2)<min_ratios_to_compute_distances
     return distance_when_not_enough_ratios
  end
  calc_kl_distance_ab(samp1,samp2,0.0,1.0)
end

function calc_distance_kl_duration_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  min_ratios_to_compute_distances = 3
  distance_when_not_enough_ratios = 10.0
  samp1 = log( summary1.stat["duration_ratio_list"] )
  samp2 = log( summary2.stat["duration_ratio_list"] )
  if length(samp1)<min_ratios_to_compute_distances || length(samp2)<min_ratios_to_compute_distances
     return distance_when_not_enough_ratios
  end
  calc_kl_distance_ab(samp1,samp2,-log(10.0),log(10.0) )
end



## eval_model
function test_model()
  global sim_param_closure = setup_sim_param_model()
  cat_phys = generate_kepler_physical_catalog(sim_param_closure)
  cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_closure)
  global summary_stat_ref_closure = calc_summary_stats_obs_demo(cat_obs,sim_param_closure)
  global cat_phys_try_closure = generate_kepler_physical_catalog(sim_param_closure)
  global cat_obs_try_closure  = observe_kepler_targets_sky_avg(cat_phys_try_closure,sim_param_closure)
  global summary_stat_try_closure  = calc_summary_stats_sim_pass_one_demo(cat_obs_try_closure,cat_phys_try_closure,sim_param_closure)
  summary_stat_try_closure   = calc_summary_stats_sim_pass_two_demo(cat_obs_try_closure,cat_phys_try_closure,summary_stat_try_closure,sim_param_closure)
  param_guess = make_vector_of_sim_param(sim_xparam_closure)
  evaluate_model_scalar_ret( param_guess)
end


