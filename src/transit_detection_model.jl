## ExoplanetsSysSim/src/transit_detection_model.jl
## (c) 2015 Eric B. Ford

# Several functions below based on https://github.com/christopherburke/KeplerPORTs/blob/master/KeplerPORTs_utils.py
# That follows the procedure outlined in Burke et al.(2015).
# However we don't currently interpolate the CDDP or mesthreshold to the relevant duration

function real_log_choose(m::Real, n::Real)
  lgamma(m+1)-lgamma(n+1)-lgamma(m-n+1.0)
end
  
function real_binom(k::Real, BigM::Real, f::Real)
    F1 = real_log_choose(BigM,k)
    F2 = k*log(f)
    F3 = (BigM-k)*log(1.0-f)
    x = exp(F1+F2+F3)
    return x
end
  
function kepler_window_function(num_transits_no_gaps::Real, duty_cycle::Real; min_transits::Real = 3.0)  # TODO: SCI DARIN:  This is where we'd use your window function data
  if num_transits_no_gaps < min_transits
     return 0.0
  else
     return max(1.0 - real_binom(min_transits,num_transits_no_gaps,duty_cycle), 0.0)
  end
end
    
# See Christiansen et al. (2015)  This assumes a linear limbdarkening coefficient of 0.6
function frac_depth_to_tps_depth(frac_depth::Real)
    const alp = 1.0874
    const bet = 1.0187
    const REALDEPTH2TPSSQUARE = 1.0  # WARNING: Waiting for this to be confirmed
    k = sqrt(frac_depth)
    tps_depth = min( (alp-bet*k) * frac_depth* REALDEPTH2TPSSQUARE, 1.0)   # NOTE: I added the max based on common sense
    return tps_depth
end

function detection_efficiency_theory(mes::Real; min_pdet_nonzero::Float64 = 0.0)
   const muoffset =  0.0
   const sig =  1.0
   const mesthresh = 7.1
   if mes > (9.0 - mesthresh - muoffset)
      return 1.0
   else
      pdet = 0.5 + 0.5*erf((mes - mesthresh - muoffset) / sqrt(2.0*sig*sig))
      pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
      return pdet
   end
end

function detection_efficiency_fressin2013(mes::Real)
   const mesmin =  6.0
   const mesmax =  16.0
   if mes <= mesmin
      return 0.0
   elseif mes >= mesmax
      return 1.0
   else
     return (mes - mesmin) / (mesmax - mesmin)
  end
end


function detection_efficiency_christiansen2015(mes::Real; mes_threshold::Real = 7.1, min_pdet_nonzero::Float64 = 0.0)
   const a =  4.65  # from code for detection_efficiency(...) at https://github.com/christopherburke/KeplerPORTs/blob/master/KeplerPORTs_utils.py
   const b =  0.98
   #const a =  4.35 # from arxiv abstract. Informal testing showed it didn't matter
   #const b =  1.05
   usemes = max(0.0,mes - 7.1 - (mes_threshold - 7.1))
   pdet = cdf(Gamma(a,b), usemes)
   pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
   return pdet
end

function detection_efficiency_dr25_simple(mes::Real; min_pdet_nonzero::Float64 = 0.0)
   const a =  30.87  # from pg 16 of https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19110-001.pdf
   const b =  0.271
   const c = 0.940
   pdet = c*cdf(Gamma(a,b), mes)
   pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
   return pdet
end

# WARNING: Hardcoded choice of transit detection efficiency here for speed and so as to not have it hardcoded in multiple places
#detection_efficiency_model = detection_efficiency_christiansen2015  
detection_efficiency_model = detection_efficiency_dr25_simple 

# Resume code original to SysSim


function interpolate_cdpp_to_duration(t::KeplerTarget, duration::Real)
   duration_in_hours = duration *24.0
   dur_idx = searchsortedlast(cdpp_durations,duration_in_hours)   # cdpp_durations is defined in constants.jl
   if dur_idx <= 0 
      cdpp = t.cdpp[1,1]
   elseif dur_idx==length(cdpp_durations) && (duration_in_hours > cdpp_durations[end]) # Should be 15 cdpp_durations.
      cdpp = t.cdpp[length(cdpp_durations),1]
   else
      w = ((duration_in_hours)-cdpp_durations[dur_idx]) / (cdpp_durations[dur_idx+1]-cdpp_durations[dur_idx])
      cdpp = w*t.cdpp[dur_idx+1,1] + (1-w)*t.cdpp[dur_idx,1]
   end
   return cdpp
end


function calc_snr_if_transit(t::KeplerTarget, depth::Real, duration::Real, cdpp::Real, sim_param::SimParam; num_transit::Real = 1)
   depth_tps = frac_depth_to_tps_depth(depth)                  # WARNING: Hardcoded this conversion 
   snr = depth_tps*sqrt(num_transit*duration*LC_rate)/cdpp     # WARNING: Assumes measurement uncertainties are uncorrelated & CDPP based on LC
end

function calc_snr_if_transit_central(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  depth = calc_transit_depth(t,s,p)
  duration_central = calc_transit_duration_central(t,s,p)
  cdpp = interpolate_cdpp_to_duration(t, duration_central)
  num_transit = calc_expected_num_transits(t,s,p,sim_param)
  calc_snr_if_transit(t,depth,duration_central,cdpp, sim_param,num_transit=num_transit)
end

function calc_prob_detect_if_transit(t::KeplerTarget, snr::Real, sim_param::SimParam; num_transit::Real = 1)
  const min_transits = 3.0                                                    # WARNING: Hard coded 3 transit minimum
  const min_pdet_nonzero = 1.0e-4                                                # TODO OPT: Consider raising threshold to prevent a plethora of planets that are very unlikely to be detected due to using 0.0 or other small value here
  wf = kepler_window_function(num_transit, t.duty_cycle, min_transits=min_transits)   # TODO SCI DETAIL: Replace statistical model with checking actual transit times for long period planets
  return wf*detection_efficiency_model(snr, min_pdet_nonzero=min_pdet_nonzero)	
end

function calc_prob_detect_if_transit(t::KeplerTarget, depth::Real, duration::Real, cdpp::Real, sim_param::SimParam; num_transit::Real = 1)
  snr = calc_snr_if_transit(t,depth,duration,cdpp, sim_param, num_transit=num_transit)
  return calc_prob_detect_if_transit(t, snr, sim_param, num_transit=num_transit)
end

function calc_prob_detect_if_transit_central(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  depth = calc_transit_depth(t,s,p)
  duration_central = calc_transit_duration_central(t,s,p)
  cdpp = interpolate_cdpp_to_duration(t, duration_central)
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_prob_detect_if_transit(t,depth,duration_central,cdpp, sim_param, num_transit=ntr)
end

function calc_prob_detect_if_transit_with_actual_b(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
  depth = calc_transit_depth(t,s,p)
  duration = calc_transit_duration(t,s,p)
  cdpp = interpolate_cdpp_to_duration(t, duration)
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_prob_detect_if_transit(t,depth,duration,cdpp, sim_param, num_transit=ntr)
end

# Compute probability of detection if we average over impact parameters b~U[0,1)
function calc_ave_prob_detect_if_transit_from_snr(t::KeplerTarget, snr_central::Real, duration_central::Real, size_ratio::Real, cdpp_central::Real, sim_param::SimParam; num_transit::Real = 1)
  const min_transits = 3.0                                                    # WARNING: Hard coded 3 transit minimum
  const min_pdet_nonzero = 1.0e-4
  wf = kepler_window_function(num_transit, t.duty_cycle, min_transits=min_transits)   
  
  detection_efficiency_central = detection_efficiency_model(snr_central, min_pdet_nonzero=min_pdet_nonzero) 
  if wf*detection_efficiency_central <= min_pdet_nonzero
     return 0.
  end 

  # Breaking integral into two sections [0,1-b_boundary) and [1-b_boundary,1], so need at least 5 points to evaluate integral via trapezoid rule
  const num_impact_param_low_b =  7                            # Number of points to evaluate integral over [0,1-b_boundary) via trapezoid rule
  const num_impact_param_high_b = 5 # (size_ratio<=0.05) ? 5 : 11  # Number of points to evaluate integral over [1-b_boudnary,1) via trapezoid rule.  If using 2*size_ratio for bondary for small planets, then keep this odd, so one point lands on 1-size_ratio.
  @assert(num_impact_param_low_b >= 5)
  @assert(num_impact_param_high_b >= 3)
  const num_impact_param = num_impact_param_low_b+num_impact_param_high_b-1 # One point is shared
  const b_boundary = (size_ratio <= 0.15) ? 2*size_ratio : min(max(0.3,size_ratio),0.5)
  b = Array{Float64}(num_impact_param)
  weight = Array{Float64}(num_impact_param)
  b[1:num_impact_param_low_b] = linspace(0.0,1-b_boundary,num_impact_param_low_b)
  b[num_impact_param_low_b:num_impact_param] = linspace(1-b_boundary,1.0,num_impact_param_high_b)
  weight[1:num_impact_param_low_b] = (1-b_boundary)/(num_impact_param_low_b-1)  # Points for first integral
  weight[1] *= 0.5                        # Lower endpoint of first integral
  weight[num_impact_param_low_b] *= 0.5   # Upper endpoint of first integral
  weight[num_impact_param_low_b] += 0.5*(b_boundary)/(num_impact_param_high_b-1) # Also lower endpoint of second integral
  weight[(num_impact_param_low_b+1):num_impact_param] = b_boundary/(num_impact_param_high_b-1)
  weight[num_impact_param] *= 0.5         # Upper endpoint of second integral
  @assert isapprox(sum(weight),1.0)

  function integrand(b::Real) 
     depth_factor = calc_depth_correction_for_grazing_transit(b,size_ratio)
     duration_factor = calc_transit_duration_factor_for_impact_parameter_b(b,size_ratio) 
     cdpp = interpolate_cdpp_to_duration(t,duration_central*duration_factor)
     snr_factor = depth_factor*sqrt(duration_factor)*(cdpp_central/cdpp)
     detection_efficiency_model(snr_central*snr_factor, min_pdet_nonzero=min_pdet_nonzero)
  end 

  ave_detection_efficiency = sum(weight .* map(integrand,b)::Vector{Float64} )    # WARNING:  Doesn't account for cdpp chaning for shorter duration transits 
  #ave_detection_efficiency = sum(weight .* map(b->detection_efficiency_model(snr_central*sqrt(calc_effective_transit_duration_factor_for_impact_parameter_b(b,size_ratio)), min_pdet_nonzero=min_pdet_nonzero),b)::Vector{Float64} )    # WARNING:  Doesn't account for cdpp chaning for shorter duration transits 

  return wf*ave_detection_efficiency
end


function calc_ave_prob_detect_if_transit(t::KeplerTarget, depth::Real, duration_central::Real, size_ratio::Real, sim_param::SimParam; num_transit::Real = 1)
  cdpp_central = interpolate_cdpp_to_duration(t, duration_central)
  snr_central = calc_snr_if_transit(t,depth,duration_central,cdpp_central, sim_param, num_transit=num_transit)
  return calc_ave_prob_detect_if_transit_from_snr(t, snr_central, duration_central, size_ratio, cdpp_central, sim_param, num_transit=num_transit)
end

function calc_ave_prob_detect_if_transit(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
  depth = calc_transit_depth(t,s,p)
  duration_central = calc_transit_duration_central(t,s,p)
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_ave_prob_detect_if_transit(t,depth,duration_central, size_ratio, sim_param, num_transit=ntr)
end

