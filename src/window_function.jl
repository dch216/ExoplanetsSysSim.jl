## ExoplanetsSysSim/src/window_function.jl
## (c) 2018 Darin Ragozzine

# Gather and prepare the window function data
module WindowFunction

export setup_window_function, get_window_function_data, get_window_function_id, eval_window_function

#using DataArrays
using DataFrames
#using CSV
using JLD
# using PyPlot	#for testing OSD interpolator. Remove when finished testing
using ApproXD	#for interpolating OSD table
using ExoplanetsSysSim.SimulationParameters


# Object to hold window function data
immutable window_function_data
  window_func_array::Array{Float64,3}      # Value of window function (window_function_id, duration_id, period_id).  Maybe rename to wf_value or data?
  wf_durations_in_hrs::Array{Float64,1}    # Boundaries for duration bins in window_func_array.  Maybe rename to durations?
  wf_periods_in_days::Array{Float64,1}     # Boundaries for periods bins in window_func_array.  Maybe rename to periods?
  sorted_quarter_strings::Array{Int64,1}   # TODO OPT: Is there a reason to keep this?  Maybe rename to quarter_strings?
  allsortedkepids::Array{Int64,1}          # value is Kepler ID.  Index is same as index to window_function_id_arr
  window_function_id_arr::Array{Int64,1}   # value is index for window_func_array.  Index is same as index to allsortedkepids
  default_wf_id::Int64                     # Index corresponding to the default window function
end


function window_function_data()
  window_function_data( Array{Float64,3}(0,0,0), Array{Float64,1}(0),Array{Float64,1}(0), Array{Int64,1}(0),Array{Int64,1}(0),Array{Int64,1}(0), 0 )
end


win_func_data = window_function_data()

function setup(sim_param::SimParam; force_reread::Bool = false)
  global win_func_data
  if haskey(sim_param,"read_window_function") && !force_reread
     return win_func_data
  end
  window_function_filename = convert(String,joinpath(Pkg.dir("ExoplanetsSysSim"), "data", convert(String,get(sim_param,"window_function","DR25topwinfuncs.jld")) ) )
  setup(window_function_filename)
  add_param_fixed(sim_param,"read_window_function",true)
  @assert( size(win_func_data.window_func_array,2) == length(win_func_data.wf_durations_in_hrs) )
  @assert( size(win_func_data.window_func_array,3) == length(win_func_data.wf_periods_in_days) )
  @assert( size(win_func_data.window_func_array,1) >= maximum(win_func_data.window_function_id_arr) )
  @assert( size(win_func_data.window_func_array,1) >= win_func_data.default_wf_id )
  return win_func_data
end



function setup(filename::String)
# Reads in the window function data collected from the Kepler Completeness Products
# see Darin Ragozzine's get/cleanDR25winfuncs.jl

  if ismatch(r".jld$",filename)
    try
      wfdata = load(filename)
      window_func_array = wfdata["window_func_array"]
      wf_durations_in_hrs = wfdata["wf_durations_in_hrs"]  # TODO OPT DETAIL: Should we convert units to days here?
      wf_periods_in_days = wfdata["wf_periods_in_days"]
      sorted_quarter_strings = wfdata["sorted_quarter_strings"]
      allsortedkepids = wfdata["allsortedkepids"]
      window_function_id_arr = wfdata["window_function_id_arr"]

      global win_func_data = window_function_data(window_func_array, wf_durations_in_hrs, wf_periods_in_days, sorted_quarter_strings, 
                                                allsortedkepids, window_function_id_arr, maximum(window_function_id_arr) )

    catch
      error(string("# Failed to read window function data > ", filename," < in jld format."))
    end
  end
 
  return win_func_data
end

setup_window_function(sim_param::SimParam; force_reread::Bool = false) = setup(sim_param, force_reread=force_reread)
setup_window_function(filename::String; force_reread::Bool = false) = setup(filename, force_reread=force_reread)

function get_window_function_data()::window_function_data
   #global win_func_data
   return win_func_data
end

function get_window_function_id(kepid::Int64; use_default_for_unknown::Bool = true)::Int64
  # takes the quarter string from the stellar catalog and determines the window function id
  # from DR25topwinfuncs.jld made by Darin Ragozzine's cleanDR25winfuncs.jl script.
  const no_win_func_available::Int64 = -1        # hardcoding this in, should match convention in window function input file

  wf_id = win_func_data.window_function_id_arr[searchsortedfirst(win_func_data.allsortedkepids,kepid)] # all Kepler kepids are in allsortedkepids

  if wf_id == no_win_func_available && use_default_for_unknown
    # if a target is observed for less than 4 quarters, then it won't have a corresponding
    # window function in this list, so throw a warning and use the last window_function_id
    # which corresponds to an "averaged" window function
    warn("Window function data is not avaialble for kepid $kepid, using default.")
    wf_id = win_func_data.default_wf_id
  end
  # TODO SCI DETAIL IMPORTANT? This does not include TPS timeouts or MESthresholds (see DR25 Completeness Products)

  return wf_id 
end


function calc_period_idx(P::Float64)::Int64
  @assert(P>zero(P))
  idx = searchsortedlast(win_func_data.wf_periods_in_days,P)
  if idx == 0
     return 1
  elseif idx<length(win_func_data.wf_periods_in_days)
     if P-win_func_data.wf_periods_in_days[idx]>win_func_data.wf_periods_in_days[idx+1]-P
        idx += 1
     end 
  end
  return idx   # TODO IMPORTANT: IMPLEMENT / TEST
end

function calc_duration_idx(D::Float64)::Int64 
  # NOTE IMPORTANT: Currently assumes we left wf data in hours, so deal with that conversion here
  @assert(D>zero(D))
  const hours_in_day = 24 
  idx = searchsortedlast(win_func_data.wf_durations_in_hrs,D*hours_in_day)
  if idx == 0
     return 1
  elseif idx<length(win_func_data.wf_durations_in_hrs)
     if D*hours_in_day-win_func_data.wf_durations_in_hrs[idx]>win_func_data.wf_durations_in_hrs[idx+1]-D*hours_in_day
        idx += 1
     end
  end
  return idx   # TODO IMPORTANT: IMPLEMENT / TEST
end


function eval_window_function(wf_idx::Int64=-1; Duration::Float64=0., Period::Float64=0.)::Float64
  D_idx = calc_duration_idx(Duration)
  P_idx = calc_period_idx(Period)
  wf = eval_window_function(wf_idx,D_idx,P_idx)
  # TODO IMPORTANT: Improve way deal with missing wf values for some durations. Interpolate?
  while wf<=zero(wf) && D_idx<length(win_func_data.wf_durations_in_hrs)  
     D_idx += 1
     wf = eval_window_function(wf_idx,D_idx,P_idx)
  end
  return wf
end

function eval_window_function(wf_idx::Int64, D_idx::Int64, P_idx::Int64)::Float64
   global win_func_data
   #@assert(1<=wf_idx<maximum(win_func_data.window_function_id_arr))
   #@assert(1<=P_idx<=length(win_func_data.wf_periods_in_days))
   #@assert(1<=D_idx<=length(win_func_data.wf_durations_in_hrs))
   return win_func_data.window_func_array[wf_idx,D_idx,P_idx]     
end

#Object for storing data necessary for OSD_interpolator
immutable OSD_data
  allosds::Array{Float64,3}
  # allosds::Array{Float32,3}  # For subset OSD files
  kepids::Array{Float64,1}
  periods_length::Int64
  durations_length::Int64
  grid::Array{Array{Float64,1},1}
end

function setup_OSD_interp(sim_param::SimParam)			#reads in 3D table of OSD values and sets up global variables to be used in interpolation
  global OSD_setup
  OSD_file = load(joinpath(Pkg.dir(), "ExoplanetsSysSim", "data", convert(String,get(sim_param,"osd_file","allosds.jld"))))
  allosds = OSD_file["allosds"]			#table of OSDs with dimensions: kepids,durations,periods
  periods = OSD_file["periods"][1,:]		#1000 period values corresponding to OSD values in the third dimension of the allosds table
  kepids = OSD_file["kepids"]			#kepids corresponding to OSD values in the first dimension of the allosds table
  OSD_file = 0 # unload OSD file to save memory  
  durations = [1.5,2.,2.5,3.,3.5,4.5,5.,6.,7.5,9.,10.5,12.,12.5,15.] #14 durations corresponding to OSD values in the first dimension of theh allosds table
  periods_length = length(allosds[1,1,:])
  durations_length = length(allosds[1,:,1])
  grid = Array{Float64,1}[]			#grid used in OSD_interpolator
  push!(grid, durations)   
  push!(grid, periods)
  global compareNoise = Float64[]		#testing variable used to make sure OSD_interpolator is producing reasonable snrs
  OSD_setup = OSD_data(allosds, kepids, periods_length, durations_length, grid)
  allosds = 0 # unload OSD table to save memory  
  return OSD_setup
end

function interp_OSD_from_table(kepid::Int64, period::Real, duration::Real)
  kepid = convert(Float64,kepid)
  meskep = OSD_setup.kepids			#we need to find the index that this planet's kepid corresponds to in allosds.jld
  kepid_index = findfirst(meskep, kepid)
  if kepid_index == 0
     kepid_index = rand(1:88807)		#if we don't find the kepid in allosds.jld, then we make a random one
  end
  olOSD = view(OSD_setup.allosds,kepid_index,:,:)    #use correct kepid index to extract 2D table from 3D OSD table
  # olOSD = convert(Array{Float64,2},olOSD)
  lint = Lininterp(olOSD, OSD_setup.grid)	#sets up linear interpolation object
  osd = ApproXD.eval2D(lint, [duration*24,period])[1]	#interpolates osd
  return osd
end

# function cdpp_vs_osd(ratio::Float64, cuantos::Int64)
# #testing function that takes ratios of cdpp_snr/osd_snr and plots a histogram to make sure the results are reasonable.
#   global compareNoise
#   push!(compareNoise,ratio)
#   if length(compareNoise) == cuantos
#     PyPlot.plt[:hist](compareNoise,100)
#     println("MES median: ",median(compareNoise)," MES mean: ",mean(compareNoise), " Standard deviation: ",std(compareNoise))
#     cuantos = 100000000
#   end
#   return cuantos
# end

end  # module WindowFunction



