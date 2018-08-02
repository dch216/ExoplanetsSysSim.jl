## ExoplanetsSysSim/src/stellar_table.jl
## (c) 2015 Eric B. Ford

module StellarTable
using ExoplanetsSysSim
using DataArrays
using DataFrames
using CSV
using JLD

#if VERSION >= v"0.5-"
#  import Compat: UTF8String, ASCIIString
#end

export setup_star_table, star_table, num_usable_in_star_table, set_star_table

df = DataFrame()
usable = Array{Int64}(0)
#data = Array{Float64}(0,0)
#colid = Dict()
        
function setup(sim_param::SimParam; force_reread::Bool = false)
  global df
  if haskey(sim_param,"read_stellar_catalog") && !force_reread
     return df
     #return data
  end
  stellar_catalog_filename = convert(String,joinpath(Pkg.dir("ExoplanetsSysSim"), "data", convert(String,get(sim_param,"stellar_catalog","q1_q17_dr24_stellar.csv")) ) )
  df = setup(stellar_catalog_filename)
  add_param_fixed(sim_param,"read_stellar_catalog",true)
  add_param_fixed(sim_param,"num_kepler_targets",num_usable_in_star_table())
  return df  
end

function setup(filename::String; force_reread::Bool = false)
  global df, usable
  if ismatch(r".jld$",filename)
  try 
    data = load(filename)
    df = data["stellar_catalog"]
    usable = data["stellar_catalog_usable"]
    Core.typeassert(df,DataFrame)
    Core.typeassert(usable,Array{Int64,1})
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in jld format."))
  end
  else
  try 
    #df = readtable(filename)
    df = CSV.read(filename,nullable=true)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in ascii format."))
  end

  #=  Hoping we can soon get rid of this mess with map below.  There's probably a more efficient way.
  has_mass = .!(ismissing.(df[:mass]) .| ismissing.(df[:mass_err1]) .| ismissing.(df[:mass_err2]))
  has_radius = .!(ismissing.(df[:radius]) .| ismissing.(df[:radius_err1]) .| ismissing.(df[:radius_err2]))
  has_dens = .!(ismissing.(df[:dens]) .| ismissing.(df[:dens_err1]) .| ismissing.(df[:dens_err2]))
  has_cdpp = .!(ismissing.(df[:rrmscdpp04p5]) .| ismissing.(df[:rrmscdpp01p5]) .| ismissing.(df[:rrmscdpp15p0])) # TODO: WARNING: Doesn't include all CDPPs
  has_rest = .!(ismissing.(df[:rrmscdpp04p5]) .| ismissing.(df[:dataspan]) .| ismissing.(df[:dutycycle]))
  is_usable = .&(has_mass, has_radius, has_dens, has_cdpp, has_rest)
  =# 

  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  # Now we read in all CDPP's, so can interpolate to transit duration
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp01p5, :rrmscdpp02p0, :rrmscdpp02p5, :rrmscdpp03p0, :rrmscdpp03p5, :rrmscdpp04p5, :rrmscdpp05p0, :rrmscdpp06p0, :rrmscdpp07p5, :rrmscdpp09p0, :rrmscdpp10p5, :rrmscdpp12p0, :rrmscdpp12p5, :rrmscdpp15p0, :cdppslplong, :cdppslpshrt, :dataspan, :dutycycle ]

  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  is_usable = [ !any(ismissing.([ df[i,j] for j in 1:size(df,2) ])) for i in 1:size(df,1) ]
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]
  end
  return df
  #global data = convert(Array{Float64,2}, df) # df[usable, symbols_to_keep] )
  #global colid = Dict(zip(names(df),[1:length(names(df))]))
  #return data
end

setup_star_table(sim_param::SimParam; force_reread::Bool = false) = setup(sim_param, force_reread=force_reread)
setup_star_table(filename::String; force_reread::Bool = false) = setup(filename, force_reread=force_reread)

function num_usable()
  global usable
  @assert typeof(usable) == Array{Int64,1}
  length(usable)
end
num_usable_in_star_table() = num_usable()

function idx(i::Integer)
  global usable
  usable[i]
end

#function col( sym::Symbol )
#  global colid
#  colid[sym]
#end

function star_table(i::Integer, sym::Symbol)
  global df, usable
  return df[i,sym]
  #return df[usable[i],sym]
end

#function star_table_data(i::Integer, sym::Symbol)
#  global data
#  return data[i,col(sym)]
#end

function star_table(i::Integer)
  global df, usable
  return df[i,:]
  #return df[usable[i],:]
end

function star_table(i::Integer, sym::Vector{Symbol})
  global df, usable
  return df[i,sym]
  #return df[usable[i],sym]
end

function star_table(i::Vector{Integer}, sym::Symbol)
  global df, usable
  return df[i,sym]
  #return df[usable[i],sym]
end

function star_table(i::Vector{Integer}, sym::Vector{Symbol})
  global df, usable
  return df[i,sym]
  #return df[usable[i],sym]
end

function set_star_table(df2::DataFrame)
  global df
  df = df2
end

function set_star_table(df2::DataFrame, usable2::Array{Int64,1})
  global df, usable
  df = df2
  usable = usable2
end

end # module StellarTable

# using ExoplanetsSysSim.StellarTable

function generate_star_from_table(sim_param::SimParam, id::Integer)  # WARNING:  To be renamed once there's a working/tested version that uses a stellar catalog with GAIA data
  mu_r = StellarTable.star_table(id,:radius)
  sig_r1 = StellarTable.star_table(id,:radius_err1)
  sig_r2 = StellarTable.star_table(id,:radius_err2)
  z = randn() 
  r = mu_r + (z>0) ?  z*sig_r1 : z*sig_r2
  m = rand(Normal(r,0.1))::Float64
  while m<0.0
    m = rand(Normal(r,0.1))::Float64
  end
  f = rand(Normal(1.0,0.1))::Float64
  while f<0.0
    f = 1.0+0.1*randn()
  end
  # ld = LimbDarkeningParamQuadratic(0.5,0.5)
  return SingleStar(r,m,f,id)
end

function generate_star_from_table(sim_param::SimParam)
  id = rand(1:StellarTable.num_usable_in_star_table())
  generate_star_from_table(sim_param, id)
end

