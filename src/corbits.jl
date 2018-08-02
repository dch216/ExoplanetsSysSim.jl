## ExoplanetsSysSim/src/corbits.jl
## (c) 2015 Eric B. Ford

if VERSION >= v"0.4.0-"
   global const LIB_CORBITS = Libdl.find_library(["libcorbits.so"],[".",joinpath(Pkg.dir(),"ExoplanetsSysSim/"),joinpath(Pkg.dir(),"CORBITS/"),"/usr/local/lib"])  # WARNING: Assumes can find libcorbits.so
else
   global const LIB_CORBITS = find_library(["libcorbits.so"],[".",joinpath(Pkg.dir(),"ExoplanetsSysSim/"),joinpath(Pkg.dir(),"CORBITS/"),"/usr/local/lib"])  # WARNING: Assumes can find libcorbits.so
end

#if !isdefined(is_corbits_initialized)
#  @windows_only println("# Warning:  Can't call external C libraries from Windows.  CORBITS won't work")
#end
#is_corbits_initialized = true

# Call CORBITS's function prob_of_transits_approx_arrays(a, r_star, r, e, Omega, omega, inc, use)
# Returns a Cdouble (aka Float64)
# For documentation see https://github.com/jbrakensiek/CORBITS
function prob_of_transits_approx(a::Vector{Cdouble},r_star::Cdouble,r::Vector{Cdouble}, e::Vector{Cdouble},
                                 Omega::Vector{Cdouble}, omega::Vector{Cdouble}, inc::Vector{Cdouble}, use::Vector{Cint})
  @assert(length(a) == length(r) == length(e) == length(Omega) == length(omega) == length(inc) >= length(use) )
  @assert(length(use) >=1 )
  num_pl = length(use)
  return ccall( (:prob_of_transits_approx_arrays, LIB_CORBITS), Cdouble,
      (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint),
       a, r_star, r, e, Omega, omega, inc, use, num_pl)
end

function test_corbits()
  #@compat is_windows() ? error("# CORBITS won't work on windows") : nothing
  a =  Cdouble[0.05, 0.15]
  r_star = convert(Cdouble,0.005)
  r = Cdouble[0.0001,0.0001]
  ecc = Cdouble[0.02, 0.1]
  Omega = Cdouble[0.0, 0.0]
  omega = Cdouble[ 0.0, 0.5]
  inc = Cdouble[pi/2, pi/2]
  use_pl = Cint[1,1]
  prob_of_transits_approx(a, r_star, r, ecc, Omega, omega, inc, use_pl)
end


