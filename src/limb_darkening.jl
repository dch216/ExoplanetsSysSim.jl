## ExoplanetsSysSim/src/limb_darkening.jl
## (c) 2015 Eric B. Ford

# Note this file is currently not used by SysSim.

abstract LimbDarkeningParamAbstract

immutable LimbDarkeningParamQuadratic <: LimbDarkeningParamAbstract
  u1::Float64
  u2::Float64
  # Demo of how to enforce constraints on parameter when constructing class
  # TODO SCI DETAIL: Repalce with sensible limits on LD params
  LimbDarkeningParamQuadratic(_u1::Real, _u2::Real ) = (!( (-2.0<=_u1<=2.0) && (-2.0<=_u2<=2.0) )) ? error(string("Invalid quadratic limb darkening parameters: ",_u1," & ", _u2)) : new([_u1,_u2])
end


