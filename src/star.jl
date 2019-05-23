## ExoplanetsSysSim/src/star.jl
## (c) 2015 Eric B. Ford

#using Distributions

@compat abstract type StarAbstract end               # Check does using StarAbstract cause a significant performance hit

immutable Star <: StarAbstract                
  radius::Float64
  mass::Float64
  flux::Float64                      # relevant once have multiple stars in one target
  ld::LimbDarkeningParamAbstract         # TODO SCI DETAIL: add limb darkening param?
  id::Int64                          # id for looking up properties in stellar catalog

end
#typealias SingleStar Star
SingleStar = Star

immutable BinaryStar <: StarAbstract
  primary::Star
  secondary::Star
  orbit::Orbit
end

immutable MultipleStar <: StarAbstract                      # Will we want to handle triple, quad systems?
  component::Vector{StarAbstract}
  orbit::Vector{Orbit}
end

flux(s::Star) = s.flux                              # Demo of how to specify function behavior that depends on the derived type
flux(s::BinaryStar) = s.primary.flux + s.secondary.flux
flux(s::MultipleStar) = sum( flux, s.component)


mass(s::Star) = s.mass
mass(s::BinaryStar) = s.primary.mass + s.secondary.mass
mass(s::MultipleStar) = sum( mass, s.component)::Float64

function generate_stars(sim_param::SimParam)
   generate_star = get_function(sim_param,"generate_star")
   num_target_stars = get_int(sim_param,"num_targets_sim_pass_one")
   star_list = Array{StarAbstract}(num_target_stars)
   for i in 1:num_target_stars
     s =  generate_star(sim_param)
     star_list[i] = s
     #star_list[i] = generate_star(sim_param)
   end
  return star_list
end

function generate_star_dumb(sim_param::SimParam) 
  r = rand(Uniform(0.8,1.3))::Float64                
  m = rand(Normal(r,0.1))::Float64
  while m<0.0
    m = rand(Normal(r,0.1))::Float64
  end
  f = rand(Normal(1.0,0.1))::Float64
  while f<0.0
    f = 1.0+0.1*randn()
  end
  ld = LimbDarkeningParamQuadratic(0.4603,0.2291)   # TODO: Once we implement limb darkening
  return SingleStar(r,m,f,ld,0) 
  #return SingleStar(r,m,f,0) 
end


function test_star_constructors(sim_param::SimParam)
  star_tmp = generate_star_dumb(sim_param)
  f1 = flux(star_tmp)
  f2 = flux(BinaryStar(star_tmp,star_tmp,Orbit(10.0,0.0,0.0,0.0,0.0,0.0)))
  f4 = flux(MultipleStar([star_tmp for i in 1:4], [Orbit(10.0,0.0,0.0,0.0,0.0,0.0) for i in 1:4]) )
  # println("# Fluxes: ", f1, " ", f2, " ", f4)
  star_list = generate_stars(sim_param)
  return true
end

