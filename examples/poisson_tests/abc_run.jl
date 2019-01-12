## ExoplanetsSysSim/examples/hsu_etal_2018/abc_run.jl
## (c) 2018 Danley C. Hsu & Eric B. Ford
# Script for producing Q1-Q16 FGK planet candidate occurrence rate estimates

include("abc_setup.jl")

using SysSimABC
using ExoplanetsSysSim
using JLD
using StatsBase

out2txt = false # Write occurrence rates & densities to text files
expandpart = true # Expand final generation for robust posteriors

println("Setting up simulation...")
@time abc_plan = setup_abc()
println("")
println("Running simulation...")
@time output = run_abc(abc_plan)
println("")
println("Running simulation (part 2)...")
@time abc_plan = setup_abc_p2(abc_plan)
@time output = run_abc(abc_plan, output)
#@time abc_plan = change_distance()
#@time output = run_abc(abc_plan, output)
println("")

save(string("test-pop-out.jld"), "output", output, "ss_true", EvalSysSimModel.get_ss_obs())

if expandpart
    println("Expanding to large generation...")
    @time theta_largegen, weights_largegen = run_abc_largegen(output, EvalSysSimModel.get_ss_obs(), output.accept_log.epsilon[end-1], npart=200)
    println("")

    save(string("test-pop-out.jld"), "output", output, "ss_true", EvalSysSimModel.get_ss_obs(), "theta_largegen", theta_largegen, "weights_largegen", weights_largegen)
end

if out2txt
    file_rate = open("rate_output.txt", "w")
    file_dens = open("dens_output.txt", "w")
end

limitP = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
limitR = get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1})
const r_dim = length(limitR)-1

if expandpart
    #weight_vec = aweights(weights_largegen)
    weight_vec = aweights(fill(1.0, length(weights_largegen)))
else
    #weight_vec = aweights(output.weights)
    weight_vec = aweights(fill(1.0, length(output.weights)))
end

for p_ind = 1:(length(limitP)-1)
    col_ind = (p_ind-1)*(r_dim+1)+1
    for r_ind = 1:(length(limitR)-1)
        dens_denom = 1.0/(log2(limitP[p_ind+1])-log2(limitP[p_ind]))/(log2(limitR[r_ind+1])-log2(limitR[r_ind]))

        bin_ind = (p_ind-1)*(r_dim+1)+r_ind+1
        if expandpart
            quant_arr = quantile(theta_largegen[bin_ind,:].*theta_largegen[col_ind,:], weight_vec, [0.1587, 0.5, 0.8413])
        else
            quant_arr = quantile(output.theta[bin_ind,:].*output.theta[col_ind,:], weight_vec, [0.1587, 0.5, 0.8413])
        end

        println("-----------------------------")
        println("Orbital Period (day) = ", string(limitP[p_ind:p_ind+1]), " / Planet Radius (R_earth) = ", string(limitR[r_ind:r_ind+1]/ExoplanetsSysSim.earth_radius))
        println("")
        println("Rate = ", string(quant_arr[2], " + ", quant_arr[3]-quant_arr[2], " - ", quant_arr[2]-quant_arr[1]))
        println("Density = ", string(quant_arr[2]*dens_denom, " + ", (quant_arr[3]-quant_arr[2])*dens_denom, " - ", (quant_arr[2]-quant_arr[1])*dens_denom))

        if out2txt
            write(file_rate, string(quant_arr[2], " + ", quant_arr[3]-quant_arr[2], " - ", quant_arr[2]-quant_arr[1], "\n"))
            write(file_dens, string(quant_arr[2]*dens_denom, " + ", (quant_arr[3]-quant_arr[2])*dens_denom, " - ", (quant_arr[2]-quant_arr[1])*dens_denom, "\n"))
        end
    end
end

if out2txt
    close(file_rate)
    close(file_dens)
end

println("-----------------------------")
println("")
println(EvalSysSimModel.get_ss_obs())
