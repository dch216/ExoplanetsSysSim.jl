using DataFrames, CSV

function interpolate(x_in::AbstractArray{T1}, y_in::AbstractArray{T2}, x::T1) where {T1<:Number, T2<:Number}
    @assert x_in[1] <= x <= x_in[2]
    ((x_in[2]-x)*y_in[1]+(x-x_in[1])*y_in[2]) / (x_in[2]-x_in[1])
end

function interp_quantile_list(x::AbstractArray{T,1},q::Real;
            quantile_list=1:length(x) ) where T<:Number
   @assert issorted(quantile_list)
   @assert length(x) == length(quantile_list)
   idx_hi = searchsortedfirst(quantile_list,q)
   idx_lo = idx_hi-1
   if idx_hi>length(x)
     return x[end]
   elseif idx_lo<1
       return x[1]
   elseif q==quantile_list[idx_hi]
     return x[idx_hi]
   end
    interpolate(view(quantile_list,idx_lo:idx_hi),view(x,idx_lo:idx_hi), q)
end

mcmc_filename = "dr25_koi_mcmc_quant.csv"
koi_filename = "q1_q17_dr25_koi_original.csv"
mcmc_df = CSV.read(mcmc_filename)
tmp_koi_cat = readlines(koi_filename)
tmp_ind = 1
num_skip = 1
while tmp_koi_cat[tmp_ind][1] == '#'
    num_skip += 1
    tmp_ind += 1
end
df = CSV.read(koi_filename, allowmissing=:all, header=num_skip, rows_for_type_detect = 6000)
df = join(df, mcmc_df, on=:kepoi_name, kind=:left)

df[:koi_depth_cat] = map(x -> x, df[:koi_depth])
df[:koi_depth_err1_cat] = map(x -> x, df[:koi_depth_err1])
df[:koi_depth_err2_cat] = map(x -> x, df[:koi_depth_err2])
df[:koi_duration_cat] = map(x -> x, df[:koi_duration])
df[:koi_duration_err1_cat] = map(x -> x, df[:koi_duration_err1])
df[:koi_duration_err2_cat] = map(x -> x, df[:koi_duration_err2])

num_quantiles = 99
quantile_list = range(1/(num_quantiles+1), 1/(num_quantiles+1), num_quantiles)
q_hi = 0.5+0.34134
q_mid = 0.5
q_lo = 0.5-0.34134

for j in 1:length(df[:kepoi_name])
    if !ismissing(df[j, :depth_mean]) && df[j, :depth_mean] >= 0
        depth_quantiles = Array{Float64, 1}()
        duration_quantiles = Array{Float64, 1}()
        for i in 1:num_quantiles
            push!(depth_quantiles, df[j, Symbol("depth_q" * string(i))]*1e6)
            push!(duration_quantiles, df[j, Symbol("duration_q" * string(i))]*24)
        end
        df[j, :koi_depth] = interp_quantile_list(depth_quantiles,q_mid,quantile_list=quantile_list)
        df[j, :koi_depth_err1] =  interp_quantile_list(depth_quantiles,q_hi,quantile_list=quantile_list) - interp_quantile_list(depth_quantiles,q_mid,quantile_list=quantile_list)
        df[j, :koi_depth_err2] =  interp_quantile_list(depth_quantiles,q_lo,quantile_list=quantile_list) - interp_quantile_list(depth_quantiles,q_mid,quantile_list=quantile_list)
        df[j, :koi_duration] =  interp_quantile_list(duration_quantiles,q_mid,quantile_list=quantile_list)
        df[j, :koi_duration_err1] = interp_quantile_list(duration_quantiles,q_hi,quantile_list=quantile_list) - interp_quantile_list(duration_quantiles,q_mid,quantile_list=quantile_list)
        df[j, :koi_duration_err2] = interp_quantile_list(duration_quantiles,q_lo,quantile_list=quantile_list)- interp_quantile_list(duration_quantiles,q_mid,quantile_list=quantile_list)
    end
end

delete!(df, :depth_mean)
delete!(df, :depth_std)
delete!(df, :duration_mean)
delete!(df, :duration_std)
delete!(df, :radr_gt1_cnt)
delete!(df, :b_gt1_cnt)
for i in 1:num_quantiles
    delete!(df, Symbol("depth_q" * string(i)))
    delete!(df, Symbol("duration_q" * string(i)))
end

CSV.write("q1_q17_dr25_koi.csv", df)
