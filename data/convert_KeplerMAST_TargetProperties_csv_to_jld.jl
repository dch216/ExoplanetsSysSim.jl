using DataFrames, CSV, JLD2

function write_dataframe_to_jld2(fn::AbstractString, d::DataFrame)
  f = jldopen(fn,"w")
  for s in names(d)
    if all(ismissing.(csv[s]))
       continue
    end
    if sum(ismissing.(csv[s]))==0
		write(f,String(s),collect(skipmissing(csv[s])))
	else
		write(f,String(s),csv[s])
	end
  end
  close(f)
end


#csv = CSV.read("KeplerMAST_TargetProperties.csv",types=[Int64,String,String,String,String,Int64,Float64,Int64,Float64,Int64])

#write_dataframe_to_jld2("KeplerMAST_TargetProperties.jld2",csv)

#  #=
csv = CSV.read("q1_q17_dr25_koi.csv",header=157)
write_dataframe_to_jld2("q1_q17_dr25_koi.jld2",csv)
#  =#

