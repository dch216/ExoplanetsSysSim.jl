using HDF5, DataFrames, CSV, JLD
using StatsBase, Polynomials, CurveFit
using PyPlot

kep_filename = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q17_dr25_stellar.csv")
gaia_filename = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "gaiadr2_keplerdr25_crossref.csv")
mast_filename = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "KeplerMAST_TargetProperties.csv")
tmass_filename = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "Gaia_Kepler_2MASS_mags.csv")
stellar_catalog_file_out = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1q17_dr25_gaia_m_2mass.jld")

kep_df = CSV.read(kep_filename)
gaia_df = CSV.read(gaia_filename, rows_for_type_detect=30000)
mast_df = CSV.read(mast_filename)
tmass_df = CSV.read(tmass_filename, rows_for_type_detect=30000)

dup_gaiaid = find(nonunique(DataFrame(x = gaia_df[:source_id])))
deleterows!(gaia_df, dup_gaiaid)

println("Total crossref target stars = ", length(gaia_df[:kepid]))

mag_diff = gaia_df[:phot_g_mean_mag].-gaia_df[:kepmag]
quant_arr = quantile(mag_diff, [0.067,0.933])   # 1.5-sigma cut
mag_match = find(x->quant_arr[1]<=x<=quant_arr[2], mag_diff)
gaia_df = gaia_df[mag_match,:]

gaia_col = [:kepid,:source_id,:parallax,:parallax_error,:astrometric_gof_al,:astrometric_excess_noise_sig,:phot_g_mean_mag,:bp_rp,:a_g_val,:priam_flags,:teff_val,:teff_percentile_lower,:teff_percentile_upper,:radius_val,:radius_percentile_lower,:radius_percentile_upper,:lum_val,:lum_percentile_lower,:lum_percentile_upper]
df = join(kep_df, gaia_df[:,gaia_col], on=:kepid)
df = join(df, mast_df, on=:kepid)
df = join(df, tmass_df, on=:source_id, makeunique=true)
kep_df = nothing
gaia_df = nothing

println("Total target stars (KOIs) matching magnitude = ", length(df[:kepid]), " (", sum(df[:nkoi]),")")

df[:teff] = df[:teff_val]
df[:teff_err1] = df[:teff_percentile_upper].-df[:teff_val]
df[:teff_err2] = df[:teff_percentile_lower].-df[:teff_val]
delete!(df, :teff_val)
delete!(df, :teff_percentile_upper)
delete!(df, :teff_percentile_lower)
for x in 1:length(df[:kepid])
    if !isnan(df[x,:radius_val])
        df[x,:radius_err1] = df[x,:radius_percentile_upper]-df[x,:radius_val]
        df[x,:radius_err2] = df[x,:radius_percentile_lower]-df[x,:radius_val]
        df[x,:radius] = df[x,:radius_val]
    end
end
delete!(df, :radius_val)
delete!(df, :radius_percentile_upper)
delete!(df, :radius_percentile_lower)

not_binary_suspect = (df[:astrometric_gof_al] .<= 20) .& (df[:astrometric_excess_noise_sig] .<= 5)
astrometry_good = []
for x in 1:length(df[:kepid])
    if !(ismissing(df[x,:priam_flags]))
        pflag = string(df[x,:priam_flags])
         if (pflag[2] == '0') & (pflag[3] == '0') # WARNING: Assumes flag had first '0' removed by crossref script
             push!(astrometry_good, true)
         else
             push!(astrometry_good, false)
         end
     else
         push!(astrometry_good, false)
     end
end
astrometry_good = astrometry_good .& (df[:parallax_error] .< 0.3*df[:parallax])
planet_search = df[:kepmag] .<= 16.

has_mass = .! (ismissing.(df[:mass]) .| ismissing.(df[:mass_err1]) .| ismissing.(df[:mass_err2]))
has_radius = .! (ismissing.(df[:radius]) .| ismissing.(df[:radius_err1]) .| ismissing.(df[:radius_err2]))# .| isnan.(df[:radius]) .| isnan.(df[:radius_err1]) .| isnan.(df[:radius_err2]))
#has_dens = .! (ismissing.(df[:dens]) .| ismissing.(df[:dens_err1]) .| ismissing.(df[:dens_err2]))
has_cdpp = .! (ismissing.(df[:rrmscdpp01p5]) .| ismissing.(df[:rrmscdpp02p0]) .| ismissing.(df[:rrmscdpp02p5]) .| ismissing.(df[:rrmscdpp03p0]) .| ismissing.(df[:rrmscdpp03p5]) .| ismissing.(df[:rrmscdpp04p5]) .| ismissing.(df[:rrmscdpp05p0]) .| ismissing.(df[:rrmscdpp06p0]) .| ismissing.(df[:rrmscdpp07p5]) .| ismissing.(df[:rrmscdpp09p0]) .| ismissing.(df[:rrmscdpp10p5]) .| ismissing.(df[:rrmscdpp12p0]) .| ismissing.(df[:rrmscdpp12p5]) .| ismissing.(df[:rrmscdpp15p0]))
has_rest = .! (ismissing.(df[:dataspan]) .| ismissing.(df[:dutycycle]))
has_limbd = .! (ismissing.(df[:limbdark_coeff1]) .| ismissing.(df[:limbdark_coeff2]) .| ismissing.(df[:limbdark_coeff3]) .| ismissing.(df[:limbdark_coeff4]))
has_tmass = .! (ismissing.(df[:ks_msigcom]) .| ismissing.(df[:ks_m]) .| ismissing.(df[:j_m]))

mast_cut =.&(df[:numLCEXqtrs].>0,df[:numLCqtrs].>4)

is_usable = .&(has_rest, has_cdpp, mast_cut, astrometry_good, has_tmass)#, not_binary_suspect)

df = df[find(is_usable),:]
println("Total stars (KOIs) with valid parameters = ", length(df[:kepid]), " (", sum(df[:nkoi]),")")

m_absg = df[:phot_g_mean_mag] - df[:a_g_val] + 5 + 5*log10.(df[:parallax]/1000.)
m_absk = df[:ks_m] + 5 + 5*log10.(df[:parallax]/1000.)

m_absk_err = sqrt.(df[:ks_msigcom].^2 .+ ((5./(df[:parallax]*log(10))).^2 .* df[:parallax_error].^2))
df[:ks_absm] = m_absk
df[:ks_absm_err] = m_absk_err
m_gcolor = (1.7 .<= df[:bp_rp])# .<= 5.0)
m_tcolor = (0.617 + 0.162) .<= (df[:j_m] .- df[:ks_m])

df[:radius] = 1.9515 - (0.352*m_absk) .+ (0.0168*m_absk.^2)
m_rad_err = abs.((-0.352 + 2*0.0168*m_absk) .* m_absk_err)
df[:radius_err1] = m_rad_err
df[:radius_err2] = -m_rad_err

df[:mass] = 0.5858 + (0.3872*m_absk) .- (0.1217*m_absk.^2) .+ (0.0106*m_absk.^3) .- (2.7262e-4*m_absk.^4)
m_mass_err = abs.((0.3872 - (2*0.1217*m_absk) .+ (3*0.0106*m_absk.^2) .- (4*2.7262e-4*m_absk.^3)).*m_absk_err)
df[:mass_err1] = m_mass_err
df[:mass_err2] = -m_mass_err

m_teff = df[:teff] .<= 4000
m_logg = df[:logg] .> 3
#m_teff = (2320 .<= df[:teff] .<= 3870)
m_gmag = (7.55 .<= m_absg)# .<= 16.33)
m_kmag = (4.8 .<= m_absk)

#FGK = find(.&(m_gcolor,m_gmag,m_tcolor,m_kmag))#, df[:radius].<2))
M_samp = find(.&(m_tcolor,m_kmag, df[:radius] .< 2, df[:mass] .> 0.5))
#FGK = find(m_color)
#FGK = find(.&(m_teff, m_logg, df[:radius].<2))

println("Total M stars (KOIs) with valid parameters = ", length(M_samp), " (", sum(df[M_samp, :nkoi]),")")

plot_samp = rand(1:length(df[:kepid]), 5000)
plt[:scatter](df[M_samp, :bp_rp], m_absg[M_samp], s=3, label = "M", color="red")
plt[:scatter](df[plot_samp, :bp_rp], m_absg[plot_samp], s=3, label = "All", color="black")

plt[:ylabel](L"$M_G$")
plt[:xlabel](L"$B_p - R_p$")
plt[:ylim](reverse(plt[:ylim]()))
plt[:legend]()
plt[:savefig]("m-dwarf_samp.png")

# See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
# TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
symbols_to_keep = [ :kepid, :source_id, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :teff, :phot_g_mean_mag, :bp_rp, :j_m, :ks_m, :lum_val, :rrmscdpp01p5, :rrmscdpp02p0, :rrmscdpp02p5, :rrmscdpp03p0, :rrmscdpp03p5, :rrmscdpp04p5, :rrmscdpp05p0, :rrmscdpp06p0, :rrmscdpp07p5, :rrmscdpp09p0, :rrmscdpp10p5, :rrmscdpp12p0, :rrmscdpp12p5, :rrmscdpp15p0, :dataspan, :dutycycle, :limbdark_coeff1, :limbdark_coeff2, :limbdark_coeff3, :limbdark_coeff4, :contam]
# delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
df = df[M_samp, symbols_to_keep]
tmp_df = DataFrame()    
for col in names(df)
    tmp_df[col] = collect(skipmissing(df[col]))
end
df = tmp_df

save(stellar_catalog_file_out,"stellar_catalog", df)
