#using ExoplanetsSysSim
using HDF5, JLD, DataFrames

filename = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q17_dr25_stellar.csv")
cks_filename = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "cks_physical_merged.csv")
stellar_catalog_file_out = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q17_dr25_cks_christiansen.jld")

  df = readtable(filename)
  cks_df = readtable(cks_filename)

  for x in 1:length(cks_df[:id_kic])
    temp_ind = findfirst(df[:kepid], cks_df[x,:id_kic])
    df[temp_ind,:teff] = cks_df[x,:iso_steff]
    df[temp_ind,:logg] = cks_df[x,:iso_slogg]
    df[temp_ind,:mass] = cks_df[x,:iso_smass]
    df[temp_ind,:mass_err1] = cks_df[x,:iso_smass_err1]
    df[temp_ind,:mass_err2] = cks_df[x,:iso_smass_err2]
    df[temp_ind,:radius] = cks_df[x,:iso_srad]
    df[temp_ind,:radius_err1] = cks_df[x,:iso_srad_err1]
    df[temp_ind,:radius_err2] = cks_df[x,:iso_srad_err2]
  end

  has_mass = ! (isna(df[:mass]) | isna(df[:mass_err1]) | isna(df[:mass_err2]))
  has_radius = ! (isna(df[:radius]) | isna(df[:radius_err1]) | isna(df[:radius_err2]))
  has_dens = ! (isna(df[:dens]) | isna(df[:dens_err1]) | isna(df[:dens_err2]))
  has_rest = ! (isna(df[:rrmscdpp04p5]) | isna(df[:dataspan]) | isna(df[:dutycycle]))
  in_Q1Q12 = []
  for x in df[:st_quarters]
    subx = string(x)
    subx = ("0"^(17-length(subx)))*subx
    indQ = search(subx, '1')
    if ((indQ < 1) | (indQ > 12))
      push!(in_Q1Q12, false)
    else
      push!(in_Q1Q12, true)
    end
  end
  is_FGK = []
  for x in 1:length(df[:teff])
    if ((df[x,:teff] > 4000.0) & (df[x,:teff] < 7000.0) & (df[x,:logg] > 4.0))
      push!(is_FGK, true)
    else
      push!(is_FGK, false)
    end
  end
  is_usable = has_radius & has_mass & has_rest & has_dens & is_FGK #& in_Q1Q12
  if contains(filename,"q1_q12_christiansen.jld")
    is_usable = is_usable & in_Q1Q12
  end
  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  # TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp04p5, :dataspan, :dutycycle ]
  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]

save(stellar_catalog_file_out,"stellar_catalog", df)
