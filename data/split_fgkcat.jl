using JLD, DataFrames

star_cat = load("q1q17_dr25_gaia_fgk.jld", "stellar_catalog")

star_catF = star_cat[(6000.<=star_cat[:teff].<7300),:]
star_catG = star_cat[(5300.<=star_cat[:teff].<6000),:]
star_catK = star_cat[(3900.<=star_cat[:teff].<5300),:]

save("q1q17_dr25_gaia_f.jld", "stellar_catalog", star_catF)
save("q1q17_dr25_gaia_g.jld", "stellar_catalog", star_catG)
save("q1q17_dr25_gaia_k.jld", "stellar_catalog", star_catK)
