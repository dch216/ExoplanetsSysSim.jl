# This directory contains data files likely to be commonly used with SysSim
For more info about the file contents, see:  
http://exoplanetarchive.ipac.caltech.edu/docs/Kepler_stellar_docs.html  
http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html

# Downloading stellar catalog files
Some csv files were too big for github, so you should download them yourself.  
The commands in download_stellar_tables.sh were used to download the table of stellar properties in CSV format.  
Reading CSV files is significantly slower than JLD files, so it's recommended to make a jld version for files to be used often.

# Downloading koi catalog files
We include the csv version that we've downloaded in the ExoplanetSySim repository  
I don't know how to script the downloading of the KOI tables from:  
http://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=koi  
When you download via the interactive table, I suggest selecting all rows and all columns to make sure you don't miss something  
See info at:  
http://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html  
http://exoplanetarchive.ipac.caltech.edu/docs/Kepler_KOI_docs.html  
http://exoplanetarchive.ipac.caltech.edu/docs/PurposeOfKOITable.html

Q: What are the contents of the JLD files?  
A: The JLD files are repackaged stellar catalog data from the different csv versions available on the Exoplanet Archive.  
   These files only keep the following columns and only include targets which have valid values in all of these columns:  
   kepid, mass, mass_err1, mass_err2, radius, radius_err1, radius_err2, dens, dens_err1, dens_err2, rrmscdpp04p5, dataspan, dutycycle

   Some JLD files have had additional filters and are labeled with the appropriate suffix:  
* christiansen = Filtered on the FGK selection in the Christiansen et al. (2015) study [Teff = 4000-7000K; log g > 4; observed at least once in Q1-Q12 (only applies to Q1-Q16 file)]
* cks = Incorporates improved stellar parameters from the California Kepler Survey
* stellar_type = Filtered on stellar type as defined by the SAG13 working group [K: 3900K <= Teff < 5300K; G: 5300K <= Teff < 6000K; F: 6000K <= Teff < 7300K]

Q: What is the difference between q1_q16_koi_cand.csv and q1_q16_koi.csv?  
A: Danley Hsu generated q1_q16_koi_cand.csv based on q1_q16_koi.csv as well as previous KOI catalogs.  The primary difference is that he has used disposition information from previous KOI catalogs that was not reflected in the raw Q16 KOI catalog to select all KOIs with the "candidate" disposition.


KeplerMAST_TargetProperties.csv contains a summary of key properties for Kepler Targets that are not included in the other catalogs. 
This information was gathered from the Kepler portion of the Mikulski Archive for Space Telescopes (MAST) at archive.stsci.edu/kepler using a CasJobs query. 
The large database of all Kepler MAST data was ingested and resummarized into several columns. A focus was on identifying whether the lightcurve of a 
particular target ("kepid") during a particular month/quarter was obtained as part of the Exoplanet search (e.g., under the "Investigation ID" of "EX") or not. 
We also gathered information on the availability of Short Cadence data. Thus the "LCEXst_quarters" is the "quarter string" for Quarters 1-17 with a 1 if
this target was observed (with LC) under the EX Investigation ID and 0 if not. (Quarter 0 is not included in the quarter strings, following the stellar catalog
convention.) "LCst_quarters" reports observations of any kind and should be identical to the quarter string from the stellar catalog. Short Cadence quarter
strings are also provided, but here the digit 0, 1, 2, or 3 represents whether 0, 1, 2, or 3 months were obtained with Short Cadence data (both under the "EX"
program and overall). We also summarize these quarter strings with single numbers indicating the number of quarters obtained under "EX" / overall in LC / SC 
mode. (Note that numSCEXqtrs and numSCqtrs are floats since fractional quarters are possible.) We also elected to store the skygroup ID for each target. 
