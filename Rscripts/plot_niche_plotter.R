  
model_dir = "SHAF_mod176.5_gppl1_ts_cld"

minimal_models = c('AUS_mod400.5_gpp_gppl1_cld',
'EQAS_mod64.5_pr',
'SA_mod496.5_gpp_gppm1s_pr_ts_cld',
'TCAM_mod480.5_gpp_gppl1_pr_ts',
'SEAS_mod500.5_gpp_gppm1_pr_ts_cld_pop',
'NHAF_mod84.5_pr_cld_pop',
'CEAS_mod216.5_gppl1_pr_cld_vp',
'BOAS_mod232.5_gppm1_pr_ts_vp',
'BONA_mod448.5_gpp_gppl1_pr',
'SHAF_mod176.5_gppl1_ts_cld')#,
# 'EUME_mod80.5_pr_cld')

png(filename = paste0(fire_dir, "/niche_plots_popcrop/all_regions_occurance", "_", "crop", "_", "pop", "_", start_date, "-", end_date, ".png"), width = 250*6, height = 500*6, res=300)
# layout(matrix(c(1,2,3,4), nrow=2, byrow=T))
par(mfrow=c(5,2), mar=c(4,5,1,1), oma=c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)


for (model_dir in minimal_models){

fire_pred_filename = paste0(fire_dir,"/",output_dir, "/", model_dir, "/", fire_pred_file)
fire_pred = NcCreateOneShot(filename = fire_pred_filename, var_name = "fire")
# fire_pred$time = fire_pred$time - 15
# fire_pred$time = as.Date("2003-1-15") + 365.2524/12*(0:155)
fire_pred = NcClipTime(fire_pred, start_date, end_date)
# fire_pred$data = fire_pred$data - 0.000
# fire_pred$data[fire_pred$data < 0.00] = 0

# fire_pred = gfedl1nc

glimits = c(fire_pred$lons[1],
            fire_pred$lons[length(fire_pred$lons)],
            fire_pred$lats[1],
            fire_pred$lats[length(fire_pred$lats)])  # get limits from predicted data

slices_per_yr_pred = 365.2524/as.numeric(mean(diff(fire_pred$time[-length(fire_pred$time)])))

lat_res = mean(diff(fire_pred$lats))*111e3
lon_res = mean(diff(fire_pred$lons))*111e3
cell_area = t(matrix(ncol = length(fire_pred$lons), data = rep(lat_res*lon_res*cos(fire_pred$lats*pi/180), length(fire_pred$lons) ), byrow = F ))

fire_obs = NcCreateOneShot(filename = fire_obs_file, var_name = "ba", glimits = glimits)
fire_obs$time = as.Date("1997-1-15") + 365.2524/12*(0:239)
fire_obs = NcClipTime(fire_obs,  start_date, end_date)
fire_obs$data[is.na(fire_pred$data)] = NA

# slices_per_yr_obs = 365.2524/as.numeric(mean(diff(fire_obs$time[-length(fire_obs$time)])))
# 
# slice_pred = apply(X = fire_pred$data, FUN = function(x){mean(x, na.rm=T)}, MARGIN = c(1,2))*slices_per_yr_pred
# slice_pred = slice_pred*cell_area
# 
# slice_obs = apply(X = fire_obs$data, FUN = function(x){mean(x, na.rm=T)}, MARGIN = c(1,2))*slices_per_yr_obs
# slice_obs = slice_obs*cell_area

region = substr(x = model_dir, start = 1, stop = gregexpr("_", text = model_dir)[[1]][1]-1)
plot.niche_xy(cropnc, popnc, fire_obs, fire_pred, region_name = region, minx=0, maxx=0.8, miny=0, maxy=3.5)
#plot.niche.occurance_xy(tsnc, prnc, fire_obs, fire_pred, region_name = region, minx=-50, maxx=35, miny=0, maxy=30)
#plot.niche.occurance_xy(cldnc, gppl1nc, fire_obs, fire_pred, region_name = region, minx=0, maxx=1, miny=0, maxy=0.4)
cat(region, "\n")
}

dev.off()

# plot.niche_xy(cldnc, gppl1nc, fire_obs, fire_pred, region_name = "AUS", xlab="Cloud Cover", ylab="Cumm. GPP") 
# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "AUS")

# plot.niche_xy(cldnc, gppl1nc, fire_obs, fire_pred, region_name = "CEAS")
# plot.niche_xy(prnc, vpnc, fire_obs, fire_pred, region_name = "CEAS")
# plot.niche_xy(prnc, cldnc, fire_obs, fire_pred, region_name = "CEAS")
# plot.niche_xy(prnc, gppl1nc, fire_obs, fire_pred, region_name = "CEAS")
# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "CEAS")

# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "BOAS")
# plot.niche_xy(gppm1nc, prnc, fire_obs, fire_pred, region_name = "BOAS")
# plot.niche_xy(gppm1nc, tsnc, fire_obs, fire_pred, region_name = "BOAS")
# plot.niche_xy(prnc, vpnc, fire_obs, fire_pred, region_name = "BOAS")

# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "EQAS")

# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "SA")
# plot.niche_xy(gppm1snc, prnc, fire_obs, fire_pred, region_name = "SA")
# plot.niche_xy(cldnc, gppm1snc, fire_obs, fire_pred, region_name = "SA")

# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "SEAS")
# plot.niche_xy(gppm1nc, prnc, fire_obs, fire_pred, region_name = "SEAS")
# plot.niche_xy(tsnc, popnc, fire_obs, fire_pred, region_name = "SEAS")
# plot.niche_xy(gppm1nc, popnc, fire_obs, fire_pred, region_name = "SEAS")

# plot.niche_xy(gppl1nc, prnc, fire_obs, fire_pred, region_name = "BONA")
# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "BONA")

# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "SHAF")
# plot.niche_xy(cldnc, gppl1nc, fire_obs, fire_pred, region_name = "SHAF")

# plot.niche_xy(cldnc, popnc, fire_obs, fire_pred, region_name = "NHAF")

# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "TCAM")
# plot.niche_xy(gppl1nc, prnc, fire_obs, fire_pred, region_name = "TCAM")


# #### 2002-12
# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "SHAF")
# plot.niche_xy(cldnc, gppl1nc, fire_obs, fire_pred, region_name = "SHAF")
# plot.niche_xy(cldnc, tsnc, fire_obs, fire_pred, region_name = "SHAF")
# plot.niche_xy(gppl1nc, tsnc, fire_obs, fire_pred, region_name = "SHAF")
