fire_dir = "/home/jaideep/codes/PureNN_fire"

library(ncdf4)

# dataset = "eval"
# datg = read.fireData_gfed(dataset = dataset, dir=paste0(fire_dir, "/",output_dir, "/", model_dir), regions=reg)

regnc = nc_open(paste0(data_dir,"/Fire_BA_GFED4.1s/ancil/basis_regions_1deg.nc"))
regions = ncvar_get(regnc, varid = "region")

pfts_modis = 1:11
pftnames_modis = c("Barren", 
                   "NLE", 
                   "BLE",  
                   "NLD", 
                   "BLD",
                   "MX",
                   "CS",
                   "OS",
                   "WS",
                   "S",
                   "GR",
                   "AG")

# dft_file = nc_open(filename = paste0(data_dir, "/forest_type/MODIS/dft_MODIS_global_12lev_agri-bar_lt0.5_1deg.nc"))
# dft = ncvar_get(dft_file, "ft")


tsnc  = NcCreateOneShot(paste0(data_dir,"/Fire/ts.2002-2015.nc"), var_name = "ts", tlim=c(start_date, end_date))
vpnc  = NcCreateOneShot(paste0(data_dir,"/Fire/vp.2002-2015.nc"), var_name = "vp", tlim=c(start_date, end_date))
cldnc  = NcCreateOneShot(paste0(data_dir,"/Fire/cld.2002-2015.nc"), var_name = "cld", tlim=c(start_date, end_date))
prnc  = NcCreateOneShot(paste0(data_dir,"/Fire/pr.2002-2015.nc"), var_name = "pr", tlim=c(start_date, end_date))
gppnc  = NcCreateOneShot(paste0(data_dir,"/Fire/gpp.2002-2015.nc"), var_name = "gpp", tlim=c(start_date, end_date))
gppl1nc  = NcCreateOneShot(paste0(data_dir,"/Fire/gppl1.2002-2015.nc"), var_name = "gppl1", tlim=c(start_date, end_date))
gppm1nc  = NcCreateOneShot(paste0(data_dir,"/Fire/gppm1.2002-2015.nc"), var_name = "gppm1", tlim=c(start_date, end_date))
gppm1snc  = NcCreateOneShot(paste0(data_dir,"/Fire/gppm1s.2002-2015.nc"), var_name = "gppm1s", tlim=c(start_date, end_date))
ltnnc  = NcCreateOneShot(paste0(data_dir,"/Fire/ltn.2002-2015.nc"), var_name = "ltn", tlim=c(start_date, end_date))
popnc  = NcCreateOneShot(paste0(data_dir,"/Fire/pop.2002-2015.nc"), var_name = "pop", tlim=c(start_date, end_date))
ftnc = NcCreateOneShot(paste0(data_dir,"/forest_type/MODIS/ftmap_modis_global_1deg_12levs.nc"), var_name = "ft")
ftnc = NcRepeat(ftnc, popnc$time)

prnc$data = log(1+prnc$data)
ltnnc$data = log(1e-5+ltnnc$data)
popnc$data[popnc$data > 1e18]  = NA
popnc$data  = log10(1+popnc$data)
cropnc = NcSelectLevel(ftnc,12)
cropnc$name = "crop"

# gfedl1nc  = NcCreateOneShot(paste0(data_dir,"/Fire/gfedl1.2002-2015.nc"), var_name = "gfedl1")
# gfedl1nc$data[gfedl1nc$data > 1e18]  = NA
# gfedl1nc  = NcClipTime(gfedl1nc,  start_date, end_date)




# model_dir = "SHAF_mod176.5_gppl1_ts_cld"
# 
# minimal_models = c('AUS_mod400.5_gpp_gppl1_cld',
#                    'EQAS_mod64.5_pr',
#                    'SA_mod496.5_gpp_gppm1s_pr_ts_cld',
#                    'TCAM_mod480.5_gpp_gppl1_pr_ts',
#                    'SEAS_mod500.5_gpp_gppm1_pr_ts_cld_pop',
#                    'NHAF_mod84.5_pr_cld_pop',
#                    'CEAS_mod216.5_gppl1_pr_cld_vp',
#                    'BOAS_mod232.5_gppm1_pr_ts_vp',
#                    'BONA_mod448.5_gpp_gppl1_pr',
#                    'SHAF_mod176.5_gppl1_ts_cld')#,
# # 'EUME_mod80.5_pr_cld')
  
fire_obs = NcCreateOneShot(filename = fire_obs_file, var_name = "ba")
fire_obs$time = as.Date("1997-1-15") + 365.2524/12*(0:239)
fire_obs = NcClipTime(fire_obs,  start_date, end_date)
# fire_obs$data[is.na(fire_pred$data)] = NA

# layout(matrix(c(1,2,3,4), nrow=2, byrow=T))
# par(mfrow=c(5,2), mar=c(4,5,1,1), oma=c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
    
for (region in names(regions_list)){
  # region = "GLOBE"
  if (region %in% c("SHAF", "SA", "AUS")){
    vars = list(ltnnc, tsnc, prnc, vpnc, cldnc, gppl1nc, gppm1snc, popnc, cropnc)
  }else{
    vars = list(ltnnc, tsnc, prnc, vpnc, cldnc, gppl1nc, gppm1nc, popnc, cropnc)
  } 
  
  png(filename = paste0(fire_dir, "/niche_plots_regionwise_3/", region, "_", start_date, "-", end_date, ".png"), width = 950*3, height = 900*3, res=300)
  par(mfrow=c(length(vars), length(vars)), mar=c(2,2,0.5,0.5),  mgp=c(0.5,0.5,0))
    
  cat("\n", region, ": ")
  for (i in 1:length(vars)){
      for(j in 1:length(vars)){
        if (i < j) plot.niche_xy(vars[[i]], vars[[j]], fire_obs, region_name = region, xaxt="n", yaxt="n")
        else if (i>j) plot.niche.occurance_xy(vars[[j]], vars[[i]], fire_obs, region_name = region, xlab="", ylab="")
        else{
          plot.niche_x(vars[[i]], fire_obs, region_name = region, xaxt="n", yaxt="n")
        }
        mtext(text = region, side = 3, line = 1, outer=T)
        cat(".")
      }
  }
  
  dev.off()
}


# png(filename = paste0(fire_dir, "/niche_plots_regionwise/", region, "_pft_", start_date, "-", end_date, ".png"), width = 1250*3, height = 900*3, res=300)
# par(mfrow=c(length(regions_list), length(2:12)), mar=c(2,2,0.5,0.5),  mgp=c(0.5,0.5,0))
# for (region in names(regions_list)){
#   cat("\n", region, ": ")
#   for (i in 2:12){
#       plot.niche_x(NcSelectLevel(ftnc,i), fire_obs, region_name = region, xaxt="n", yaxt="n")
#       if (i==2) mtext(text = pftnames_modis[i], side = 3, line = 1)
#       cat(".")
#   }
#   dev.off()
# }

  # #### 2002-12
# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "SHAF")
# plot.niche_xy(cldnc, gppl1nc, fire_obs, fire_pred, region_name = "SHAF")
# plot.niche_xy(cldnc, tsnc, fire_obs, fire_pred, region_name = "SHAF")
# plot.niche_xy(gppl1nc, tsnc, fire_obs, fire_pred, region_name = "SHAF")
