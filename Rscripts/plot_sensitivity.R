rm(list = ls())

## CREATE timeseries	
library(ncdf4)

fire_dir = "/home/jaideep/codes/PureNN_fire"
output_dir = "output_globe_sensitivity_3"
model_dir = "AUS_mod440.5_gpp_gppl1_ts_cld_vp"

library(sf)

shp = st_read("D:/Water LUE Project/data/continent_shapefile/continent shapefile/continent.shp")


# cat("model_dir = ", model_dir, "\n")
# cat("output_dir = ", output_dir, "\n")

source(paste0(fire_dir, "/Rscripts/utils.R"))

mha_per_m2 = 0.0001/1e6

# for (mod in c("gfed_xcf")){ #}, "xdxl", "xlmois", "xts", "xpop", "xrh", "xwsp", "xrh_lmois")){
#
#   for (i in 1:10){
    # model = paste0(mod, "_", i)

    setwd(paste0(fire_dir,"/",output_dir, "/", model_dir ))
    system("mkdir -p figures", ignore.stderr = T)

    new = NcCreateOneShot(filename = "fire.vpd+.nc", var_name = "fire")
    orig = NcCreateOneShot(filename = "fire.nc", var_name = "fire")
    new$data = new$data*12
    orig$data = orig$data*12
    tseffect = new
    tseffect$data = (new$data - orig$data)/(0.001+orig$data)*100
    # tseffect$data[tseffect$data > 100]=100
    
    glimits = c(tseffect$lons[1],
                tseffect$lons[length(tseffect$lons)],
                tseffect$lats[1],
                tseffect$lats[length(tseffect$lats)])  # get limits from predicted data

    
    lat_res = mean(diff(tseffect$lats))*111e3
    lon_res = mean(diff(tseffect$lons))*111e3
    cell_area = t(matrix(ncol = length(tseffect$lons), data = rep(lat_res*lon_res*cos(tseffect$lats*pi/180), length(tseffect$lons) ), byrow = F ))

    lim = 20 #max(abs(tseffect$data), na.rm=T)
    cols = createPalette(c("blue4", "skyblue", "cyan", "#f0f0f0", "yellow", "orange", "red"),c(-lim, -lim/2, -0.1*lim, 0, 0.1*lim, lim/2, lim)*1000, n = 1000) #gfed
    cols = createPalette(c("purple", "blue", "cyan", "white", "greenyellow", "yellow","orange", "red"),c(-0.4, -0.3, -0.1, 0, 0.1, 0.18, 0.25,0.4)*1000, n = 1000) #gfed
    tseffect$data[tseffect$data > lim] = lim
    tseffect$data[tseffect$data < -lim] = -lim
    
    png(filename = paste0("figures/vpd_effect_AUS.png"),res = 300,width = 800*3,height = 405*3*1.5) # 520 for sasplus, india, 460 for SAS
    layout(matrix(c(1,1,
                    1,1,
                    1,1,
                    2,2), ncol=2,byrow = T))  # vertical
    par(mar=c(4,5,4,1), oma=c(1,1,1,1), cex.lab=1.5, cex.axis=1.5)

    image(tseffect$lon, tseffect$lat, tseffect$data, col = cols, zlim = c(-lim,lim), xlab="Longitude",ylab = "Latitude")
    mtext(side=3, line=1, text = sprintf("Total increase = %.2f", sum(tseffect$data*cell_area, na.rm=T)*mha_per_m2))
    plot(shp$geometry, add=T, lwd=0.1)  
    
    image(x=seq(-lim, lim, length.out=1000), y=1, z=t(matrix(data = seq(-lim, lim, length.out=1000), nrow=1)), yaxt="n", xlab="% Increase in burned fraction", col=cols)
    
    
    dev.off()

    
    
    #   }
#   cat("\n\n\n")
# }

    
dft = NcCreateOneShot("C:/Users/Jaideep/Desktop/modis/dft_MODIS_global_12lev_agri-bar_lt0.5_1deg.nc", var_name = "ft")

for (i in 0:12){
  cat(i, " ", 
      length(tseffect$data[dft$data==i]), " ", 
      sum(orig$data[dft$data==i]*cell_area[dft$data==i], na.rm=T)*mha_per_m2, " ", 
      sum(new$data[dft$data==i]*cell_area[dft$data==i], na.rm=T)*mha_per_m2, " ", 
      sum(tseffect$data[dft$data==i]*cell_area[dft$data==i], na.rm=T)*mha_per_m2, "\n")
}


### NEW SENSITIVITY MAPS - 23/12/2019

rm(list = ls())

## CREATE timeseries	
library(ncdf4)

fire_dir = "/home/jaideep/codes/PureNN_fire"
output_dir = "merged_models"
model_dir = "merged_sens_4"

# fire_dir = "/home/jaideep/codes/PureNN_fire"  
# output_dir = "output_globe"
# model_dir = "AUS_mod1464.6.2_gppm1s_gpp_gppl1_ts_cld_vp"

setwd(paste0(fire_dir,"/Rscripts"))
source("utils.R")

# library(sf)
# shp = st_read("D:/Water LUE Project/data/continent_shapefile/continent shapefile/continent.shp")

firenc_new = NcCreateOneShot(paste0(fire_dir,"/",output_dir, "/", model_dir, "/fire.vp+1pc.2002-1-1-2015-12-31.nc"), var_name = "fire")
firenc_orig = NcCreateOneShot(paste0(fire_dir,"/",output_dir, "/", model_dir, "/fire.2002-1-1-2015-12-31.nc"), var_name = "fire")

firenc_orig_yearly = firenc_orig
firenc_orig_yearly$data = apply(firenc_orig_yearly$data, MARGIN = c(1,2), FUN = function(x){mean(x,na.rm=T)})*12
cols = createPalette(c("black", "blue4", "blue", "skyblue", "cyan","mediumspringgreen","yellow","orange", "red","brown"),c(0,0.2,0.5,1,2,5,10,20,50,100)*1000, n = 1000) #gfed
plot.netcdf(firenc_orig_yearly, col=cols, zlim = c(0,1))
  
firenc_new_yearly = firenc_new
firenc_new_yearly$data = apply(firenc_new_yearly$data, MARGIN = c(1,2), FUN = function(x){mean(x,na.rm=T)})*12
cols = createPalette(c("black", "blue4", "blue", "skyblue", "cyan","mediumspringgreen","yellow","orange", "red","brown"),c(0,0.2,0.5,1,2,5,10,20,50,100)*1000, n = 1000) #gfed
plot.netcdf(firenc_new_yearly, col=cols, zlim = c(0,1))


lat_res = 111e3
lon_res = 111e3
cell_area = t(matrix(ncol = length(firenc_orig$lons), data = rep(lat_res*lon_res*cos(firenc_orig$lats*pi/180), length(firenc_orig$lons) ), byrow = F ))

mha_per_m2 = 0.0001/1e6
global_fire = sum(as.numeric(firenc_orig_yearly$data*cell_area), na.rm=T)*mha_per_m2 
global_fire_new = sum(as.numeric(firenc_new_yearly$data*cell_area), na.rm=T)*mha_per_m2 


lim = 0.25 #max(abs(tseffect$data), na.rm=T)
cols = createPalette(c("blue4", "skyblue", "cyan", "#f0f0f0", "yellow", "orange", "red"),c(-lim, -lim/2, -0.1*lim, 0, 0.1*lim, lim/2, lim)*1000, n = 1000) #gfed
diff_fire_yearly = firenc_new_yearly
diff_fire_yearly$data = firenc_new_yearly$data - firenc_orig_yearly$data
plot.netcdf(dat = diff_fire_yearly, zlim = c(-lim,lim), col = cols, ilev = 1, itime = 4,preserve_layout = T)

  

lim = 100 #max(abs(tseffect$data), na.rm=T)
cols = createPalette(c("blue4", "skyblue", "cyan", "#f0f0f0", "yellow", "orange", "red"),c(-lim, -lim/2, -0.1*lim, 0, 0.1*lim, lim/2, lim)*1000, n = 1000) #gfed
diff_fire_yearly = firenc_new_yearly
diff_fire_yearly$data = as.integer(firenc_orig_yearly$data > 0.01)*(firenc_new_yearly$data - firenc_orig_yearly$data)/(firenc_orig_yearly$data)*100
plot.netcdf(dat = diff_fire_yearly, zlim = c(-lim,lim), col = cols, ilev = 1, itime = NA,preserve_layout = T)



diff_fire = firenc_orig
diff_fire$data = (firenc_new$data - firenc_orig$data)
diff_fire_ymonmean = diff_fire
diff_fire_ymonmean$data = diff_fire_ymonmean$data[,,1:12]
diff_fire_ymonmean$month = 1:12

lim = 0.25 #max(abs(tseffect$data), na.rm=T)
# lim = 50
cols = createPalette(c("blue4", "skyblue", "cyan", "#f0f0f0", "yellow", "orange", "red"),c(-lim, -lim/2, -0.05*lim, 0, 0.1*lim, 0.2*lim, lim)*1000, n = 1000) #gfed
# cols = createPalette(c("black","blue4", "skyblue", "cyan", "#f0f0f0", "yellow", "orange", "red", "brown4"),c(-20*lim, -lim, -lim/2, -0.1*lim, 0, 0.1*lim, lim/2, lim, 20*lim)*1000, n = 1000) #gfed
for (i in c(3,5,8,12)){
  diff_month = apply(X=diff_fire$data[,,diff_fire$month == i], MARGIN = c(1,2), FUN = mean)  
  diff_month_pc = diff_month / (0.01 + apply(X=firenc_orig$data[,,firenc_orig$month == i], MARGIN = c(1,2), FUN = mean)) *100
  diff_fire_ymonmean$data[,,i] = diff_month
  plot.netcdf(dat = diff_fire_ymonmean, zlim = c(-lim,lim), col = cols, itime = i, preserve_layout = T)
  title(main=i)
}


# 
# lim = 20 #max(abs(tseffect$data), na.rm=T)
# cols = createPalette(c("black","blue4", "skyblue", "cyan", "#f0f0f0", "yellow", "orange", "red", "brown4"),c(-10*lim, -lim, -lim/2, -0.1*lim, 0, 0.1*lim, lim/2, lim, 10*lim)*1000, n = 1000) #gfed
# layout(matrix(c(1,1,
#                 1,1,
#                 1,1,
#                 2,2), ncol=2,byrow = T))  # vertical
# par(mar=c(4,5,4,1), oma=c(1,1,1,1), cex.lab=1.5, cex.axis=1.5)
# 
# plot.netcdf(dat = diff_fire, zlim = c(-lim,lim), col = cols, ilev = 1, itime = 4,preserve_layout = T)
