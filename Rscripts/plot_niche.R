# rm(list=ls())

# library(chron)
# rm(list = ls())
#### PREDICTED FIRES - CALIBRATION ####

fire_dir = "/home/jaideep/codes/PureNN_fire"
output_dir = "output_globe_vGCB"
model_dir = "BOAS_mod744.7_ltn_gppm1_pr_ts_vp"
data_dir= "/media/jaideep/San/Data"

# for (model_dir in list.files(path = paste0(fire_dir,"/",output_dir), no.. = T, pattern = "mod")){
cat(model_dir, "\n")

# region_name = strsplit(model_dir, split = "_")[[1]][1]
regions_list = list(BONA = c(1), TCAM = c(2,3), TENA=c(2), CEAM=c(3), SA = c(4,5), NHAF=c(8), SHAF = c(9), EUME = c(6,7), AF=c(8,9), NHAF = c(8), SHAF = c(9), BOAS = c(10), CEAS= c(11), SEAS=c(12), EQAS = c(13), AUS = c(14), GLOBE = 1:14)
# reg = get(region_name, regions_list)

fire_obs_file = paste0(data_dir, "/Fire_BA_GFED4.1s/nc/GFED_4.1s_1deg.1997-2016.nc")  # Need absolute path here
fire_pred_file = "fire.2002-1-1-2015-12-31.nc"
start_date  = "2002-1-1"
end_date    = "2015-12-31"

varnames = list(cld="Cloud cover", gpp="GPP", gppl1 = "Canopy GPP", ts="Temperature", pr="Precipitation", pop="Pop. density", vp="Vapour pres", gppm1="Litter GPP (N)", gppm1s="Litter GPP (S)", crop="Crop frac", ltn="Lightning", ft_12="Crop frac")

# Get model_dir from command line
args = commandArgs(trailingOnly = T)
l = (strsplit(x = args, split = "="))
opt = unlist(lapply(l, "[[", 1))
spec = unlist(lapply(l, "[[", 2))

findOpt = function(o){
  ! (is.null(spec[opt==o]) | length(which(opt == o)) == 0)
}

if (findOpt("model_dir")) model_dir = spec[opt=="model_dir"]
if (findOpt("output_dir")) output_dir = spec[opt=="output_dir"]

source(paste0(fire_dir,"/Rscripts/utils.R"))


# ts_slice = apply(X = tsnc$data, MARGIN = c(1,2), FUN = mean)
# vp_slice = apply(X = vpnc$data, MARGIN = c(1,2), FUN = mean)
# cld_slice = apply(X = cldnc$data, MARGIN = c(1,2), FUN = mean)
# pop_slice = apply(X = popnc$data, MARGIN = c(1,2), FUN = mean)
# gppl1_slice = apply(X = gppl1nc$data, MARGIN = c(1,2), FUN = mean)
# gppm1_slice = apply(X = gppm1nc$data, MARGIN = c(1,2), FUN = mean)
# gfedl1_slice = apply(X = gfedl1nc$data, MARGIN = c(1,2), FUN = mean)


# fire_pred_filename = paste0(fire_dir,"/",output_dir, "/", model_dir, "/", fire_pred_file)
# fire_pred = NcCreateOneShot(filename = fire_pred_filename, var_name = "fire")
# # fire_pred$time = fire_pred$time - 15
# # fire_pred$time = as.Date("2003-1-15") + 365.2524/12*(0:155)
# fire_pred = NcClipTime(fire_pred, start_date, end_date)
# # fire_pred$data = fire_pred$data - 0.000
# # fire_pred$data[fire_pred$data < 0.00] = 0
# 
# # fire_pred = gfedl1nc
# 
# glimits = c(fire_pred$lons[1],
#             fire_pred$lons[length(fire_pred$lons)],
#             fire_pred$lats[1],
#             fire_pred$lats[length(fire_pred$lats)])  # get limits from predicted data
# 
# slices_per_yr_pred = 365.2524/as.numeric(mean(diff(fire_pred$time[-length(fire_pred$time)])))
# 
# lat_res = mean(diff(fire_pred$lats))*111e3
# lon_res = mean(diff(fire_pred$lons))*111e3
# cell_area = t(matrix(ncol = length(fire_pred$lons), data = rep(lat_res*lon_res*cos(fire_pred$lats*pi/180), length(fire_pred$lons) ), byrow = F ))
# 
# fire_obs = NcCreateOneShot(filename = fire_obs_file, var_name = "ba", glimits = glimits)
# fire_obs$time = as.Date("1997-1-15") + 365.2524/12*(0:239)
# fire_obs = NcClipTime(fire_obs,  start_date, end_date)
# fire_obs$data[is.na(fire_pred$data)] = NA
# 
# slices_per_yr_obs = 365.2524/as.numeric(mean(diff(fire_obs$time[-length(fire_obs$time)])))
# 
# slice_pred = apply(X = fire_pred$data, FUN = function(x){mean(x, na.rm=T)}, MARGIN = c(1,2))*slices_per_yr_pred
# slice_pred = slice_pred*cell_area
# 
# slice_obs = apply(X = fire_obs$data, FUN = function(x){mean(x, na.rm=T)}, MARGIN = c(1,2))*slices_per_yr_obs
# slice_obs = slice_obs*cell_area
# 


# xnc=gfedl1nc
# ync=vpnc
# znc = fire_obs
# zpnc = fire_pred



plot.niche_xy = function(xnc, ync, znc, zpnc=NULL, region_name, ft = 1:12, bins = c(50,50), tx = function(x){x}, ty=function(y){y}, minx=NA, maxx=NA, miny=NA, maxy=NA, ...){
  region_indices = get(region_name, regions_list)
  # r1 = region_indices[1]
  # r2 = region_indices[length(region_indices)]
  reg = regions
  reg[!(reg %in% region_indices)] = NA
  
  xnc$data[is.na(reg)] = NA
  ync$data[is.na(reg)] = NA
  znc$data[is.na(reg)] = NA
  if (!is.null(zpnc)) zpnc$data[is.na(reg)] = NA
  
  xnc$data[xnc$data>1e18]=NA
  ync$data[ync$data>1e18]=NA
  
  x = as.numeric(tx(xnc$data))
  y = as.numeric(ty(ync$data))
  z = as.numeric(znc$data)
  if (!is.null(zpnc)) zp= as.numeric(zpnc$data)
  
  if (is.na(minx)) minx = min(x[!is.na(z)], na.rm=T)
  if (is.na(maxx)) maxx = max(x[!is.na(z)], na.rm=T)
  if (is.na(miny)) miny = min(y[!is.na(z)], na.rm=T)
  if (is.na(maxy)) maxy = max(y[!is.na(z)], na.rm=T)
  
  breaks_x = seq(minx, maxx, length.out=bins[1])
  breaks_y = seq(miny, maxy, length.out=bins[2])
  
  xclass = cut(x, breaks = breaks_x)
  yclass = cut(y, breaks = breaks_y)
  mat      = tapply(X = z,   INDEX = list(xclass, yclass), FUN = function(x)(mean(x,na.rm=T)))
  if (!is.null(zpnc)) mat_pred = tapply(X = zp,  INDEX = list(xclass, yclass), FUN = function(x)(mean(x,na.rm=T)))
  
  # t_x  = apply(X = xnc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_y  = apply(X = ync$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_z  = apply(X = znc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_zp  = apply(X = zpnc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # 
  # t_x_yr = tapply(X = t_x, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
  # t_y_yr = tapply(X = t_y, INDEX = strftime(ync$time, "%Y"), FUN = mean)
  # t_z_yr = tapply(X = t_z, INDEX = strftime(znc$time, "%Y"), FUN = mean)
  # t_zp_yr = tapply(X = t_zp, INDEX = strftime(zpnc$time, "%Y"), FUN = mean)
  # 
  # mod_x = lm(t_x_yr~seq(1,length(t_x_yr)))
  # mod_y = lm(t_y_yr~seq(1,length(t_y_yr)))
  
  # png(filename = paste0(fire_dir, "/niche_plots_popcrop/", region_name, "_", xnc$name, "_", ync$name, "_", start_date, "-", end_date, ".png"), width = 500*3, height = 500*3, res=300)
  # 
  # layout(matrix(c(1,2,3,4), nrow=2, byrow=T))
  # par(mar=c(4,5,1,1), oma=c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
  
  cols_niche = createPalette(c("black", "blue4", "blue", "mediumspringgreen", "yellow", "red", "brown", "white"), c(0,0.002, 0.005, 0.02, 0.05, 0.2, 0.5, 1), 10000)
  image(mat, zlim=c(0,0.25), x = breaks_x, y=breaks_y, col=cols_niche, ...) #xlab=as.character(varnames[xnc$name]), ylab=as.character(varnames[ync$name]), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )

  # image(mat_pred, zlim=c(0,0.25), x = breaks_x, y=breaks_y, col=cols_niche, xlab=as.character(varnames[xnc$name]), ylab=as.character(varnames[ync$name]), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )

  # image(mat, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(t_x_yr), max(t_x_yr)), ylim=c(min(t_y_yr), max(t_y_yr)), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  # 
  # image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(t_x_yr), max(t_x_yr)), ylim=c(min(t_y_yr), max(t_y_yr)), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )

  #mtext(region_name, side = 3, outer = TRUE, cex=1.3)
  # dev.off()  
#   image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, ...)
# 
  # lons = c(27.5, -0.5  )
  # lats = c(7.5,  +9.5  )
#   for (i in 1:length(lons)){
#   
#   tx1 = xnc$data[which(xnc$lons == lons[i]), which(xnc$lats == lats[i]),]
#   ty1 = ync$data[which(ync$lons == lons[i]), which(ync$lats == lats[i]),]
#   tz1 = znc$data[which(znc$lons == lons[i]), which(znc$lats == lats[i]),]
#   tzp1 = zpnc$data[which(zpnc$lons == lons[i]), which(zpnc$lats == lats[i]),]
#   #   
# #   segments(tx1[-length(tx1)], ty1[-length(ty1)],
# #            tx1[-1], ty1[-1],
# #            col=colorRampPalette(c("black", "white"))(12*12), lwd=1.5 )
# # 
# #   
#   tx1_yr = tapply(X = tx1, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
#   ty1_yr = tapply(X = ty1, INDEX = strftime(ync$time, "%Y"), FUN = mean)
#   tz1_yr = tapply(X = tz1, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
#   tzp1_yr = tapply(X = tzp1, INDEX = strftime(ync$time, "%Y"), FUN = mean)
  # 
# #  image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(tx1_yr), max(tx1_yr)), ylim=c(min(ty1_yr), max(ty1_yr)), ...)
#   segments(tx1_yr[-length(tx1_yr)], ty1_yr[-length(ty1_yr)],
#            tx1_yr[-1], ty1_yr[-1],
#            col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
# 
#   
#   }  
}

# plot.niche_xy(prnc, cldnc, fire_obs, fire_pred, region_name = "NHAF", xlab="Precipitation", ylab="Cloud cover")
# plot.niche_xy(prnc, popnc, fire_obs, fire_pred, region_name = "NHAF", xlab="Precipitation", ylab="Population density")
# plot.niche_xy(tsnc, prnc, fire_obs, fire_pred, region_name = "NHAF", xlab="Temperature", ylab="Precipitation")
# plot.niche_xy(tsnc, gppm1nc, fire_obs, fire_pred, region_name = "NHAF", xlab="Temperature", ylab="Yearly GPP")
# plot.niche_xy(prnc, gppm1nc, fire_obs, fire_pred, region_name = "NHAF", xlab="Precipitation", ylab="Yearly GPP")
# plot.niche_xy(popnc, gppm1nc, fire_obs, fire_pred, region_name = "NHAF", xlab="Population density", ylab="Yearly GPP")
# 
# plot.niche_xy(tsnc, cldnc, fire_obs, fire_pred, region_name = "SHAF", xlab="Temperature", ylab="Cloud cover")
# plot.niche_xy(gppl1nc, cldnc, fire_obs, fire_pred, region_name = "SHAF", xlab="Cumm. GPP", ylab="Cloud Cover")
# plot.niche_xy(gppl1nc, prnc, fire_obs, fire_pred, region_name = "SHAF", xlab="Cumm. GPP", ylab="Precipitation")
# 

plot.niche_x = function(xnc, znc, zpnc=NULL, region_name, ft = 1:12, bins = c(50), tx = function(x){x}, minx=NA, maxx=NA, ...){
  region_indices = get(region_name, regions_list)
  # r1 = region_indices[1]
  # r2 = region_indices[length(region_indices)]
  reg = regions
  reg[!(reg %in% region_indices)] = NA
  
  xnc$data[is.na(reg)] = NA
  znc$data[is.na(reg)] = NA
  if (!is.null(zpnc)) zpnc$data[is.na(reg)] = NA
  
  xnc$data[xnc$data>1e18]=NA

  x = as.numeric(tx(xnc$data))
  z = as.numeric(znc$data)
  if (!is.null(zpnc)) zp= as.numeric(zpnc$data)
  
  if (is.na(minx)) minx = min(x[!is.na(z)], na.rm=T)
  if (is.na(maxx)) maxx = max(x[!is.na(z)], na.rm=T)

  breaks_x = seq(minx, maxx, length.out=bins[1])

  xclass = cut(x, breaks = breaks_x)
  mat      = tapply(X = z,   INDEX = list(xclass), FUN = function(x)(mean(x,na.rm=T)))
  if (!is.null(zpnc)) mat_pred = tapply(X = zp,  INDEX = list(xclass), FUN = function(x)(mean(x,na.rm=T)))
  
  # t_x  = apply(X = xnc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_y  = apply(X = ync$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_z  = apply(X = znc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_zp  = apply(X = zpnc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # 
  # t_x_yr = tapply(X = t_x, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
  # t_y_yr = tapply(X = t_y, INDEX = strftime(ync$time, "%Y"), FUN = mean)
  # t_z_yr = tapply(X = t_z, INDEX = strftime(znc$time, "%Y"), FUN = mean)
  # t_zp_yr = tapply(X = t_zp, INDEX = strftime(zpnc$time, "%Y"), FUN = mean)
  # 
  # mod_x = lm(t_x_yr~seq(1,length(t_x_yr)))
  # mod_y = lm(t_y_yr~seq(1,length(t_y_yr)))
  
  # png(filename = paste0(fire_dir, "/niche_plots_popcrop/", region_name, "_", xnc$name, "_", ync$name, "_", start_date, "-", end_date, ".png"), width = 500*3, height = 500*3, res=300)
  # 
  # layout(matrix(c(1,2,3,4), nrow=2, byrow=T))
  # par(mar=c(4,5,1,1), oma=c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
  
  # cols_niche = createPalette(c("black", "blue4", "blue", "mediumspringgreen", "yellow", "red", "brown", "white"), c(0,0.002, 0.005, 0.02, 0.05, 0.2, 0.5, 1), 10000)
  plot(y = mat, x = mids(breaks_x), col="red", type="l", lwd=3, xlab=as.character(varnames[xnc$name]), ylab="Mean BA", ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  
  # image(mat_pred, zlim=c(0,0.25), x = breaks_x, y=breaks_y, col=cols_niche, xlab=as.character(varnames[xnc$name]), ylab=as.character(varnames[ync$name]), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  
  # image(mat, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(t_x_yr), max(t_x_yr)), ylim=c(min(t_y_yr), max(t_y_yr)), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  # 
  # image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(t_x_yr), max(t_x_yr)), ylim=c(min(t_y_yr), max(t_y_yr)), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  
  #mtext(region_name, side = 3, outer = TRUE, cex=1.3)
  # dev.off()  
  #   image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, ...)
  # 
  # lons = c(27.5, -0.5  )
  # lats = c(7.5,  +9.5  )
  #   for (i in 1:length(lons)){
  #   
  #   tx1 = xnc$data[which(xnc$lons == lons[i]), which(xnc$lats == lats[i]),]
  #   ty1 = ync$data[which(ync$lons == lons[i]), which(ync$lats == lats[i]),]
  #   tz1 = znc$data[which(znc$lons == lons[i]), which(znc$lats == lats[i]),]
  #   tzp1 = zpnc$data[which(zpnc$lons == lons[i]), which(zpnc$lats == lats[i]),]
  #   #   
  # #   segments(tx1[-length(tx1)], ty1[-length(ty1)],
  # #            tx1[-1], ty1[-1],
  # #            col=colorRampPalette(c("black", "white"))(12*12), lwd=1.5 )
  # # 
  # #   
  #   tx1_yr = tapply(X = tx1, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
  #   ty1_yr = tapply(X = ty1, INDEX = strftime(ync$time, "%Y"), FUN = mean)
  #   tz1_yr = tapply(X = tz1, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
  #   tzp1_yr = tapply(X = tzp1, INDEX = strftime(ync$time, "%Y"), FUN = mean)
  # 
  # #  image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(tx1_yr), max(tx1_yr)), ylim=c(min(ty1_yr), max(ty1_yr)), ...)
  #   segments(tx1_yr[-length(tx1_yr)], ty1_yr[-length(ty1_yr)],
  #            tx1_yr[-1], ty1_yr[-1],
  #            col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  # 
  #   
  #   }  
  
}

plot.niche.occurance_xy = function(xnc, ync, znc, zpnc=NULL, region_name, ft = 1:12, bins = c(50,50), tx = function(x){x}, ty=function(y){y}, minx=NA, maxx=NA, miny=NA, maxy=NA, ...){
  region_indices = get(region_name, regions_list)
  # r1 = region_indices[1]
  # r2 = region_indices[length(region_indices)]
  reg = regions
  reg[!(reg %in% region_indices)] = NA
  
  xnc$data[is.na(reg)] = NA
  ync$data[is.na(reg)] = NA
  znc$data[is.na(reg)] = NA
  # if (!is.null(zpnc)) zpnc$data[is.na(reg)] = NA
  
  xnc$data[xnc$data>1e18]=NA
  ync$data[ync$data>1e18]=NA
  
  x = as.numeric(tx(xnc$data))
  y = as.numeric(ty(ync$data))
  z = as.numeric(znc$data)
  # if (!is.null(zpnc)) zp= as.numeric(zpnc$data)
  
  if (is.na(minx)) minx = min(x[!is.na(z)], na.rm=T)
  if (is.na(maxx)) maxx = max(x[!is.na(z)], na.rm=T)
  if (is.na(miny)) miny = min(y[!is.na(z)], na.rm=T)
  if (is.na(maxy)) maxy = max(y[!is.na(z)], na.rm=T)
  
  breaks_x = seq(minx, maxx, length.out=bins[1])
  breaks_y = seq(miny, maxy, length.out=bins[2])
  
  xclass = cut(x, breaks = breaks_x)
  yclass = cut(y, breaks = breaks_y)
  mat      = tapply(X = z,   INDEX = list(xclass, yclass), FUN = function(x)(length(which(!is.na(x) ))))
  mat = mat/length(which(!is.na(z)))
  # mat = mat/3e5
  # cat(length(which(!is.na(z) )))
  # mat_pred = tapply(X = zp,  INDEX = list(xclass, yclass), FUN = function(x)(mean(x,na.rm=T)))
  
  # t_x  = apply(X = xnc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_y  = apply(X = ync$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_z  = apply(X = znc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # t_zp  = apply(X = zpnc$data, MARGIN = 3, FUN = function(x) {mean(x, na.rm=T)})
  # 
  # t_x_yr = tapply(X = t_x, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
  # t_y_yr = tapply(X = t_y, INDEX = strftime(ync$time, "%Y"), FUN = mean)
  # t_z_yr = tapply(X = t_z, INDEX = strftime(znc$time, "%Y"), FUN = mean)
  # t_zp_yr = tapply(X = t_zp, INDEX = strftime(zpnc$time, "%Y"), FUN = mean)
  # 
  # mod_x = lm(t_x_yr~seq(1,length(t_x_yr)))
  # mod_y = lm(t_y_yr~seq(1,length(t_y_yr)))
  
  # png(filename = paste0(fire_dir, "/niche_plots_popcrop/", region_name, "_", xnc$name, "_", ync$name, "_", start_date, "-", end_date, ".png"), width = 500*3, height = 500*3, res=300)
  # 
  # layout(matrix(c(1,2,3,4), nrow=2, byrow=T))
  # par(mar=c(4,5,1,1), oma=c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
  
  cols_niche = createPalette(c("black", "cyan", "white"), c(0,200,1000), 10000)
  cols_niche = createPalette(c("white", "grey85", "blue", "cyan", "white"), c(0,10,200,700,2000), 10000)
  image(mat, zlim=c(0,0.01), x = breaks_x, y=breaks_y, col=cols_niche, ...) #xlab=as.character(varnames[xnc$name]), ylab=as.character(varnames[ync$name]), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  
  # image(mat_pred, zlim=c(0,0.25), x = breaks_x, y=breaks_y, col=cols_niche, xlab=as.character(varnames[xnc$name]), ylab=as.character(varnames[ync$name]), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  
  # image(mat, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(t_x_yr), max(t_x_yr)), ylim=c(min(t_y_yr), max(t_y_yr)), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  # 
  # image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(t_x_yr), max(t_x_yr)), ylim=c(min(t_y_yr), max(t_y_yr)), ...)
  # segments(fitted(mod_x)[-length(fitted(mod_x))], fitted(mod_y)[-length(fitted(mod_y))],
  #          fitted(mod_x)[-1], fitted(mod_y)[-1],
  #          col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  
  #mtext(region_name, side = 3, outer = TRUE, cex=1.3)
  # dev.off()  
  #   image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, ...)
  # 
  # lons = c(27.5, -0.5  )
  # lats = c(7.5,  +9.5  )
  #   for (i in 1:length(lons)){
  #   
  #   tx1 = xnc$data[which(xnc$lons == lons[i]), which(xnc$lats == lats[i]),]
  #   ty1 = ync$data[which(ync$lons == lons[i]), which(ync$lats == lats[i]),]
  #   tz1 = znc$data[which(znc$lons == lons[i]), which(znc$lats == lats[i]),]
  #   tzp1 = zpnc$data[which(zpnc$lons == lons[i]), which(zpnc$lats == lats[i]),]
  #   #   
  # #   segments(tx1[-length(tx1)], ty1[-length(ty1)],
  # #            tx1[-1], ty1[-1],
  # #            col=colorRampPalette(c("black", "white"))(12*12), lwd=1.5 )
  # # 
  # #   
  #   tx1_yr = tapply(X = tx1, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
  #   ty1_yr = tapply(X = ty1, INDEX = strftime(ync$time, "%Y"), FUN = mean)
  #   tz1_yr = tapply(X = tz1, INDEX = strftime(xnc$time, "%Y"), FUN = mean)
  #   tzp1_yr = tapply(X = tzp1, INDEX = strftime(ync$time, "%Y"), FUN = mean)
  # 
  # #  image(mat_pred, zlim=c(0,0.2), x = breaks_x, y=breaks_y, col=cols_niche, xlim=c(min(tx1_yr), max(tx1_yr)), ylim=c(min(ty1_yr), max(ty1_yr)), ...)
  #   segments(tx1_yr[-length(tx1_yr)], ty1_yr[-length(ty1_yr)],
  #            tx1_yr[-1], ty1_yr[-1],
  #            col=colorRampPalette(c("black", "magenta"))(12), lwd=3 )
  # 
  #   
  #   }  
}

