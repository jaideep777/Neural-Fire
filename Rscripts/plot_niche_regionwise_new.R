var1nc = cropnc
var2nc = gppl1nc

minx = min(var1nc$data[var1nc$data<1e18], na.rm=T)
maxx = max(var1nc$data[var1nc$data<1e18], na.rm=T)
miny = min(var2nc$data[var1nc$data<1e18], na.rm=T)
maxy = max(var2nc$data[var1nc$data<1e18], na.rm=T)

png(filename = paste0(fire_dir, "/plot_niche_varwise_1/", var1nc$name, "_", var2nc$name, "_", start_date, "-", end_date, ".png"), width = 650*3, height = 700*3, res=300)

par(mfrow = c(4,3), mar=c(1,1,1,1), oma=c(4,4,1,1), cex.lab=1.5, cex.axis=1.5)
regions_now = c("GLOBE", "BONA", "BOAS", "TCAM", "CEAS", "EUME", "SEAS", "EQAS", "NHAF", "SA", "AUS", "SHAF")
for (iregion in 1:length(regions_now)){
  region = regions_now[iregion]
  cat(region, "\n")
  xaxt1 = "n"
  yaxt1 = "n"
  if ((iregion-1) %% 3 == 0) yaxt1="s"
  if (iregion-1 >= 9) xaxt1="s"
  plot.niche.occurance_xy(var1nc, var2nc, fire_obs, region_name = region, xlab="", ylab="", xaxt=xaxt1, yaxt=yaxt1, minx=minx, maxx=maxx, miny=miny, maxy=maxy, main=region)
  # plot.niche_xy(var1nc, var2nc, fire_obs, region_name = region, xlab="", ylab="", xaxt=xaxt1, yaxt=yaxt1, minx=minx, maxx=maxx, miny=miny, maxy=maxy, main=region)
}
mtext(text=varnames[[var1nc$name]], side = 1, outer = T, line=2)
mtext(text=varnames[[var2nc$name]], side = 2, outer = T, line=2)

dev.off()