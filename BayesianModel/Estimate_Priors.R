### Estimate Prior cutoffs from plot dataset

library(glue)
library(terra)
library(ggplot2)
library(sf)
library(pool)

sppDb <- dbPool(
  drv = RPostgres::Postgres(),
  dbname = "spp_feas",
  host = Sys.getenv("BCGOV_HOST"),
  port = 5432, 
  user = Sys.getenv("BCGOV_USR"),
  password = Sys.getenv("BCGOV_PWD")
)


QRY <- glue_sql("select plotnum, spp, ST_X(ST_Transform(geometry, 4326)) AS lon, 
              ST_Y(ST_Transform(geometry, 4326)) AS lat 
              from plotdata", .con = sppDb)
dat <- dbGetQuery(sppDb, QRY)
setDT(dat)
dat <- unique(dat)
dat[,id := seq_along(plotnum)]
setkey(dat, id)
dem <- rast("../Common_Files/NA_DEM_250.tif")
temp <- vect(dat)
crs(temp) <- "epsg:4326"

bgc <- vect("../Common_Files/WNA_BGCv13/WNA_v13_4326.gpkg")
bc_bgc <- bgc[bgc$Region == "BC",]
bgc_association <- terra::intersect(temp, bgc)
bgc_ass <- as.data.frame(bgc_association) |> 
  as.data.table()
fwrite(bgc_ass, "Plot_data_bgc_all.csv")

bgc_coast <- bgc_ass[grep("CWH.*|CDF.*|MH.*", BGC),]

elev <- extract(dem, temp)
dat[,elev := elev$northamerica_elevation_cec_2023]

spp_use <- c("Fd", "Pl", "Lw", "Py", "At", "Cw", "Hw", "Bl", "Ba", "Ss", "Sx", "Sb", "Hm", "Yc", "Pa", "Pw", "Ra", "Act")
dat_spp <- dat[spp %in% spp_use,]

clim <- downscale(dat_spp, return_refperiod = TRUE, vars = c("DD5_an","CMD_an","VPD_sm", "CMI_an"))
setkey(clim, id)

tree_clim <- clim[dat_spp]
tree_clim[bgc_ass, BGC := i.BGC, on = "id"]
tree_clim[,coast_int := "Interior"]
tree_clim[grep("CWH.*|CDF.*|MH.*", BGC), coast_int := "Coast"]

coast_spp <- tree_clim[coast_int == "Coast" & spp %in% c("Yc","Ss","Ba"),]
int_spp <- tree_clim[coast_int == "Interior" & !(spp %in% c("Yc","Ss","Ba")),]
comb_spp <- rbind(coast_spp, int_spp)

vpd_cutoffs <- tree_clim[,.(VPD75 = quantile(VPD_sm, 0.75, na.rm = T), 
                            VPD50 = quantile(VPD_sm, 0.5, na.rm = T)), by = spp]

dd5_cutoffs <- tree_clim[,.(DD502 = quantile(DD5_an, 0.02, na.rm = T), 
                            DD550 = quantile(DD5_an, 0.5, na.rm = T)), by = spp]

setorder(vpd_cutoffs, VPD75)
tree_clim[,spp := factor(spp, levels = vpd_cutoffs$spp)]

setorder(dd5_cutoffs, DD550)
tree_clim[,spp := factor(spp, levels = dd5_cutoffs$spp)]

ggplot(tree_clim, aes(x = spp, y = VPD_sm)) +
  geom_violin(quantile.linetype = 1, quantiles = c(0.2,0.5,0.98)) +
  geom_point(data = vpd_cutoffs, aes(x = spp, y = VPD75))

ggplot(tree_clim, aes(x = spp, y = CMI_an)) +
  geom_violin(quantile.linetype = 1, quantiles = c(0.2,0.5,0.98))

ggsave("DD5_Cutoffs_byspp.png", width = 12, height = 8)

ggplot(tree_clim, aes(x = DD5_an)) +
  geom_histogram() +
  facet_wrap(~spp, scales = "free")

ggplot(tree_clim, aes(x = spp, y = VPD_sm, col = coast_int)) +
  geom_violin(quantile.linetype = 1, quantiles = c(0.1,0.5,0.99))


p[vpd_cutoffs, V50_median := i.VPD75, on = c(species = "spp")]

tree_clim <- na.omit(tree_clim)
tree_stats <- tree_clim[,.(VPD_Width = quantile(VPD_sm, 0.99) - quantile(VPD_sm, 0.95)), by = spp]
p[tree_stats, w_V_median := i.VPD_Width, on = c(species = "spp")]
