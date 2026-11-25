### Bayesian Sprint Model\
## Kiri Daust

library(data.table)
library(climr)
library(cmdstanr)
library(DBI)
library(dplyr)
library(pool)

sibec.anal <- readRDS("data/SIBEC_Analysis_data_2026.rds")

sibec.tst <- sibec.anal
# Step 1: Calculate non-zero SI50 count per TreeSpp
valid_species <- sibec.tst %>%
  filter(SI50 > 0.5) %>%
  count(TreeSpp) %>%
  filter(n >= 20) %>%
  pull(TreeSpp)

# Step 2: Filter main dataset to retain only those species
sibec.tst <- sibec.tst %>%
  filter(TreeSpp %in% valid_species)

final_dat <- as.data.table(sibec.tst)
final_dat <- final_dat[,.(PlotNumber, BGC, SI50, TreeSpp, rSMR, SNR, DD5_an, CMD_an, VPD_sm, PAS_an, DD18_an, Eref_an)]

##only keep 0.5 SI for BG units and only for 4C
tmp <- final_dat[grepl("BG.*",BGC) & rSMR %in% c(4) & SNR %in% c("C"),]

##just use circumzonal data
fit_dat <- final_dat[SI50 > 0.5 & rSMR %in% c(3,4,5,6) & SNR %in% c("B", "C", "D", "E"),]
fit_dat <- rbind(fit_dat, tmp)

# vars <- c("DD5_an", "VPD_sm", "CMD_an")
# spp <- "Sx"
# 
# pairs(
#   as.data.frame(
#     fit_dat[TreeSpp == spp, ..vars]
#   ),
#   main = paste(spp, "— observed climate relationships"),
#   pch = 16,
#   cex = 0.4,
#   col = adjustcolor("steelblue", alpha.f = 0.25)
# )

###piecewise with DD5:CMD interaction
source("BayesianModel/site_index_piecewise_v2.R")
p <- data.table::fread("BayesianModel/piecewise_priors.csv")
setnames(fit_dat, c("plotnum", "bgc", "site_index", "species", "moisture", "nutrient","DD5","CMD","VPD", "PAS", "DD18", "Eref"))

##make sure we have the same species in the data and in the priors
fit_dat <- fit_dat[species %in% p$species,]
p <- p[species %in% unique(fit_dat$species),]

##prep the data
prepared_piecewise <- prepare_si_data_piecewise_v2(fit_dat,p, gamma_CMD_mean=0, gamma_CMD_mean_sd=.25, gamma_CMD_tau_sd=.20)
##fit!
fit_piecewise <- fit_si_piecewise_v2(prepared_piecewise)

plot_si_curves_piecewise_v2(fit_piecewise,variable="DD5",species=c("Pl","Fd","Sx", "Py"),VPD=0.6,CMD=200,
                         moisture=4,nutrient="C",observations=TRUE, xlim = c(0,2500))
plot_si_curves_piecewise_v2(fit_piecewise,variable="CMD",species=c("Pl","Fd","Sx", "Py"), DD5 = 1200,VPD=0.2,
                         moisture=4,nutrient="C")
for (temperature in c(500,1500,2500)) {
  plot_si_curves_piecewise_v2(fit_piecewise,variable='CMD',
                              species=c('Fd','Lw'),DD5=temperature,VPD=.6,moisture=4,nutrient='C')
}


###piecewise fit
source("BayesianModel/fitting_functions.R")
p <- data.table::fread("BayesianModel/piecewise_priors.csv")
setnames(fit_dat, c("plotnum", "bgc", "site_index", "species", "moisture", "nutrient","DD5","CMD","VPD", "PAS", "DD18", "Eref"))

##plot some priors
check_si_piecewise_plateau_priors(p, species = c("Fd","Sx"))

##make sure we have the same species in the data and in the priors
fit_dat <- fit_dat[species %in% p$species,]
p <- p[species %in% unique(fit_dat$species),]

##prep the data
prepared_piecewise <- prepare_si_data_piecewise(fit_dat,p)
##fit!
fit_piecewise <- fit_si_piecewise(prepared_piecewise)

plot_si_curves_piecewise(fit_piecewise,variable="DD5",species=c("Pl","Fd","Sx", "Py"),VPD=0.6,CMD=200,
                         moisture=4,nutrient="C",observations=TRUE, xlim = c(0,2500))
plot_si_curves_piecewise(fit_piecewise,variable="VPD",species=c("Pl","Fd","Sx", "Py"), DD5 = 1200,CMD=200,
                         moisture=4,nutrient="C")
plot_si_edatopic_piecewise(fit_piecewise,"Sx",DD5=1300,VPD=.6,CMD=100)

#####################################################################################
##plot residuals vs climate variables

plot_si_residuals(
  fit_piecewise, fit_dat,
  variable = "DD5"
)
plot_si_residuals(
  fit_piecewise, dat_circzonal,
  variable = "VPD"
)
plot_si_residuals(
  fit_piecewise, dat_circzonal,
  variable = "CMD"
)
plot_si_residuals(
  fit_piecewise, fit_dat,
  variable = "DD18"
)
plot_si_residuals(
  fit_piecewise, dat_circzonal,
  variable = "Eref"
)


## check if species are lost in zones where they should not exist
bgc <- vect("../Common_Files/WNA_BGCv13/BGC_BC_v13_Jan21.gdb/")
bgc_tmp <- bgc[bgc$MAP_LABEL %in% c("IDFdk3","SBSdk", "BGxh1", "CWHvm1"),]
bgc_pts <- spatSample(bgc_tmp, 200, strata = "MAP_LABEL")
elev <- extract(dem, bgc_pts)

points <- project(bgc_pts, "epsg:4326") |> as.data.frame(geom = "XY") 
setDT(points)
points$elev <- elev$northamerica_elevation_cec_2023
points[,id := seq_along(elev)]
points <- points[,.(id, BGC = MAP_LABEL, x, y, elev)]
points <- na.omit(points)
setnames(points, c("id","BGC", "lon", "lat", "elev"))

clim <- downscale(points, return_refperiod = TRUE, vars = c("DD5_an","CMD_an", "VPD_sm"))
clim[points, BGC := i.BGC, on = "id"]

val_dat <- clim[,.(id, BGC, rSMR = 4, rSNR = "C", DD5_an, CMD_an, VPD_sm)]
setnames(val_dat, c("id", "BGC", "moisture", "nutrient", "DD5", "CMD", "VPD"))
val_dat <- unique(val_dat)

preds <- predict_si_validation_piecewise(fit_piecewise, val_dat)

ggplot(temp_clim, aes(x = BGC, y = VPD_sm)) +
  geom_violin()

preds <- preds[species %in% c("Ss","Pl","Fd","Py", "Cw"),]
ggplot(preds, aes(x = species, y = SI_median)) +
  geom_violin() +
  facet_wrap(~BGC)
