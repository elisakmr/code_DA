




lsm_latlon <- get(load(file = file.path(WD.tmp, paste0(year, "_mapbiomas_lsm",i,".RData"))))
lsm_utm <- get(load(file = file.path(WD.tmp, paste0(year, "_mapbiomasbis_lsm",i,".RData"))))

sum((lsm_latlon %>% filter(metric=='pland'))$value)
sum((lsm_utm %>% filter(metric=='pland'))$value)


head(lsm_latlon)
head(lsm_utm)
