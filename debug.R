




lsm_latlon <- get(load(file = file.path(WD.tmp, paste0(year, "_mapbiomas_lsm",i,".RData"))))
lsm_utm <- get(load(file = file.path(WD.tmp, paste0(year, "_mapbiomasbis_lsm",i,".RData"))))

sum((lsm_latlon %>% filter(metric=='pland'))$value)
sum((lsm_utm %>% filter(metric=='pland'))$value)


head(lsm_latlon)
head(lsm_utm)


list_lsm <- list()
list_buf <- list()

for (i in 1:200){
  
  WD.tmp = list.files(path=dir_data, pattern = paste0("Basin_",i,"_"), full.names = TRUE)
  
  list_buf[[i]] <- get(load(file=file.path(WD.tmp, paste0(year, "_mapbiomas_lsm_buf",tampon, "_",i,".RData")))) 
  

  
  
  
}

which(is.na(list_buf))
df_lsm <- do.call("rbind",list_lsm)
df_buf <- do.call("rbind",list_buf)
dim(df_buf)
dim(df_lsm)

which(list_buf=="nonascente")


ras3=mapbiomas[[3]]
ras30=mapbiomas[[30]]

plot(ras3)
plot(ras30)

lsm3<-calculate_lsm(ras3,level="class",metric = c("ca","pland","np"),directions =8)
lsm30 <- calculate_lsm(ras30,level="class",metric = c("ca","pland","np"),directions =8)

are3 <- lsm_l_area_mn(mapbiomas)
are30 <- lsm_l_area_mn(ras30)
lsm3<-calculate_lsm(mapbiomas,level="class",metric = c("ca","pland","np"),directions =8, progress = TRUE)
df=list_lsm[[i]]
d1 <- df %>% filter (layer==2)
d4 <- df %>% filter (layer==14)


###########reminder#########################
unlist(str_split(WD.tmp,'_', n=4))[4]
############################################


for (i in 1:2031){
  
  WD.tmp = list.files(path=dir_data, pattern = paste0("Basin_",i,"_"), full.names = TRUE)
  
  basin_id <- unlist(str_split(WD.tmp,'_', n=4))[4]
  
  lsm <- get(load(file=file.path(WD.tmp, paste0("mapbiomas_lsm",i,".RData")))) 
  lsm_buf100 <- get(load(file=file.path(WD.tmp, paste0("mapbiomas_lsm_buf100_",i,".RData")))) 
  lsm_buf500 <- get(load(file=file.path(WD.tmp, paste0("mapbiomas_lsm_buf500_",i,".RData")))) 
  hydrometric <- get(load(file=file.path(WD.tmp, paste0("hydrometric_",i,".RData")))) 
  mapbiomas <- rast(file.path(WD.tmp, paste0("mapbiomas_",i, ".tif")))
  rios <- st_read(file.path(WD.tmp,paste0("riossimple_",i,".shp")))
  nascente <- st_read(file.path(WD.tmp,paste0("nascente_",i,".shp")))
  reserv <- st_read(file.path(WD.tmp,paste0("reserv_",i,".shp")))
  
  st_write(rios,file.path(WD.tmp,paste0("riossimple_",basin_id,".shp")))
  st_write(nascente,file.path(WD.tmp,paste0("nascente_",basin_id,".shp")))
  st_write(reserv,file.path(WD.tmp,paste0("reserv_",basin_id,".shp")))
  terra::writeRaster(mapbiomas, file.path(WD.tmp, paste0("mapbiomas_", basin_id, ".tif")), overwrite = TRUE)
  save(hydrometric, file = file.path(WD.tmp, paste0("hydrometric_",basin_id,".RData")))
  save(lsm, file = file.path(WD.tmp, paste0("mapbiomas_lsm_",basin_id,".RData")))
  save(lsm_buf100, file = file.path(WD.tmp, paste0("mapbiomas_lsm_buf100_",basin_id,".RData")))
  save(lsm_buf500, file = file.path(WD.tmp, paste0("mapbiomas_lsm_buf500_",basin_id,".RData")))
  
}
 

basin_select <- st_read(file.path(dir_data, "hydroshed", "basin_select.shp"))
hydrobasin <- st_read(file.path(dir_data, "hydroshed", "subbas_3bas.shp"))
length(which(basin_select$HYBAS_ID%in%hydrobasin$HYBAS_ID))
length(which(hydrobasin$HYBAS_ID%in%basin_select$HYBAS_ID))
length(which(duplicated(basin_select$HYBAS_ID)))


vec_lsm <- vector()
j=0
for (i in basin_select$HYBAS_ID[21:1703]){
  j=j+1
  vec_lsm[j] = list.files(path=dir_data, pattern = paste0("_",i), full.names = TRUE)[1]
  
}

which(is.na(vec_lsm))
i=20
id <- 6090469940 
ID.name = paste0("Basin_",i,"_",id)
#check if directory exist
  dir.create(file.path(dir_data, ID.name))
  WD.new = paste(dir_data, ID.name,sep="/")
  
  mapbiomas_crop <- terra::crop(mapbiomas_full,basin_select[i,])
  mapbiomas_mask <- terra::mask(mapbiomas_crop,vect(basin_select[i,]))
  mapbiomas_mask_utm <- terra::project(mapbiomas_mask,"+init=epsg:32721",method="ngb")
  
  WD.tmp = list.files(path=dir_data, pattern = paste0("Basin_",i,"_"), full.names = TRUE)
  
  terra::writeRaster(mapbiomas_mask_utm, file.path(WD.tmp, paste0("mapbiomas_", 6090469940, ".tif")), overwrite = TRUE)
  
