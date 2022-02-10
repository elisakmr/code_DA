




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
  
  
  basin_rios
  df_bv
  
  for (i in df_bv$ID){
    df_bv[which(i==df_bv$ID),'d.rios']<-basin_rios[which(i==basin_rios$HYBAS_ID),]$d.rios
  }
  
  head(df_bv)

  var(unlist(df_bv %>% filter(maj=='pasture') %>% dplyr::select(d.rios)  ))
    
  
  var(unlist(df_bv %>% filter(maj=='soybean') %>% dplyr::select(value)  ))
  mean(unlist(df_rv %>% filter(maj=='pasture') %>% dplyr::select(value)  ))
  
  df_rv

  st_write(rios.buf,file.path(WD.tmp,paste0("riosbuf_",basin_id,".shp")))
  terra::writeRaster(mapbiomas_buf[[35]], file.path(WD.tmp, paste0("mapbiomasbuf_", basin_id, ".tif")), overwrite = TRUE)

  
  for(basin_id in basin_select$HYBAS_ID[880:890])   { 
    
    # directory
    WD.tmp = list.files(path=dir_data, pattern = paste0("_",basin_id), full.names = TRUE)[1]
    
    # loading landuse data
    mapbiomas <- rast(file.path(WD.tmp, paste0("mapbiomas_", basin_id, ".tif")))
    rios <- st_read(file.path(WD.tmp,paste0("riossimple_",basin_id,".shp"))) 
    

      rios.buf = st_buffer(rios, dist = tampon)
      mapbiomas_buf = terra::mask(mapbiomas,vect(rios.buf))
      
      #computing landscape metrics at basin scale (MapBiomas)
      mapbiomas.basin.lsm = calculate_lsm(mapbiomas_buf,level="class",metric = c("ca","pland","np"),directions =8)
      
      mapbiomas.basin.df = as.data.frame(mapbiomas.basin.lsm)
      

    save(mapbiomas.basin.df, file = file.path(WD.tmp, paste0("lsm_rip",tampon, "_",basin_id,".RData")))
    
  }
  
st_simplify(rios)
l<-which(basin_id==hydrobasin$HYBAS_ID)
crop_list <- st_intersects(hydrobasin[l,],hydronetwork)
crop_shp <- hydronetwork[unlist(crop_list),]
crop.utm <- st_transform(crop_shp, crs="EPSG:32721") 

rios.buf = st_buffer(crop.utm[-633,], dist = tampon)
length(st_is_valid(hydronetwork))
st_write(crop.utm, file.path(WD.tmp, paste0("riossimple_",basin_id,".shp")), delete_layer = TRUE)


#########################################   buffers to be plotted
basin_id = 6090573160 # 6090499060
WD.tmp = list.files(path=dir_data, pattern = paste0("_",basin_id), full.names = TRUE)[1]
rios <- st_read(file.path(WD.tmp,paste0("riossimple_",basin_id,".shp"))) 
rios.buf = st_buffer(rios, dist = 500)
st_write(rios.buf, file.path(dir_data, 'processed', 'water', paste0("riosbuf_", basin_id,".shp")), delete_layer = TRUE)

##############################################

WD.tmp = list.files(path=dir_data, pattern = paste0("_",basin_id), full.names = TRUE)[1]

# loading landuse data
mapbiomas <- rast(file.path(WD.tmp, paste0("mapbiomas_", basin_id, ".tif")))
rios <- st_read(file.path(WD.tmp,paste0("riossimple_",basin_id,".shp"))) 
rios.buf = st_buffer(rios, dist = tampon)
area1 <- sum(st_area(rios.buf))
rios.buf.union <- st_union(rios.buf)
area.uni <- st_area(rios.buf.union)
plot(rios.buf)
plot(rios.buf.union)
st_write(rios.buf.union, file.path(dir_data, 'processed', 'water', paste0("riosbufuni_", basin_id,".shp")))
mapbiomas_buf1 = terra::mask(mapbiomas,vect(rios.buf))
mapbiomas_buf2 = terra::mask(mapbiomas,vect(rios.buf.union))
plot(mapbiomas_buf1[[1]])
plot(mapbiomas_buf2[[2]])

############# crop past maximum ################

max((df_lsm_relatif %>% filter(past.rel>60))$soy.rel)

mean((df_basin_foret %>% filter(maj=='pasture'))$value)
mean((df_rios_foret %>% filter(maj=='soybean'))$value)



