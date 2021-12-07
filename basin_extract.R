

library(plyr)
library(dplyr)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(maps)
library(nnet)
library(latticeExtra)
library(sp)
library(rgdal)
library(raster)
library(sf)
library(stars)
library(exactextractr)
library(tictoc)
library(terra)
library(landscapemetrics)

#FUNCTIONS
##################################################################################
i<-4
MapBiomas.im <- rast(file.path(dir_data, "mapbiomas", "biomas1985_2019_mg.tif"))
Modis.im <- rast(file.path(dir_data, "pangaea", "panga_braz84_mt.tif"))
basin.sf <- st_read(file.path(dir_data, "hydroshed","hydrobasins_level9", "hydrobasin_mg_lev9.shp"))
rios.sf <- st_read(file.path(dir_data, "water_surface", "rios_simples", "rios_simple84.shp"))
nascentes.sf <- st_read(file.path(dir_data, "water_surface", "nascentes", "nascente84.shp"))
reservatorios.sf <- st_read(file.path(dir_data, "water_surface", "lago", "input_lago_reservatorio_icv.shp"))
  
extractBasin2 = function(i,MapBiomas.im,Modis.im,basin.sf,rios.sf,nascentes.sf,reservatorios.sf){
  #Map biomas land cover areas

  basin.tmp = basin.sf[i,]
  basin.tmp.utm = st_transform(basin.sf[i,], crs="+init=epsg:32721")
  
  ID.name = paste0("Basin_",i,"_",basin.tmp$HYBAS_ID)
  #check if directory exist
  if (!file.exists(ID.name)){
    dir.create(file.path(dir_data, ID.name))
    WD.new = paste(dir_data, ID.name,sep="/")
  }
  
  #crop and mask map biomas over basin 
  MapBiomas.basin.name = paste0(WD.new,"/MapBiomas_basin_",ID.name,".tif")
  MapBiomas.basin.im = terra::crop(MapBiomas.im,basin.tmp)
  MapBiomas.basin.im = terra::project(MapBiomas.basin.im,as.character(crs(basin.tmp.utm)),method="ngb")
  MapBiomas.basin.im = terra::mask(MapBiomas.basin.im,vect(basin.tmp.utm))
  terra::writeRaster(MapBiomas.basin.im,MapBiomas.basin.name,filetype="GTiff",datatype = "INT1U",overwrite=T)
  
  #crop and mask map modis over basin 
  Modis.basin.name = paste0(WD.new,"/Modis_basin_",ID.name,".tif")
  Modis.basin.im = terra::crop(Modis.im,basin.tmp)
  Modis.basin.im = terra::project(Modis.basin.im,as.character(crs(basin.tmp.utm)),method="ngb")
  Modis.basin.im = terra::mask(Modis.basin.im,vect(basin.tmp.utm))
  terra::writeRaster(Modis.basin.im,Modis.basin.name,filetype="GTiff",datatype = "INT1U",overwrite=T)
  
  #intersection with basin
  sf::sf_use_s2(FALSE)
  Nascentes.basin = st_intersection(nascentes.sf,basin.sf[i,])
  Reservatorios.basin = st_intersection(reservatorios.sf,basin.sf[i,])
  #Rios.basin = st_intersection(Rios.poly,basin.sf[i,])
  
  if (nrow(Nascentes.basin) > 0){
    #project to UTM
    Nascentes.basin.utm = st_transform(Nascentes.basin,st_crs(MapBiomas.basin.im))
    Basin.nascentes.name = paste0(WD.new,"/Nascentes_basin_",ID.name,".shp")
    st_write(Nascentes.basin.utm,Basin.nascentes.name,driver = "ESRI Shapefile")
  }
  
  if (nrow(Reservatorios.basin) > 0){
    Reservatorios.basin.utm = st_transform(Reservatorios.basin,st_crs(MapBiomas.basin.im))
    Basin.reservatorios.name = paste0(WD.new,"/Reservatorios_basin_",ID.name,".shp")
    st_write(Reservatorios.basin.utm,Basin.reservatorios.name,driver = "ESRI Shapefile")
  }
  
  # if (nrow(Rios.basin) > 0){
  #   Rios.basin.utm = st_transform(Rios.basin,st_crs(MapBiomas.basin.im))
  #   Basin.rios.name = paste0(WD.new,"/Rios_basin_",ID.name,".shp")
  #   st_write(Rios.basin.utm,Basin.rios.name,driver = "ESRI Shapefile")
  # }
}  


#function
processMetricsBasin2 = function(i,MapBiomas.im,Modis.im,basin.sf,rios.sf,nascentes.sf,reservatorios.sf){
  #Map biomas land cover areas
  #crop and mask map biomas over basin
  MapBiomas.basin.im = terra::crop(MapBiomas.im,basin.sf[i,])
  MapBiomas.basin.im = raster::mask(raster(MapBiomas.basin.im),basin.sf[i,])
  plot(MapBiomas.basin.im)
  plot(st_geometry(basin.sf[i,]),add=T)
  # 
  #project to utm
  MapBiomas.basin.im.utm = terra::project(rast(MapBiomas.basin.im),"+init=epsg:32721",method="ngb")
  basin.utm.sf = st_transform(basin.sf[i,],st_crs(MapBiomas.basin.im.utm))
  BasinAreaHa = st_area(basin.utm.sf)/10000
  
  #compute landscape metrics at basin scale (MapBiomas)
  MapBiomas.basin.lsm = calculate_lsm(raster(MapBiomas.basin.im.utm),level="class",metric = c("ca","pland","np"),directions =8)
  MapBiomas.basin.lsm.df = as.data.frame(MapBiomas.basin.lsm)
  MapBiomas.basin.lsm.df$ID = i
  
  #compute landscape metrics at basin scale (MODIS-INPE)
  Modis.basin.im = terra::crop(Modis.im,basin.sf[i,])
  Modis.basin.im = raster::mask(raster(Modis.basin.im),basin.sf[i,])
  plot(Modis.basin.im)
  plot(st_geometry(basin.sf[i,]),add=T)
  
  #project to utm
  Modis.basin.im.utm = terra::project(rast(Modis.basin.im),"+init=epsg:32721",method="ngb")
  Modis.basin.lsm = calculate_lsm(raster(Modis.basin.im.utm),level="class",metric = c("ca"),directions =8)
  Modis.basin.lsm.df = as.data.frame(Modis.basin.lsm)
  Modis.basin.lsm.df$ID = i
  
  #intersection with basin
  sf::sf_use_s2(FALSE)
  Nascentes.basin = st_intersection(Nascentes.poly,basin.sf[i,])
  Reservatorios.basin = st_intersection(Reservatorios.poly,basin.sf[i,])
  Rios.basin = st_intersection(Rios.poly,basin.sf[i,])
  
  #plot
  plot(st_geometry(Nascentes.basin),add=T,col="red",lwd=2)
  plot(st_geometry(Reservatorios.basin),add=T,col="red")
  plot(st_geometry(Rios.basin),add=T,col="blue")
  # 
  #project to UTM
  Nascentes.basin.utm = st_transform(Nascentes.basin,st_crs(MapBiomas.basin.im.utm))
  Rios.basin.utm = st_transform(Rios.basin,st_crs(MapBiomas.basin.im.utm))
  Reservatorios.basin.utm = st_transform(Reservatorios.basin,st_crs(MapBiomas.basin.im.utm))
  
  #number of sources in AOI
  Basin.nb.nascentes = nrow(Nascentes.basin.utm)
  #cumulated length of rivers in AOI
  Basin.length.rios = sum(st_length(Rios.basin.utm))
  
  #compute buffers (100 m)
  Basin.nascentes.utm.buf100 = st_buffer(Nascentes.basin.utm,dist = 100,)
  Basin.rios.utm.buf100 = st_buffer(Rios.basin.utm,dist = 100,)
  Basin.reservatorios.utm.buf100 = st_buffer(Reservatorios.basin.utm,dist = 100,)
  plot(MapBiomas.basin.im.utm)
  plot(st_geometry(basin.utm.sf),add=T,lwd=3)
  plot(st_geometry(Basin.nascentes.utm.buf100),border="red",lwd=2,add=T)
  plot(st_geometry(Basin.rios.utm.buf100),add=T,border="green")
  plot(st_geometry(Basin.reservatorios.utm.buf100),add=T,border="blue")
  # 
  #compute deforestation in buffer (based on MapBiomas but could also be based on PRODES and CRA2008)
  Buf.im = setValues(MapBiomas.basin.im.utm,numeric(ncell(MapBiomas.basin.im.utm)))
  Basin.nascentes.utm.buf100.im = terra::rasterize(vect(Basin.nascentes.utm.buf100),Buf.im,background=0)
  Basin.rios.utm.buf100.im = terra::rasterize(vect(Basin.rios.utm.buf100),Buf.im,background=0)
  Basin.reservatorios.utm.buf100.im = terra::rasterize(vect(Basin.reservatorios.utm.buf100),Buf.im,background=0)
  Basin.hydro.Buf100.im = ((Basin.nascentes.utm.buf100.im + Basin.rios.utm.buf100.im + Basin.reservatorios.utm.buf100.im) > 0)
  # AOI.hydro.Buf100.im = ((AOI.nascentes.utm.buf100.im + AOI.rios.utm.buf100.im) > 0)
  
  plot(Basin.hydro.Buf100.im)
  # 
  Basin.hydro.Buf100.MapBiomas.im = Basin.hydro.Buf100.im * MapBiomas.basin.im.utm
  Basin.hydro.Buf100.MapBiomas.im[Basin.hydro.Buf100.MapBiomas.im == 0]=NA
  MapBiomas.basin.lsm.buf = calculate_lsm(raster(Basin.hydro.Buf100.MapBiomas.im),level="class",metric = c("ca","pland","np"),directions =8)
  MapBiomas.basin.lsm.buf.df = as.data.frame(MapBiomas.basin.lsm.buf)
  MapBiomas.basin.lsm.buf.df$ID = i
  
  #set list with final metrics for basin
  Basin.metrics.ls = list("BasinID" = basin.sf[i,]$ID,
                          "BasinArea" = BasinAreaHa,
                          "BasinMapBiomas" = MapBiomas.basin.im.utm,
                          "BasinModis" = Modis.basin.im.utm,
                          "BasinRios" = Rios.basin.utm,
                          "BasinNascentes" = Nascentes.basin.utm,
                          "BasinReservatorios" = Reservatorios.basin.utm,
                          "BasinMapBiomas2019" = MapBiomas.basin.lsm.df,
                          "BasinBufMapBiomas2019" = MapBiomas.basin.lsm.buf.df,
                          "BasinModis2017" = Modis.basin.lsm.df,
                          "NascentesNB" = Basin.nb.nascentes,
                          "RiosLength" = Basin.length.rios)
  
  saveRDS(Basin.metrics.ls,file = paste0("Basin.metrics.",i,".rds"))
  # write.table("",paste0("Basin.",i,".txt"))
  return(Basin.metrics.ls)
}

#function
processMetricsBasin3 = function(i){
  
  #set working directory
  WD.tmp = list.files(path=dir_data, pattern = paste0("Basin_",i,"_"), full.names = TRUE)
  
  # loading files
  MapBiomas.im = rast(file.path(dir_data, "mapbiomas", "biomas1985_2019_mg.tif"))
  basin.sf = st_read(file.path(dir_data, "hydroshed","hydrobasins_level9", "hydrobasin_mg_lev9.shp"))
  Nascentes.poly = st_read(file.path(dir_data, "water_surface", "nascentes", "nascente84.shp"))
  Reservatorios.poly = st_read(file.path(dir_data, "water_surface", "lago", "input_lago_reservatorio_icv.shp"))
  Rios.poly = st_read(file.path(dir_data, "water_surface", "rios_simples", "rios_simple84.shp"))
  
  #Map biomas land cover areas
  #crop and mask map biomas over basin
  MapBiomas.basin.im = terra::crop(MapBiomas.im,basin.sf[i,])
  MapBiomas.basin.im = raster::mask(raster(MapBiomas.basin.im),basin.sf[i,])
  plot(MapBiomas.basin.im)
  plot(st_geometry(basin.sf[i,]),add=T)
  # 
  #project to utm
  MapBiomas.basin.im.utm = terra::project(rast(MapBiomas.basin.im),"+init=epsg:32721",method="ngb")
  basin.utm.sf = st_transform(basin.sf[i,],st_crs(MapBiomas.basin.im.utm))
  
  BasinAreaHa = st_area(basin.utm.sf)/10000
  
  #compute landscape metrics at basin scale (MapBiomas)
  MapBiomas.basin.lsm = calculate_lsm(raster(MapBiomas.basin.im.utm),level="class",metric = c("ca","pland","np"),directions =8)
  MapBiomas.basin.lsm.df = as.data.frame(MapBiomas.basin.lsm)
  MapBiomas.basin.lsm.df$ID = i
  
  # #compute landscape metrics at basin scale (MODIS-INPE)
  # Modis.basin.im = terra::crop(Modis.im,basin.sf[i,])
  # Modis.basin.im = raster::mask(raster(Modis.basin.im),basin.sf[i,])
  # plot(Modis.basin.im)
  # plot(st_geometry(basin.sf[i,]),add=T)
  # 
  # #project to utm
  # Modis.basin.im.utm = terra::project(rast(Modis.basin.im),"+init=epsg:32721",method="ngb")
  Modis.basin.lsm = calculate_lsm(raster(Modis.basin.im.utm),level="class",metric = c("ca"),directions =8)
  Modis.basin.lsm.df = as.data.frame(Modis.basin.lsm)
  Modis.basin.lsm.df$ID = i
  
  #intersection with basin
  sf::sf_use_s2(FALSE)
  Nascentes.basin = st_intersection(Nascentes.poly,basin.sf[i,])
  Reservatorios.basin = st_intersection(Reservatorios.poly,basin.sf[i,])
  Rios.basin = st_intersection(Rios.poly,basin.sf[i,])
  
  #plot
  plot(st_geometry(Nascentes.basin),add=T,col="red",lwd=2)
  plot(st_geometry(Reservatorios.basin),add=T,col="red")
  plot(st_geometry(Rios.basin),add=T,col="blue")
  # 
  #project to UTM
  Nascentes.basin.utm = st_transform(Nascentes.basin,st_crs(MapBiomas.basin.im.utm))
  Rios.basin.utm = st_transform(Rios.basin,st_crs(MapBiomas.basin.im.utm))
  Reservatorios.basin.utm = st_transform(Reservatorios.basin,st_crs(MapBiomas.basin.im.utm))
  
  #number of sources in AOI
  Basin.nb.nascentes = nrow(Nascentes.basin.utm)
  #cumulated length of rivers in AOI
  Basin.length.rios = sum(st_length(Rios.basin.utm))
  
  #compute buffers (100 m)
  Basin.nascentes.utm.buf100 = st_buffer(Nascentes.basin.utm,dist = 100,)
  Basin.rios.utm.buf100 = st_buffer(Rios.basin.utm,dist = 100,)
  Basin.reservatorios.utm.buf100 = st_buffer(Reservatorios.basin.utm,dist = 100,)
  plot(MapBiomas.basin.im.utm)
  plot(st_geometry(basin.utm.sf),add=T,lwd=3)
  plot(st_geometry(Basin.nascentes.utm.buf100),border="red",lwd=2,add=T)
  plot(st_geometry(Basin.rios.utm.buf100),add=T,border="green")
  plot(st_geometry(Basin.reservatorios.utm.buf100),add=T,border="blue")
  # 
  #compute deforestation in buffer (based on MapBiomas but could also be based on PRODES and CRA2008)
  Buf.im = setValues(MapBiomas.basin.im.utm,numeric(ncell(MapBiomas.basin.im.utm)))
  Basin.nascentes.utm.buf100.im = terra::rasterize(vect(Basin.nascentes.utm.buf100),Buf.im,background=0)
  Basin.rios.utm.buf100.im = terra::rasterize(vect(Basin.rios.utm.buf100),Buf.im,background=0)
  Basin.reservatorios.utm.buf100.im = terra::rasterize(vect(Basin.reservatorios.utm.buf100),Buf.im,background=0)
  Basin.hydro.Buf100.im = ((Basin.nascentes.utm.buf100.im + Basin.rios.utm.buf100.im + Basin.reservatorios.utm.buf100.im) > 0)
  # AOI.hydro.Buf100.im = ((AOI.nascentes.utm.buf100.im + AOI.rios.utm.buf100.im) > 0)
  
  plot(Basin.hydro.Buf100.im)
  # 
  Basin.hydro.Buf100.MapBiomas.im = Basin.hydro.Buf100.im * MapBiomas.basin.im.utm
  Basin.hydro.Buf100.MapBiomas.im[Basin.hydro.Buf100.MapBiomas.im == 0]=NA
  MapBiomas.basin.lsm.buf = calculate_lsm(raster(Basin.hydro.Buf100.MapBiomas.im),level="class",metric = c("ca","pland","np"),directions =8)
  MapBiomas.basin.lsm.buf.df = as.data.frame(MapBiomas.basin.lsm.buf)
  MapBiomas.basin.lsm.buf.df$ID = i
  
  #set list with final metrics for basin
  Basin.metrics.ls = list("BasinID" = basin.sf[i,]$ID,
                          "BasinArea" = BasinAreaHa,
                          "BasinMapBiomas" = MapBiomas.basin.im.utm,
                          "BasinModis" = Modis.basin.im.utm,
                          "BasinRios" = Rios.basin.utm,
                          "BasinNascentes" = Nascentes.basin.utm,
                          "BasinReservatorios" = Reservatorios.basin.utm,
                          "BasinMapBiomas2019" = MapBiomas.basin.lsm.df,
                          "BasinBufMapBiomas2019" = MapBiomas.basin.lsm.buf.df,
                          "BasinModis2017" = Modis.basin.lsm.df,
                          "NascentesNB" = Basin.nb.nascentes,
                          "RiosLength" = Basin.length.rios)
  
  saveRDS(Basin.metrics.ls,file = paste0("Basin.metrics.",i,".rds"))
  # write.table("",paste0("Basin.",i,".txt"))
  return(Basin.metrics.ls)
}




##################################################################################
WD = "C:/Users/damie/Documents/CHOVE-CHUVA/Data/MapBiomas"
setwd(WD)

#open mapbiomas
List.MapBiomas = list.files(".",pattern = ".tif", recursive = F)
MapBiomas85_19_mt <- rast(List.MapBiomas)
MapBiomas19.im = biomas85_19_mt[[nlyr(MapBiomas85_19_mt)]]

#open MODIS INPE data
List.modis.Files = list.files('C:/Users/damie/Documents/CHOVE-CHUVA/Data/MODIS-INPE/', pattern = '.tif', recursive=T, include.dirs=FALSE,full.names = T)
List.modis.Files # 17 Fichiers / 2001-2017 

# TRANSFORMER LISTE EN SPAT.RASTER 
List.modis.im = rast(List.modis.Files)
List.modis.im
modis17_mt = List.modis.im[[nlyr(List.modis.im)]]
# Afficher une image
# plot(List.modis.im[[nlyr(List.modis.im)]],main="Image INPE de 2017 sur le Mato Grosso") 


# open hydrobasins from hydroshed 
hydrobasin <- st_read("C:/Users/damie/Documents/SHP/HydroSheds", layer= "hydrobasin_mg")
NumBasin = nrow(hydrobasin)  

#add ID
hydrobasin$ID = c(1:NumBasin)


#export data for each basin
####################################################################################

strt = Sys.time()
#step 2: parallel processing foreach
#do parallel loop to extract proportions
test.ids = sample(1:NumBasin,10) 
foreach(i = test.ids, .export=c("extractBasin2"), .packages=c("raster","terra","sf","parallel","doSNOW"),.inorder=T) %do% {
  cat("\n Basin = ", i)
  extractBasin2(i,MapBiomas19.im,modis17_mt,hydrobasin,Rios.poly,Nascentes.poly,Reservatorios.poly)
}
print(Sys.time() - strt)


#export data for each basin
####################################################################################

#step 1: prepare cluster
strt = Sys.time()
n.cores <- detectCores()-2
if (NumBasin < n.cores) { n.cores = NumBasin}
cl<-makeCluster(n.cores)
registerDoSNOW(cl)
cat("\t\t\t Number of cores = ",n.cores, " for number of stations = ", NumBasin ,"\n")
flush.console()

strt = Sys.time()
#step 2: parallel processing foreach
#do parallel loop to extract proportions
Metrics.basin.all.ls.par = foreach(i = 337:NumBasin, .export=c("processMetricsBasin2"), .packages=c("raster","landscapemetrics","terra","sf","parallel","doSNOW"),.inorder=T) %do% {
  #tic()
  cat("\n Basin = ", i)
  processMetricsBasin2(i)
  #toc()
}
print(Sys.time() - strt)
#step 3: Close cluster
stopCluster(cl)










#list all possible landscape metrics
as.data.frame(list_lsm())


#open shapefiles of rios, nascentes and reservatorios
Nascentes.poly = st_read("C:/Users/damie/Documents/CHOVE-CHUVA/Data/Hydro/es_mt_nascentes.shp", layer = "es_mt_nascentes")
Rios.poly = st_read("C:/Users/damie/Documents/CHOVE-CHUVA/Data/Hydro/rios_simple84_ok.shp", layer = "rios_simple84_ok")
Reservatorios.poly = st_read("C:/Users/damie/Documents/CHOVE-CHUVA/Data/Hydro/output_hidro_poly.shp", layer = "output_hidro_poly")

sf::sf_use_s2(FALSE)
Nascentes.poly = st_transform(Nascentes.poly,st_crs(hydrobasin))
Reservatorios.poly = st_transform(Reservatorios.poly,st_crs(hydrobasin))
Rios.poly = st_transform(Rios.poly,st_crs(hydrobasin))


#step 1: prepare cluster
strt = Sys.time()
n.cores <- detectCores()-2
if (NumBasin < n.cores) { n.cores = NumBasin}
cl<-makeCluster(n.cores)
registerDoSNOW(cl)
cat("\t\t\t Number of cores = ",n.cores, " for number of stations = ", NumBasin ,"\n")
flush.console()

strt = Sys.time()
#step 2: parallel processing foreach
#do parallel loop to extract proportions
Metrics.basin.all.ls.par = foreach(i = 337:NumBasin, .export=c("processMetricsBasin2"), .packages=c("raster","landscapemetrics","terra","sf","parallel","doSNOW"),.inorder=T) %do% {
  #tic()
  cat("\n Basin = ", i)
  processMetricsBasin2(i,MapBiomas19.im,modis17_mt,hydrobasin,Rios.poly,Nascentes.poly,Reservatorios.poly)
  #toc()
}
print(Sys.time() - strt)
#step 3: Close cluster
stopCluster(cl)

#3 mins



Basin.test = readRDS("Basin.metrics.13.rds")

