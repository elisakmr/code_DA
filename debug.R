




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
