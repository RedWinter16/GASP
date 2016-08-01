#****This function will correctly set the workspace for the GASP scripts. It will set the correct proxy, load libraries, find working directory and create 
#****a new folder with todays date and then load in all other functions.

GASP_SetWorkspace<-function()
{
  #clears current values and functions from workspace
  rm(list=ls())
  
  #set proxy for MIT
  Sys.setenv(http_proxy="llproxy.llan.ll.mit.edu:8080")
  
  #Load libraries
  library(dplyr)
  library(reshape2)
  library(magrittr)
  library(ggmap)
  library(rgbif)
  library(spocc)
  library(scrubr)
  
  #Set working directory, this will need to be able to change for all users but it does not right now.
  wd<-'C:/Users/ch26466/Documents/GASP/R'
  setwd(wd)
  
  #create new data directory with todays date and return
  today_date<-Sys.Date()
  data_dir<-file.path(wd, paste(today_date, "data"))
  dir.create(data_dir, showWarnings = FALSE)
  
  #load all needed functions from the working directory
  source(file.path(wd, 'GASP_Main.R'))
  source(file.path(wd, 'GASP_CheckLabreps.R'))
  source(file.path(wd, 'GASP_RMzeroTaxa.R'))
  source(file.path(wd, 'GASP_CheckSitereps.R'))
  source(file.path(wd, 'GASP_CombineLabreps.R'))
  source(file.path(wd, 'GASP_getdata.R'))
  source(file.path(wd, 'GASP_Output.R'))
  source(file.path(wd, 'GASP_Geolocate.R'))
  source(file.path(wd, 'GASP_GeoMap.R'))
  
  return(data_dir)
}