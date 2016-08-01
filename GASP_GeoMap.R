GASP_GeoMap<-function(data.dir)
{

  
  files<-list.files(data.dir, pattern = 'via GBIF.csv')
  USAmap<-ggmap(get_stamenmap(bbox = c( -125.0011, 24.9493, -66.9326, 49.5904), zoom = 5, maptype= 'toner-lite'))
  
  for( i in 1:length(files))
  {
      filename<-files[i]
      df<-read.csv(file.path(data.dir, filename))
      title<-data.frame(strsplit(filename, '/'))
      title<-as.character(title[[nrow(title),1]])
      
      geomap<-USAmap + geom_point(aes(x=longitude, y=latitude), colour = 'orange', size = 2, alpha = 0.1,  data=df)+
      geom_density2d(aes(x=longitude, y=latitude), data=df, col="orange")+
      labs(title = title)
      filename<-strsplit(filename, ".csv")
      ggsave(geomap, filename = file.path(data.dir, paste(filename, '.jpeg')))
      print(geomap)
  }
  
}