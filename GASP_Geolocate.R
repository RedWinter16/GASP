GASP_Geolocate<-function(df, site, filename)
  {
      ####This section will remove taxa without genus and species  
      df<-filter(df, Genus != "")
      df$Species<-gsub("$", "*", df$Species)
      ####trnl only
      df<-filter(df, Primer.Set == "trnl")
    
      #query_genus<-filter(df, df$Species == "*")
      #query_genus<-paste(query_genus$Genus, query_genus$Species)
      #query_genus<-unique(query_genus)
      
      query_species<-filter(df, df$Species != '*')
      query_species<-query_species$Species
      query_species<-unique(query_species)
    
      #map("state")
    
      #if(length(query_genus) > 1)
      #{
      #  Gout<-occ(query = query_genus, from =  'gbif', limit = 100, has_coords = TRUE, gbifopts = list(continent='north_america', hasGeospatialIssue=FALSE))
      #  Gout<-na.omit(occ2df(Gout))
      #  Gout<-Gout %>% 
      #    dedup() %>%   
      #    coord_impossible(lon = 'longitude', lat = 'latitude') %>%
      #    coord_incomplete(lon = 'longitude', lat = 'latitude') %>%
      #    coord_unlikely(lon = 'longitude', lat = 'latitude')
      #}
      
      if(length(query_species) > 1)
          {
            Sout<-occ(query = query_species, from =  'gbif' , has_coords = TRUE, gbifopts = list(country = 'US', hasGeospatialIssue=FALSE))
            Sout<-na.omit(occ2df(Sout))
            Sout<-Sout %>% 
              dedup() %>%   
              coord_impossible(lon = 'longitude', lat = 'latitude') %>%
              coord_incomplete(lon = 'longitude', lat = 'latitude') %>%
              coord_unlikely(lon = 'longitude', lat = 'latitude')
            result<-Sout#rbind.data.frame(Gout, Sout, make.row.names = TRUE)
            write.csv(result, paste(filename, site, 'Geoloctaion of species data via GBIF.csv' ))
          }

}
    