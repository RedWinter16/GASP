GASP_Geolocate<-function(df, site, filename)
  {
  library(rgbif)
  library(spocc)
  library(scrubr)
  library(mapr)
  
  ####This section will remove taxa without genus and species  
  df<-filter(df, Genus != "")
  df$Species<-gsub("$", "*", df$Species)

  query_genus<-filter(df, df$Species == "*")
  query_genus<-paste(query_genus$Genus, query_genus$Species)
  query_genus<-unique(query_genus)
  
  query_species<-filter(df, df$Species != '*')
  query_species<-query_species$Species
  query_species<-unique(query_species)

  #map("state")
 for (i in c(query_genus, query_species))
  {
    df<-i
    if(length(df) > 0)
        {
        out<-occ(query = df, from =  'gbif', limit = 1, gbifopts = list(hasGeospatialIssue=FALSE))
        out<-na.omit(occ2df(out))
        out<-out %>% dedup()
        result_counties<-data.frame(matrix(1, nrow = 0, ncol = 5))
        colnames(result_counties)<-c("record_id", "total", "county_name", "state", "species")
        specimen<-as.character(df[1])
        out<-bison(specimen, count = 100, what = 'counties')
          if(is.null(out$counties) == FALSE)
          {
                out<-as.data.frame(out$counties)
                out<-mutate(out, species = specimen)
                result_counties<-rbind.data.frame(result_counties, out)
          }
          if(length(df)>1)
          {
              for(i in 2:length(df))
              {
                    specimen<-as.character(df[i])
                    out<-bison(specimen, count = 100, what = "counties")
                    if(is.null(out$counties) == FALSE)
                          {
                                out<-as.data.frame(out$counties)
                                out<-mutate(out, species = specimen)
                                result_counties<-rbind.data.frame(result_counties, out)
                          }
              }
          }
            write.csv(result_counties, paste(filename, site, "species.csv"), row.names = FALSE)
    }
 }
  }
