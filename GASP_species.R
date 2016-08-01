GASP_species<-function(df, filename)
  
  {
  library(reshape2)
  library(rbison)
  library(spocc)
  library(mapr)
  
  filename<-strsplit(filename, ".csv")
  df_MA<-select(df, Primer.Set:ConsensusLineage, contains("MA"))
  df_MA<-GASP_RMzeroTaxa(df_MA)
  MA<-length(df_MA[,1])
  
  df_SC<-select(df, Primer.Set:ConsensusLineage, contains("SC"))
  df_SC<-GASP_RMzeroTaxa(df_SC)
  SC<-length(df_SC[,1])
  
  df_FL<-select(df, Primer.Set:ConsensusLineage, contains("FL"))
  df_FL<-GASP_RMzeroTaxa(df_FL)
  FL<-length(df_FL[,1])
  
  result<-data.frame(matrix(0, nrow = 3, ncol = 5))

  colnames(result)<-c('ALL', '18S', 'ITS', 'trnl', 'CO1')
  result[1,1]<-MA
  result[2,1]<-SC
  result[3,1]<-FL
  
  
  MA<-filter(df_MA, Primer.Set == "18S")
  MA<-GASP_RMzeroTaxa(MA)
  MA<-length(MA[,1])
  result[1,2]<-MA
  SC<-filter(df_SC, Primer.Set == "18S")
  SC<-GASP_RMzeroTaxa(SC)
  SC<-length(SC[,1])
  result[2,2]<-SC
  FL<-filter(df_FL, Primer.Set == "18S")
  FL<-GASP_RMzeroTaxa(FL)
  FL<-length(FL[,1])
  result[3,2]<-FL
  
  MA<-filter(df_MA, Primer.Set == "ITS")
  MA<-GASP_RMzeroTaxa(MA)
  MA<-length(MA[,1])
  result[1,3]<-MA
  SC<-filter(df_SC, Primer.Set == "ITS")
  SC<-GASP_RMzeroTaxa(SC)
  SC<-length(SC[,1])
  result[2,3]<-SC
  FL<-filter(df_FL, Primer.Set == "ITS")
  FL<-GASP_RMzeroTaxa(FL)
  FL<-length(FL[,1])
  result[3,3]<-FL
  
  MA<-filter(df_MA, Primer.Set == "trnl")
  MA<-GASP_RMzeroTaxa(MA)
  MA<-length(MA[,1])
  result[1,4]<-MA
  SC<-filter(df_SC, Primer.Set == "trnl")
  SC<-GASP_RMzeroTaxa(SC)
  SC<-length(SC[,1])
  result[2,4]<-SC
  FL<-filter(df_FL, Primer.Set == "trnl")
  FL<-GASP_RMzeroTaxa(FL)
  FL<-length(FL[,1])
  result[3,4]<-FL
  
  MA<-filter(df_MA, Primer.Set == "CO1")
  MA<-GASP_RMzeroTaxa(MA)
  MA<-length(MA[,1])
  result[1,5]<-MA
  SC<-filter(df_SC, Primer.Set == "CO1")
  SC<-GASP_RMzeroTaxa(SC)
  SC<-length(SC[,1])
  result[2,5]<-SC
  FL<-filter(df_FL, Primer.Set == "CO1")
  FL<-GASP_RMzeroTaxa(FL)
  FL<-length(FL[,1])
  result[3,5]<-FL
  
  row.names(result)<-c("MA Site", "SC Site", "FL Site")
  write.csv(result, paste(filename, "number of taxa per site.csv"))
  
  #map("state")
  #text(-120, 49.5, paste(filename), cex = 0.6)
  #query(df_MA, 'MA', filename, 'red', data_dir)
  #query(df_SC, 'SC', filename, 'blue', data_dir)
  #query(df_FL,'FL', filename, 'green', data_dir)

}






query<-function(df, site, filename, color, data_dir)
{
  
  df<-df[,c(1,4)]
  df<-filter(df, df[,1] != c("18S", "CO1"))
  df<-df[,-1]
  df<-colsplit(df, c("s:", "; s_"), c("Consensus", "Species"))
 
  df<-unique(df[,2])
  df<-df[-1]
  
  df<-gsub("_sp.*$", "", df)
  df<-gsub(" sp .*$", "", df)
  df<-gsub(" sp.*$", "", df)
  df<-gsub(" cf .*$", "", df)
  df<-gsub("_", " ", df)
  df<-gsub(";.*$", "", df)
  df<-gsub("uncultured", "", df)
  df<-gsub("^\\s+|\\s+&", "", df)
  df<-unique(df)
  
  
  if(length(df) > 0)
      {
      #out<-occ(query = df, from =  'bison', limit = 500, bisonopts = what = 'counties')
      #out<-na.omit(occ2df(out))
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
