GASP_Output<-function(df, filename)
{
  filename<-strsplit(filename, ".csv")
  df_MA<-select(df, Primer.Set:ConsensusLineage, contains("MA"))
  df_MA<-GASP_RMzeroTaxa(df_MA)
  MA<-length(df_MA[,1])
  split_consensus(df_MA, 'MA', filename)
  
  df_SC<-select(df, Primer.Set:ConsensusLineage, contains("SC"))
  df_SC<-GASP_RMzeroTaxa(df_SC)
  SC<-length(df_SC[,1])
  split_consensus(df_SC, 'SC', filename)
  
  df_FL<-select(df, Primer.Set:ConsensusLineage, contains("FL"))
  df_FL<-GASP_RMzeroTaxa(df_FL)
  FL<-length(df_FL[,1])
  split_consensus(df_FL, 'FL', filename)
  
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
  
  
}

split_consensus<-function(df, site, filename)
{

  df<-filter(df, df[,1] != "18S")
  df<-filter(df, df[,1] != "CO1")
  df<-df[,c(1:3, 5:ncol(df), 4)]
  
  taxa<-colsplit(df[,ncol(df)],"g:", c("Consensus", "Genus"))
  taxa<-colsplit(taxa[,"Genus"], "s:", c("Genus", "Species"))

  df<-cbind.data.frame(df, taxa)

    for(i in 0:1)
    {
      i<-ncol(df)-i
      df[, i]<-gsub("_sp.*$", "", df[,i])
      df[, i]<-gsub(" sp .*$", "", df[, i])
      df[, i]<-gsub(" sp.*$", "", df[, i])
      df[, i]<-gsub(" cf .*$", "", df[, i])
      df[, i]<-gsub("_", " ", df[, i])
      df[, i]<-gsub(";.*$", "", df[, i])
      df[, i]<-gsub("uncultured", "", df[, i])
      df[, i]<-gsub(",", "", df[, i])
      df[, i]<-gsub("^\\s+|\\s+&", "", df[, i])
    }
  
  GASP_Geolocate(df, site, filename)
  write.csv(df, paste(filename, site, "Taxa Only.csv" ), row.names = FALSE)

}

