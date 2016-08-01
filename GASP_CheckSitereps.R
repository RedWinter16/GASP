GASP_CheckSitereps<-function(df, replicates, collectVec, siteVec)
{
  #####Modify rows that do not meet the minimum requiremnt in a number of replicates
  #vec indicates sets of sites
  #replicates is the minmum of the laboratory replicates that must meet the requiements for relative abundance

  columns<-colnames(df)
  sel<-which(match(siteVec, collectVec)  > 0)
  sel2<-siteVec[sel]
  sel<-which(match(columns, sel2) > 0)
  
      if(length(collectVec) > 0 && length(siteVec) > 0 && length(sel) > replicates)
      {
          for (i in 1:nrow(df))
          {
                if (length(which(df[i, sel] > 0)) < replicates)
                {
                      df[i, sel]<-0 #set to zero
                }
          }
      
      }
  return(df)
  
}