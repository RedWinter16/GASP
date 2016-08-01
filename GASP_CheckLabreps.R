GASP_CheckLabreps<-function(df, minimum, replicates,index)
{
  #####Modify rows that do not meet the minimum requiremnt in a number of replicates
  #index indicates sets of 3/4/10/etc.lab replicates
  #minimum is the minimum relative abundance that a entry must have
  #replicates is the minmum of the laboratory replicates that must meet the requiements for relative abundance

  #create vector based on the index so we can look at a subset of the main data
  vec<-seq(5, ncol(df), index)

      for (i in 1:length(vec))
      {
        if (vec[i] <= (ncol(df)-(index-1)))
        {
                for (j in 1:nrow(df))
                {
                      Cells2check<-df[j,(vec[i]):(vec[i]+(index-1))]
                      if (length(which(Cells2check >= minimum)) < replicates)
                      {
                            df[j,(vec[i]):(vec[i]+(index-1))]<-0 #set the cell to zero
                      }
                      
                }
        }
      }
  return(df)
  }
