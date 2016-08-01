# Inputs are a data frame of Jonah data split into slides and swabs, the type (slides or swabs), and the answers (Y/N) to rm_controls and rm_collect1

GASP_CombineLabreps<- function(df, index)
{
  #####Combine lab replicates, this is very similar to the check lab replicates function but with no minimum requirements
  #index indicates sets of 3/4/10/etc.lab replicates
  
  #create vector based on the index so we can look at a subset of the main data
  df_new<-df[,1:4]
  vec<-seq(5, ncol(df), index)
  colname<-colnames(df)
  colname<-colname[vec]
  colname<-gsub('[0-9]+', '', colname)
  
  for (i in 1:length(vec))
  {
      if(vec[i] < ncol(df))
      {
          Columns2combine<-df[,(vec[i]):(vec[i]+(index-1))]
          coltemp<-apply(Columns2combine, 1, mean)
          df_new<-cbind.data.frame(df_new, coltemp)
      }
  }
  
  colnames(df_new)[5:ncol(df_new)]<-colname
  return(df_new)
  
}
