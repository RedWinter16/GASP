GASP_RMzeroTaxa<-function(df)
{
  #####Remove taxa rows that become zero across all columns, these are now empty
  checkcol<-ncol(df)
  if(checkcol > 5)
      {
        df<-mutate(df, Eliminate_if_0 = apply(df[,5:ncol(df)],1, sum))
      }  
  i<-1
      while (i <= nrow(df))
      {
        ifelse (df[[i, ncol(df)]] == 0, df<-df[-i,], i<-i+1)
            
      }
  if(checkcol > 5)
      {
      df<-df[, -ncol(df)]
      }
  return(df)
  
}