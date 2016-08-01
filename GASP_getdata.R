GASP_getdata<-function(filename, data_dir)
{
  ####Read in  the file containing Jonah data
  df<-read.csv(filename)
  
  ####Move the last column over so Meta data is in the first four columns
  df<-df[, c(1:3, ncol(df), 4:(ncol(df)-1))]
  tempdf<-df[,5:ncol(df)]
  df<-df[,1:4]
  ####Calculate the relative abundance
  #First find the total reads for the sample
  sumrow<-apply(tempdf, 2, sum)
  #Second loop through all taxa for the sample and divide bu the total reads per sample
  for (i in 1:nrow(tempdf))
        {
            for (j in 1:ncol(tempdf))
              {
              tempdf[i,j]<-tempdf[[i,j]]/sumrow[[j]]
              }
  }
  df<-cbind.data.frame(df, tempdf)
  ####Write out the relative abudace into a new file
  #find system date
  date<-Sys.Date()
  #create file name
  df_save<-file.path(data_dir, paste(date, "Relative Abundance.csv"))
  #write out as .csv into data folder
  write.csv(df, df_save, row.names = FALSE)
  return(df)
}