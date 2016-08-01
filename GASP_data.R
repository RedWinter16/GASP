library(dplyr)
GASP.data.cleanup<-function(filename, cutoff_reads, replicates, replicates_sites)
{
  #Read in  the file containing Jonah data
  df<-read.csv(filename)
  
  #These vectors create list of samples and help to group the columns by slides, swabs, collect, and site

  #collects
  Site1_collect1<-c("LL1A", "LL1B", "LL1C")
  Site2_collect1<-c("LL2A", "LL2B", "LL2C", "LL2D", "LL2E")
  
  Site1_collect2<-c("LL1D", "LL1E", "LL1F", "LL1G")
  Site2_collect2<-c("LL2F", "LL2G", "LL2H", "LL2I")
  Site3_collect2<- c("LL3A", "LL3B", "LL3C", "LL3D", "LL3E")
  
  Site2_collect23<-c("LL2N", "LL2O", "LL2P")
  Site2_collect3<-c("LL2J", "LL2K", "LL2L", "LL2M")
  
  #controls
  Site2_controls<-c("LL2A", "LL2Q", "LL2R", "LL2S")
  Site3_controls<-c("LL3A")
  
  #swabs
  Site1_swabs<-c("LL1C", "LL1G")
  Site2_swabs<-c("LL2E", "LL2I", "LL2M")
  Site3_swabs<-c("LL3E")
  
  #slides
  Site1_slides<-c("LL1A", "LL1B", "LL1D", "LL1E", "LL1F")
  Site2_slides<-c("LL2A", "LL2B", "LL2C", "LL2D", "LL2F", "LL2G", "LL2H","LL2J", "LL2K", "LL2L","LL2N", "LL2O", "LL2P")
  Site3_slides<-c("LL3A", "LL3B", "LL3C", "LL3D")
 
  #Remove controls 
  #Prompt user if they would like to remove control samples from the data
  rm_controls<-'Y' #readline(prompt="Remove controls? (Y/N) ")
  if (rm_controls == 'Y')
  {
    for(i in 1:length(Site2_controls)){               #there are controls in the site 2 samples
      df<-select(df, -contains(Site2_controls[i]))
      Site2_collect1<-setdiff(Site2_collect1, Site2_controls) #if you remove the controls from the data frame also remove from our sample list
      Site2_slides<-setdiff(Site2_slides, Site2_controls)
    }
    for(i in 1:length(Site3_controls)){               #there are controls in the site 3 samples
      df<- select(df, -contains(Site3_controls[i]))
      Site3_collect2<-setdiff(Site3_collect2, Site3_controls)  #if you remove the controls from the data frame also remove from our sample list
      Site3_slides<-setdiff(Site3_slides, Site3_controls)
    }
  }

  
  #Remove collect 1
  rm_collect1<-'Y' #readline(prompt = "Remove collect 1? (Y/N) ")
  if (rm_collect1 == 'Y')
  {
    for(i in 1:length(Site1_collect1)){
      df<-select(df, -contains(Site1_collect1[i]))
    }
    Site1_swabs<-setdiff(Site1_swabs, Site1_collect1)
    Site1_slides<-setdiff(Site1_slides, Site1_collect1)
    
    for(i in 1:length(Site2_collect1)){
      df<-select(df, -contains(Site2_collect1[i]))
    }
    Site2_swabs<-setdiff(Site2_swabs, Site2_collect1)
    Site2_slides<-setdiff(Site2_slides, Site2_collect1)
  }
  
  #Move the last column over so Meta data is in the first four columns
  df<-df[, c(1:3, ncol(df), 4:(ncol(df)-1))]
  
  #Separate out swabs
  swabsdf<-df[,1:4]
  for (i in 1:length(Site1_swabs)){
    columns<-select(df, contains(Site1_swabs[i]))
    swabsdf<-cbind.data.frame(swabsdf, columns)
  }
  for (i in 1:length(Site2_swabs)){
    columns<-select(df, contains(Site2_swabs[i]))
    swabsdf<-cbind.data.frame(swabsdf, columns)
  }
  for (i in 1:length(Site3_swabs)){
    columns<-select(df, contains(Site3_swabs[i]))
    swabsdf<-cbind.data.frame(swabsdf, columns)
  }

  #Separate out slides
  slidesdf<-df[,1:4]
  for (i in 1:length(Site1_slides)){
    columns<-select(df, contains(Site1_slides[i]))
    slidesdf<-cbind.data.frame(slidesdf, columns)
  }
  for (i in 1:length(Site2_slides)){
    columns<-select(df, contains(Site2_slides[i]))
    slidesdf<-cbind.data.frame(slidesdf, columns)
  }
  for (i in 1:length(Site3_slides)){
    columns<-select(df, contains(Site3_slides[i]))
    slidesdf<-cbind.data.frame(slidesdf, columns)
  }

  #Modify rows that do not have greater than a certain number of reads in 2 of the three replicates for swabs
  #cutoff_reads<-as.numeric(readline(prompt = "Please enter the minimum number of reads for the swabs data: "))
  #replicates<-as.numeric(readline(prompt = paste("Please enter the minimum number of replicates that must have a minimum of", cutoff_reads, " for swabs: ")))
  vec<-seq(5, ncol(swabsdf), 3)
  for (i in 1:length(vec))
    {
    for (j in 1:nrow(swabsdf))
      {
      if (length(which(swabsdf[j,(vec[i]):(vec[i]+2)] >= cutoff_reads)) < replicates)
        {
        swabsdf[j,(vec[i]):(vec[i]+2)]<-0
        }
      }
    }
  
  #Modify rows that do not have >3 reads in 2 of the three replicates for slides
  #cutoff_reads<-as.numeric(readline(prompt = "Please enter the minimum number of reads for the slides data: "))
  #replicates<-as.numeric(readline(prompt = paste("Please enter the minimum number of replicates that must have a minimum of", cutoff_reads, " for slides: ")))
  vec<-seq(5, ncol(slidesdf)-8, 3)
  for (i in 1:length(vec))
  {
    for (j in 1:nrow(slidesdf))
    {
      if (length(which(slidesdf[j,(vec[i]):(vec[i]+2)] >= cutoff_reads)) < replicates)
      {
        slidesdf[j,(vec[i]):(vec[i]+2)]<-0
      }
    }
  }
  #Separte the last 8 columns becuase these are weird replicates of 4 not 3, still need 2/4 to have >= reads
  vec<-seq(ncol(slidesdf)-7, ncol(slidesdf), 4)
  for (i in 1:length(vec))
  {
    for (j in 1:nrow(slidesdf))
    {
      if (length(which(slidesdf[j,(vec[i]):(vec[i]+3)] >=cutoff_reads)) < (replicates+1))
      {
        slidesdf[j,(vec[i]):(vec[i]+3)]<-0
      }
    }
  }
  
  slidesdf<-mutate(slidesdf, Eliminate_if_0 = apply(slidesdf[,5:ncol(slidesdf)],1, sum))
  
  i<-1
  while (i <= nrow(slidesdf))
  {
    if ( slidesdf[i, ncol(slidesdf)] == 0)
    {
      slidesdf<-slidesdf[-i,]
    } else
    {
      i<-i+1
    }
    
  }
  
  slidesdf<-slidesdf[,-ncol(slidesdf)]
  
  swabsdf<-mutate(swabsdf, Eliminate_if_0 = apply(swabsdf[,5:ncol(swabsdf)],1, sum))
  
  i<-1
  while (i <= nrow(swabsdf))
  {
    if ( swabsdf[i, ncol(swabsdf)] == 0)
    {
      swabsdf<-swabsdf[-i,]
    } else
    {
      i<-i+1
    }
    
  }
  
  swabsdf<-swabsdf[,-ncol(swabsdf)]

  #Write out swabs  
  write.csv(swabsdf, file.path("C:/Users/ch26466/Documents/GASP/JonahVentures/CP Second Swing/Swabs/Lab Replicates/", paste0("Swabs_All Primers_Lab_Replicates_", replicates, "_of_3_GT", cutoff_reads,".csv")), row.names = FALSE)  
  
  #Write out slides  
  write.csv(slidesdf, file.path("C:/Users/ch26466/Documents/GASP/JonahVentures/CP Second Swing/Slides/Lab Replicates/", paste0("Slides_All Primers_Lab_Replicates_", replicates, "_of_3_GT", cutoff_reads,".csv")), row.names = FALSE) 
    
  
  #now work on site replicates. Remove data that does not exist in a minimum of the slides for a collect.
  swabsdf<-(sum_LL_sample_reads(swabsdf, 'swab', rm_controls, rm_collect1))
  slidesdf<-(sum_LL_sample_reads(slidesdf, 'slide', rm_controls, rm_collect1))
  
  #replicates_sites<-as.numeric(readline(prompt = paste("Please enter the minimum number of Site replicates_sites that must have a minimum of", cutoff_reads, " for slides: ")))
  ###site1, collect 2
  #swabs
  columns<-colnames(swabsdf)
  sel<-which(match(Site1_collect2, Site1_swabs)  > 0)
  sel2<-Site1_collect2[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(swabsdf))
  {
    if (length(which(swabsdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      swabsdf[i, sel]<-0
    }
  }
  #slides
  columns<-colnames(slidesdf)
  sel<-which(match(Site1_collect2, Site1_slides)  > 0)
  sel2<-Site1_collect2[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(slidesdf))
  {
    if (length(which(slidesdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      slidesdf[i, sel]<-0
    }
  }
  
  ###Site2, collect 2
  #swabs
  columns<-colnames(swabsdf)
  sel<-which(match(Site2_collect2, Site2_swabs)  > 0)
  sel2<-Site2_collect2[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(swabsdf))
  {
    if (length(which(swabsdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      swabsdf[i, sel]<-0
    }
  }
  #slides
  columns<-colnames(slidesdf)
  sel<-which(match(Site2_collect2, Site2_slides)  > 0)
  sel2<-Site2_collect2[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(slidesdf))
  {
    if (length(which(slidesdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      slidesdf[i, sel]<-0
    }
  }
  ###Site2, collect 23
  #slides
  columns<-colnames(slidesdf)
  sel<-which(match(Site2_collect23, Site2_slides)  > 0)
  sel2<-Site2_collect23[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(slidesdf))
  {
    if (length(which(slidesdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      slidesdf[i, sel]<-0
    }
  }
  ###Site2, collect 3
  #swabs
  columns<-colnames(swabsdf)
  sel<-which(match(Site2_collect3, Site2_swabs)  > 0)
  sel2<-Site2_collect3[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(swabsdf))
  {
    if (length(which(swabsdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      swabsdf[i, sel]<-0
    }
  }
  #slides
  columns<-colnames(slidesdf)
  sel<-which(match(Site2_collect3, Site2_slides)  > 0)
  sel2<-Site2_collect3[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(slidesdf))
  {
    if (length(which(slidesdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      slidesdf[i, sel]<-0
    }
  }
  
  ###Site3 collect 2
  #swabs
  columns<-colnames(swabsdf)
  sel<-which(match(Site3_collect2, Site3_swabs)  > 0)
  sel2<-Site3_collect2[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(swabsdf))
  {
    if (length(which(swabsdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      swabsdf[i, sel]<-0
    }
  }
  #slides
  columns<-colnames(slidesdf)
  sel<-which(match(Site3_collect2, Site3_slides)  > 0)
  sel2<-Site3_collect2[sel]
  sel<-which(match(columns, sel2) >0)
  
  for (i in 1:nrow(slidesdf))
  {
    if (length(which(slidesdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
    {
      slidesdf[i, sel]<-0
    }
  }
  
  ###collect 1
  if (rm_collect1 == "N")
  {
    ###Site 1 collect 1
    #swabs
    columns<-colnames(swabsdf)
    sel<-which(match(Site1_collect1, Site1_swabs)  > 0)
    sel2<-Site1_collect1[sel]
    sel<-which(match(columns, sel2) >0)
    
    for (i in 1:nrow(swabsdf))
    {
      if (length(which(swabsdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
      {
        swabsdf[i, sel]<-0
      }
    }
    #slides
    columns<-colnames(slidesdf)
    sel<-which(match(Site1_collect1, Site1_slides)  > 0)
    sel2<-Site1_collect1[sel]
    sel<-which(match(columns, sel2) >0)
    
    for (i in 1:nrow(slidesdf))
    {
      if (length(which(slidesdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
      {
        slidesdf[i, sel]<-0
      }
    }
    
    ###Site 2 collect 1
    #swabs
    columns<-colnames(swabsdf)
    sel<-which(match(Site2_collect1, Site2_swabs)  > 0)
    sel2<-Site2_collect1[sel]
    sel<-which(match(columns, sel2) >0)
    
    for (i in 1:nrow(swabsdf))
    {
      if (length(which(swabsdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
      {
        swabsdf[i, sel]<-0
      }
    }
    #slides
    columns<-colnames(slidesdf)
    sel<-which(match(Site2_collect1, Site2_slides)  > 0)
    sel2<-Site2_collect1[sel]
    sel<-which(match(columns, sel2) >0)
    
    for (i in 1:nrow(slidesdf))
    {
      if (length(which(slidesdf[i, sel] >= cutoff_reads)) < replicates_sites && length(sel) >= 3)
      {
        slidesdf[i, sel]<-0
      }
    }
  }
  
  slidesdf<-mutate(slidesdf, Eliminate_if_0 = apply(slidesdf[,5:ncol(slidesdf)],1, sum))
  
  i<-1
  while (i <= nrow(slidesdf))
  {
    if ( slidesdf[i, ncol(slidesdf)] == 0)
    {
      slidesdf<-slidesdf[-i,]
    } else
    {
      i<-i+1
    }
    
  }
  
  
  
  slidesdf<-slidesdf[,-ncol(slidesdf)]
  
  #Write out swabs  
  #write.csv(swabsdf, file.path("C:/Users/ch26466/Documents/GASP/JonahVentures/CP Second Swing/Swabs/Site Replicates/", paste0("Swabs_All Primers_Site_Replicates_", replicates, "_of_3_GT", cutoff_reads,".csv")), row.names = FALSE)  
  
  #Write out slides  
  write.csv(slidesdf, file.path("C:/Users/ch26466/Documents/GASP/JonahVentures/CP Second Swing/Slides/Site Replicates/", paste0("Slides_All Primers_", "Lab_Replicates_", replicates, "_of_3_Site_Replicates_", replicates_sites, "_of_3_GT", cutoff_reads,".csv")), row.names = FALSE)   
    

}
###########################################################################################################################################################



# Inputs are a data frame of Jonah data split into slides and swabs, the type (slides or swabs), and the answers (Y/N) to rm_controls and rm_collect1

sum_LL_sample_reads<- function(df, type, rm_controls, rm_collect1)
{
  columns_init<-ncol(df)
  if(rm_controls == "N" && type == 'slide')
  {
    df<-mutate(df, LL2A = rowMeans(select(df, contains("2A")))) #slide control
    df<-mutate(df, LL2Q = rowMeans(select(df, contains("2Q")))) #slide control
    df<-mutate(df, LL2R = rowMeans(select(df, contains("2R")))) #slide control
    df<-mutate(df, LL2S = rowMeans(select(df, contains("2S")))) #slide control
    df<-mutate(df, LL3A = rowMeans(select(df, contains("3A")))) #slide control
    
  } else if(type == "slide" && rm_collect1 == 'Y')
  {
    df<-mutate(df, LL1D = rowMeans(select(df, contains("1D")))) #slide
    df<-mutate(df, LL1E = rowMeans(select(df, contains("1E")))) #slide
    df<-mutate(df, LL1F = rowMeans(select(df, contains("1F")))) #slide
    df<-mutate(df, LL2F = rowMeans(select(df, contains("2F")))) #slide
    df<-mutate(df, LL2G = rowMeans(select(df, contains("2G")))) #slide
    df<-mutate(df, LL2H = rowMeans(select(df, contains("2H")))) #slide
    df<-mutate(df, LL2J = rowMeans(select(df, contains("2J")))) #slide
    df<-mutate(df, LL2K = rowMeans(select(df, contains("2K")))) #slide
    df<-mutate(df, LL2L = rowMeans(select(df, contains("2L")))) #slide
    df<-mutate(df, LL2N = rowMeans(select(df, contains("2N")))) #slide
    df<-mutate(df, LL2O = rowMeans(select(df, contains("2O")))) #slide
    df<-mutate(df, LL2P = rowMeans(select(df, contains("2P")))) #slide
    df<-mutate(df, LL3B = rowMeans(select(df, contains("3B")))) #slide
    df<-mutate(df, LL3C = rowMeans(select(df, contains("3C")))) #slide
    df<-mutate(df, LL3D = rowMeans(select(df, contains("3D")))) #slide
    
  } else if(type == "slide" && rm_collect1 == 'N')
  {
    df<-mutate(df, LL1A = rowMeans(select(df, contains("1A")))) #slide collect 1
    df<-mutate(df, LL1B = rowMeans(select(df, contains("1B")))) #slide collect 1
    df<-mutate(df, LL2B = rowMeans(select(df, contains("2B")))) #slide collect 1
    df<-mutate(df, LL2C = rowMeans(select(df, contains("2C")))) #slide collect 1
    df<-mutate(df, LL2D = rowMeans(select(df, contains("2D")))) #slide collect 1
    
    
  } else if(type == 'swab' && rm_collect1 == 'Y')
  {
    df<-mutate(df, LL1G = rowMeans(select(df, contains("1G")))) #swab
    df<-mutate(df, LL2I = rowMeans(select(df, contains("2I")))) #swab
    df<-mutate(df, LL2M = rowMeans(select(df, contains("2M")))) #swab
    df<-mutate(df, LL3E = rowMeans(select(df, contains("3E")))) #swab
    
  } else if(type == 'swab' && rm_collect1 == 'N')
  {
    df<-mutate(df, LL1C = rowMeans(select(df, contains("1C")))) #swab collect 1
    df<-mutate(df, LL2E = rowMeans(select(df, contains("2E")))) #swab collect 1
  }
  
  df<-df[, c(1:4, (columns_init+1):ncol(df))]
  return(df)
}
