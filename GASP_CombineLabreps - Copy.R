# Inputs are a data frame of Jonah data split into slides and swabs, the type (slides or swabs), and the answers (Y/N) to rm_controls and rm_collect1

GASP_CombineLabreps<- function(df, type, rm_controls, rm_collect1)
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

