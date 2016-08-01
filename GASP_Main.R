#****The main script. This will take the input Jonah data file, munge it into a table of relative abundance split by slides and swabs samples/sites/collects.
#****This is heavily dependent on the format of the input table, especially the column headers and the order of columns.

GASP_Main<-function(filename, minimums, replicates)
{
  ####Set up workspace using function call and return data directory
  source('~/GASP/R/GASP_SetWorkspace.R')
  data_dir<-GASP_SetWorkspace()
  
  
  ####Call function to convert Jonah data to realtive abundance and import the file into workspace
  relativeAB<-GASP_getdata(filename, data_dir)
  
  
  ############################################################################################################################################################
  #These vectors create list of samples and help to group the columns. THIS controls every part of analysis downstream.
  
  #controls
  Controls<-c("FLA", "MAA", "MAQ", "MAR", "MAS")
  
  #collect 1
  
  Collect1<-c("SCA", "SCB", "SCC", "MAA", "MAB", "MAC", "MAD", "MAE")
  
  #collects
  Site1_collect1<-c("SCA", "SCB", "SCC")
  Site2_collect1<-c("MAA", "MAB", "MAC", "MAD", "MAE")
  
  Site1_collect2<-c("SCD", "SCE", "SCF", "SCG")
  Site2_collect2<-c("MAF", "MAG", "MAH", "MAI")
  Site3_collect2<- c("FLA", "FLB", "FLC", "FLD", "FLE")
  
  Site2_collect23<-c("MAN", "MAO", "MAP")
  Site2_collect3<-c("MAJ", "MAK", "MAL", "MAM")
  
  #swabs
  Site1_swabs<-c("SCC", "SCG")
  Site2_swabs<-c("MAE", "MAI", "MAM")
  Site3_swabs<-c("FLE")
  
  #slides
  Site1_slides<-c("SCA", "SCB", "SCD", "SCE", "SCF")
  Site2_slides<-c("MAA", "MAB", "MAC", "MAD", "MAF", "MAG", "MAH","MAJ", "MAK", "MAL","MAN", "MAO", "MAP")
  Site3_slides<-c("FLA", "FLB", "FLC", "FLD")     
  
  ############################################################################################################################################################
  
  
  
  ####Remove Data from anaylsis 
  
  #Prompt user if they would like to remove control samples from the data, or comment out if running several iterations with the same settings
  rm_controls<-'Y'#readline(prompt="Remove controls? (Y/N) ")
      if (rm_controls == 'Y')
      {
            for(i in 1:length(Controls))
            {             
                  relativeAB<-select(relativeAB, -(contains(Controls[[i]])))
            }
            #if you remove the controls from the data frame also remove from our sample lists
            #collects
            Collect1<-setdiff(Collect1, Controls)
            Site1_collect1<-setdiff(Site1_collect1, Controls)
            Site2_collect1<-setdiff(Site2_collect1, Controls)
            Site1_collect2<-setdiff(Site1_collect2, Controls)
            Site2_collect2<-setdiff(Site2_collect2, Controls)
            Site3_collect2<-setdiff(Site3_collect2, Controls)
            Site2_collect23<-setdiff(Site2_collect23, Controls)
            Site2_collect3<-setdiff(Site2_collect3, Controls)
            #swabs
            Site1_swabs<-setdiff(Site1_swabs, Controls)
            Site2_swabs<-setdiff(Site2_swabs, Controls)
            Site3_swabs<-setdiff(Site3_swabs, Controls)
            #slides
            Site1_slides<-setdiff(Site1_slides, Controls)
            Site2_slides<-setdiff(Site2_slides, Controls)
            Site3_slides<-setdiff(Site3_slides, Controls)
      }
  
  ####Remove Collect 1 THIS IS ONLY FOR FIRST SET OF DATA since there is only 1 collect 1 and it was not performed correctly
  #Prompt user if they would like to remove control samples from the data, or comment out if running several iterations with the same settings
  rm_collect1<-"Y"#readline(prompt = "Remove collect 1? (Y/N) ")
      if (rm_collect1 == 'Y')
      {
            for(i in 1:length(Collect1))
            {   
              relativeAB<-select(relativeAB, -(contains(Collect1[[i]])))
            }
            #if you remove the Collect1 from the data frame also remove from our sample lists too
            #collects
            Site1_collect1<-NULL
            Site2_collect1<-NULL
            #swabs
            Site1_swabs<-setdiff(Site1_swabs, Collect1)
            Site2_swabs<-setdiff(Site2_swabs, Collect1)
            Site3_swabs<-setdiff(Site3_swabs, Collect1)
            #slides
            Site1_slides<-setdiff(Site1_slides, Collect1)
            Site2_slides<-setdiff(Site2_slides, Collect1)
            Site3_slides<-setdiff(Site3_slides, Collect1)
      }
  
  
  
  ####Separate swabs and slides into separate tables, the data is different enough that they should be separated.
  
  #start by create vectors of each, do all removing of data before this step or you will need to re-edit.
  
  slides<-c(Site1_slides,Site2_slides,Site3_slides)
  swabs<-c(Site1_swabs,Site2_swabs,Site3_swabs)
  
  #Separate out swabs from rest of the data
  swabsdf<-relativeAB[,1:4]
      for (i in 1:length(swabs))
        {
        column<-select(relativeAB, contains(swabs[[i]]))
        swabsdf<-cbind.data.frame(swabsdf, column)
        }
  #Separate out slides
  slidesdf<-relativeAB[,1:4]
      for (i in 1:length(slides))
        {
        column<-select(relativeAB, contains(slides[[i]]))
        slidesdf<-cbind.data.frame(slidesdf, column)
        }

  rm(relativeAB)
  
  
  ####Check Lab Replicates
  
  #This piece of code will need to be modified. Call the function GASP_CheckLabreps but change the data fed to it and the index to match the number of 
  #lab replicates in the sample. For example, if columns 8+9 are 10 lab replicates feed the program df[,c(1:4,8,9)] and index = 10.
  # this is also set up to loop through a number of iterations of the minimum and repliactes
  print(Sys.time())
      for(i in 1:length(minimums))
      {
            for(j in 1:length(replicates))
            {
              index<-3
              rep<-as.integer(replicates[[j]])
              
              #swabs
              tempdf<-GASP_CheckLabreps(swabsdf, minimums[[i]], replicates[[j]], index)
              tempdf<-GASP_RMzeroTaxa(tempdf)
              newfile<-file.path(data_dir, paste("Swabs Lab Replicates ", replicates[[j]],"of", index, "GT", minimums[[i]],".csv"))
              #####write.csv(tempdf, newfile, row.names = FALSE)
              #GASP_bison(tempdf, newfile)
              
              #Begin by calling function to take the average relative abundace across lab replicates
              #The combinelabreps relies on an index similar to the checklabreps function
              tempdf<-GASP_CombineLabreps(tempdf, index)
              tempdf<-GASP_RMzeroTaxa(tempdf)
              
              GASP_Output(tempdf, newfile)
              print(Sys.time())
              
              #######################################################################################               
              #######################################################################################               
              ####Check Site Replicates##############################################################
              
              #!!!!this entire piece of code will also need modification or elimination.
        
              #Begin by calling function to take the average relative abundace across lab replicates
              #The combinelabreps relies on an index similar to the checklabreps function
              #tempdf<-GASP_CombineLabreps(tempdf, index)
              #Now compare site replicates, using the checksitereps function
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site1_collect1,  Site1_swabs)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site2_collect1,  Site2_swabs)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site1_collect2,  Site1_swabs)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site2_collect2,  Site2_swabs)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site3_collect2,  Site3_swabs)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site2_collect23,  Site2_swabs)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site2_collect3,  Site2_swabs)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site3_collect2,  Site3_swabs)
              
              #####Remove taxa rows that become zero across all columns, these are now empty
              #tempdf<-GASP_RMzeroTaxa(tempdf)
              
              #Save file
              #newfile<-file.path(data_dir, paste("Lab Replicates Swabs", replicates[[j]],"of", index, "GT", minimums[[i]], "Site Replicates", replicates[[j]],"of", index, ".csv"))
              #write.csv(tempdf, newfile, row.names = FALSE)
              #GASP_bison(tempdf, newfile)
              #rm(tempdf)
              #######################################################################################                                 
              #######################################################################################           
              #######################################################################################               
              
              
              #slides (the last 8 columns have a different number of replicates so we must treat them separetly)
              tempdf_A<-GASP_CheckLabreps(slidesdf[,1:(ncol(slidesdf)-8)], minimums[[i]], replicates[[j]], index)

              index<-4
              #the last columns of slides
              tempdf_B<-GASP_CheckLabreps(slidesdf[, c(1:4, (ncol(slidesdf)-7):ncol(slidesdf))], minimums[[i]], replicates[[j]], index)
              
              #recombine the 2 parts of slides df, this is the most annoying part but there is nothing else I can do
              tempdf<-cbind.data.frame(tempdf_A, tempdf_B[,-c(1:4)])
              tempdf<-GASP_RMzeroTaxa(tempdf)
      
              index<-3
              newfile<-file.path(data_dir, paste("Slides Lab Replicates ", replicates[[j]],"of", index, "GT", minimums[[i]],".csv"))
              ####write.csv(tempdf, newfile, row.names = FALSE)
              #GASP_bison(tempdf, newfile)
              
              #The combinelabreps relies on an index similar to the checklabreps function
              index<-3
              tempdf_A<-GASP_CombineLabreps(tempdf_A, index)
              
              index<-4
              tempdf_B<-GASP_CombineLabreps(tempdf_B, index)
              
              
              #recombine the 2 parts of slides df, this is the most annoying part but there is nothing else I can do
              tempdf<-cbind.data.frame(tempdf_A, tempdf_B[,-c(1:4)])
              tempdf<-GASP_RMzeroTaxa(tempdf)
              
              
              rm(tempdf_A)
              rm(tempdf_B)
              
              GASP_Output(tempdf, newfile)
              print(Sys.time())
              #######################################################################################               
              #######################################################################################               
              ####Check Site Replicates##############################################################
              
              #!!!!this entire piece of code will also need modification or elimination.
              
              #call function to take the average relative abundace across lab replicates
              #The combinelabreps relies on an index similar to the checklabreps function
              #index<-3
              #tempdf_A<-GASP_CombineLabreps(tempdf_A, index)
              
              #index<-4
              #tempdf_B<-GASP_CombineLabreps(tempdf_B, index)
              
              
              #recombine the 2 parts of slides df, this is the most annoying part but there is nothing else I can do
              #tempdf<-cbind.data.frame(tempdf_A, tempdf_B[,-c(1:4)])
              #tempdf<-GASP_RMzeroTaxa(tempdf)
              
              
              #rm(tempdf_A)
              #rm(tempdf_B)
              
              #Now compare site replicates, using the checksitereps function
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site1_collect1,  Site1_slides)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site2_collect1,  Site2_slides)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site1_collect2,  Site1_slides)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site2_collect2,  Site2_slides)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site3_collect2,  Site3_slides)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site2_collect23,  Site2_slides)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site2_collect3,  Site2_slides)
              #tempdf<-GASP_CheckSitereps(tempdf, rep,  Site3_collect2,  Site3_slides)
              
              #####Remove taxa rows that become zero across all columns, these are now empty
              #tempdf<-GASP_RMzeroTaxa(tempdf)
              
              #index<-3
              #Save file
              #newfile<-file.path(data_dir, paste("Lab Replicates Slides", replicates[[j]],"of", index, "GT", minimums[[i]], "Site Replicates", replicates[[j]],"of", index, ".csv"))
              #write.csv(tempdf, newfile, row.names = FALSE)
              #GASP_bison(tempdf, newfile)
              rm(tempdf)
              #######################################################################################                                 
              #######################################################################################           
              #######################################################################################            
              
            }
      }
  GASP_GeoMap(data_dir)
}
