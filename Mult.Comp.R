Mult.Comp<-function(filename, vec_cutoff, vec_labreps, vec_sitereps)
{
    for (i in 1:length(vec_cutoff))
    {   
      for(j in 1:length(vec_labreps))
      {
        for (k in 1:length(vec_sitereps))
        {
          GASP.data.cleanup(filename, vec_cutoff[i], vec_labreps[j], vec_sitereps[k])
        }
      }
    }
}