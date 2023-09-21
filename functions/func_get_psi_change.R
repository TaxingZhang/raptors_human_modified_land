get_psi_change<-function(s) {
  load(all_sp[s])
  re <- ranef(fm)
  modes <- colSums(bup(re, stat="mode"))/sampleSize(fm)
  mean_psi_change<-mean(modes[20])-mean(modes[1])
  return(mean_psi_change)
}