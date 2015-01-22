#this is a function used in Lovell and McKay (2015) that computes bootstrapped estimates
#of the basic stats function available in the Hierfstat package
#it takes four arguments
#####
#fs is a fstat input object that can be created from adegenet via genind2hierfstat
#pop is the grouping variable used to calculate stats
#nboot is the number of bootstrap replicates requested
#rarify is the population size to rarify to. Must be larger than the smallest population size
#####
boot.basicstats<-function(fs,pop="pop",nboot=10,rarify=5){
  rare_in<-split(fs,as.character(fs[,pop]))
  fs.rare<-ldply(lapply(rare_in,function(x) x[sample(rownames(x),size=rarify,replace=F),]),data.frame)[,-1]
  bs.out<-data.frame()
  while(length(rownames(bs.out))<nboot){
    boot.in<-fs.rare[sample(1:length(fs.rare[,1]), replace=T),]
    boot.order<-boot.in[order(boot.in$pop),]
    bs<-basic.stats(boot.order)
    bs.out<-rbind(bs.out,bs$overall)
  }
  colnames(bs.out)<-names(bs$overall)
  return(bs.out)
}
