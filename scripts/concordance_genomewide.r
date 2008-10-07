# Examples of use: 
# sample.pp = c(0,0,0,0,0,0,0,.2,.3,.4,.1)
# mat4.5 = dpPosteriorWeights(alpha=1,N=6000,n=10,prob.clade(9,4))
# genome.pp = mat4.5 %*% sample.pp
# distribution.summary(genome.pp)
# xx = 0:6000/6000
# plot(xx, genome.pp,type="h")




UB = function(s){
# Number of unrooted binary trees with s taxa
 s = floor(s)
 if (s<4){return(1)}
 else{    return( prod(2*(4:s)-5))}
}

logUB = function(s){
# Logarithm of the number of unrooted binary trees with s taxa
 s = floor(s)
 if (s<4){ return(0)}
 else{ 
  return( sum(log(2*(4:s)-5)))
 }
}

prob.clade = function(Ntax,Nclade){
 # returns the prior probability of a clade.
 # Ntax = total # of taxa, Nclade= # of taxa in one part of the clade.
 exp(logUB(Nclade+1)+logUB(Ntax-Nclade+1)-logUB(Ntax))
}

lA = function(alpha,n,nonly=TRUE){
# calculates log(A_n(alpha)/A_n(1)) = log( alpha*...*(alpha+n-1) / n! )
# if nonly is TRUE. Otherwise, it returns the aforementioned value
# for all integers between 1 and n.
 if (nonly){
  res=sum(log(1+(alpha-1)/1:n))
 } else{
  res=c(0,cumsum(log(1+(alpha-1)/1:n)))
 }
 return(res)
}

dpPosteriorWeights = function(alpha,N,n,pclade=0,Ntax=0,Nclade=0){
# N = total number of genes in the genome
# n = number of genes in the sample.
# pclade = prior probability of the clade (or feature of interest). If
# not provided, then Ntax and Nclade must be provided.
# Ntax = total number of taxa
# Nclade = # of taxa on one side of the bipartition 
# output: (N+1) x (n+1) matrix of weights.
# The vector of genome-wide posterior probabilities
# P{k genes out of N have the clade|Data} (k varying between 0 and N)
# will then be obtained by multiplying the matrix of weights by the
# vector of sample-wide posterior probabilities
# P{j genes out of n have the clade|Data} (j varying between 0 and n)

 if (pclade == 0 && (Ntax == 0 || Nclade == 0)){
  cat("Ntax or Nclade is missing, as well as the clade probability (pclade)\n")
  return(NULL)
 } else {
  if (pclade ==0){ pclade = prob.clade(Ntax,Nclade)}
  wts1 = sapply(alpha*pclade+(0:n) , lA, n=N-n,nonly=F)
  wts2 = sapply(alpha*(1-pclade)+(0:n), lA, n=N-n,nonly=F)
  wts3 = lA(alpha+n,N-n,nonly=T)
  # these are (N-n+1)x(n+1) matrices
  wts4 = matrix(0,N+1,n+1)
  for (j in 0:n){
   wts4[1+j+0:(N-n),j+1] = exp(wts1[,j+1] + wts2[(N-n):0+1,n-j+1] - wts3)
  }
  return(wts4)
 }
}

genomewide = function(samplewide,alpha,N,n,Ntax,Nclade,conf.level=.95){
# samplewide = vector of posterior probabilities for the
#              sample-wide concordance factors. 
# alpha: prior level of discordance
# N = total number of genes in the genome
# n = number of genes in the sample.
# Ntax = total number of taxa
# Nclade = # of taxa in one side of the bipartition
# output = vector of posterior probabilities for the
#          genome-wide concordance factors. 

 mat = dpPosteriorWeights(alpha,N,n,pclade=0,Ntax,Nclade)
 genomePP = mat %*% samplewide
 plot((0:N)/N,genomePP,type="h",xlab="genome-wide concordance factor",
      ylab="posterior probability")
 distribution.summary(genomePP,conf.level=conf.level)
}


distribution.summary = function(probs, values=0, conf.level=.95){
# returns mean, standard deviation and 95% confidence limits
# of a probability distribution. probs = vector of probabilities.
# values: if 0 (default), the values of the random variable are
# assumed to be 0,1/N,2/N,...,1 where (N+1) is the length of the
# vector "probs" (adapted to concordance factors).
# Otherwise, "values" needs to be a vector of same length as "probs".  
 probs = probs/sum(probs) 
 if (length(values)==1 && values==0) {
  Ngenes = length(probs)-1
  values = 0:Ngenes/Ngenes
 }
 N = length(values);
 a = (1-conf.level)/2;
 u = cumsum(probs);
 Spinf = max(which(u <= a),1)
 Spinf = values[Spinf]
 Spsup = min(which(u >= 1-a),N)
 Spsup = values[Spsup]
 mu = weighted.mean(values,probs);
 va = weighted.mean(values^2,probs) - mu^2;

 return(list(mean=mu,std.dev=sqrt(va), lower=Spinf,upper=Spsup));
}
