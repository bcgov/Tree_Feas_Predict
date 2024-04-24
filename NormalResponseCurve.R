
##########Create Normal Distribution Graphic

#allocate list size to store means
meanOfSampleMeansVector <- numeric(1000)
#for 1000 iterations create 40 exponential random variable with variance of 0.2 units
for (i in 1:1000 ){ 
  sample <- rexp(n=40,0.2) 
  #get mean of sample
  meanOfSample <- mean(sample) 
  #set the mean in list 
  meanOfSampleMeansVector[i] <- meanOfSample
}
propDensity=dnorm(meanOfSampleMeansVector,mean(meanOfSampleMeansVector),sd(meanOfSampleMeansVector))
require(ggplot2)
qplot(meanOfSampleMeansVector,propDensity,geom="line")+
  xlab("Environmental Gradient")+ylab("Density of Tree Species Occurence (Suitability)")
  #ggtitle("Sample Means of Exponential Distribution")
