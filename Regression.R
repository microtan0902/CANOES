source("CANOES.R")

gc <- read.table("gc.txt")$V2
canoes.reads <- read.table("canoes.reads.txt")
# rename the columns of canoes.reads
sample.names <- paste("S", seq(1:26), sep="")
names(canoes.reads) <- c("chromosome", "start", "end", sample.names)
# create a vector of consecutive target ids
target <- seq(1, nrow(canoes.reads))
# combine the data into one data frame
canoes.reads <- cbind(target, gc, canoes.reads)

head(canoes.reads[,'S1'])
case_sample <- canoes.reads[,'S1']
mean(case_sample)
gc
plot(gc,case_sample)
mean_rc <- mean(case_sample)

# We assume that RD of a given sample follows negative binomial, 
# then if we want to calcuate the likelihood of the RD give CN, we need to get the parameters of the NB model.
# the NB model has two parameters, mu and sigma, where 他们受两个因素的影响，即平均RD和GC含量。 

#D ~ Negative binomial (mu, s)
#mu ~ f (GC, mean depth)
#s ~ f (GC, mean depth)
library(mgcv)
fit <- gam(mu ~ s(mean_rc) + s(gc), family=Gamma(link=log), data=case_sample)
# we don't want variance less than Poisson
# we take maximum of genome-wide estimate, method of moments estimate
# and Poisson variance
v.estimate <- pmax(predict(fit, counts, type="response"), counts$var, 
                   counts$mean * 1.01)