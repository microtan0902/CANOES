# read in the data
gc <- read.table("gc.txt")$V2
canoes.reads <- read.table("canoes.reads.txt")
# rename the columns of canoes.reads
sample.names <- paste("S", seq(1:26), sep="")
names(canoes.reads) <- c("chromosome", "start", "end", sample.names)
# create a vector of consecutive target ids
target <- seq(1, nrow(canoes.reads))
# combine the data into one data frame
canoes.reads <- cbind(target, gc, canoes.reads)

######################
#cova test
system.time(covariance <- cor(canoes.reads[, sample.names], canoes.reads[, sample.names]))
system.time(cor_fast <- cora(as.matrix(canoes.reads[, sample.names])))
system.time(cor_fast2 <- fastCor(canoes.reads[, sample.names],nSplit = 1, upperTri = FALSE, optBLAS = FALSE, verbose = TRUE))
all.equal(covariance, cor_fast)
all.equal(covariance,cor_fast2)

x <- matrnorm(100, 40)
s1 <- cov(x) 
s2 <- cova(x)
all.equal(s1, s2)
x <- NULL

