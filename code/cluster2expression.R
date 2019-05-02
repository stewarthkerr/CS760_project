library(tidyr)

seed         <- 7331
num_reps     <- 20
cluster_file <- "../output/Clustering.RDATA.RData"
output_file  <- "../data/Arab.Meristem/arabidopsis.meristem.modules.expression.tsv"
output       <- file(output_file, open="wt")

args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}

printOutput <- function(MEs, s_rep) {
    rMEs <- MEs[s_rep,]
    write("",output, append=TRUE) #Prefixed newline
    for( i in 1:nrow(rMEs) ) {
        expression_str <- paste(c(s_rep[i],rMEs[i,]), collapse='\t')
        write(expression_str,output, append=TRUE)
    }
}

cluster <- get(load(cluster_file))
MEs     <- cluster$MEs
nbits   <- nrow(MEs)/2
sample_range <- 1:(2^nbits)
set.seed(seed) #To have reproducible behavior
samples <- sample(sample_range, num_reps, replace = FALSE)
samples_bits <- sapply(samples, function(s) { as.integer(intToBits(s)) })[1:nbits,] + 1 #Add 1 for indexing with R
lut   <- matrix(rownames(MEs), 2, nbits)
header <- paste(sprintf("\"%s\"",c("Development",colnames(MEs))),collapse="\t")
write(header,output)
for( i in 1:ncol(samples_bits) )
{
    s <- samples_bits[,i]
    #browser()
    s_rep <- diag(lut[s,])
    printOutput(MEs, s_rep)
}
close(output)
