#Packages
library(bnstruct)
library(dplyr)
library(reshape)
library(tidyr)
library(bnlearn)

#Load in data
df <- read.delim("./data/DREAM4_InSilico_Size10/insilico_size10_1/insilico_size10_1_timeseries.tsv", comment.char="#")
gold <- read.delim("./data/DREAM4_Challenge2_GoldStandards/Size 10/DREAM4_GoldStandard_InSilico_Size10_1.tsv", comment.char = "#",
                   stringsAsFactors = FALSE)
#df <- read.delim("./data/DREAM4_InSilico_Size100/insilico_size100_1/insilico_size100_1_timeseries.tsv", comment.char="#")

#Get number of timepoints
n_timepoints <- 2

#Get gene names/numbers
gene_names <- colnames(df)[-1]
n_genes <- length(gene_names)

#Format data as 2 timepoint data

#Get integer for each time point
df$int_time <- (df$Time / 50) + 1
#Get unique values for times (there are 21)
int_times <- unique(df$int_time)

#For each time point grab data and data from each n_timepoint. 
#Cbind rows and then append to a datafarme with rbind
final_df = data.frame()

for (t in int_times) {
  
  if (t <= max(int_times)-(n_timepoints-1)) {
    
    old_names <- colnames(df %>% select(-c(Time, int_time)))
    new_names <- paste(old_names,paste0(".",1), sep='')
    data_t <- df %>% filter(int_time == t) %>% select(-c(Time, int_time)) %>% rename_at(vars(old_names), function(x) new_names)
    
    for (i in 1:(n_timepoints-1)) {
      # append .i to i+1 column names
      old_names <- colnames(df %>% select(-c(Time, int_time)))
      new_names <- paste(old_names,paste0(".",i+1), sep='')
      
      data_t_plus <- df %>% filter(int_time == t + i) %>% 
        select(-c(Time,int_time)) %>% rename_at(vars(old_names), function(x) new_names)
      
      data_t <- cbind(data_t, data_t_plus)
      
    }
    
    final_df <- rbind(final_df, data_t)
    
  }
}

#Creates a BN dataset
bn_df <- BNDataset(data = final_df,
                   variables = colnames(final_df),
                   discreteness = rep(rep('c',n_genes),n_timepoints),
                   num.time.steps = n_timepoints,
                   node.sizes = rep(rep(2,n_genes),n_timepoints))

#Creates a BN object with pre-trained structure and learn the parameter
dbn <- BN(dataset = bn_df)
struct <- read.csv(file="./code/GlobalMIT/structure.csv", header=FALSE, sep=",")
dbn@dag[1:10, 11:20] = as.matrix(struct)
dbn <- learn.params(dbn, bn_df, ess=1)
show(dbn)
plot(dbn)

