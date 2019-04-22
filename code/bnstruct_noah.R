# packages
library(bnstruct)
library(dplyr)
library(reshape)
library(tidyr)

# setwd (specific to me)
setwd("~/GithubRepositories/CS760_project")

# Load in data
df <- read.delim("./data/DREAM4_InSilico_Size10/insilico_size10_1/insilico_size10_1_timeseries.tsv", comment.char="#")

# Get number of timepoints
timepoints <- df[1] %>% unique()
n_timepoints <- 2

# Get gene names/numbers
gene_names <- colnames(df)[-1]
n_genes <- length(gene_names)


# Format data as 2 timepoint data

# get integer for each time point
df$int_time <- (df$Time / 50) + 1
# get unique values for times (there are 21)
int_times <- unique(df$int_time)

# for each time point grab data and data from t+1.  cbind rows and then append to a datafarme
final_df = data.frame()

for (t in int_times) {

  if (t < max(int_times)) {
    
    time_pair <- cbind(df %>% filter(int_time == t) %>% select(-c(Time, int_time)), df %>% filter(int_time == t + 1) %>% select(-c(Time,int_time)))
  
    final_df <- rbind(final_df, time_pair)
    
  }
}

#Creates a BN dataset
bn_df <- BNDataset(data = final_df,
                   variables = rep(gene_names,n_timepoints),
                   discreteness = rep(rep('c',n_genes),n_timepoints),
                   num.time.steps = n_timepoints,
                   node.sizes = rep(rep(2,n_genes),n_timepoints))

#Examine the dataset
#show(bn_df)

#Attempt to learn network
dbn <- learn.dynamic.network(bn_df, num.time.steps = n_timepoints)
show(dbn)

################################################################################################################################################
################################################################################################################################################

library(bnstruct)
library(dplyr)
library(reshape)
library(tidyr)
df <- read.delim("../data/DREAM4_InSilico_Size10/insilico_size10_1/insilico_size10_1_timeseries.tsv", comment.char="#")

#Get number of timepoints
timepoints <- df[1] %>% unique()
n_timepoints <- nrow(timepoints)

#Get gene names/numbers
gene_names <- colnames(df)[-1]
n_genes <- length(gene_names)

#Tranpose the data for BNDataset
df_unstack <- function(df,start){
  #Height of the stacks
  timepoints <- df[1] %>% unique()
  n_timepoints <- nrow(timepoints)
  
  #Pulls the first n_timepoints rows
  unstack_df <- df[start:(start+n_timepoints-1),]
  
  return(unstack_df)
}

df_widen <- function(unstack_df){
  #Transposes the dataset to get it in the right form for BN
  wide_unstack_df <- data.frame()
  for (tp in 1:n_timepoints){
    wide_unstack_df = c(wide_unstack_df, unstack_df[tp,-1])
  }
  
  return(as.data.frame(wide_unstack_df))
}

wide_unstack_df_combine <- function(df){
  timepoints <- df[1] %>% unique()
  n_timepoints = nrow(timepoints)
  n_obs = nrow(df)/n_timepoints
  
  final_df = data.frame()
  for (i in 1:n_obs){
    start = ((i-1)*n_timepoints)+1
    unstack_df = df_unstack(df,start)
    wide_unstack_df = df_widen(unstack_df)
    final_df = rbind(final_df, wide_unstack_df)
  }
  
  return(final_df)
}

bn_ready_df = wide_unstack_df_combine(df)

#Creates a BN dataset
bn_df <- BNDataset(data = bn_ready_df,
                   variables = rep(gene_names,n_timepoints),
                   discreteness = rep(rep('c',n_genes),n_timepoints),
                   num.time.steps = n_timepoints,
                   node.sizes = rep(rep(2,n_genes),n_timepoints))

#Examine the dataset
#show(bn_df)

#Attempt to learn network
dbn <- learn.dynamic.network(bn_df, num.time.steps = n_timepoints)
show(dbn)