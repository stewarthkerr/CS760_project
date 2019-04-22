#Packages
library(bnstruct)
library(dplyr)
library(reshape)
library(tidyr)

#Setwd (specific to me)
setwd("~/GithubRepositories/CS760_project")

#Load in data
df <- read.delim("./data/DREAM4_InSilico_Size10/insilico_size10_1/insilico_size10_1_timeseries.tsv", comment.char="#")
#df <- read.delim("./data/DREAM4_InSilico_Size100/insilico_size100_1/insilico_size100_1_timeseries.tsv", comment.char="#")

#Get number of timepoints
n_timepoints <- 6

#Get gene names/numbers
gene_names <- colnames(df)[-1]
n_genes <- length(gene_names)


#Format data as 2 timepoint data

#Get integer for each time point
df$int_time <- (df$Time / 50) + 1
#Get unique values for times (there are 21)
int_times <- unique(df$int_time)

#For each time point grab data and data from each n_timepoint. cbind rows and then append to a datafarme
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


#Attempt to learn network
dbn <- learn.dynamic.network(bn_df, num.time.steps = n_timepoints)
show(dbn)
plot(dbn)

###Loop for plots
################################################################################################################################################
################################################################################################################################################

pdf(paste0("bnstructtimepointsplot.pdf"))
for(n in 2:15){
  #Setwd (specific to me)
  setwd("~/GithubRepositories/CS760_project")
  
  #Load in data
  df <- read.delim("./data/DREAM4_InSilico_Size10/insilico_size10_1/insilico_size10_1_timeseries.tsv", comment.char="#")
  #df <- read.delim("./data/DREAM4_InSilico_Size100/insilico_size100_1/insilico_size100_1_timeseries.tsv", comment.char="#")
print(n)
#Get number of timepoints
n_timepoints <- n

#Get gene names/numbers
gene_names <- colnames(df)[-1]
n_genes <- length(gene_names)


#Format data as 2 timepoint data

#Get integer for each time point
df$int_time <- (df$Time / 50) + 1
#Get unique values for times (there are 21)
int_times <- unique(df$int_time)

#For each time point grab data and data from each n_timepoint. cbind rows and then append to a datafarme
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


#Attempt to learn network
dbn <- learn.dynamic.network(bn_df, num.time.steps = n_timepoints)
#show(dbn)
plot(dbn)
title(sub = paste0(n, " # of timepoints"))

}
dev.off()
dev.off()

###Stuart's Code
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