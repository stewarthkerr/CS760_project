### Load package, data, and define functions
#####################################################################################################################
#####################################################################################################################
library(bnstruct)
library(dplyr)
library(reshape)
library(tidyr)
#library(graph)
#library(Rgraphviz) required to plot DBN

#Read in the RData for Arab
load("../output/Clustering.RDATA.RData")
#Read in gold standard
gold <- read.delim("../data/Arab.Meristem/arabidopsis.meristem.modules.interactions.tsv", header=FALSE) %>%
  filter(V1 != 0, V2 != 0) %>%
  mutate(V1 = paste("C",V1,sep=""), V2 = paste("C",V2,sep=""))

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
  timepoints <- df[1] %>% unique()
  n_timepoints <- nrow(timepoints)
  
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


###Data input and formatting
#####################################################################################################################
#####################################################################################################################
df_clustered = cbind(net$colors, data_no0) %>%
  filter(`net$colors` > 0) %>%
  group_by(`net$colors`)

#Average across clusters
df_clust_avg = summarise_all(df_clustered, mean)

#1 row per cluster per replication
tempdf1 = select(df_clust_avg, cluster = `net$colors`, contains("_1")) %>%
  rename_at(.vars = vars(ends_with("_1")),
            .funs = funs(sub("_1", "", .)))
tempdf2 = select(df_clust_avg, cluster = `net$colors`, contains("_2")) %>%
  rename_at(.vars = vars(ends_with("_2")),
            .funs = funs(sub("_2", "", .)))
#Grab clusters
clusters <- as.character(pull(tempdf1, cluster))
clusters <- paste("C",clusters,sep="")
#Grab timepoints
n_timepoints <- ncol(tempdf1)-1
#Tranpose data
tempdf1 <- t(tempdf1)[2:11,]
colnames(tempdf1) <- clusters
tempdf1 <- as.data.frame(tempdf1) %>%
  rename(clusters) %>%
  mutate(rep = 1, timepoint = rownames(tempdf1))
tempdf2 <- t(tempdf2)[2:11,]
colnames(tempdf2) <- clusters
tempdf2 <- as.data.frame(tempdf2) %>%
  rename(clusters) %>%
  mutate(rep = 2, timepoint = rownames(tempdf2))

arab_clust_avg = rbind(tempdf1,tempdf2) 

arab_clust_avg_reshaped = data.frame()
col.from = colnames(arab_clust_avg)
col.to1 = paste(colnames(arab_clust_avg),"T1",sep="_") #Needed for timepoint names
col.to2 = paste(colnames(arab_clust_avg),"T2",sep="_") #Needed for timepoint names
for (i in seq(from = 1, to = n_timepoints-1)){
  
  d1 = filter(arab_clust_avg, timepoint == paste("M",i,sep="")) %>%
    rename_at(vars(col.from), ~col.to1)
  d2 = filter(arab_clust_avg, timepoint == paste("M",i+1,sep="")) %>%
    rename_at(vars(col.from), ~col.to2)
  d1_1 = filter(d1, rep_T1 == 1)
  d1_2 = filter(d1, rep_T1 == 2)
  d2_1 = filter(d2, rep_T2 == 1)
  d2_2 = filter(d2, rep_T2 == 2)
  all_comb = rbind(cbind(d1_1,d2_1),
                   cbind(d1_1,d2_2),
                   cbind(d1_2,d2_1),
                   cbind(d1_2,d2_2))
  arab_clust_avg_reshaped = rbind(arab_clust_avg_reshaped,all_comb)
}

arab_clust_avg_reshaped = unite_(arab_clust_avg_reshaped,"dp", c("rep_T1", "timepoint_T1","rep_T2","timepoint_T2"), sep = "_")


###Learn the DBN
#####################################################################################################################
#####################################################################################################################
#Temporary
var_names <- colnames(select(arab_clust_avg_reshaped, -dp))
n_clusters = length(clusters)
n_timepoints = 2

#Creates a BN dataset
bn_df <- BNDataset(data = select(arab_clust_avg_reshaped, -dp),
                   variables = var_names,
                   discreteness = rep('c',n_clusters),
                   num.time.steps = n_timepoints,
                   node.sizes = rep(3,n_clusters))

#Examine the dataset
#show(bn_df)

#Attempt to learn network
dbn <- learn.dynamic.network(bn_df, num.time.steps = 2)
#show(dbn)
#cpts(dbn)
#x = dag(dbn)
#plot(dbn)


###Assessing BN
#####################################################################################################################
#####################################################################################################################
###Collapse Network
#####################################################################################################################
#####################################################################################################################
edges <- dbn@dag
add_elementwise <- function(x) Reduce("+", x)

#Collapse matrix across timepoints
collapse <- list();
k=1
for(i in seq(1,n_clusters*n_timepoints,n_clusters))
{
  for(j in seq(1,n_clusters*n_timepoints,n_clusters))
  {
    collapse[[k]] <- edges[i:(i+n_clusters-1),j:(j+n_clusters-1)]
    k = k+1
  }
}

collapsed_edges <- add_elementwise(collapse)
collapsed_edges <- collapsed_edges > 0

edge_names <- c()
for(row in 1:nrow(collapsed_edges)) {
  for(col in 1:ncol(collapsed_edges)) {
    if (collapsed_edges[row, col]){
      if (clusters[row] != clusters[col]){
        edge_names <- c(edge_names,clusters[row],clusters[col])
      }
    }
  }
}


for(diag in 1:nrow(collapsed_edges)){
  collapsed_edges[diag,diag] = FALSE
}

#if (verbose == TRUE) {
#cond_plot <- graph_from_adjacency_matrix(collapsed_edges,mode = "directed")
#plot(cond_plot)
#}

### Get precision/recall for predicted edges, evaluated on gold standard
#####################################################################################################################
#####################################################################################################################

format_edge <- matrix(edge_names,
                      ncol = 2, byrow = TRUE,
                      dimnames = list(NULL, c("from", "to")))
edge_preds <- data.frame(format_edge)
names(edge_preds) <- c("V1","V2")

n_edges <- nrow(gold)

#Get True Positives
pos <- gold %>% filter(V3 == 1)
tp <- nrow(merge(pos, edge_preds, by = c("V1", "V2")))

#Get False Negatives
fn <- nrow(pos) - tp

#Get False Positives
fp <- nrow(edge_preds) - tp

#Get Precision
precision <- tp / (tp + fp)

#Get Recall
recall <- tp / (tp + fn)

#Get Accuracy
accuracy <- (n_edges - fn - fp) / n_edges
