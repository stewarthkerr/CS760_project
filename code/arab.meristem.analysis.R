### Load package, data, and define functions
#####################################################################################################################
#####################################################################################################################
library(bnstruct)
library(dplyr)
library(gtools)
library(tidyr)
#library(graph)
#library(Rgraphviz) required to plot DBN

expression = "average"
n_dbn_timepoints = 2
n_reps = 2
d_file="../output/Clustering.RDATA.RData"
g_file="../data/Arab.Meristem/arabidopsis.meristem.modules.interactions.t0.02.tsv"

args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}

save_file = sprintf("../output/arab.meristem.analysis.%s.%d.RData", expression, n_dbn_timepoints)
#Read in the RData for Arab
load(d_file)
#Read in gold standard
gold <- read.delim(g_file, header=FALSE) %>%
  filter(V1 != 0, V2 != 0) %>%
  mutate(V1 = sub("ME","C",V1), V2 = sub("ME","C",V2))

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
if( expression == "average" ) {
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
    mutate(rep = 1, timepoint = rownames(tempdf1))
    tempdf2 <- t(tempdf2)[2:11,]
    colnames(tempdf2) <- clusters
    tempdf2 <- as.data.frame(tempdf2) %>%
    mutate(rep = 2, timepoint = rownames(tempdf2))

    arab_clust_avg = rbind(tempdf1,tempdf2) 
    n_clusters = length(clusters)
} else { 
    #1 row per cluster per replication
    mes <- data.frame(t(net$MEs))
    mes <- mes[2:nrow(mes),]
    clusters = rownames(mes)
    n_clusters = nrow(mes)
    rownames(mes) = sub('ME','',rownames(mes))
    mes <- mes[order(as.numeric(rownames(mes))),]
    rownames(mes) = sprintf('C%s', rownames(mes))
    tempdf1 = select(mes, contains("_1")) %>%
        rename_at(.vars = vars(ends_with("_1")),
                  .funs = funs(sub("_1", "", .)))
    tempdf2 = select(mes, contains("_2")) %>%
        rename_at(.vars = vars(ends_with("_2")),
                  .funs = funs(sub("_2", "", .)))
    #Grab timepoints
    n_timepoints <- ncol(tempdf1)
    tempdf1 <- as.data.frame(t(tempdf1))
    tempdf1 <- mutate(tempdf1, rep = 1, timepoint = rownames(tempdf1))
    tempdf2 <- as.data.frame(t(tempdf2))
    tempdf2 <- mutate(tempdf2, rep = 2, timepoint = rownames(tempdf2))

    arab_clust_avg = rbind(tempdf1,tempdf2) 
}

arab_clust_avg_reshaped = data.frame()
col.from = colnames(arab_clust_avg)
for (i in 1:n_dbn_timepoints) {
    eval(parse(text=sprintf("col.to%d = paste(colnames(arab_clust_avg),\"T%d\",sep=\"_\")", i, i)))
}
#Generate rep suffixes for permutations
suffixes = c()
for( k in 1:n_reps ) {
    suffixes = cbind(suffixes, sprintf("%d", k))
}
perms = permutations(n=n_reps,r=n_dbn_timepoints,v=suffixes,repeats.allowed=TRUE)
for( k in 1:n_dbn_timepoints ) {
    perms[,k] = sprintf("d%d_%s", k, perms[,k])
}
perms <- as.data.frame(perms)
mutate_str <- sprintf("all_comb_str <- perms %%>%% mutate(bind_str=sprintf(\"cbind(%s)\",%s))", paste0(rep("%s",n_dbn_timepoints), collapse=','), paste0(colnames(perms), collapse=','))
eval(parse(text=(mutate_str)))
rbind_str <- sprintf("rbind(%s)",paste0(all_comb_str$bind_str,collapse=','))
    

for (i in seq(from = 1, to = n_timepoints-(n_dbn_timepoints-1))){
  for (j in 1:n_dbn_timepoints) {
      eval(parse(text=sprintf("d%d = filter(as.data.frame(arab_clust_avg), timepoint == paste(\"M\",i+%d,sep=\"\")) %%>%% rename_at(vars(col.from), ~col.to%d)", j,j-1,j)))
      for( k in 1:n_reps ) {
        eval(parse(text=sprintf("d%d_%d = filter(d%d, rep_T%d == %d)", j, k, j, j, k)))
      }
  }
  eval(parse(text=sprintf("all_comb <- %s", rbind_str)))
  #all_comb = c()
  arab_clust_avg_reshaped = rbind(arab_clust_avg_reshaped,all_comb)
}

rep_str = paste0(sprintf("\"rep_T%d\",\"timepoint_T%d\"", 1:n_dbn_timepoints, 1:n_dbn_timepoints),collapse=',')
eval(parse(text=sprintf("arab_clust_avg_reshaped = unite_(arab_clust_avg_reshaped,\"dp\", c(%s), sep = \"_\")", rep_str)))


###Learn the DBN
#####################################################################################################################
#####################################################################################################################
#Temporary
var_names <- colnames(select(arab_clust_avg_reshaped, -dp))

#Creates a BN dataset
bn_df <- BNDataset(data = select(arab_clust_avg_reshaped, -dp),
                   variables = var_names,
                   discreteness = rep('c',n_clusters),
                   num.time.steps = n_dbn_timepoints,
                   node.sizes = rep(3,n_clusters))

#Examine the dataset
#show(bn_df)

#Attempt to learn network
save.image(file=save_file)
dbn <- learn.dynamic.network(bn_df, num.time.steps = n_dbn_timepoints)
save.image(file=save_file)
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
for(i in seq(1,n_clusters*n_dbn_timepoints,n_clusters))
{
  for(j in seq(1,n_clusters*n_dbn_timepoints,n_clusters))
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
save.image(file=save_file)
