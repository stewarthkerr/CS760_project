#Packages
library(bnstruct)
library(dplyr)
library(reshape)
library(tidyr)
library(bnlearn)
library(igraph)

#Setwd (specific to me)
setwd("~/GithubRepositories/CS760_project")
#Load in data
d_path <- "./data/DREAM4_InSilico_Size100/insilico_size100_1/insilico_size100_1_timeseries.tsv"
g_path <- "./data/DREAM4_Challenge2_GoldStandards/Size 100/DREAM4_GoldStandard_InSilico_Size100_1.tsv"
d <- read.delim(d_path, comment.char = "#")
g <- read.delim(g_path, header = FALSE, stringsAsFactors = FALSE)

#### connect to server

#library(remoter)
#remoter::client(addr = "DESKTOP-M7AFUVD", port = 55655, password = "noahtaft")
#c2s(d)
#c2s(g)

bnstruct_analysis <- function(df, gold, n_timepoints, 
                              node_size = 2, scoring_function = "BDeu", fitting_algorithm = "mmhc",
                              verbose = FALSE) {
  

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
                   node.sizes = rep(rep(node_size,n_genes),n_timepoints))


#Attempt to learn network
dbn <- learn.dynamic.network(bn_df, num.time.steps = n_timepoints,
                             scoring.func = scoring_function, algo = fitting_algorithm)
if (verbose == TRUE){
show(dbn)
plot(dbn)
}

###Collapse Network
#####################################################################################################################
#####################################################################################################################
edges <- dbn@dag
add_elementwise <- function(x) Reduce("+", x)

#Collapse matrix across timepoints
collapse <- list();
k=1
for(i in seq(1,n_genes*n_timepoints,n_genes))
{
  for(j in seq(1,n_genes*n_timepoints,n_genes))
  {
    collapse[[k]] <- edges[i:(i+n_genes-1),j:(j+n_genes-1)]
    k = k+1
  }
}

collapsed_edges <- add_elementwise(collapse)
collapsed_edges <- collapsed_edges > 0

edge_names <- c()
for(row in 1:nrow(collapsed_edges)) {
  for(col in 1:ncol(collapsed_edges)) {
    if (collapsed_edges[row, col]){
      if (gene_names[row] != gene_names[col]){
        edge_names <- c(edge_names,gene_names[row],gene_names[col])
      }
    }
  }
}


for(diag in 1:nrow(collapsed_edges)){
  collapsed_edges[diag,diag] = FALSE
}

cond_plot <- graph_from_adjacency_matrix(collapsed_edges,mode = "directed")
if (verbose == TRUE) {
plot(cond_plot)
}

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

output <- list()

output[[1]] = c(precision, recall, accuracy)
output[[2]] = dbn
output[[3]] = collapsed_edges
output[[4]] = cond_plot

return(output)

}

#############################################################################
bnstruct_analysis(d, g, n_timepoints = 2, node_size = 2, verbose = FALSE)

acc_list = list()
dbn_list = list()
collapse_list = list()
cond_list = list()

k = 1
for (n_time in 2:4){
  for (n_size in 2:3){
    print(paste(k, n_time, n_size, sep = " "))
    out = bnstruct_analysis(d, g, n_timepoints = n_time, node_size = n_size,
                                        verbose = FALSE)
    acc_list[[k]] = out[[1]]
    dbn_list[[k]] = out[[2]]
    collapse_list[[k]] = out[[3]]
    cond_list[[k]] = out[[4]]
    k = k + 1
  }}

#s2c(acc_list)
