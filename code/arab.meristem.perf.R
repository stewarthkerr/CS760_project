library(dplyr)

save_file = "arab.meristem.analysis.me.2.RData"
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}
load(save_file)

#Read in gold standard
gold <- read.delim(g_file, header=FALSE) %>%
  filter(V1 != 0, V2 != 0) %>%
  mutate(V1 = sub("ME","C",V1), V2 = sub("ME","C",V2))

if ( expression != "average" ) {
    clusters <- rownames(mes)
}

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

for(row in 1:nrow(collapsed_edges)) {
  for(col in 1:ncol(collapsed_edges)) {
    if (collapsed_edges[row, col]){
        collapsed_edges[col,row] = TRUE #Convert to undirected, as gold is undirected
    }
  }
}

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
