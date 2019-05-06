#Setwd (specific to me)
setwd("~/GithubRepositories/CS760_project")
#Load in data
d_path <- "./data/DREAM4_InSilico_Size100/insilico_size100_1/insilico_size100_1_timeseries.tsv"
g_path <- "./data/DREAM4_Challenge2_GoldStandards/Size 100/DREAM4_GoldStandard_InSilico_Size100_1.tsv"
d <- read.delim(d_path, comment.char = "#")
g <- read.delim(g_path, header = FALSE, stringsAsFactors = FALSE)


library(dplyr)
g_true_edges <- g %>% filter(V3 == 1)
table(g_true_edges[,2])
# Choose G10 because it has the most number of incoming edges in the gold standard

### must load in bnstruct_analysis_output.RData to get collapse_list[[]]
### This is the model with three time points and three bins.
collapsed_edges <- collapse_list[[4]]

gene_names <- colnames(d)[-1]

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

format_edge <- matrix(edge_names,
                      ncol = 2, byrow = TRUE,
                      dimnames = list(NULL, c("from", "to")))
edge_preds <- data.frame(format_edge)
names(edge_preds) <- c("V1","V2")

# Lasso
library(glmnet)

d_lasso <- d[-1]
d_lasso_y <- as.matrix(d_lasso[10])
d_lasso_x <- as.matrix(d_lasso[-10])

cv_lasso <- cv.glmnet(d_lasso_x, d_lasso_y, alpha= 1)
b <- coef(cv_lasso)

b<- as.matrix(b)
b <- data.frame(b)

write.csv(b, "./code/lasso_coefficients.csv")
