cluster <- get(load('../output/Clustering.RDATA.RData'))
colors  <- cluster$colors
colnames(colors) <- c('id','module')
write.csv(colors, "../output/Clustering.csv")
