library (ape)
#library(glmnet)

##read data and make the pairwise distance table for multiple gene trees
#read.tree("6148treefile.pro.nwk") -> trees
species <- "bat"
read.tree(paste(species,".nuc.nwk",sep="")) -> trees
read.table("6148seq_length.txt") -> seq_length
#seq_length <- seq_length/3 #for protein
trees_all <- names(trees)
length_cut <- 10 ##paramaters
pairwise_dist <- c()
species_no <- 3#length(species)
count <- 0
mean_branch_length <- c()
media_branch_length <- c()
#trun_nu <- 5

for (i in trees_all){
  cophenetic.phylo(trees[[i]]) -> tree_dis
  if ((dim(tree_dis)[1] >= species_no)&(seq_length[i,] >= length_cut)){
    count <- count + 1
    edge_length <- trees[[i]]$edge.length
    tree_edge <- sort(trees[[i]]$edge.length)
    #for (t in 1:trun_nu){
    #  tree_edge <- tree_edge[-1]
    #}
    mean(tree_edge) -> mean_length
    median(trees[[i]]$edge.length) -> media_length
    gene_mean_branch_length <- cbind(i,mean_length)
    gene_media_branch_length <- cbind(i, media_length)
    mean_branch_length <- rbind(mean_branch_length, gene_mean_branch_length)
    media_branch_length <- rbind(media_branch_length, gene_media_branch_length)
  }
}

head(mean_branch_length);dim(mean_branch_length)
mean_branch_length <- data.frame(mean_branch_length)
colnames(mean_branch_length) <- c("0gene", "mean_branch_length_nuc")
write.table(mean_branch_length, file = paste(species,"_mean_branch.nuc.txt",sep=""), quote = F, row.names = F, sep = "\t")

head(media_branch_length);dim(media_branch_length)
media_branch_length <- data.frame(media_branch_length)
colnames(media_branch_length) <- c("0gene", "media_branch_length_nuc")
write.table(media_branch_length, file = paste(species,"_median_branch.nuc.txt",sep=""), quote = F, row.names = F, sep = "\t")
