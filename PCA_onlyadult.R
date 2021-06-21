library(cluster)
library(Biobase)
library(qvalue)
NO_REUSE = F

# try to reuse earlier-loaded data if possible
if (file.exists("77samples.genes.counts.matrix.RData") && ! NO_REUSE) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    load("77samples.genes.counts.matrix.RData")
} else {
    print('Reading matrix file.')
    primary_data = read.table("../../77samples.genes.counts.matrix", header=T, com='', sep="\t", row.names=1, check.names=F)
    primary_data = as.matrix(primary_data)
}
source("/Users/danielle/Desktop/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/R/heatmap.3.R")
source("/Users/danielle/Desktop/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
source("/Users/danielle/Desktop/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/R/pairs3.R")
data = primary_data
samples_data = read.table("/Users/danielle/Desktop/2018_04_LSP/samples/samples_onlyadult.txt", header=F, check.names=F, fill=T)
samples_data = samples_data[samples_data[,2] != '',]
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])
data = data[, colnames(data) %in% samples_data[,2], drop=F ]
nsamples = length(sample_types)
larval_colors = c('#FFCF1AFF','#FFCF1AFF','#FFCF1AFF','#595959FF','#595959FF','#595959FF')
adult_colors= c('#4b984b','#595959FF','#1f497dFF','#9474BDFF','#b12334FF','#DD7E00')
larval_symbols = c(0,4,5)
female_symbols = c(1,1,1,1,1,1)
male_symbols = c(2,2,2,2,2,2)
sample_colors = c(adult_colors,adult_colors)
female_symbols = c(1,1,1,1,1,1)
male_symbols = c(2,2,2,2,2,2)
sample_symbols = c(female_symbols,male_symbols)
print(sample_colors)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
    samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
    sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
# reorder according to sample type.
tmp_sample_reordering = order(sample_factoring)
data = data[,tmp_sample_reordering,drop=F]
sample_factoring = sample_factoring[tmp_sample_reordering]
initial_matrix = data # store before doing various data transformations
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
write.table(data, file="77samples.genes.counts.matrix.log2.dat", quote=F, sep='	');
pdf("77samples.genes.counts.matrix.log2.principal_components.pdf")
data = as.matrix(data)
# Z-scale and center the genes across all the samples for PCA
prin_comp_data = initial_matrix
prin_comp_data = log2(prin_comp_data+1)
prin_comp_data = scale(prin_comp_data)
prin_comp_data = t(scale(t(prin_comp_data), center=TRUE, scale=F)) # just center trans expr level, retain original effect size.
pca = prcomp(t(prin_comp_data), center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1));
write.table(pca$rotation, file="77samples.genes.counts.matrix.log2.PCA.loadings", quote=F, sep="	")
write.table(pca$x, file="77samples.genes.counts.matrix.log2.PCA.scores", quote=F, sep="	")
print(sample_colors)
for (i in 1:(max(3,2)-1)) {
    xrange = range(pca$x[,i])
    yrange = range(pca$x[,i+1])
    samples_want = rownames(pca$x) %in% sample_type_list[[sample_types[1]]]
    pc_i_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i]*100)
    pc_i_1_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i+1]*100)
    plot(pca$x[samples_want,i], pca$x[samples_want,i+1], xlab=paste('PC',i, pc_i_pct_var), ylab=paste('PC',i+1, pc_i_1_pct_var), xlim=xrange, ylim=yrange, col=sample_colors[1])
    for (j in 2:nsamples) {
        samples_want = rownames(pca$x) %in% sample_type_list[[sample_types[j]]]
        points(pca$x[samples_want,i], pca$x[samples_want,i+1], col=sample_colors[j], pch=sample_symbols[j])
    }
    plot.new()
    legend('topleft', as.vector(sample_types), col=sample_colors, pch=sample_symbols, ncol=2)
}

par(def.par)
pcscore_mat_vals = pca$rotation[,1:3]
pcscore_mat = matrix_to_color_assignments(pcscore_mat_vals, col=colorpanel(256,'purple','black','yellow'), by='row')
colnames(pcscore_mat) = paste('PC', 1:ncol(pcscore_mat))
## generate heatmaps for PC extreme vals
uniq_genes = c()
for (i in 1:3) {
    ## get genes with extreme vals
    print(paste('range', range(pca$rotation[,i])))
    ordered_gene_indices = order(pca$rotation[,i])
    num_genes = length(ordered_gene_indices)
    extreme_ordered_gene_indices = unique(c(1:100, (num_genes-100):num_genes))
    print('extreme ordered gene indices')
    print(extreme_ordered_gene_indices)
    selected_gene_indices = ordered_gene_indices[extreme_ordered_gene_indices]
    print('selected gene indices');print(selected_gene_indices);
    print('PC scores:');print(pca$rotation[selected_gene_indices,i])
    selected_genes_matrix = prin_comp_data[selected_gene_indices,]
    pc_color_bar_vals = pcscore_mat_vals[selected_gene_indices,i]
    print(pc_color_bar_vals)
    pc_color_bar = as.matrix(pcscore_mat[selected_gene_indices,i])
uniq_genes = unique(c(uniq_genes, rownames(selected_genes_matrix)))
    write.table(selected_genes_matrix, file=paste("77samples.genes.counts.matrix.log2", '.PC',i,'_extreme',100,'.matrix', sep=''), quote=F, sep="	")
    heatmap.3(selected_genes_matrix, col=greenred(256), scale='none', density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, margins=c(10,10), cex.main=0.75, RowSideColors=pc_color_bar, cexRow=0.5, main=paste('heatmap for', 100, ' extreme of PC', i), ColSideColors=sampleAnnotations)
}
all_selected_genes_matrix = prin_comp_data[uniq_genes,]
write.table(all_selected_genes_matrix, file=paste("77samples.genes.counts.matrix.log2", '.PC_all','_extreme',100,'.matrix', sep=''), quote=F, sep="	")
heatmap.3(all_selected_genes_matrix, col=greenred(256), scale='none', density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, margins=c(10,10), cex.main=0.75, cexRow=0.5, main=paste('heatmap for ALL selected ', 100, ' extreme of all PCs'), ColSideColors=sampleAnnotations)
dev.off()
gene_cor = NULL
