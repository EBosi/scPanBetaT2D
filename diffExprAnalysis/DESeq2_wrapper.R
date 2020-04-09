# import gene expr matrix and factors
sample.raws <- read.csv("/home/bosi/analisi/beta_raw.csv")
myFact <- read.csv("/home/bosi/analisi/beta_factors.csv")

###############
# DESeq2
BiocManager::install("DESeq2")
library("DESeq2")
library("BiocParallel")

# to increase speed I used 20 proc
register(MulticoreParam(20))

rownames(myFact) <- rownames(sample.raws)
dds <- DESeqDataSetFromMatrix(countData = t(sample.raws),
                              colData = myFact,
                              design= ~ dataset + diabetes)

dds <- DESeq(dds,
             sfType="poscounts",
             parallel=TRUE,
             useT=TRUE,
             minmu=1e-6,
             minReplicatesForReplace=Inf)

deseq_res <- results(dds)
deseq_res_nona <- deseq_res[!is.na(deseq_res$padj),]
dim(deseq_res_nona[deseq_res_nona$padj<0.05,])
up_reg <- deseq_res_nona[(deseq_res_nona$padj<0.05)&(deseq_res_nona$log2FoldChange>=1),]
down_reg <- deseq_res_nona[(deseq_res_nona$padj<0.05)&(deseq_res_nona$log2FoldChange<=-1),]

dim(up_reg)
dim(down_reg)

deseq_res_df <- data.frame(genes=c(rownames(down_reg),rownames(up_reg)),up_or_down=c(rep("down",16),rep("up",210)))
write.csv(deseq_res_df,"diffExprAnalysis/deseq_res.csv")

deseq_res[!is.na(deseq_res$pvalue),]$log2FoldChange

deseq_res_df_ALL <- data.frame(genes=rownames(deseq_res[!is.na(deseq_res$pvalue),]),
                               foldChange=deseq_res[!is.na(deseq_res$pvalue),]$log2FoldChange,
                               pval=deseq_res[!is.na(deseq_res$pvalue),]$pvalue,
                               FDR=deseq_res[!is.na(deseq_res$pvalue),]$padj)
write.csv(deseq_res_df_ALL,"diffExprAnalysis/deseq_res_ALL.csv")

# volcano plot
plot(deseq_res_nona$log2FoldChange, -log10(deseq_res_nona$padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6, xlim=c(-4,4),col=deseq_res_nona$color)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(0.05), col="brown")

