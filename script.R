#! /usr/bin/env Rscript

#Package upload
library(magrittr)
library(gprofiler2)
library(ggplot2)
library(dplyr)
library(R.devices)
library(data.table)
ensembl()
library(fgsea)
library(BiocParallel)

#setwd("/project/workfile")
loglevel <- 10

#GENCODE version 36 was chosed as the reference genome and has been indexed with the default k-mer values.
GRCh38.salmon.index <- "/path/salmon_index"
GRCh38.gtf <- "/path/GRCh38/gencode.V36.annotation.gtf.gz"

#There is a template for metadata.tsv
pl <- metadata("metadata.tsv") %>% pipeline() %>%
     options(salmon.index=GRCh38.salmon.index, gtf=GRCh38.gtf,nthreads=4,njobs=4)

BiocParallel::register(BiocParallel::MulticoreParam(pl$option$nthreads))


add.prefix <- function (ids, REND)
{
    sapply(ids,
           function (id)
           {
               if (is.na(id) || "" == id || is.null(id))
                   NA
               else
                   file.path("output", "rnaseq", id)
           })
}
pl$metadata$fastq1 <- add.prefix(pl$metadata$fastq1, "R1")
pl$metadata$fastq2 <- add.prefix(pl$metadata$fastq2, "R2")


pl <- pl  %>% salmon() %>%
    salmonGenes()

filename <- file.path("output", "salmon", "salmon.tpm.tsv")
if (! file.exists(filename)) {
    isoforms <- getResultsByClass(pl, "salmon_isoforms")
    write.table(isoforms$txi$abundance,
                filename,
                sep="\t", quote=FALSE)
}
filename <- file.path("output", "salmon", "salmon.genes.tpm.tsv")
if (! file.exists(filename)) {
    genes <- getResultsByClass(pl, "salmon_genes")
    write.table(genes$txi$abundance,
                filename,
                sep="\t", quote=FALSE)
}


FCFilter <- function (res, threshold=-1)
{
    res <- res[complete.cases(res),]
    res[abs(res$log2FoldChange) > threshold,]
}

de.write <- function (pl, res.name, filename)
{
    
    res <- getResultsByName(pl, res.name)
    res <- res[order(res$padj),]

    mariadb <- ensembl()
    enrich <- geneDescription(mariadb, row.names(res))

    res <- cbind(res, enrich)
    
    write.table(res, filename, sep="\t", quote = FALSE)

    return(res)
}
    
####Differential gene expression analysis for four diseases

#### type 2 diabetes (T2D), non-diabetic (ND)
print("===== T2D vs ND")
filename.t2d <- "output/de-t2d.tsv"
if (! file.exists(filename.t2d)) {
    pl <- pl %>% 
        rucdr::filter((condition %in% c("ND_T2D", "T2D")) &
                      (! is.na(quant.sf.fn)), unfiltered=T) %>%
        deseq2(design=~batch+condition) %>%
        deseq2Results("T2D", "ND_T2D", condition="condition",
                      independentFiltering=FALSE, name="T2D")
    resT2D <- de.write(pl, "T2D", filename.t2d)
} else {
    resT2D <- read.table(filename.t2d, header = TRUE,
                         sep="\t", quote="")
}

#### multiple sclerosis (MS)
print("===== MS vs control")
filename.ms <- "output/de-ms.tsv"
if (! file.exists(filename.ms)) {
    pl <- pl %>%
        rucdr::filter((condition %in% c("CTL_MS", "MS")) &
                      (! is.na(quant.sf.fn)), unfiltered=T) %>%
        deseq2(design=~condition) %>%
        deseq2Results("MS", "CTL_MS", condition="condition",
                      independentFiltering=FALSE, name="MS")
    resMS <- de.write(pl, "MS", filename.ms)
} else {
    resMS <- read.table(filename.ms, header = TRUE,
                         sep="\t",quote="")
}
                      
#### type 1 diabetes (T1D), non-diabetic (ND)
print("===== T1D vs ND")
filename.t1d <- "output/de-t1d.tsv"
if (! file.exists(filename.t1d)) {
    pl <- pl %>%
        rucdr::filter((condition %in% c("ND_T1D", "T1D")) &
                      (! is.na(quant.sf.fn)), unfiltered=T) %>%
        deseq2(design=~batch+condition) %>%
        deseq2Results("T1D", "ND_T1D", condition="condition",
                      independentFiltering=FALSE, name="T1D")
    resT1D <- de.write(pl, "T1D", filename.t1d)
} else {
    resT1D <- read.table(filename.t1d, header = TRUE,
                         sep="\t", quote="")
}

#### Alzheimer's diasese
print("===== AD vs control")
filename.ad <- "output/de-ad.tsv"
if (! file.exists(filename.ad)) {
    pl <- pl %>%
        rucdr::filter((condition %in% c("CTL_AD", "AD")) &
                      (! is.na(quant.sf.fn)), unfiltered=T) %>%
        deseq2(design=~condition) %>%
        deseq2Results("AD", "CTL_AD", condition="condition",
                      independentFiltering=FALSE, name="AD")
    resAD <- de.write(pl, "AD", filename.ad)
} else {
    resAD <- read.table(filename.ad, header = TRUE,
                         sep="\t", quote="")
}

dir.create("output/rrho/", showWarnings=FALSE)

#### RRHO analysis for pairwise diesease

res2rrho <- function (res)
{
    results <- res$log2FoldChange
    names(results) <- row.names(res)
    return(results)
}


colfunc <- colorRampPalette(c("#eb3434", "#eb9334", "#ebeb34", "#49eb34", "#34eba5", "#34b4eb", "#3446eb"))
mycol = colfunc(1000)

call.rrho <- function(res1, res2, labels=c("res1","res2"),
                      alternative="two.sided")
{
    mariadb <- ensembl()

    fn1 <- gsub(pattern=" ", replacement="-", x=labels[1])
    fn2 <- gsub(pattern=" ", replacement="-", x=labels[2])
    fn <- file.path("output", "rrho", paste0(fn1, "_", fn2))
    res1.res2 <- RRHO(res2rrho(res1), res2rrho(res2), stepsize=500)
    ggsave(paste0(fn, ".tiff"),
           plot=plot(res1.res2,
                     labels=labels, colors = c(mycol, rev(mycol)),repel.force=100,show.legend=TRUE,base_size=20),dpi=300, units="cm",width=17, height=15)
    tiff(paste0(fn, "_3D.tiff"))
    plot3d(res1.res2,labels=labels)
    dev.off()

    lapply(c("upup", "downdown", "updown", "downup"),
          function (direction)
           {
               logging("call.rrho(): Getting genes")
               genes <- getEnrichment(res1.res2,
                                      direction)

               enrich <- geneDescription(mariadb, genes)
               print(paste("Writting", fn, "for direction", direction))
               print(head(enrich))
               write.table(enrich,
                           paste0(fn, "_", direction, ".tsv"),
                           sep="\t", quote=FALSE)
           })

    return(res1.res2)
}


results <- list(
    list(resT1D, "T1D vs ND"),
    list(resT2D, "T2D vs ND"),
    list(resMS, "MS vs CTL"),
    list(resAD, "AD vs CTL"))

nresults <- length(results)

cls <- c(
    apply(expand.grid(1:nresults, 1:nresults), 1,
          function(row)
          {
              i <- row[1]
              j <- row[2]
              if (i < j)
              {
                  x <- results[[i]]
                  y <- results[[j]]
                  function ()
                      call.rrho(x[[1]],y[[1]],
                                labels=c(x[[2]], y[[2]]))
              } else
              {
                     function () {NULL}
              }
          }),
    function ()
        call.rrho(resT1D, resAD, labels=c("T1D vs CTL", "AD vs CTL")),
    function ()
	call.rrho(resT1D, resT2D, labels=c("T1D vs CTL", "T2D vs CTL")),
    function ()
        call.rrho(resT1D, resMS, labels=c("T1D vs CTL", "MS vs CTL")),
    function ()
        call.rrho(resT2D, resT1D, labels=c("T2D vs CTL", "T1D vs CTL")),
    function ()
        call.rrho(resT2D, resAD,labels=c("T2D vs CTL", "AD vs CTL")),
    function ()
        call.rrho(resT2D, resMS, labels=c("T2D vs CTL", "MS vs CTL")),
    function ()
	call.rrho(resAD, resT2D, labels=c("AD vs CTL", "T2D vs CTL")),
    function ()
        call.rrho(resAD, resMS, labels=c("AD vs CTL", "MS vs CTL")),
    function ()
	call.rrho(resAD, resT1D, labels=c("AD vs CTL", "T1D vs CTL")),
    function ()
	call.rrho(resMS, resT1D, labels=c("MS vs CTL", "T1D vs CTL")),
    function ()
	call.rrho(resMS, resT2D, labels=c("MS vs CTL", "T2D vs CTL")),
    function ()
	call.rrho(resMS, resAD, labels=c("MS vs CTL", "AD vs CTL"))
)


print(paste("Length", length(cls)))
rrho <- mclapply(cls,
                 function (x) x(),
                 mc.cores=pl$option$nthreads * pl$option$njobs
                 )

print (rrho)

call.gprofiler <- function (genes.fn, method="top-5-bygroup", top.n=5)
{
    output.fn <- sub(pattern = "(.*)\\.tsv$",
                     replacement = paste0("\\1-gprofiler-", method,
                                         ".png"), genes.fn)
    result.fn <- sub(pattern = "(.*)\\.tsv$",
                     replacement = paste0("\\1-gprofiler-", method,
                                          ".tsv"), genes.fn)
    genes <- read.table(genes.fn, header = TRUE, sep="\t", quote="")
    if (nrow(genes) == 0)
    {
        logging(message=paste0("call.gprofiler(): ", genes.fn, " no genes"))
        return(NULL)
    }
    sources <- c("KEGG",  "REAC")
    gprof <- gost(row.names(genes), organism = "hsapiens", sources=sources)
    if (is.null(gprof$result) || nrow(gprof$result) == 0)
    {
        logging(message=paste0("call.gprofiler(): ", genes.fn, " no results"))
        return(NULL)
    }

    write.table(gprof$result[,2:(ncol(gprof$result)-1)],
                result.fn,
                sep="\t", quote=FALSE)

    res <- if ( "top-5-bygroup" == method)
           {
               gprof$result %>%
                   group_by(source) %>%
                   top_n(n = -top.n, wt = p_value)
           } else {
               gprof$result %>%
                   top_n(n = -top.n, wt = p_value)
           }

    term.ids <- res$term_id
    gprof.plot <- gostplot(gprof, interactive=FALSE)
    R.devices::suppressGraphics(
    {

        gg <- publish_gostplot(gprof.plot,
                               highlight_terms=term.ids)
        w <- max(sapply(res$term_name, nchar))
        ggsave(output.fn, height = 7 + nrow(res)/5, width = 3+w/10)
    })
}


files <- list.files(path="output/rrho/", pattern = "(down|up)(down|up)\\.tsv$",
                    full.names=TRUE)
mclapply(c(lapply(files, function (x) function() call.gprofiler(x)),
            lapply(files, function (x) function() call.gprofiler(x, method="top-30", top.n=30))),
         function (x) x(),
         mc.cores=pl$option$nthreads * pl$option$njobs)


#### TODO: gene set enrichment analysis (GSEA)
library(fgsea)
ensembl2entrez <- read.table("/srv/genomic_data/ensembl/GRCh38/hum_gene2ensembl.gz", 
                             header=T, sep="\t",quote = "", comment.char = "")
ensembl2entrez <- unique(ensembl2entrez[, c("Ensembl_gene_identifier", "GeneID")])
colnames(ensembl2entrez) <- c("ensembl_id", "entrez")

pathways.reactome <- gmtPathways("/home/xiaoyan/05.reference/geneset/c2.cp.reactome.v7.4.entrez.gmt")
pathways.kegg <- gmtPathways("/home/xiaoyan/02.project/01.AS/02.workfile/06.pipe/kegg.168.gmt")

getRank <- function(de.genes)
{
    df <- data.frame(ensembl_id=rownames(de.genes),
                     stat=de.genes$stat)
    df <- merge(df,ensembl2entrez,by="ensembl_id")
    df <- df[ ! duplicated(df$entrez),]
    ranks <- df$stat
    names(ranks) <- df$entrez
    return(ranks)
}


sortPathways <- function (fgseaRes) 
{ 
    fgseaRes$leadingEdge <- vapply(fgseaRes$leadingEdge, 
                                   paste, collapse = ", ", character(1L))
    fgseaRes[order(fgseaRes$padj,abs(1/fgseaRes$NES)),] 
}

kablePathways <- function (fgseaRes.sorted, alpha=0.05)
{
    fgseaRes.sorted[fgseaRes.sorted$padj < alpha,
                    c("pathway", "pval",    "padj",    "NES",     "size")] %>%
        kable(booktabs = TRUE, longtable = TRUE) %>%
        kable_styling(font_size = 7, latex_options = "striped") %>%
        column_spec(c(1), width = "50em")

}

ggplotPathways <- function (fgseaRes.sorted, name="Hallmark", n=42, 
                            alpha=0.05)
{
    fgseaRes.sorted$pathway <- strtrim(fgseaRes.sorted$pathway,45)
    ggplot(head(fgseaRes.sorted,n), aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj < alpha)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title=paste0("Top ", n, " ",
                          name, " pathways (fGSEA)")) + 
        theme_minimal() +
        theme(plot.title = element_text(hjust=1))
 
}

####GSEA for T1D
ranks <- getRank(resT1D)
fgseaRes.reactome <- sortPathways(fgsea(pathways=pathways.reactome, stats=ranks,
                                        minSize = 5, maxSize = 500,
                                        nperm=50000, nproc=4))
fgseaRes.kegg <- sortPathways(fgsea(pathways=pathways.kegg, stats=ranks,
                                    minSize = 5, maxSize = 500,
                                    nperm=50000, nproc=4))

write.table(fgseaRes.reactome,
            "output/t1d-fgsea-reactome.tsv",
            sep="\t", quote=FALSE)
write.table(fgseaRes.kegg,
            "output/t1d-fgsea-kegg.tsv",
            sep="\t", quote=FALSE)

gg <- ggplotPathways(fgseaRes.reactome)
ggsave("output/t1d-fgsea-reactome.png",
       plot=gg)
gg <- ggplotPathways(fgseaRes.kegg, name="KEGG",)
ggsave("output/t1d-fgsea-kegg.png", 
       plot=gg)

####GSEA for T2D
ranks <- getRank(resT2D)
fgseaRes.reactome <- sortPathways(fgsea(pathways=pathways.reactome, stats=ranks,
                                        minSize = 5, maxSize = 300,
                                        nperm=50000, nproc=4))
fgseaRes.kegg <- sortPathways(fgsea(pathways=pathways.kegg, stats=ranks,
                                    minSize = 5, maxSize = 300,
                                    nperm=50000, nproc=4))

write.table(fgseaRes.reactome,
            "output/t2d-fgsea-reactome.tsv",
            sep="\t", quote=FALSE)
write.table(fgseaRes.kegg,
            "output/t2d-fgsea-kegg.tsv",
            sep="\t", quote=FALSE)

gg <- ggplotPathways(fgseaRes.reactome)
ggsave("output/t2d-fgsea-reactome.png",
       plot=gg)
gg <- ggplotPathways(fgseaRes.kegg, name="KEGG",)
ggsave("output/t2d-fgsea-kegg.png", 
       plot=gg)

##GSEA for AD
ranks <- getRank(resAD)
fgseaRes.reactome <- sortPathways(fgsea(pathways=pathways.reactome, stats=ranks,
                                        minSize = 5, maxSize = 500,
                                        nperm=50000, nproc=4))
fgseaRes.kegg <- sortPathways(fgsea(pathways=pathways.kegg, stats=ranks,
                                    minSize = 5, maxSize = 500,
                                    nperm=50000, nproc=4))

write.table(fgseaRes.reactome,
            "output/ad-fgsea-reactome.tsv",
            sep="\t", quote=FALSE)
write.table(fgseaRes.kegg,
            "output/ad-fgsea-kegg.tsv",
            sep="\t", quote=FALSE)

gg <- ggplotPathways(fgseaRes.reactome)
ggsave("output/ad-fgsea-reactome.png",
       plot=gg)
gg <- ggplotPathways(fgseaRes.kegg, name="KEGG",)
ggsave("output/ad-fgsea-kegg.png",
       plot=gg)

##GSEA for MS
ranks <- getRank(resMS)
fgseaRes.reactome <- sortPathways(fgsea(pathways=pathways.reactome, stats=ranks,
                                        minSize = 5, maxSize = 500,
                                        nperm=50000, nproc=4))
fgseaRes.kegg <- sortPathways(fgsea(pathways=pathways.kegg, stats=ranks,
                                    minSize = 5, maxSize = 500,
                                    nperm=50000, nproc=4))

write.table(fgseaRes.reactome,
            "output/ms-fgsea-reactome.tsv",
            sep="\t", quote=FALSE)
write.table(fgseaRes.kegg,
            "output/ms-fgsea-kegg.tsv",
           sep="\t", quote=FALSE)

gg <- ggplotPathways(fgseaRes.reactome)
ggsave("output/ms-fgsea-reactome.png",
       plot=gg)
gg <- ggplotPathways(fgseaRes.kegg, name="KEGG",)
ggsave("output/ms-fgsea-kegg.png",
       plot=gg)

