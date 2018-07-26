#' Filters DADA2 taxonomy table for low bootstrap values.
#'
#' This function takes the output taxonomy table from the DADA2::assignTaxonomy and filters it based on bootstrap value.
#'
#' @param tax (Required) A list of two matrices. As returned from DADA2::assignTaxonomy(seqtab, database, outputBootstraps =  TRUE).
#' @param minboot (Optional, default=75) A Number between 0-100. Minimal accepted bootstrap value.
#' @return A list of two matrices (modified tax object).
#' @examples
#' tax1 <- FilterTaxonomy(tax)
#' tax2 <- FilterTaxonomy(tax, minboot=80)
#' @export
FilterTaxonomy <- function(tax, minboot=75) {

  TAX <- as.matrix(tax$tax)
  BOOT <- as.matrix(tax$boot)
  for (row in 1:nrow(TAX)) {
    for (col in 1:ncol(TAX)) {
      if (BOOT[row,col] < minboot) {
        TAX[row,col] <- NA
      }
    }
  }
  tax$tax <- TAX
  return(tax)
}












#' Fills unassigned taxonomic levels with useful information.
#'
#' This function assigns the higher taxonomic ranks to unassigned taxoomic levels, using the output of DADA2::assignTaxonomy.
#'
#' @param tax (Required) A list of two matrices. As returned from DADA2::assignTaxonomy(seqtab, database, outputBootstraps =  TRUE).
#' @param string (Optional, default="_x") A string. Will be attached at the end of the name comming from the higher taxonomic level.
#' @return A list of two matrices (modified tax object).
#' @examples
#' tax1 <- CompleteTaxonomy(tax)
#' tax2 <- CompleteTaxonomy(tax, string="X")
#' @export
CompleteTaxonomy <- function(tax, string="_x") {

  TAX <- as.matrix(tax$tax)
  for (row in 1:nrow(TAX)) {
    for (col in 1:ncol(TAX)) {
      if (is.na(TAX[row,col])) {
        TAX[row,col] <- paste(TAX[row,(col-1)],string,sep="")
      }
    }
  }
  tax$tax <- TAX
  return(tax)
}










#' Create phyloseq object from DADA2 output.
#'
#' This function takes the output from DADA2 package and returns a complete phyloseq object.
#'
#' @param seqtab (Required) A sequence table matrix as returned from dada2::makeSequenceTable()
#' @param tax (Required) A list of two matrices. As returned from DADA2::assignTaxonomy(seqtab, database, outputBootstraps =  TRUE).
#' @param metadata (Required) A data frame.
#' @return A S4 class Phyloseq object.
#' @examples
#' ps <- Dada2ToPhyloseq(seqtab, tax, metadata)
#' @export
Dada2ToPhyloseq <- function(seqtab, tax, metadata) {
  require(phyloseq)
  OTU <- otu_table(as.matrix(seqtab), taxa_are_rows = FALSE)
  TAX <- tax_table(as.matrix(tax$tax))
  SAM <- sample_data(as.data.frame(metadata))
  ps <- phyloseq(OTU, TAX, SAM)
  return(ps)
}










#' Writes a phyloseq object as standard OTU-like table
#'
#' This function takes a full phyloseq object and merges the count matrix with the taxonomy table. It also provides the option to give sample names.
#'
#' @param PHYLOSEQ (Required) An class S4 class phyloseq object
#' @param HEADER (Required) A character string. Representing a column in the sample dataset.
#' @param FILENAME (Optional, defalt="./ASV.csv") Character vector, filename.
#' @return A S4 class Phyloseq object.
#' @examples
#' ps <- PhyloseqToCSV(ps, "SORTED_type", "./phyloseq.csv")
#' @export
PhyloseqToCSV <- function(PHYLOSEQ, HEADER, FILENAME="./ASV.csv") {
  require(phyloseq)
  tax <- as.matrix(PHYLOSEQ@tax_table@.Data)
  count <- t(as.matrix(PHYLOSEQ@otu_table@.Data))
  metadata <- as.data.frame(PHYLOSEQ@sam_data)
  colnames(count) <- metadata[[HEADER]]
  ASV <- merge(tax,count, by=0)
  write.csv(ASV, FILENAME)
  return(ASV)
}












#' Extract Taxonomy table of a phyloseq object

#'
#' @param PHYLOSEQ (Required) An class S4 class phyloseq object.
#' @return A matrix.
#' @export

GetTax <- function(ps) {
  require(phyloseq)
  return(as.matrix(ps@tax_table@.Data))
}














#' Extract Sample metadata table of a phyloseq object

#' @param PHYLOSEQ (Required) An class S4 class phyloseq object.
#' @return A data frame
#' @export

GetSample <- function(ps) {

  require(phyloseq)

  return(as.data.frame(ps_16S@sam_data@.Data))
}
















#' Makes a DESeq2 analysis of a phyloseq object
#'
#' This function makes a DESeq2 stastistical analysis and writes the results to a redable table.
#'
#' @param PHYLOSEQ (Required) An class S4 class phyloseq object
#' @param VARIABLE (Required) A ~ variable character string. Representing a column in the sample dataset.
#' @param speciesAsNames (Optional) Boolean (Default = TRUE). If true, species names will be the tables row name.
#' @return A DEseq2 results table with significant expression level values.
#' @examples
#' ps <- DEseqTable(ps, ~ SORTED_type)
#' @export

deseqTable <- function(ps, design, speciesAsNames=TRUE) {

  require(phyloseq)
  require(DESeq2)

  tmp <- phyloseq_to_deseq2(ps, design)
  tmp <- DESeq(tmp, test="Wald", fitType="parametric")
  res = results(tmp, cooksCutoff = FALSE)
  alpha = 0.01
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PHYLOSEQ)[rownames(sigtab), ], "matrix"))
  if (speciesAsNames == TRUE) {
    sigtab$Sequence <- rownames(sigtab)
    rownames(sigtab) <- sigtab$Species
  }
  return(sigtab)
}
















#' Filtering of relevant species
#'
#' This function filters away taxa only present in few samples.
#'
#'
#' @param ps (Required) An class S4 class phyloseq object.
#' @param method (Optional, default="logRaw") "logRaw" (log transformation) or "vst" (variance stabilising DESeq2 transformation)
#' @param cutoff (Optional, default=0.8) A numeric value between 0:1.
#' @param mincount (Optional, default=10) An integer, bigger than 0. Abundance value as been counted as presence. Only effective if method="logRaw".
#' @param dsign (Required if method="rsv") A ~ variable character string. Representing a column in the sample dataset.
#' @return A DEseq2 results table with significant expression level values.
#' @examples
#' ps <- abundancePlot(ps_18S, cutoff=0.6)
#' @export

relevantSpecies <- function(ps, method="logRaw", cutoff=0.75, mincount=10, design) {

  require(phyloseq)

  if (method=="logRaw") {
    ps <- filter_taxa(ps, function(x) sum(x>mincount)>=((length(x)*cutoff)),TRUE)
    ps <- transform_sample_counts(ps, function(x) log(x))
  }

  if (method=="vst") {
    ps <- getVst(ps, design)
    ps <- filter_taxa(ps, function(x) sum(x>0)>=((length(x)*cutoff)),TRUE)
  }

  return(ps)
}















#' Boxplot with abundance of relevant species
#'
#' This function make boxplots of filtered phyloseq objects.
#'
#'
#' @param ps (Required) An class S4 class phyloseq object.
#' @param x (Optional, default="Species") A character string. Representing a taxonomic level.
#' @param fill (Optional, default="Genus") A character string. Representing a taxonomic level.
#' @return A DEseq2 results table with significant expression level values.
#' @examples
#' ps <- abundancePlot(ps_18S, x="Genus")
#' @export

abundancePlot <- function(ps, x="Species", fill="Genus") {

  require(phyloseq)
  require(ggplot2)

  graph_theme<- theme(
    #axis.text.y = element_text(size=20),
    axis.text.x = element_text(face='plain',size=10, angle=-90, hjust=0),
    #axis.title.y = element_text(size=24),
    #axis.title.x = element_text(face='plain',size=24, vjust=-0.1, hjust=0.6),
    #legend.text = element_text(size = 20),
    #legend.title = element_text(size= 24),
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
  plot <- ps %>%
    psmelt() %>%
    ggplot(aes_string(x=x, y="Abundance", fill=fill)) +
    geom_boxplot() +
    geom_point() +
    #labs(x="Genus", y="Log(Gene Abundance)") +
    theme_bw() +
    graph_theme
  return(plot)
}










#' Variance stabilizing transformation of phyloseq object using DESeq2
#'
#' @param ps (Required) An S4 class phyloseq object.
#' @param design (Required) A ~ variable character string. Representing a column in the sample dataset.
#' @return An S4 class phyloseq object, with transformed counts.
#' @examples
#' ps_18S_1 <- getVST(ps_18S, ~MONTH+SORTED_genus)
#' @export

getVst <- function(ps, design) {

  require(phyloseq)
  require(DESeq2)

  dds <- phyloseq_to_deseq2(ps, design)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- getVarianceStabilizedData(dds)
  otu_table(ps) <- otu_table(dds, taxa_are_rows = TRUE)
  return(ps)
}








#' Filters Phyloseq object by taxa with significantly different abundance
#'
#' @param ps (Required) An S4 class phyloseq object.
#' @param design (Required) A ~ variable character string. Representing a column in the sample dataset.
#' @param number (Optional, default=50) An integer. Biggest number of taxa to keep in the filtered dataset.
#' @param alpha (Optional, default=0.01) A value. Threshold for p-value.
#' @param vst (Optional, default=FALSE) Boolean operator. If TRUE, the dataset gets variance stabalising transformation.
#' @return An S4 class phyloseq object, filtered
#' @examples
#' ps_18S_1 <- getSignificant(ps_18S, ~MONTH+SORTED_genus)
#' @export

getSignificant <- function(ps, design, number = 50, alpha = 0.01, vst=FALSE) {

  require(phyloseq)
  require(DESeq2)

  dds <- phyloseq_to_deseq2(ps, design) # Import data to DESeq2
  dds <- DESeq(dds) # Run DESeq2 analysis
  res <- results(dds)
  res <- res[order(res$padj, na.last = NA), ] #Sort result table in order, most significant first.
  keepOTUs <- rownames(res[res$padj > alpha, ])[1:number] # Keep only the top 50 most significant taxa.

  if (vst==FALSE) {
    ps1 <- prune_taxa(keepOTUs, ps) # Remove nonsignificant taxa from phyloseq
  }

  if (vst==TRUE) {
    psVst <- getVst(ps, design) #Make variance stabalizing transformation
    ps1 <- prune_taxa(keepOTUs, psVst) # Remove nonsignificant taxa from phyloseq
  }

  return(ps1)

}



