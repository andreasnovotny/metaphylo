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

###################
###################

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

###################
###################

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


###################
###################


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

###################
###################

#' Extract Taxonomy table of a phyloseq object

#'
#' @param PHYLOSEQ (Required) An class S4 class phyloseq object.
#' @return A matrix.
#' @export

GetTax <- function(ps) {
  require(phyloseq)
  return(as.matrix(ps@tax_table@.Data))
}

###################
###################

#' Extract Sample metadata table of a phyloseq object

#' @param PHYLOSEQ (Required) An class S4 class phyloseq object.
#' @return A data frame
#' @export

GetSample <- function(ps) {
  require(phyloseq)
  return(as.data.frame(ps_16S@sam_data@.Data))
}
