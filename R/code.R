#' Filters DADA2 taxonomy table for low bootstrap values.
#'
#' This function takes the output taxonomy table from the DADA2::assignTaxonomy and filters it based on bootstrap value.
#'
#' @param tax (Required) A list of two matrices. As returned from DADA2::assignTaxonomy(seqtab, database, outputBootstraps =  TRUE).
#' @param minboot (Optional, default=75) A Number between 0-100. Minimal accepted bootstrap value.
#' @return A list of two matrices (modified tax object).
#' @examples
#' data("tax_16S")
#' tax_16S <- filterTaxonomy(tax_16S, minboot = 75)
#'
#' @export
filterTaxonomy <- function(tax, minboot=75) {

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
#' data("tax_16S")
#' tax_16S <- filterTaxonomy(tax_16S, minboot = 75)
#' tax_16S <- completeTaxonomy(tax_16S, string = "_x")
#'
#' @export
completeTaxonomy <- function(tax, string="_x") {

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
#' data("metadata_16S")
#' data("seqtab_16S")
#' data("tax_16S")
#' ps_16S <- dada2ToPhyloseq(seqtab_16S, tax_16S, metadata_16S)
#'
#' @export
dada2ToPhyloseq <- function(seqtab, tax, metadata) {
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
#' @param ps (Required) An class S4 class phyloseq object
#' @param header (Required) A character string. Representing a column in the sample dataset.
#' @param filename (Optional, defalt="./ASV.csv") Character vector, filename.
#' @return A S4 class Phyloseq object.
#' @examples
#' data(ps_16S)
#' phyloseqToCSV(ps_16S, header="ID_friendly", filename="~/Desktop/RSV_16S.csv")
#' @export
phyloseqToCSV <- function(ps, header, filename="./ASV.csv") {
  require(phyloseq)
  tax <- as.matrix(ps@tax_table@.Data)
  count <- t(as.matrix(ps@otu_table@.Data))
  metadata <- as.data.frame(ps@sam_data)
  colnames(count) <- metadata[[header]]
  ASV <- merge(tax,count, by=0)
  write.csv(ASV, filename)
  return(ASV)
}


#' Writes a phyloseq object as matrix, specifying rownames as prey and colnames as species.
#'
#' This function takes a full phyloseq object and merges the count matrix with the taxonomy table. It also provides the option to give sample names.
#'
#' @param ps (Required) An class S4 class phyloseq object
#' @param colname (Required) A character string. Representing a column in the sample dataset.
#' @param rowname (Optional, defalt="Species") A character string. Representing a column in the sample dataset.
#' @return A S4 class Phyloseq object.
#' @examples
#' data(ps_18S)
#' phyloseqToMatrix(ps_18S, colnames="ID_friendly")
#' @export

phyloseqToMatrix <- function(ps, colname="NAME", rowname="Species") {
  require(phyloseq)
  tax <- as.matrix(ps@tax_table@.Data)
  tax_levels <- length(rank_names(ps))+1
  count <- t(as.matrix(ps@otu_table@.Data))
  metadata <- as.data.frame(ps@sam_data)
  metadata$NAME <- rownames(metadata)
  colnames(count) <- metadata[[colname]]
  ASV <- as.data.frame(merge(tax,count, by=0))
  rownames(ASV) <- ASV[[rowname]]
  ASV[1:tax_levels] <- list(NULL)
  ASV <- as.matrix(ASV)
  return(ASV)
}





#' Makes a DESeq2 analysis of a phyloseq object
#'
#' This function makes a DESeq2 stastistical analysis and writes the results to a redable table.
#'
#' @param ps (Required) An class S4 class phyloseq object
#' @param design (Required) A ~ variable character string. Representing a column in the sample dataset.
#' @param speciesAsNames (Optional) Boolean (Default = TRUE). If true, species names will be the tables row name.
#' @return A DEseq2 results table with significant expression level values.
#' @examples
#' data(ps_16S)
#' ps_16S <- subset_samples(ps_16S, MONTH == "aug")
#' ps_16S <- subset_samples(ps_16S, SORTED_type == "Rotifer")
#' ddstab_16S <- deseqTable(ps_16S, ~ SORTED_type)
#' @export

deseqTable <- function(ps, design, speciesAsNames=TRUE) {

  require(phyloseq)
  require(DESeq2)

  tmp <- phyloseq_to_deseq2(ps, design)
  tmp <- DESeq(tmp, test="Wald", fitType="parametric")
  res <- results(tmp, cooksCutoff = FALSE)
  alpha <- 0.01
  sigtab <- res[which(res$padj < alpha), ]
  sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  if (speciesAsNames == TRUE) {
    sigtab$Sequence <- rownames(sigtab)
    rownames(sigtab) <- sigtab$Species
  }
  return(sigtab)
}












#' Filtering of prevalent species.
#'
#' This function filters away taxa only present in few samples.
#'
#'
#' @param ps (Required) A class S4 class phyloseq object.
#' @param prevalence (Optional, default=0.75) A numeric value between 0:1. Least fraction of required presence.
#' @param minpercent (Optional, default=0.1) Numeric value between 0:100. A percentage of total library size. Lowest number to count as presence.
#' @return A class S4 class phyloseq object.
#' @examples
#' data(ps_16S)
#' ps_16S <- subset_samples(ps_16S, SORTED_genus == "Synchaeta")
#' ps_16S <- prevalentSpecies(ps_16S, prevalence=0.75, minpercent=1)
#' @export
prevalentSpecies  <- function(ps, prevalence=0.75, minpercent=0.1) {

  require(phyloseq)
  require(magrittr)

  ps <- ps %>%
    transform_sample_counts(function(x) (x/sum(x))*100) %>%
    filter_taxa(function(x) sum(x>minpercent)>=((length(x)*prevalence)),TRUE) %>%
    taxa_names() %>%
    prune_taxa(ps)

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
#' data(ps_16S)
#' ps_16S <- subset_samples(ps_16S, SORTED_genus == "Synchaeta")
#' ps_16S <- relevantSpecies(ps_16S, method = "vst", cutoff = 0.75, design = ~ ID_friendly)
#' plot <- abundancePlot(ps_16S, x = "Species", fill = "Genus")
#'
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
    geom_point(aes(shape="FILTER_poresize")) +
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
#' data(ps_16S)
#' ps_16S_vst <- getVST(ps_16S, ~MONTH+SORTED_genus)
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
#' data(ps_16S)
#' ps_16S <- subset_samples(ps_16S, SORTED_genus == "Synchaeta")
#' ps_16S_sig <- getSignificant(ps_16S, ~MONTH)
#' @export

getSignificant <- function(ps, design, number = 50, alpha = 0.01, vst=FALSE) {

  require(phyloseq)
  require(DESeq2)

  dds <- phyloseq_to_deseq2(ps, design) # Import data to DESeq2
  dds <- DESeq(dds) # Run DESeq2 analysis
  res <- results(dds)
  res <- res[order(res$padj, na.last = NA), ] #Sort result table in order, most significant first.
  keepOTUs <- rownames(res[res$padj <= alpha, ])[1:number] # Keep only the top 50 most significant taxa.

  if (vst==FALSE) {
    ps1 <- prune_taxa(keepOTUs, ps) # Remove nonsignificant taxa from phyloseq
  }

  if (vst==TRUE) {
    psVst <- getVst(ps, design) #Make variance stabalizing transformation
    ps1 <- prune_taxa(keepOTUs, psVst) # Remove nonsignificant taxa from phyloseq
  }

  return(ps1)

}


#' DESeq2 analysis: p.adj named vector.
#'
#' @param ps (Required) An S4 class phyloseq object.
#' @param design (Required) A ~ variable character string. Representing a column in the sample dataset.
#' @return A Named vector of p.adj values per species.
#' @examples
#' data(ps_18S)


#' @export

pVector <- function(ps, design){

  require(magrittr)
  require(phyloseq)
  require(DESeq2)

  sigtab <- phyloseq_to_deseq2(ps, design) %>%
    DESeq() %>%
    results()

  sigtab <- cbind(as(sigtab, "data.frame"),
                  as(tax_table(ps)[rownames(sigtab), ], "matrix"))

  p.vector <- sigtab$padj
  names(p.vector) <- sigtab$Species

  return(p.vector)
}








#' Not exported. Fit colors to matrix.
#'
#' @param psmatrix (Required) A matrix (see phyloseqToMatrix).
#' @param pvector (Required) A named vector with of p-values per taxa (see pVector)
#' @param sigcol (Optional, default="Red") Character string. R-base color to represent significans.
#' @param nonsigcol (Optional, default="Green") Character string. R-base color to represent significans.
#' @param medsigcol (Optional, default="Blue") Character string. R-base color to represent significans.
#' @return A vector with colrs.


pColors <- function(psmatrix, pvector,
                    sigcol="Red",
                    nonsigcol="Green",
                    medsigcol="Blue"){

  pcolors <- vector()
  for (x in row.names(psmatrix)) {
    if (is.na(pvector[[x]])) {
      pvector[[x]]<-1
    }
    if (pvector[[x]]<=0.01) {
      pcolors <- c(pcolors, sigcol)
    } else {
      if (pvector[[x]]<=0.05) {
        pcolors <- c(pcolors, medsigcol)
      } else {
        pcolors <- c(pcolors, nonsigcol)
      }
    }
  }
  return(pcolors)
}






#' A wrapper function for bipartite::plotweb
#'
#' @param psmatrix (Required) A matrix (see phyloseqToMatrix).
#' @param pvector (Required) A named vector with of p-values per taxa (see pVector)
#' @param sigcol (Optional) Character string. R-base color to represent significans.
#' @param nonsigcol (Optional, default="Green") Character string. R-base color to represent significans.
#' @param medsigcol (Optional, default="Blue") Character string. R-base color to represent significans.
#' @param ... (Opltional) See bipartite::plotweb for more options.
#' @return A bipartite::plotweb graph.
#' @examples
#' data(ps_18S)
#'
#' @export

foodWeb <- function(psmatrix, pvector,
                    col.high = "blue4",bor.col.high = "blue4",
                    bor.col.low = "white",
                    col.interaction = c("darkgray", "dimgrey"),
                    bor.col.interaction = c("dimgrey", "darkgrey"),
                    sigcol="Red", nonsigcol="Green", medsigcol="Blue",
                    ...){

  require(bipartite)

  pcolors <- pColors(psmatrix, pvector,
                     sigcol=sigcol,
                     nonsigcol=nonsigcol,
                     medsigcol=medsigcol)

  plotweb(psmatrix, method="normal",
          col.low = pcolors, bor.col.low = bor.col.low,
          col.high = col.high, bor.col.high = bor.col.high,
          col.interaction = col.interaction,
          bor.col.interaction = bor.col.interaction,
          ...)
}



#' Sort a phyloseq object by presence of prevalent taxa in another phyloseq object
#'
#' A wrapper function for phyloseq::prune_taxa
#'
#' @param ps_1 (Required) A phyloseq object to be filtered.
#' @param ps_2 (Required) A phyloseq object to filter by.
#' @param minpercent (Optional, default=1) Character string. R-base color to represent significans.
#' @return A filtered phyloseq object.
#' @examples
#' data(ps_18S)
#'
#' @export

sortByWater <- function(ps_1, ps_2, minpercent=1) {

  require(phyloseq)
  require(magrittr)

  ps_1 <- ps_2 %>%
    transform_sample_counts(function(x) x/sum(x)*100) %>%
    filter_taxa(function(x) sum(x)> minpercent, TRUE) %>%
    taxa_names() %>%
    prune_taxa(ps_1)

  return(ps_1)
}




#' Construct Summary table for barplots
#'
#' This function takes a phyloseq as unput (or possibly a tidy data-frame) and transforms it into a summary table.
#'
#' @param data (Required) A phyloseq object or a data frame.
#' @param x (Optional, default="Genus") Character string referring to the main explainatory variable.
#' @param y (Optional, default="Abundance") Character string reffering to the dependent variable. For phyloseq objects, this is always "Abundance".
#' @param variables (optional, default=c("Family", "Class", "Order")) A vector of strings reffereng to any aditional explainatory variable for plotting/colouring/splitting/grouping ect.
#' @return A formated data.frame. The aim of the data frame is to use for plotting in ggplot (or the wrapper metaphylo::barChart)
#' @examples
#' data(ps_18S)
#' ps_18S %>%
#'  subset_samples(SORTED_genus=="Synchaeta") %>%
#'  subset_samples(MONTH=="aug") %>%
#'  subset_taxa(Class!="Rotifera") %>%
#'  transform_sample_counts(function(x) (x/sum(x))*100) %>%
#'  filter_taxa(function(x) sum(x>1)>=((length(x)*0.5)),TRUE) %>%
#'  sumTable(x = "Species", variables = c("Class", "Order"))
#' @export

sumTable <- function(data, x="Genus", y="Abundance", variables=c("Family", "Class", "Order")) {

  require(phyloseq)
  require(magrittr)

  # 1. Converts the phyloseq object to a dataframe, and deleats empty rows
  if (class(data)=="phyloseq") {
    data <- data %>%
      filter_taxa(function(x) sum(x)>0, TRUE) %>%
      psmelt()
  }

  # 2. Construct the summary table
  BY <- list(x = data[[x]])
  for (entry in variables)  {
    BY[[entry]] <- data[[entry]]
  }

  SumTable <- aggregate(data[[y]],
                        by =BY,
                        FUN = function(x) c(mean = mean(x),
                                            sd = sd(x),
                                            n = length(x)))

  SumTable <- do.call(data.frame, SumTable)
  SumTable$se <- SumTable$x.sd / sqrt(SumTable$x.n)
  colnames(SumTable) <- c(x, variables,"mean", "sd", "n", "se")

  # Sort the summary table in nice order
  SumTable <- SumTable[with(SumTable, order(SumTable[[variables[1]]],-mean)),]
  SumTable$my.order <- seq(from=1, to=length(SumTable[,1]), by=1)


  return(SumTable)
}





#' barCharts out of a summaryTable using ggplot2
#'
#' This function takes a phyloseq as unput (or possibly a tidy data-frame) and transforms it into a summary table.
#'
#' @param sumtable (Required) A dataframe (preferably the output of metaphylo::sumTable function)
#' @param x (Optional, default="Genus") Character string referring to the main explainatory variable.
#' @param fill (optional, default="Order") A string reffereng to any aditional explainatory variable for filling.
#' @param errorbars (Optional, default=TRUE) Boolean operator. If True, standard error is plotted as errorbars.
#' @param ... Any aditional variables passed on to geom_bar(...)
#' @return A ggplot object. The aim of the data frame is to use for plotting in ggplot (or the wrapper metaphylo::barChart)
#' @examples
#' data(ps_18S)
#' ps_18S %>%
#'  subset_samples(SORTED_genus=="Synchaeta") %>%
#'  subset_samples(MONTH=="aug") %>%
#'  subset_taxa(Class!="Rotifera") %>%
#'  transform_sample_counts(function(x) (x/sum(x))*100) %>%
#'  filter_taxa(function(x) sum(x>1)>=((length(x)*0.5)),TRUE) %>%
#'  sumTable(x = "Species", variables = c("Class", "Order")) %>%
#'  barChart(x = "Family", fill = "Class")
#' @export

barChart <- function(sumtable, x="Genus", fill="Order", errorbars=TRUE, ...) {

  require(ggplot2)

  plot <- ggplot(sumtable) +
    geom_bar(aes(x = reorder(sumtable[[x]],sumtable$my.order), y = mean,fill=sumtable[[fill]]), stat = "identity", ...) +
    xlab(x) + labs(fill=fill) +
    theme(axis.text.x = element_text(angle = 90))

  if (errorbars==TRUE) {
    plot <- plot +
      geom_errorbar(aes(x = reorder(sumtable[[x]],sumtable$my.order),
                        ymax = sumtable$mean+sumtable$se,
                        ymin = sumtable$mean-sumtable$se,
                        position = position_dodge(width = 0.9),
                        width = 0.25)
      )
  }

  ?geom_errorbar

  return(plot)
}


#' Plot BarCharts from phyloseq object
#'
#' This is a wrapper function around the two functions sumTable and barChart for phyloseq objects.
#'
#' @param ps (Required) Phyloseq object
#' @param x (Optional, default="Genus") Character string referring to the main explainatory variable.
#' @param variables (optional, default=c("Order")) A vector of strings reffereng to any aditional explainatory variables. First string will be used for "fill".
#' @param errorbars (Optional, default=TRUE) Boolean operator. If True, standard error is plotted as errorbars.
#' @param ... Any aditional variables passed on to geom_bar(...)
#' @return A ggplot object.
#' @examples
#' data(ps_18S)
#' ps_18S %>%
#'  subset_samples(SORTED_genus=="Synchaeta") %>%
#'  subset_samples(MONTH=="aug") %>%
#'  subset_taxa(Class!="Rotifera") %>%
#'  transform_sample_counts(function(x) (x/sum(x))*100) %>%
#'  filter_taxa(function(x) sum(x>1)>=((length(x)*0.5)),TRUE) %>%
#'  phyloseqBarChart(x="Genus", variables=c("Order"))
#' @export

phyloseqBarChart <- function(ps, x="Genus", y= "Abudance", variables=c("Order"), errorbars=TRUE, ...) {

  plot <- ps %>%
    sumTable(x = x, y= y, variables=variables) %>%
    barChart(x = x, fill=variables[1], errorbars=errorbars, ...)
  return(plot)
}







#' Subset phyloseq object based on sample data
#'
#' Instead of using phyloseq::subset_samples.
#'
#' @param physeq (Required) Phyloseq object
#' @param variable (Required) Character string referring to variable name in Sample data.
#' @param value (Required) Character string, or vector of character string.
#' @return Philtered phyloseq
#' @examples
#' data(ps_18S)
#' ps_18S %>%
#'  select_samples("SORTED_genus","Keratella")
#' ps_18S %>%
#'  select_samples("SORTED_genus",c("Keratella", "Synchaeta"))
#' @export

select_samples <- function (physeq, variable, value) {

  if (is.null(sample_data(physeq))) {
    cat("Nothing subset. No sample_data in physeq.\n")
    return(physeq)

  } else {
    oldDF <- as(sample_data(physeq), "data.frame")
    newDF <- oldDF[which(oldDF[[variable]] %in% value),]

    if (class(physeq) == "sample_data") {
      return(sample_data(newDF))

    } else {
      sample_data(physeq) <- sample_data(newDF)
      return(physeq)
    }
  }
}





#' Filtering of prevalent species, multiple sample groups.
#'
#' This function filters away taxa only present in few samples.
#'
#'
#' @param ps (Required) A class S4 class phyloseq object.
#' @param variable (Required) Character string referring to variable name in Sample data. Subseting variable.
#' @param prevalence (Optional, default=0.75) A numeric value between 0:1. Least fraction of required presence.
#' @param minpercent (Optional, default=0.1) Numeric value between 0:100. A percentage of total library size. Lowest number to count as presence.
#' @return A class S4 class phyloseq object, filtered.
#' @examples
#' data(ps_16S)
#' ps_16S <- subset_samples(ps_16S, SORTED_genus == "Synchaeta")
#' ps_16S <- prevalentSpecies(ps_16S, prevalence=0.75, minpercent=1)
#' @export

filterPrevalentSpecies <- function(ps, variable, prevalence=0.6, minpercent=0) {

  require(metaphylo)
  require(phyloseq)
  require(magrittr)

  uniques <- ps %>%
    get_variable(variable) %>%
    base::unique()

  c <- 0
  for (x in uniques) {
    c <- c+1
    single.ps <- ps %>%
      select_samples(variable, x) %>%
      prevalentSpecies(prevalence=prevalence, minpercent=minpercent)
    if (c==1) {
      merged.ps <- single.ps
    } else {
      merged.ps <- merge_phyloseq(merged.ps, single.ps)
    }
  }
  filtered.ps <- merged.ps %>%
    taxa_names() %>%
    prune_taxa(ps)

  print(paste("Prevalence filtering based on",length(uniques),"sample groups:"))
  print(paste(uniques))
  return(filtered.ps)
}








#' LARGE WRAPPER: Summarizing RRA and FOC for phyloseq objet
#'
#' This function filters away taxa only present in few samples.
#' Constructs Relative read count tables and frequency of occurance tables
#'
#'
#' @param ps (Required) A class S4 class phyloseq object.
#' @param variable (Required) Variable name in Sample data. Subseting variable.
#' @param prevalence (Optional, default=0.6) A numeric value between 0:1. Least fraction of required presence.
#' @param rowname (Required) Only evaluated if PS_output=FALSE. This has to be the taxonomic rank with the highest available resolution. (Refer to phyloseq::tax_glom)
#' @return A list of two phyloseq objects, OR A list of two matrices. (See PS_output)
#' @examples
#' data(ps_18S)
#' ps_18S <- subset_samples(ps_18S, SORTED_genus == "Synchaeta") %>%
#' tax_glom("Family") %>%
#' summarizedPhyloseq(variable="MONTH")
#'
#' @export

summarizedPhyloseq <- function(ps = synchaeta_18S, variable = MONTH, rowname=Family, prevalence=0.6) {

  require(phyloseq)
  require(magrittr)
  require(metagMisc)

  # Passing values for tidy evaluation
  variable__ <- deparse(substitute(variable))
  rowname__ <- deparse(substitute(rowname))
  variable_ <- enquo(variable)
  rowname_ <- enquo(rowname)

  # First rare species are discarded
  filtered_ps <- ps %>%
    transform_sample_counts(function(x) x/sum(x)) %>% # Sanple counts are ransformed into
    filterPrevalentSpecies(variable = variable__, prevalence=prevalence, minpercent=0)

  # Next results are summarized acording to relative read abundance (mean).
  RRA_mx <- filtered_ps %>%
    psmelt() %>%
    group_by(!!variable_, !!rowname_) %>%
    summarise(Abundance = mean(Abundance)) %>%
    spread(key = variable__, value = Abundance) %>%
    column_to_rownames(rowname__)

  # The concept of frequency of occurance is defined.
  foc <- function(vector) {
    FOC <- vector %>%
      map_dbl(function(x) ifelse(x>0, 1, 0)) %>%
      mean()
    return(FOC)
  }

  # Phyloseq is summarized acording to frequency of occurance.
  FOC_mx <- filtered_ps %>%
    psmelt() %>%
    group_by(!!variable_, !!rowname_) %>%
    summarise(Abundance = foc(Abundance)) %>%
    spread(key = variable__, value = Abundance) %>%
    column_to_rownames(rowname__)

  # Outputs are combined to list object.
  result <- list(frequency_of_occurence = FOC_mx, relative_read_abundance = RRA_mx)



  # Finaly add a taxonomic for "others", i.e taxa that has been filtered out.
  x1 <- t(base::as.matrix(1-colSums(result$relative_read_abundance)))
  rownames(x1) <- "OTHERS"
  result$relative_read_abundance <- rbind(result$relative_read_abundance, x1)

  x2 <- matrix(1:1, ncol=ncol(x1))
  rownames(x2) <- "OTHERS"
  colnames(x2) <- colnames(x1)
  result$frequency_of_occurence <- rbind(result$frequency_of_occurence, x2)

  return(result)
}





#' LARGE WRAPPER: Pipartite plot of FOC and RRA data
#'
#' Produces bipartite plot out of summarized phyloseqobject (output of summarizedPhyloseq function)
#' Thickness of interaction represents Relative read abundance. Color of interaction represents Frequency of occurance.
#'
#'
#' @param dualmat (Required) The output of summarizedPhyloseq (i.e. a list of two matrixes.)
#' @param col.order (Optional) Vector of strings, representing column order.
#' @param row.order (Optional) Vector of strings, representing row order.
#' @param foc.color (Optional, c("gray35", "gray40", "gray45", "gray50", "gray90")).
#' @param ... (Optional). Refer to bipartite::plotweb options.
#' @return A bipartite::plotweb plot.
#' @examples
#' data(ps_18S)
#' ps_18S <- subset_samples(ps_18S, SORTED_genus == "Synchaeta") %>%
#' tax_glom("Family") %>%
#' summarizedPhyloseq(variable="MONTH")
#'
#' @export

colorWeb <- function(dualmat,
                     foc.color = c("gray35", "gray40", "gray45", "gray50", "gray90"),
                     method = "normal",
                     text.rot = -90,
                     labsize =0.5,
                     bor.col.interaction = NULL,
                     ...) {

  require(bipartite)

  col.interaction <- c()
  c <- 0
  for (row in 1:nrow(mat$frequency_of_occurence)) {
    for (col in 1:ncol(mat$relative_read_abundance)) {
      c <- c+1
      if (mat$frequency_of_occurence[row,col]>=0.9) {
        col.interaction[c] <- foc.color[1]
      } else {
        if (mat$frequency_of_occurence[row,col]>=0.8) {
          col.interaction[c] <- foc.color[2]
        } else {
          if (mat$frequency_of_occurence[row,col]>=0.7) {
            col.interaction[c] <- foc.color[3]
          } else {
            if (mat$frequency_of_occurence[row,col]>=0.6) {
              col.interaction[c] <- foc.color[4]
            }else{
              col.interaction[c] <- foc.color[5]
            }
          }
        }
      }
    }
  } # End of colmat

  plotweb(dualmat$relative_read_abundance,
          method = method,
          text.rot = text.rot,
          labsize = labsize,
          col.interaction = col.interaction,
          bor.col.interaction = NULL)

} # End of color-web
