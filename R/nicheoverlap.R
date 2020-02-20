#' Schoeners Index of Nieche Overlap
#'
#' Calculates the percentage overlap of diet between two species in a food web matrix.
#'
#' @param mat (Required) A food web matrix, with consurs as columns and prey as rows. Relative proportion.
#' @param species1 (Required) Character string, Colname of Species 1
#' @param species2 (Required) Character string, Colname of Species 2
#' @return Double, percentage of Diet Overlap.
#'
#' @export


nicheOverlap <- function(mat, species1, species2) {

  require(rlang)
  require(dplyr)

  # Error message 1: Normalized input
  if(any(colSums(mat)<0.99|colSums(mat)>1)) stop(
    'Require input normalized to proportion, columns has to sum up to 1')
  if(any(mat < 0 | mat > 1) ) stop(
    'Require input normalized to proportion, values between 0 and 1')

  species1 <- sym(species1)
  species2 <- sym(species2)

  alfa <- mat %>%
    as.data.frame() %>%
    rownames_to_column(var="Species") %>%
    mutate(diff=abs(!!species1-!!species2)) %>%
    pull(var=diff) %>%
    sum()*-0.5+1

  return(alfa*100)
} # End of nicheOverlap function




#' Schoeners Pairwise Index of Nieche Overlap
#'
#' Calculates the percentage overlap of diet between ALL pairs of species in a food web matrix.
#'
#' @param mat (Required) A food web matrix, with consurs as columns and prey as rows. Relative proportion.
#' @return A matrix, of type distance matrix, percentage of Diet Overlap.
#'
#' @export


pairwiseNicheOverlap <- function(mat) {

  # 2. Create a matrix for the pairwise comparaisons based on colnames in matrix
  overlap_mat <- matrix(nrow = ncol(mat) , ncol = ncol(mat))
  colnames(overlap_mat) <- colnames(mat)
  rownames(overlap_mat) <- colnames(mat)


  # 3. Fill the matrix with paired nieche overlap indexes
  for (x in 1:ncol(mat)) {
    for (y in 1:ncol(mat)) {

      j <- colnames(mat)[x]
      k <- colnames(mat)[y]
      a <- nicheOverlap(mat=mat, species1=j, species2 =k)

      overlap_mat[j,k] <- a
    } # End of x loop
  } # End of y loop


  # Error message 3: 100% Overlap between the same species
  for (x in colnames(overlap_mat)) {
    if (overlap_mat[x,x]!=100) stop('Output: Same species do not sum up to 100')
  }# End of error 3


  return(overlap_mat)
} # End of pairwiseNicheOverlap loop


#' Schoeners Index of Nieche Overlap, Per Sample.
#'
#' Calculates the percentage overlap of diet between ALL possible pairs of samples in a phyloseq object,
#'
#' @param ps (Required) A S4 type Phyloseq object. OBS currently tax-glomed to Family level.
#' @param sample_group (Optional, defaut="SORTED_genus). A character string, grouping variable.
#' @return A dataframe, with
#'
#' @export

nicheOverlapPerSample <- function(ps, sample_group="SORTED_genus") {


  twoSpeciesNiecheOverlap <-function(mat1, mat2) {

    b <- c()
    for (x in 2:ncol(mat1)) {
      for (y in 2:ncol(mat2)) {
        a <- select(mat1, Family, colnames(mat1)[x]) %>%
          left_join(select(mat2, Family, colnames(mat2)[y]), by="Family") %>%
          column_to_rownames("Family") %>%
          nicheOverlap(species1 = colnames(mat1)[x], species2 = colnames(mat2)[y])
        b<- c(b,a)
      }
    }
    return(b)
  }





  group_var <- sym(sample_group)

  groups <- ps %>%
    get_variable(variable) %>%
    unique()

  new_ps <- ps %>%
    transform_sample_counts(function(x)x/sum(x)) %>%
    psmelt()


  # Construct list of sample matrixes, on group basis
  matlist <- list()

  for (z in 1:length(groups)) {

    matlist[[groups[z]]]<- new_ps %>%
      filter(!!group_var == groups[z]) %>%
      select(Sample, Abundance, Family) %>%
      spread(key=Sample, value=Abundance)

  }



  nichdata <- data.frame()
  for (m in 1:length(groups)) {
    for (n in 1:length(groups)) {

      if(n==m) next

      mat1 <- matlist[[groups[m]]]
      mat2 <- matlist[[groups[n]]]

      testgroup <- paste(c(groups[m], "VS", groups[n]), collapse = "")
      values <- twoSpeciesNiecheOverlap(mat1, mat2)

      data <- data.frame(NicheOverlap =as.double(values),
                         Testgroup = as.character(rep(testgroup, length(values))))
      nichdata <- bind_rows(nichdata, data)
    }
  }
  return(nichdata)
}
