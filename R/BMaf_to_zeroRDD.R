#' Computing zero-valued Root Distance differences from allele frequencies
#'
#' This function computes zero-valued Root Distance (RD) differences from allele frequencies.
#'
#'\code{BMaf_to_zeroRDD()} calculate the smallest number (noted as \eqn{H}) of the observed Root Distance (RD) differences that must be zero in order for it to agree to a tree-topology (If we observed infinite data from the model without error, we would have exactly \eqn{H} of the estimated RD differences equal to zero. In a finite sample, however, these will not be identically zero, but rather will be close to zero.)
#'
#' @param trans_mat_allele_freq A \eqn{(P+1) \times L} matrix, containing allele frequencies of \eqn{P} taxa and an outgroup for \eqn{L} loci.
#' @param use Specify which part of data is used to compute the covariance matrix. Details can be checked in \code{stats::cov()} in R. The options are "\code{complete.obs}", "\code{pairwise.complete.obs}", "\code{everything}", "\code{all.obs}", and "\code{na.or.complete}".
#'
#' @return An array where zero valued RD differences have zero values, and the rest have a value of \eqn{1}. Note that the array is \eqn{P \times P \times P}, and for each \eqn{i}, \code{array_zero_ID[i,,]} is symmetric.
#' @noRd
#'
#' @import stats
#'
#' @examples
#'
#' # load example data from rapidphylo package
#' data("Human_Allele_Frequencies")
#' mat_allele_freq <- Human_Allele_Frequencies
#' # perform logistic transformation
#' mat_allele_freq[mat_allele_freq==1] <- 0.99
#' mat_allele_freq[mat_allele_freq==0] <- 0.01
#' trans_mat_allele_freq <- log(mat_allele_freq/(1-mat_allele_freq))
#' # convert type of object into data frame
#' trans_mat_allele_freq <- as.data.frame(trans_mat_allele_freq)
#' outgroup <- 'Han'
#' names<-row.names(trans_mat_allele_freq)
#' # use the population names as the row names of your transformed allele frequency matrix")
#' if (is.character(outgroup)){
#'   index<-which(names==outgroup)
#'   }else {
#'     index<-outgroup
#'     }
#' trans_mat_allele_freq<-rbind(trans_mat_allele_freq[-index, ],trans_mat_allele_freq[index,])
#' label<-row.names(trans_mat_allele_freq)
#' # run BMaf_to_zeroRDD function
#' array_zero_ID<-BMaf_to_zeroRDD(trans_mat_allele_freq,use="pairwise.complete.obs")
#' array_zero_ID[1:14,1:14,1]
#'
BMaf_to_zeroRDD<-function(trans_mat_allele_freq,
                          use=c("complete.obs","pairwise.complete.obs","everything","all.obs","na.or.complete")){

  P <- nrow(trans_mat_allele_freq) - 1

  H <- P*(P-1)*(P-2)/6
  ## H is the smallest number of  RD differences that have to be zero in order for it
  ## to agree to a tree-topology

  use<-match.arg(use)
  mat_normalized <- t(trans_mat_allele_freq)
  mat_sigma <- cov(mat_normalized,use = use)
  for (x in 1:(P+1))
    for (y in 1:(P+1))
      if (mat_sigma[x,y] < 0)
        mat_sigma[x,y] <- 0

  #######################################################################

  raw_array_ID <- array(dim=c(P,P,P),0)

  n_list_ID <- P*(P-1)*(P-2)/2

  list_ID <- matrix(nrow=n_list_ID,ncol=4)
  ## list_ID lists the relevant absolute IDs from array_ID

  array_zero_ID <- array(dim=dim(raw_array_ID),1)

  i_list_ID <- 0


  for (i in 1:P)
    for (j in 1:P)
      for (k in 1:P)
      {
        raw_array_ID[i,j,k] <- abs(mat_sigma[i,j] - mat_sigma[i,k])

        if ( ((i != j) && (i != k)) && (j < k))
        {
          i_list_ID <- i_list_ID + 1

          list_ID[i_list_ID,] <- c(raw_array_ID[i,j,k],i,j,k)
        }


        if ( ((i < j) && (i < k)) && (j < k))
        {

          maxrd <- which.max(c(mat_sigma[i,j],mat_sigma[i,k],mat_sigma[j,k]))

          if (maxrd == 1)
          {
            array_zero_ID[k,i,j] <- 0
            array_zero_ID[k,j,i] <- 0

          }
          if (maxrd == 2)
          {
            array_zero_ID[j,i,k] <- 0
            array_zero_ID[j,k,i] <- 0

          }
          if (maxrd == 3)
          {
            array_zero_ID[i,j,k] <- 0
            array_zero_ID[i,k,j] <- 0

          }

        }

      }

  ## array_zero_ID
  ## This is an array where zero valued RD differences have a value zero, and the rest have a value of 1
  ## Note that the array is P X P X P, and for each i, array_zero_ID[i,,] is symmetric

  return(array_zero_ID)

}





