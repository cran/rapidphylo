#' Estimating tree-topology from allele frequency data
#'
#' \code{RDM()} estimates a tree-topology from allele frequencies.
#'
#' The input matrix is the observed values of the frequencies at tips \eqn{1, 2, ..., P, P+1}. A logit transformation is performed on the allele frequency data, so that the observed values are approximately normal (The logit transformation of r refers to \eqn{\log\frac{r}{1-r}}). The transformed matrix is converted into a data frame for further analyses.
#'
#' @param outgroup A variable that can be either the population name or a numerical row number of the outgroup data.
#' @param use Specify which part of data is used to compute the covariance matrix. The options are "\code{complete.obs}", "\code{pairwise.complete.obs}", "\code{everything}", "\code{all.obs}", and "\code{na.or.complete}". See \code{stats::cov} for more details.
#' @param mat_allele_freq A \eqn{(P+1) \times L} matrix containing the allele frequencies, where there are \eqn{P} taxa, plus one outgroup, and \eqn{L} loci.
#'
#' @return An estimated tree-topology in Newick format.
#' @export
#'
#' @references Jing Peng, Haseena Rajeevan, Laura Kubatko, and Arindam RoyChoudhury (2021) \emph{A fast likelihood approach for estimation of large phylogenies from continuous trait data}. Molecular Phylogenetics and Evolution 161 107142.
#'
#' @examples
#' # A dataset "Human_Allele_Frequencies" is loaded with the package;
#' # it has allele frequencies in 44,000 sites for
#' # 5 human populations and one outgroup human population.
#'
#' # check data dimension
#' dim(Human_Allele_Frequencies)
#'
#' # run RDM function
#' rd_tre <- RDM(Human_Allele_Frequencies, outgroup = "San", use = "pairwise.complete.obs")
#'
#' # result visualization
#' plot(rd_tre, use.edge.length = FALSE, cex = 0.5)
#'
RDM<-function(mat_allele_freq,outgroup,
              use=c("complete.obs","pairwise.complete.obs","everything","all.obs","na.or.complete")){

  mat_allele_freq[mat_allele_freq==1]<-0.99
  mat_allele_freq[mat_allele_freq==0]<-0.01
  trans_mat_allele_freq<-log(mat_allele_freq/(1-mat_allele_freq))
  trans_mat_allele_freq <- as.data.frame(trans_mat_allele_freq)

  use<-match.arg(use)
  names<-row.names(trans_mat_allele_freq)
  if (is.character(outgroup)){
    index<-which(names==outgroup)
  }else {
    index<-outgroup
  }
  trans_mat_allele_freq<-rbind(trans_mat_allele_freq[-index, ],trans_mat_allele_freq[index,])
  label<-row.names(trans_mat_allele_freq)
  array_zero_ID<-BMaf_to_zeroRDD(trans_mat_allele_freq,use=use)
  base_tree<-zeroRDD_to_splits(array_zero_ID)
  rd_tre<-topology_to_newick(base_tree,label=label)
  return(rd_tre)
}
