#' Compiling component subtrees into an object of class ‘multiPhylo’
#'
#' \code{multiphylo()} is used in function \code{topology_to_newick()}. It complies component subtrees and transform them into an object of class ‘multiPhylo’.
#'
#' @param base_tree A matrix, where each row represents a split, but without repetitions. In a given row, a k-way split is represented by assigning each set of taxa that is descendent of a split as an identifying number from 1 to k, and numbering the position corresponding to each taxa that is part of such a set by that identifying number. The other positions are set to zero.
#' @param label A character vector containing the names of the taxa.
#'
#' @import stats
#' @import phangorn
#' @import ape
#'
#' @return An object of class ‘multiPhylo’.
#' @noRd
#'
#' @examples NA
#'
multiphylo<-function(base_tree,label){

  outgroup<-label[length(label)]

  for (i in 1:dim(base_tree)[1]){
    ind3<-which(base_tree[i,]==0)
    ind4<-which(base_tree[i,]==1)
    ind5<-which(base_tree[i,]==2)
    if (length(ind3)==0){
      grp1<-label[ind4]
      if (length(ind4)>1){
        grp1<-paste("(",paste(grp1,collapse = ","),")",sep = "")
      }
      grp2<-label[ind5]
      if (length(ind5)>1){
        grp2<-paste("(",paste(grp2,collapse = ","),")",sep = "")
      }
      grp<-paste("(",outgroup,",",grp1,",",grp2,");",sep="")
    }else{
      grp1<-label[ind4]
      if (length(ind4)>1){
        grp1<-paste("(",paste(grp1,collapse = ","),")",sep = "")
      }
      grp2<-label[ind5]
      if (length(ind5)>1){
        grp2<-paste("(",paste(grp2,collapse = ","),")",sep = "")
      }
      grp3<-label[ind3]
      if (length(ind3)>1){
        grp3<-paste("(",paste(grp3,collapse = ","),")",sep = "")
      }
      grp<-paste("(",outgroup,",",grp1,",",grp2,",",grp3,");",sep="")
    }
    tree<-read.tree(text = grp)
    mp<-list(tree)
    class(mp)<- "multiPhylo"
    if (i==1) mp2<-mp else mp2<-c(mp2,mp)
  }
  return(mp2)
}
