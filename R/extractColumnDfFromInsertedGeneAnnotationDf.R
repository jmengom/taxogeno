#' Extract Semicolon-Separated Column into 1:N Relationship
#'
#' Converts a column containing semicolon-separated values into a normalized
#' 1:N relationship data frame. Each unique value gets its own row linked
#' to the original gene ID.
#'
#' @param insertedGeneAnnotationDf Data frame containing gene annotations.
#'   Must have a "geneid" column.
#' @param colName Character string. Name of the column to expand.
#' @param newColName Character string. Name for the expanded column.
#'   Defaults to \code{colName}.
#'
#' @return A data frame with two columns:
#'   \item{geneid}{The gene identifier (repeated for each value)}
#'   \item{newColName}{The extracted values (one per row)}
#'
#' @examples
#' \dontrun{
#' # Input: geneid | enzyme
#' #        1      | EC1;EC2;EC3
#' # Output: geneid | ec
#' #         1      | EC1
#' #         1      | EC2
#' #         1      | EC3
#' enzyme_df <- extractColumnDfFromInsertedGeneAnnotationDf(annot_df, "enzyme", "ec")
#' }
#'
#' @export
extractColumnDfFromInsertedGeneAnnotationDf<-function(insertedGeneAnnotationDf,colName,newColName=colName){
    if(length(colName)>1){
        stop("length(colName)>1")
    }
    if(length(newColName)>1){
        stop("length(newColName)>1")
    }
    colVecList<-strsplit(insertedGeneAnnotationDf[,colName],";")

    geneidVec<-rep(insertedGeneAnnotationDf[,"geneid"],lengths(colVecList))
    colVec<-unlist(colVecList) ; rm(colVecList)
    outputDf<-data.frame(geneidVec,colVec,stringsAsFactors=FALSE) ; rm(geneidVec,colVec)
    colnames(outputDf)<-c("geneid",newColName)
    
    outputDf
}
