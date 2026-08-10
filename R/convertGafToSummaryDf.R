#' Convert GAF Data Frame to Summary Data Frame
#'
#' Aggregates a GAF data frame to count annotations per term and adds
#' the proteome ID as a reference.
#'
#' @param gafDf A GAF-format data frame.
#' @param proteomeId Integer. The proteome ID to associate with the summary.
#'
#' @return A data frame with columns:
#'   \item{annotkwid}{Annotation term ID}
#'   \item{annotkwcount}{Count of annotations for this term}
#'   \item{proteomeid}{The proteome ID}
#'
#' @examples
#' \dontrun{
#' summary_df <- convertGafToSummaryDf(gaf_df, proteome_id = 1)
#' }
#'
#' @export
convertGafToSummaryDf<-function(gafDf,proteomeId){
    if(length(proteomeId)>1){
        stop("length(proteomeId)>1")
    }
    summaryDf<-aggregate(dbobjectid~annotkwid,data=gafDf,length)
    summaryDf[,"proteomeid"]<-rep(proteomeId,nrow(summaryDf))
    colnames(summaryDf)[colnames(summaryDf)=="dbobjectid"]<-"annotkwcount"   
    summaryDf
}
