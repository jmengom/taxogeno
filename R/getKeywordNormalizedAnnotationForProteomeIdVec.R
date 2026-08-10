#' Get Keyword Normalized Annotation for Proteomes
#'
#' Retrieves normalized keyword annotation vectors for a set of proteomes.
#' Uses a Cartesian product of proteomes and all available keywords.
#'
#' @param dbConn A DBI database connection object.
#' @param proteomeIdVec Integer vector. Proteome IDs to retrieve.
#'
#' @return A data frame with columns:
#'   \item{proteomeid}{Proteome identifier}
#'   \item{annotkwid}{Keyword ID}
#'   \item{annotkwrelval}{Normalized value (count/genecount)}
#'
#' @importFrom DBI dbGetQuery
#'
#' @examples
#' \dontrun{
#' norm_df <- getKeywordNormalizedAnnotationForProteomeIdVec(conn, c(1, 2, 3))
#' }
#'
#' @export
getKeywordNormalizedAnnotationForProteomeIdVec<-function(dbConn,proteomeIdVec){
    proteomeIdKeywordIdCartesianDf<-dbGetQuery(
        dbConn,
        "select proteomeid as proteomeid,
                 genecount as genecount,
                   keyword as keyword,
                 keywordid as annotkwid
           from taxogeno.proteome,
                uniprot.uniprot_keyword
          where proteomeid in ($1)",
        params=list(proteomeIdVec))

    proteomeIdkeywordCountDf<-dbGetQuery(
        dbConn,
        "select tg.proteomeid      as proteomeid,
                tgk.keyword        as keyword,
                count(tgk.keyword) as annotkwcount
         from taxogeno.gene tg
         inner join taxogeno.gene_keyword_rel tgk
	 on tg.geneid=tgk.geneid and tg.proteomeid in ($1)
         where tg.proteomeid in ($1)
         group by tg.proteomeid,tgk.keyword",
        params=list(proteomeIdVec))

    summaryDf<-merge(proteomeIdKeywordIdCartesianDf,proteomeIdkeywordCountDf,by=c("proteomeid","keyword"),all.x=TRUE)
    summaryDf[,"annotkwrelval"]<-summaryDf[,"annotkwcount"]/summaryDf[,"genecount"]
    summaryDf[,c("proteomeid","annotkwid","annotkwrelval")]
}

