#' Get Normalized Annotation for Proteomes
#'
#' Retrieves normalized annotation vectors (GO Slim or Keywords) for
#' a set of proteomes from the database. The normalized value is
#' calculated as count / genecount.
#'
#' @param dbConn A DBI database connection object.
#' @param annotationType Character string. Must be "generated_goslim" or "keyword".
#' @param proteomeIdVec Integer vector. Proteome IDs to retrieve.
#'
#' @return A data frame with columns:
#'   \item{proteomeid}{Proteome identifier}
#'   \item{annotkwid}{Annotation term ID}
#'   \item{annotkwcount}{Raw annotation count}
#'   \item{annotkwrelval}{Normalized value (count/genecount)}
#'
#' @importFrom DBI dbGetQuery
#'
#' @examples
#' \dontrun{
#' norm_df <- getNormalizedAnnotationForProteomeIdVec(
#'   conn,
#'   "generated_goslim",
#'   c(1, 2, 3)
#' )
#' }
#'
#' @export
getNormalizedAnnotationForProteomeIdVec<-function(dbConn,annotationType,proteomeIdVec){
    if(length(annotationType)>1){
        stop("length(annotationType)>1")
    }

    normalizedAnnotationDf<-data.frame(proteomeid=character(),annotkwid=character(),annotkwrelval=numeric())
    if(annotationType=="generated_goslim"){
        normalizedAnnotationDf<-dbGetQuery(
           dbConn,
           "SELECT
              taxogeno.proteome.proteomeid as proteomeid,
              annotkwid                     as annotkwid,
              annotkwcount                  as annotkwcount,
              annotkwcount::float/genecount as annotkwrelval

            FROM taxogeno.generated_goslim_summary
            INNER JOIN taxogeno.proteome
              ON taxogeno.generated_goslim_summary.proteomeid=taxogeno.proteome.proteomeid
              AND taxogeno.generated_goslim_summary.proteomeid IN ($1)

            ORDER BY taxogeno.proteome.proteomeid, annotkwid",
           params=list(proteomeIdVec)
        )
    } else if(annotationType=="keyword"){
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

        proteomeIdKeywordCountDf<-dbGetQuery(
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

        summaryDf<-merge(proteomeIdKeywordIdCartesianDf,proteomeIdKeywordCountDf,by=c("proteomeid","keyword"),all.x=TRUE)
        summaryDf[,"annotkwrelval"]<-summaryDf[,"annotkwcount"]/summaryDf[,"genecount"]
        
        normalizedAnnotationDf<-summaryDf[,c("proteomeid","annotkwid","annotkwcount","annotkwrelval")]
        rm(summaryDf)
    }
    normalizedAnnotationDf
}

