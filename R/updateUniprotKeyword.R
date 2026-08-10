#' Update UniProt Keywords
#'
#' Downloads the latest UniProt keywords in both OBO and TSV formats
#' and updates the database tables.
#'
#' @param dbConn A DBI database connection object.
#'
#' @return NULL (invisible). The function updates the database directly.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Downloads \code{keyword.obo} from UniProt
#'   \item Downloads \code{keyword.tsv} from UniProt
#'   \item Reads the TSV file
#'   \item Updates \code{uniprot.uniprot_keyword_category}
#'   \item Updates \code{uniprot.uniprot_keyword}
#' }
#'
#' @importFrom DBI dbExecute dbWriteTable SQL
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' updateUniprotKeyword(conn)
#' }
#'
#' @export
updateUniprotKeyword<-function(dbConn){
    ## semiconstants
    uniprotKeywordOboUrl<-"'http://www.uniprot.org/keywords/?query=*&format=obo'"
    uniprotKeywordOboFilePath<-"keyword.obo"
    system2("wget", args=c(uniprotKeywordOboUrl,"-O",uniprotKeywordOboFilePath), wait=TRUE)

    uniprotKeywordTsvUrl<-"'http://www.uniprot.org/keywords/?query=*&format=tab'"
    uniprotKeywordColNames<-c("keywordid","keyword","keyworddescription","keywordcategory")
    uniprotKeywordTsvFilePath<-"uniprot_keyword.tsv"
    system2("wget", args=c(uniprotKeywordTsvUrl,"-O",uniprotKeywordTsvFilePath), wait=TRUE)
    uniprotKeywordDf<-read.table(uniprotKeywordTsvFilePath,
                                 header=TRUE,
                                 row.names=NULL,
                                 col.names=uniprotKeywordColNames,
                                 quote="",
                                 sep="\t",
                                 stringsAsFactors=FALSE)
    unlink(uniprotKeywordTsvFilePath)
    dbExecute(dbConn,"DELETE FROM uniprot.uniprot_keyword")
    dbExecute(dbConn,"DELETE FROM uniprot.uniprot_keyword_category")
    dbWriteTable(dbConn,SQL("uniprot.uniprot_keyword_category"), unique(uniprotKeywordDf[,"keywordcategory",drop=FALSE]), append=TRUE,row.names=FALSE)
    dbWriteTable(dbConn,SQL("uniprot.uniprot_keyword"), uniprotKeywordDf,append=TRUE,row.names=FALSE)
}
