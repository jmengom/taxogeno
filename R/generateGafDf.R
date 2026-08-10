#' Generate Gene Association Format (GAF) Data Frame
#'
#' Creates a GAF data frame from ontology annotations, following the
#' Gene Ontology Consortium's GAF 2.0 specification.
#'
#' @param dbConn A DBI database connection object (used for future extensions).
#' @param ontologyAnnotDf Data frame with columns:
#'   \item{geneid}{Gene identifier}
#'   \item{annotkwid}{Annotation term ID}
#'   \item{aspect}{GO aspect (C, F, or P)}
#'
#' @return A data frame in GAF 2.0 format with 17 columns:
#'   \item{db}{Database source ("taxogeno")}
#'   \item{dbobjectid}{Gene ID}
#'   \item{dbobjectsymbol}{Gene symbol (same as ID)}
#'   \item{qualifier}{Qualifier (NA)}
#'   \item{annotkwid}{GO term ID}
#'   \item{dbreference}{Reference (NA)}
#'   \item{evidencecode}{Evidence code ("IEA")}
#'   \item{with}{With/from (NA)}
#'   \item{aspect}{GO aspect}
#'   \item{dbobjectname}{Object name (NA)}
#'   \item{dbobjectsynonym}{Synonyms (NA)}
#'   \item{dbobjecttype}{Object type ("protein")}
#'   \item{taxon}{Taxon (NA)}
#'   \item{date}{Current timestamp in ISO format}
#'   \item{assignedby}{Source ("sma3s")}
#'   \item{annotationextension}{Extensions (NA)}
#'   \item{geneproductformid}{Gene product form ID (NA)}
#'
#' @examples
#' \dontrun{
#' annot_df <- data.frame(
#'   geneid = c(1, 2),
#'   annotkwid = c("GO:0003674", "GO:0005575"),
#'   aspect = c("F", "C")
#' )
#' gaf_df <- generateGafDf(conn, annot_df)
#' }
#'
#' @export
generateGafDf<- function(dbConn,ontologyAnnotDf) {
    gafDf<-data.frame(
        "db"                 =rep("taxogeno", nrow(ontologyAnnotDf)),
        "dbobjectid"         =ontologyAnnotDf[,"geneid"],
        "dbobjectsymbol"     =ontologyAnnotDf[,"geneid"],
        "qualifier"          =rep(NA,nrow(ontologyAnnotDf)),
        "annotkwid"         =ontologyAnnotDf[,"annotkwid"],
        "dbreference"        =rep(NA,nrow(ontologyAnnotDf)),
        "evidencecode"       =rep("IEA",nrow(ontologyAnnotDf)),
        "with"               =rep(NA,nrow(ontologyAnnotDf)),
        "aspect"             =ontologyAnnotDf[,"aspect"],
        "dbobjectname"       =rep(NA,nrow(ontologyAnnotDf)),
        "dbobjectsynonym"    =rep(NA,nrow(ontologyAnnotDf)),
        "dbobjecttype"       =rep("protein",nrow(ontologyAnnotDf)),
        "taxon"              =rep(NA,nrow(ontologyAnnotDf)),
        "date"               =rep( strftime( as.POSIXlt(Sys.time()), "%Y-%m-%dT%H:%M:%S%z"), nrow(ontologyAnnotDf) ),
        "assignedby"         =rep("sma3s",nrow(ontologyAnnotDf)),
        "annotationextension"=rep(NA,nrow(ontologyAnnotDf)),
        "geneproductformid"  =rep(NA,nrow(ontologyAnnotDf))
    )
    gafDf
}
