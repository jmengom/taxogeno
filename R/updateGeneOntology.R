#' Update Gene Ontology Tables
#'
#' Downloads the latest Gene Ontology OWL files and parses them into
#' the database using PostgreSQL's XMLTABLE functionality.
#'
#' @param dbConn A DBI database connection object.
#'
#' @return NULL (invisible). The function updates the database directly.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Downloads \code{go.owl} from Gene Ontology
#'   \item Downloads \code{gislim_generic.owl} from Gene Ontology
#'   \item Inserts OWL files into \code{gene_ontology.gene_ontology_xml}
#'   \item Parses OWL XML using \code{XMLTABLE} to extract:
#'     \itemize{
#'       \item GO ID (\code{goid})
#'       \item GO Label (\code{golabel})
#'       \item GO Aspect (\code{goaspect})
#'     }
#'   \item Updates \code{gene_ontology.gene_ontology}
#'   \item Updates \code{gene_ontology.goslim_generic}
#' }
#'
#' @details
#' The parsing uses XMLTABLE with appropriate XML namespaces:
#' \itemize{
#'   \item \code{rdf}: \url{http://www.w3.org/1999/02/22-rdf-syntax-ns#}
#'   \item \code{owl}: \url{http://www.w3.org/2002/07/owl#}
#'   \item \code{oboInOwl}: \url{http://www.geneontology.org/formats/oboInOwl#}
#'   \item \code{rdfs}: \url{http://www.w3.org/2000/01/rdf-schema#}
#' }
#'
#' @importFrom DBI dbExecute
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' updateGeneOntology(conn)
#' }
#'
#' @export
updateGeneOntology<-function(dbConn){
    insertGoOwl<-function(dbConn,owlFileName){
        dbExecute(dbConn,"DELETE from gene_ontology.gene_ontology_xml WHERE filename=$1",params=list(owlFileName))
        owlContents <- rawToChar(readBin(owlFileName, "raw", file.info(owlFileName)$size))
        dbExecute(dbConn,"")
        dbExecute(dbConn,
                  "insert into gene_ontology.gene_ontology_xml (filename, xmldata) values($1,XMLPARSE(DOCUMENT $2))",
                  params=list(owlFileName,owlContents))
        rm(owlContents)
    }
    deleteTablifiedGoOwl<-function(dbConn,tableName){
        dbExecute(dbConn,sprintf("DELETE from gene_ontology.%s",tableName))
    }
    insertTablifiedGoOwl<-function(dbConn,tableName,owlFileName){
        dbExecute(dbConn,
                  sprintf("insert into gene_ontology.%s(goid,golabel,goaspect)
                       select
                         gene_ontology_tablified.goid     as goid,
                         gene_ontology_tablified.golabel  as golabel,
                         gene_ontology_tablified.goaspect as goaspect
                       --
                        from
                          -- table with xml source column inside
                          gene_ontology.gene_ontology_xml,
                          -- generated table from xml source
                          xmltable(
                            -- namespaces
                            xmlnamespaces(
                              'http://www.w3.org/1999/02/22-rdf-syntax-ns#'   as \"rdf\",
                              'http://www.w3.org/2002/07/owl#'                as \"owl\",
                              'http://www.geneontology.org/formats/oboInOwl#' as \"oboInOwl\",
                              'http://www.w3.org/2000/01/rdf-schema#'         as \"rdfs\"
                            ),
                            -- path
                            '/rdf:RDF/owl:Class'
                            -- xml source
                            passing
                              gene_ontology.gene_ontology_xml.xmldata
                            -- columns
                            columns
                              goid     text PATH 'oboInOwl:id/text()',
                              golabel  text PATH 'rdfs:label/text()',
                              goaspect text PATH 'oboInOwl:hasOBONamespace/text()'
                          ) gene_ontology_tablified -- generated table alias
                       --
                        where
                         gene_ontology.gene_ontology_xml.filename=$1
                        and
                         gene_ontology_tablified.goid is not null",
                       tableName),
                  params=list(owlFileName))
    }

    ## Descargar los owl
    
    goOwlUrl<-"http://current.geneontology.org/ontology/go.owl"
    goOwlFileName<-"go.owl"
    system2("wget", args=c(goOwlUrl,"-O",goOwlFileName), wait=TRUE)

    goSlimGenericOwlUrl<-"http://current.geneontology.org/ontology/subsets/goslim_generic.owl"
    goSlimGenericOwlFileName<-"goslim_generic.owl"
    system2("wget", args=c(goSlimGenericOwlUrl,"-O",goSlimGenericOwlFileName), wait=TRUE)

    ## ## ## ## ##

    insertGoOwl(dbConn,"go.owl")
    insertGoOwl(dbConn,"goslim_generic.owl")
    
    deleteTablifiedGoOwl(dbConn,"goslim_generic")
    deleteTablifiedGoOwl(dbConn,"gene_ontology")

    insertTablifiedGoOwl(dbConn,"gene_ontology","go.owl")
    insertTablifiedGoOwl(dbConn,"goslim_generic","goslim_generic.owl")
}
