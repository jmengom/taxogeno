#' Map GO Terms to GO Slim Using OWLTools
#'
#' Runs the OWLTools map2slim command to map specific GO terms to generic
#' GO Slim terms. This is a crucial step for creating a simplified,
#' high-level view of functional annotations.
#'
#' @param gafDf A GAF-format data frame (from \code{\link{generateGafDf}}).
#' @param owl Character string. Path to the OWL ontology file.
#' @param subset Character string. Optional subset name (e.g., "goslim_generic").
#'
#' @return A GAF-format data frame with mapped GO Slim terms.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Writes the GAF data to a temporary file
#'   \item Executes OWLTools with the map2slim command
#'   \item Reads the output from the temporary file
#'   \item Cleans up temporary files
#' }
#'
#' @examples
#' \dontrun{
#' mapped_gaf <- owltoolsMap2SlimGafDf(
#'   gaf_df,
#'   owl = "goslim_generic.owl",
#'   subset = "goslim_generic"
#' )
#' }
#'
#' @export
owltoolsMap2SlimGafDf<-function(gafDf,owl,subset=NA){
    gafNamesVec<-c("db","dbobjectid","dbobjectsymbol","qualifier",
                   "annotkwid","dbreference","evidencecode","with",
                   "aspect","dbobjectname","bobjectsynonym","dbobjecttype",
                   "proteomeid","date","assignedby","annotationextension",
                   "geneproductformid")
    inputGafFilePath<- tempfile()
    write.table(gafDf, file = inputGafFilePath, append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "", dec = ".", row.names = FALSE,
                col.names = FALSE)
    
    outputGafFilePath<- tempfile()

    owltoolsArgs<-if(!is.na(subset))
                       c(owl,
                        "--gaf",inputGafFilePath,
                        "--map2slim",
                        "--subset",subset,
                        "--write-gaf",outputGafFilePath)
                   else
                       c(owl,
                        "--gaf",inputGafFilePath,
                        "--map2slim",
                        "--write-gaf",outputGafFilePath)
    
    map2slimGafStdOut<-system2("./bin/owltools",args=owltoolsArgs,stdout=TRUE)
    #unlink(inputGafFilePath)
    map2slimGafDf<-read.csv(text=map2slimGafStdOut,header=FALSE,quote="",stringsAsFactors=FALSE,comment.char='!',sep="\t")
    outputGafDf<-read.csv(outputGafFilePath,header=FALSE,quote="",stringsAsFactors=FALSE, col.names=gafNamesVec,comment.char='!',sep="\t")
    #unlink(outputGafFilePath)
    outputGafDf
}
