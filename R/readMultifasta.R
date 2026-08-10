#' Read Multi-FASTA File Line by Line
#'
#' Reads a FASTA file containing multiple protein sequences using a streaming
#' approach to avoid loading the entire file into memory. This is essential
#' for processing large proteome files that may be several gigabytes.
#'
#' @param multifastaPath Character string. Path to the FASTA file.
#'
#' @return A data frame with two columns:
#'   \item{fastaheader}{The FASTA header (identifier after ">")}
#'   \item{aasequence}{The amino acid sequence (concatenated from all lines)}
#'
#' @details
#' The function processes the file line by line using a state machine:
#' \enumerate{
#'   \item Looks for lines starting with ">" (FASTA header)
#'   \item Collects subsequent lines as sequence until the next header
#'   \item Handles sequences that span multiple lines
#' }
#'
#' @examples
#' \dontrun{
#' fasta_df <- readMultifasta("path/to/proteome.faa")
#' nrow(fasta_df)  # Number of protein sequences
#' }
#'
#' @export
readMultifasta<-function(multifastaPath){   
    ## Leer fichero línea a línea para no petar cargando el fichero completo en RAM
    ## https://stackoverflow.com/questions/12626637/read-a-text-file-in-r-line-by-line#35761217

    ##inicializar variables que almacenan FASTA uno a uno
    fastaHeader=""
    fastaSequence=""

    ## booleanos control de flujo
    startSequenceBody=FALSE
    insideSequenceBody=FALSE

    fastaHeaderVec<-character()
    fastaSequenceVec<-character()
    fileConn = file(multifastaPath, "r")
    while ( TRUE ) {
        oneLine = readLines(fileConn, n = 1)
        if ( length(oneLine) == 0 ) {
            if(insideSequenceBody){
                ##cerrando cosa anterior
                fastaHeaderVec<-c(fastaHeaderVec,fastaHeader)
                fastaSequenceVec<-c(fastaSequenceVec,fastaSequence)
                
                ##limpiar variables de último elemento
                fastaHeader=""
                fastaSequence=""

                ##indicar que ya no se está leyendo el cuerpo de la secuencia
                insideSequenceBody = FALSE
            }
            break
        }

        ## Cosas de FASTA
        if(startsWith(oneLine, ">")){
            ## Encontrar una cabecera FASTA tiene que interrumpir la secuencia
            if(insideSequenceBody){
                ##cerrando cosa anterior

                fastaHeaderVec<-c(fastaHeaderVec,fastaHeader)
                fastaSequenceVec<-c(fastaSequenceVec,fastaSequence)
                
                ##limpiar variables de último elemento
                fastaHeader=""
                fastaSequence=""

                ##indicar que ya no se está leyendo el cuerpo de la secuencia
                insideSequenceBody = FALSE
            }

            ## empezando cosa nueva
            
            ## Guardar id 
            fastaHeader = gsub("^>([^\\s]*)", "\\1", oneLine)
            ##En la siguiente vuelta debe comenzar el cuerpo de la secuencia
            startSequenceBody = TRUE
            next
        }
        if(startSequenceBody){
            ##Cambio de booleano
            startSequenceBody=FALSE
            insideSequenceBody=TRUE
        }
        if(insideSequenceBody){
            ##Concatenar secuencia, por si esta está en varias líneas
            fastaSequence<-paste(fastaSequence, oneLine, sep='')
        }
    }
    close(fileConn)
    data.frame(fastaheader=fastaHeaderVec,aasequence=fastaSequenceVec)
}
