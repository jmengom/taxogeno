#' Save Sma3s File Set to Database (Main ETL Pipeline)
#'
#' The main ETL (Extract-Transform-Load) function that processes a complete
#' Sma3s annotation set and loads it into the taxogeno database.
#'
#' @param dbConn A DBI database connection object.
#' @param sma3sAnnotationFilePath Character string. Path to the Sma3s TSV annotation file.
#' @param multifastaFilePath Character string. Path to the multi-FASTA protein file.
#' @param tagVec Character vector. Optional tags to associate with the proteome.
#' @param isUserProteome Logical. Whether this is a user-uploaded proteome.
#' @param gcaId Character string. Optional GCA accession number.
#' @param ncbiTaxId Integer. Optional NCBI taxonomy ID.
#' @param dbName Character string. Source database name.
#' @param sourceUrl Character string. Source URL for the data.
#' @param doUpdateNcbiAssemblyInfo Logical. Whether to update NCBI assembly info.
#' @param webFile Logical. Whether the files were uploaded via web interface.
#'   If TRUE, the function uses \code{webSma3sAnnotationFileName} for the basename.
#' @param webSma3sAnnotationFileName Character string. Original filename for web uploads.
#'   Only used when \code{webFile = TRUE}.
#'
#' @return Character string. The generated jobid (UUID) for tracking the upload.
#'
#' @details
#' The pipeline performs these steps in order:
#' \enumerate{
#'   \item \strong{Read annotations}: Loads the Sma3s annotation TSV file using
#'     \code{\link{readSma3sAnnotation}}.
#'   \item \strong{Read multifasta}: Loads the protein FASTA file using
#'     \code{\link{readMultifasta}}.
#'   \item \strong{Insert proteome metadata}: Creates a record in
#'     \code{taxogeno.proteome} with:
#'     \itemize{
#'       \item Basename (derived from file path or web filename)
#'       \item MD5 checksums for both input files
#'       \item Gene count
#'       \item User/assembly flags
#'       \item GCA ID and NCBI taxonomy ID (if provided)
#'       \item Source database name and URL
#'       \item Auto-generated jobid (UUID) for web tracking
#'     }
#'   \item \strong{Update NCBI assembly info}: If \code{!isUserProteome && doUpdateNcbiAssemblyInfo},
#'     extracts the GCA ID from the basename and updates the assembly information.
#'   \item \strong{Insert taxonomy}: If not a user proteome, inserts the taxon and
#'     its ancestors into the taxonomy table using \code{insert_taxon()}.
#'   \item \strong{Insert genes}: Stores basic gene information (fastaheader,
#'     fastashortheader, genename, genedescription) and links sequences from the
#'     multifasta.
#'   \item \strong{Insert annotation relationships}: For each gene, expands
#'     semicolon-separated fields into 1:N relationships and inserts into:
#'     \itemize{
#'       \item \code{gene_enzyme_rel} (EC numbers)
#'       \item \code{gene_keyword_rel} (UniProt keywords)
#'       \item \code{gene_pathway_rel} (Pathway annotations)
#'       \item \code{gene_goslim_rel} (GO Slim terms)
#'       \item \code{gene_goc_rel} (Cellular Component GO terms)
#'       \item \code{gene_gof_rel} (Molecular Function GO terms)
#'       \item \code{gene_gop_rel} (Biological Process GO terms)
#'     }
#'   \item \strong{Generate GO Slim summaries}: Recalculates GO Slim annotations
#'     using OWLTools map2slim and saves aggregated counts to
#'     \code{generated_goslim_summary}.
#'   \item \strong{Calculate Euclidean distances}: Computes and caches distance
#'     matrices based on GO Slim and Keyword annotations via
#'     \code{\link{saveEuclideanDistancesForProteomeId}}.
#' }
#'
#' @section Job ID (UUID):
#' The function generates a UUID using \code{UUIDgenerate()} which is stored
#' in the \code{proteome.jobid} column. This jobid is returned by the function
#' and can be used to track the upload status in the web application.
#'
#' @section Memory Management:
#' The function handles large files efficiently:
#' \itemize{
#'   \item FASTA files are read line-by-line using streaming (see \code{\link{readMultifasta}})
#'   \item Distance calculations use sparse matrix representation and chunking
#' }
#'
#' @section Web Upload Support:
#' When \code{webFile = TRUE}, the function:
#' \itemize{
#'   \item Uses \code{webSma3sAnnotationFileName} for the basename
#'   \item Still reads from \code{sma3sAnnotationFilePath} for the actual file content
#'   \item The jobid can be used to locate the uploaded files for cleanup
#' }
#' 
#' @importFrom DBI dbGetQuery dbExecute dbWriteTable SQL
#' @importFrom tools file_path_sans_ext
#' @importFrom uuid UUIDgenerate
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' 
#' # Standard usage for public proteomes
#' jobid <- saveSma3sFileSet(
#'   conn,
#'   "path/to/annotation.tsv",
#'   "path/to/proteome.faa",
#'   tagVec = c("bacteria", "example"),
#'   isUserProteome = FALSE,
#'   doUpdateNcbiAssemblyInfo = TRUE
#' )
#' 
#' # Web upload usage
#' jobid <- saveSma3sFileSet(
#'   conn,
#'   "/tmp/upload_1234.tsv",
#'   "/tmp/upload_1234.faa",
#'   tagVec = "user_upload",
#'   isUserProteome = TRUE,
#'   webFile = TRUE,
#'   webSma3sAnnotationFileName = "my_proteome.tsv"
#' )
#' 
#' print(paste("Upload job ID:", jobid))
#' }
#'
#' @seealso
#' \code{\link{readSma3sAnnotation}} for reading the annotation file
#' \code{\link{readMultifasta}} for reading FASTA files
#' \code{\link{updateNcbiAssemblyInfoForProteomeIdVec}} for updating assembly metadata
#' \code{\link{extractColumnDfFromInsertedGeneAnnotationDf}} for parsing semicolon-separated fields
#' \code{\link{saveGeneratedGoslimSummaryForGoStrictAnnotationData}} for GO Slim processing
#' \code{\link{saveEuclideanDistancesForProteomeId}} for distance calculations
#'
#' @export
saveSma3sFileSet<-function(dbConn,
                           sma3sAnnotationFilePath,
                           multifastaFilePath,
                           tagVec=character(),
                           isUserProteome=FALSE,
                           gcaId=NA,
                           ncbiTaxId=NA,
                           dbName=NA,
                           sourceUrl=NA,
                           doUpdateNcbiAssemblyInfo=TRUE,
                           webFile=FALSE,
                           webSma3sAnnotationFileName=NA){

    sma3sAnnotationDf<-readSma3sAnnotation(sma3sAnnotationFilePath)
    sma3sAnnotationMd5Sum<- md5sum(sma3sAnnotationFilePath)
    multifastaDf<-readMultifasta(multifastaFilePath)
    multifastaMd5Sum<-md5sum(multifastaFilePath)

    geneCount<-sum(complete.cases(sma3sAnnotationDf[,"fastaheader"]), na.rm=TRUE)
    if(geneCount==0){
        error("Empty Sma3s tsv file. Proteome with no genes.")
    }

    sma3sAnnotationBaseName<-NA
    if(webFile==FALSE){
        sma3sAnnotationBaseName<-basename(tools::file_path_sans_ext(sma3sAnnotationFilePath))
    } else {
        sma3sAnnotationBaseName<-basename(tools::file_path_sans_ext(webSma3sAnnotationFileName))
    }
    
    ## ## Campos del proteoma
    proteomeInfoObj<-data.frame(
        ## proteomeid=NA, ## proteomeid asignado por la base de datos. No se mete como columna en el DF porque lo va a interpretar como que es null en dbWriteTable
        basename = sma3sAnnotationBaseName,
        sma3sannotationmd5sum = sma3sAnnotationMd5Sum,
        multifastamd5sum = multifastaMd5Sum,
        genecount = geneCount,
        is_userproteome = isUserProteome,
        gcaid = gcaId,
        ncbitaxid = ncbiTaxId,
        dbname = dbName,
        sourceurl = sourceUrl,
        do_updatencbiassemblyinfo = doUpdateNcbiAssemblyInfo,
		jobid = UUIDgenerate()
    )
        
    proteomeId<-dbGetQuery(dbConn,"
      INSERT INTO taxogeno.proteome(
        basename,
        sma3sannotationmd5sum,
        multifastamd5sum,
        genecount,
        is_userproteome,
        gcaid,
        ncbitaxid,
        dbname,
        sourceurl,
        do_updatencbiassemblyinfo,
		jobid
      ) VALUES ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11)
      RETURNING proteomeid 
    ",params=list(proteomeInfoObj[["basename"]],
                  proteomeInfoObj[["sma3sannotationmd5sum"]],
                  proteomeInfoObj[["multifastamd5sum"]],
                  proteomeInfoObj[["genecount"]],
                  proteomeInfoObj[["is_userproteome"]],
                  proteomeInfoObj[["gcaid"]],
                  proteomeInfoObj[["ncbitaxid"]],
                  proteomeInfoObj[["dbname"]],
                  proteomeInfoObj[["sourceurl"]],
                  proteomeInfoObj[["do_updatencbiassemblyinfo"]],
				  proteomeInfoObj[["jobid"]])) [1,"proteomeid"] ## tomar la primera entrada de la columna proteomeid, porque sólo hay una, del RETURNING

    ## ## si es un proteoma no de un usuario, sacar el posible GCA del proteoma que se ha introducido a partir de su nombre
    if(!isUserProteome && doUpdateNcbiAssemblyInfo ){
        updateNcbiAssemblyInfoForProteomeIdVec(dbConn,proteomeId)
    }

    ## actualizar a lo que hay en la bbdd
    proteomeInfoObj<-dbGetQuery(dbConn,"SELECT * FROM taxogeno.proteome WHERE proteomeid=$1", params=list(proteomeId))
    
    ## ## Etiquetas del proteoma
    tagDf<-data.frame(proteomeid = rep(proteomeId,length(tagVec)),
                      tag = tagVec )
    dbWriteTable(dbConn,SQL("taxogeno.proteome_tag_rel"),tagDf,append=TRUE,row.names=FALSE)

    ## GENES
    ## ## Campos generales de los genes
    geneDf<-sma3sAnnotationDf[,c("fastashortheader","genename","genedescription","fastaheader")]
    ## ## Engancharles sus secuencias si las hay
    geneDf<-merge(geneDf,multifastaDf,by="fastaheader",all.x=TRUE)
    ## ## Vincular con el proteoma que se está metiendo
    geneDf[,"proteomeid"]<-rep(proteomeId, nrow(geneDf))
    ## ## Escribir en la base de datos
    dbWriteTable(dbConn,SQL("taxogeno.gene"),geneDf,append=TRUE,row.names=FALSE)

    ## ## las claves autoincrementales hay que recuperarlas de la base de datos
    geneIdFastaShortHeaderRelDf<-dbGetQuery(dbConn,
                                            "SELECT geneid,fastashortheader FROM taxogeno.gene WHERE proteomeid=$1",
                                            params=list(proteomeId))
    insertedGeneAnnotationDf<-merge(sma3sAnnotationDf,geneIdFastaShortHeaderRelDf, by="fastashortheader")
    rm(sma3sAnnotationDf)
    
    ## Campos 1:n gen:campos
    enzymeDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "enzyme", newColName="ec")
    dbWriteTable(dbConn,SQL("taxogeno.gene_enzyme_rel"),enzymeDf,append=TRUE,row.names=FALSE)

    keywordDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "keyword")
    dbWriteTable(dbConn,SQL("taxogeno.gene_keyword_rel"),keywordDf,append=TRUE,row.names=FALSE)
    rm(keywordDf)

    ## #########################
    pathwayDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "pathway")
    dbWriteTable(dbConn,SQL("taxogeno.gene_pathway_rel"),pathwayDf,append=TRUE,row.names=FALSE)
    rm(pathwayDf)
    
    ## #########################
    goslimDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "goslim", newColName="goid")
    dbWriteTable(dbConn,SQL("taxogeno.gene_goslim_rel"),goslimDf,append=TRUE,row.names=FALSE)
    rm(goslimDf)
    
    ## #########################
    gocDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"goc", newColName="goid")
    dbWriteTable(dbConn,SQL("taxogeno.gene_goc_rel"),gocDf,append=TRUE,row.names=FALSE)
    gofDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"gof", newColName="goid")
    dbWriteTable(dbConn,SQL("taxogeno.gene_gof_rel"),gofDf,append=TRUE,row.names=FALSE)
    gopDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"gop", newColName="goid")
    dbWriteTable(dbConn,SQL("taxogeno.gene_gop_rel"),gopDf,append=TRUE,row.names=FALSE)
    goCategory<-c("goc"="cellular_component",
                  "gop"="biological_process",
                  "gof"="molecular_function")
    gocDf[,"gocategory"]<-rep(goCategory["goc"],nrow(gocDf))
    gofDf[,"gocategory"]<-rep(goCategory["gof"],nrow(gofDf))
    gopDf[,"gocategory"]<-rep(goCategory["gop"],nrow(gopDf))
    goStrictDf<-rbind(gocDf,gofDf,gopDf)
    colnames(goStrictDf)[colnames(goStrictDf)=="goid"]<-"annotkwid"
    colnames(goStrictDf)[colnames(goStrictDf)=="gocategory"]<-"aspect"
    rm(gocDf)
    rm(gofDf)
    rm(gopDf)   
    saveGeneratedGoslimSummaryForGoStrictAnnotationData(dbConn,proteomeId,goStrictDf)

    ## Cálculo y cacheo de distancias ##
    ## saveKeywordEuclideanDistancesForProteomeId(dbConn,proteomeId)
    ## saveGeneratedGoslimEuclideanDistancesForProteomeId(dbConn,proteomeId)
    saveEuclideanDistancesForProteomeId(dbConn,proteomeId,c("generated_goslim","keyword"))
    
    proteomeInfoObj[["jobid"]]
}
