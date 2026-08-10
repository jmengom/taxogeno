proteomeGoAnnotationSummary <- function(conn, proteomeid){
    proteomeGoAnnotationSummaryResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.generated_goslim_summary.proteomeid,
        taxogeno.generated_goslim_summary.goid,
        taxogeno.generated_goslim_summary.annotkwcount,
        gene_ontology.gene_ontology.golabel,
        gene_ontology.gene_ontology.goaspect

      FROM taxogeno.generated_goslim_summary
        INNER JOIN gene_ontology.gene_ontology on (taxogeno.generated_goslim_summary.goid = gene_ontology.gene_ontology.goid)
      WHERE geneid= $1
    ")
    dbBind(proteomeGoAnnotationSummaryResult, list(proteomeid))
    proteomeGoAnnotationSummaryDF <- dbFetch(proteomeGoAnnotationSummaryResult)
    dbClearResult(proteomeGoAnnotationSummaryResult)
    rm(proteomeGoAnnotationSummaryResult)
    
    proteomeGoAnnotationSummaryDF
}

proteomeList<-function(conn){
    proteomeListResult <- dbSendQuery(conn, "
      SELECT
        proteomeid,
        basename,
        gcaid,
        prot.ncbitaxid as ncbitaxid,
        scientific_name
      FROM taxogeno.proteome prot
      LEFT JOIN taxogeno.taxonomy tax
      ON prot.ncbitaxid=tax.ncbitaxid
      WHERE is_userproteome = FALSE
    ")
    proteomeListDF<- dbFetch(proteomeListResult)
    dbClearResult(proteomeListResult)
    rm(proteomeListResult)

    proteomeListDF
}
userProteomeList<-function(conn,jobId){
    proteomeListDf <- dbGetQuery(conn, "
      SELECT
        proteomeid,
        basename,
        gcaid,
        prot.ncbitaxid as ncbitaxid,
        scientific_name
      FROM taxogeno.proteome prot
      LEFT JOIN taxogeno.taxonomy tax
      ON prot.ncbitaxid=tax.ncbitaxid
      WHERE is_userproteome = TRUE
      AND jobid=$1
    ",params=list(jobId))
    proteomeListDf
}

proteomeInfo <- function (conn, proteomeid){
    proteomeInfoResult <- dbSendQuery(conn, "
      SELECT
        prot.proteomeid,
        prot.basename,
        prot.creationtimestamp,
        prot.is_userproteome,
        prot.genecount,
        prot.gcaid,
        prot.ncbitaxid,
        prot.dbname,
        prot.sourceurl,

        tax.scientific_name,
        tax.node_rank

      FROM taxogeno.proteome prot
        LEFT JOIN taxogeno.taxonomy tax ON (prot.ncbitaxid=tax.ncbitaxid)
      WHERE proteomeid= $1
    ")
    dbBind(proteomeInfoResult, list(proteomeid))
    proteomeInfoDF<- dbFetch(proteomeInfoResult)
    dbClearResult(proteomeInfoResult)
    rm(proteomeInfoResult)

    proteomeInfoDF[1,]
}

proteomeIdFromJobId <- function (conn, jobid){
    proteomeIdResult <- dbSendQuery(conn, "
      SELECT
        prot.proteomeid as proteomeid,
      FROM taxogeno.proteome prot
      WHERE jobid= $1
    ")
    dbBind(proteomeIdResult, list(jobid))
    proteomeIdDF<- dbFetch(proteomeIdResult)
    dbClearResult(proteomeIdResult)
    rm(proteomeIdResult)

    proteomeInfoDF[1,"proteomeid"]
}

geneList <- function(conn, proteomeid){
    geneListResult <- dbSendQuery(conn,"
      SELECT
        geneid,
        genename,
        genedescription,
        fastaheader
      FROM taxogeno.gene WHERE proteomeid= $1
    ")
    dbBind(geneListResult, list(proteomeid))
    geneListDF <- dbFetch(geneListResult)
    dbClearResult(geneListResult)
    rm(geneListResult)
    
    geneListDF
}
geneListColumns<-list("geneid","genename","genedescription","fastaheader","aasequence")

geneInfo <- function(conn, geneid){
    geneInfoResult <- dbSendQuery(conn,"
      SELECT
        geneid,
        genename,
        genedescription,
        fastaheader,
        aasequence
      FROM taxogeno.gene WHERE geneid= $1
    ")
    dbBind(geneInfoResult, list(geneid))
    geneInfoDF <- dbFetch(geneInfoResult)
    dbClearResult(geneInfoResult)
    rm(geneInfoResult)
    
    geneInfoDF[1:1,]
}
geneInfoColumns<-list("geneid","genename","genedescription","fastaheader","aasequence")

#########################################
geneGoAnnotationList <- function(conn, geneid){
    geneGoAnnotationListResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.gene_goc_rel.geneid,
        taxogeno.gene_goc_rel.goid,
        gene_ontology.gene_ontology.golabel,
        gene_ontology.gene_ontology.goaspect

      FROM taxogeno.gene_goc_rel
        INNER JOIN gene_ontology.gene_ontology on (taxogeno.gene_goc_rel.goid = gene_ontology.gene_ontology.goid)
      WHERE geneid= $1
UNION
      SELECT
        taxogeno.gene_gof_rel.geneid,
        taxogeno.gene_gof_rel.goid,
        gene_ontology.gene_ontology.golabel,
        gene_ontology.gene_ontology.goaspect

      FROM taxogeno.gene_gof_rel
        INNER JOIN gene_ontology.gene_ontology on (taxogeno.gene_gof_rel.goid = gene_ontology.gene_ontology.goid)
      WHERE geneid= $1
UNION
      SELECT
        taxogeno.gene_gop_rel.geneid,
        taxogeno.gene_gop_rel.goid,
        gene_ontology.gene_ontology.golabel,
        gene_ontology.gene_ontology.goaspect

      FROM taxogeno.gene_gop_rel
        INNER JOIN gene_ontology.gene_ontology on (taxogeno.gene_gop_rel.goid = gene_ontology.gene_ontology.goid)
      WHERE geneid= $1

    ")
    dbBind(geneGoAnnotationListResult, list(geneid))
    geneGoAnnotationListDF <- dbFetch(geneGoAnnotationListResult)
    dbClearResult(geneGoAnnotationListResult)
    rm(geneGoAnnotationListResult)
    
    geneGoAnnotationListDF
}

#############################################
geneKeywordAnnotationList <- function(conn, geneid){
    geneKeywordAnnotationListResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.gene_keyword_rel.geneid,
        taxogeno.gene_keyword_rel.keyword,
        uniprot.uniprot_keyword.keywordid,
        uniprot.uniprot_keyword.keywordcategory,
        uniprot.uniprot_keyword.keyworddescription

      FROM taxogeno.gene_keyword_rel
        INNER JOIN uniprot.uniprot_keyword on (taxogeno.gene_keyword_rel.keyword = uniprot.uniprot_keyword.keyword)
      WHERE geneid= $1
    ")
    dbBind(geneKeywordAnnotationListResult, list(geneid))
    geneKeywordAnnotationListDF <- dbFetch(geneKeywordAnnotationListResult)
    dbClearResult(geneKeywordAnnotationListResult)
    rm(geneKeywordAnnotationListResult)
    
    geneKeywordAnnotationListDF
}

#############################################
geneEnzymeAnnotationList <- function(conn, geneid){
    geneEnzymeAnnotationListResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.gene_enzyme_rel.geneid,
        taxogeno.gene_enzyme_rel.ec
      FROM taxogeno.gene_enzyme_rel
      WHERE geneid= $1
    ")
    dbBind(geneEnzymeAnnotationListResult, list(geneid))
    geneEnzymeAnnotationListDF <- dbFetch(geneEnzymeAnnotationListResult)
    dbClearResult(geneEnzymeAnnotationListResult)
    rm(geneEnzymeAnnotationListResult)
    
    geneEnzymeAnnotationListDF
}

#############################################
genePathwayAnnotationList <- function(conn, geneid){
    genePathwayAnnotationListResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.gene_pathway_rel.geneid,
        taxogeno.gene_pathway_rel.pathway
      FROM taxogeno.gene_pathway_rel
      WHERE geneid= $1
    ")
    dbBind(genePathwayAnnotationListResult, list(geneid))
    genePathwayAnnotationListDF <- dbFetch(genePathwayAnnotationListResult)
    dbClearResult(genePathwayAnnotationListResult)
    rm(genePathwayAnnotationListResult)
    
    genePathwayAnnotationListDF
}

getGoslimAnnotationForProteomeIdVector<- function(conn, proteomeIdVector, aspectVector){

    proteomeIdArrayString<-paste('{',
                                 paste(proteomeIdVector,collapse=','),
                                 '}',
                                 sep='')
    aspectArrayString<-paste('{',
                             paste(aspectVector,collapse=','),
                             '}',
                             sep='')
    goslimAnnotationResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.generated_goslim_summary.proteomeid as proteomeid,
        taxogeno.generated_goslim_summary.annotkwid  as goid,
        gene_ontology.goslim_generic.goaspect  as goaspect,
        gene_ontology.goslim_generic.golabel   as golabel,

        taxogeno.generated_goslim_summary.annotkwcount as gocount,
        taxogeno.proteome.genecount                    as genecount,
        taxogeno.generated_goslim_summary.annotkwcount::float/taxogeno.proteome.genecount::float as relval
      FROM taxogeno.generated_goslim_summary
      INNER JOIN taxogeno.proteome
        ON (taxogeno.generated_goslim_summary.proteomeid=taxogeno.proteome.proteomeid)
      INNER JOIN gene_ontology.goslim_generic
        ON taxogeno.generated_goslim_summary.annotkwid = gene_ontology.goslim_generic.goid
      WHERE taxogeno.generated_goslim_summary.proteomeid = ANY($1)
      AND taxogeno.generated_goslim_summary.annotkwcount > 0
      AND gene_ontology.goslim_generic.golabel NOT IN ('molecular_function','cellular_component','biological_process')
      AND gene_ontology.goslim_generic.goaspect = ANY($2)
    ")
    dbBind(goslimAnnotationResult, list(proteomeIdArrayString, aspectArrayString))
    goslimAnnotationDF <- dbFetch(goslimAnnotationResult)
    dbClearResult(goslimAnnotationResult)
    rm(goslimAnnotationResult)
    
    goslimAnnotationDF
}

getTaxonomyNestedStructure<-function(conn, ncbitaxid){
    taxonResult <- dbSendQuery(conn,"
      SELECT taxogeno.taxonomy_jsonb_children($1) as json_result
    ")
    dbBind(taxonResult, list(ncbitaxid))
    taxonDF <- dbFetch(taxonResult)
    dbClearResult(taxonResult)
    rm(taxonResult)
    taxonList<-fromJSON(taxonDF$json_result[1], simplifyVector=FALSE)
    rm(taxonDF)
    
    childrenTreeTraversal<-function(nodeListObj, initial=FALSE){
        
        childrenObjList<-lapply(
            nodeListObj[["children"]],
            childrenTreeTraversal
        )

        namesChildrenObjList<- function (childrenObjList) {
            lapply(
                childrenObjList,
                function(childObj){
                    objListLabel<-paste(
                        "ncbitaxid:",attributes(childObj)[["ncbitaxid"]],"; ",
                        attributes(childObj)[["node_rank"]]," ",
                        attributes(childObj)[["scientific_name"]],
                        sep=""
                    ) ## end paste
                    
                    ## return:
                    objListLabel
                } ## end function
            ) ## end lapply
        } ## end function

        names(childrenObjList)<-namesChildrenObjList(childrenObjList)
        
        objList<-structure(childrenObjList,
                           ncbitaxid=nodeListObj[["ncbitaxid"]],
                           parent_ncbitaxid=nodeListObj[["parent_ncbitaxid"]],
                           scientific_name=nodeListObj[["scientific_name"]],
                           node_rank=nodeListObj[["node_rank"]])

        if (initial==TRUE){
            wrappedObjList<-list(objList)
            names(wrappedObjList) <- namesChildrenObjList(wrappedObjList)
            wrappedObjList
        } else {
            objList
        }
    }
    childrenTreeTraversal(taxonList,TRUE)
}

getTaxonomySelectedNodes <- function(tree, ancestry=NULL, vec=list()){
    if (is.list(tree)){
        for (i in 1:length(tree)){
            anc <- c(ancestry, names(tree)[i])
            vec <- getTaxonomySelectedNodes(tree[[i]], anc, vec)
        }    
    }
    a <- attr(tree, "stselected", TRUE)
    if (!is.null(a) && a == TRUE){
                                        # Get the element name
        len_anc <- length(ancestry)
        el <- ancestry[len_anc]
        vec[length(vec)+1] <- el

        ## Save some attributes
        lapply(names(attributes(tree)),function(attribute){
            if(grepl("^st",attribute)
               || attribute == "ncbitaxid"
               || attribute == "parent_ncbitaxid"
               || attribute == "scientific_name"
               || attribute == "node_rank"){
                attr(vec[[length(vec)]], attribute) <<- attr(tree,attribute)
            }
        })
    }
    return(vec)
}

getReferenceProteomeIdsForTaxonomySelectedNodes<-function(conn, tree){
    ncbiTaxIdVector<-as.integer(
        sapply(
            getTaxonomySelectedNodes(tree),
            function(node){attr(node,"ncbitaxid")}
        )
    )
    refProteomeListResult <- dbSendQuery(conn,"
      SELECT
        proteomeid
      FROM taxogeno.proteome
      WHERE ncbitaxid in ($1)
      AND is_userproteome = FALSE
    ")
    dbBind(refProteomeListResult, list(ncbiTaxIdVector))
    refProteomeListDF <- dbFetch(refProteomeListResult)
    dbClearResult(refProteomeListResult)
    rm(refProteomeListResult)
    
    refProteomeListDF$proteomeid
}

getEuclideanInterproteomeDistances <- function(conn,
                                               annotationType,
                                               proteomeId,
                                               proteomeIdVector){
    ## distanceListResult <- dbSendQuery(conn,sprintf("
    ##    SELECT
    ##      proteomeid_least,
    ##      (SELECT tax.scientific_name
    ##       FROM taxogeno.proteome AS prot
    ##       INNER JOIN taxogeno.taxonomy AS tax
    ##       ON (prot.ncbitaxid=tax.ncbitaxid)
    ##       WHERE prot.proteomeid=eucdisttab.proteomeid_least) as ref_scientific_name,

    ##      proteomeid_greatest,
    ##      (SELECT tax.scientific_name
    ##       FROM taxogeno.proteome AS prot
    ##       INNER JOIN taxogeno.taxonomy AS tax
    ##       ON (prot.ncbitaxid=tax.ncbitaxid)
    ##       WHERE prot.proteomeid=eucdisttab.proteomeid_greatest) as comp_scientific_name,
    ##      distance
    ##    FROM taxogeno.euclidean_distances_%s AS eucdisttab
    ##    WHERE 
    ##            (proteomeid_least = %2$s and proteomeid_greatest in (%3$s))
    ##            or
    ##            (proteomeid_least in (%3$s) and proteomeid_greatest = %2$s)
    ## ",annotationType,proteomeId, paste(proteomeIdVector,collapse=",")))


    distanceListResult <- dbSendQuery(conn,sprintf("
       SELECT
         proteomeid_least,
         (SELECT basename
          FROM taxogeno.proteome AS prot
          WHERE prot.proteomeid=eucdisttab.proteomeid_least) as ref_scientific_name,

         proteomeid_greatest,
         (SELECT basename
          FROM taxogeno.proteome AS prot
          WHERE prot.proteomeid=eucdisttab.proteomeid_greatest) as comp_scientific_name,
         distance
       FROM taxogeno.euclidean_distances_%s AS eucdisttab
       WHERE 
               (proteomeid_least = %2$s and proteomeid_greatest in (%3$s))
               or
               (proteomeid_least in (%3$s) and proteomeid_greatest = %2$s)
    ",annotationType,proteomeId, paste(proteomeIdVector,collapse=",")))
    
    
    distanceListDF <- dbFetch(distanceListResult)
    dbClearResult(distanceListResult)
    rm(distanceListResult)
    
    distanceListDF
##########
}
