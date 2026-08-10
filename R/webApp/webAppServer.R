taxogenoWebAppServer <- function(input, output) {
    options(shiny.maxRequestSize=128*1024^2) 
    conn <- dbConnect(RPostgres::Postgres(),
                      host = '/tmp/',
                      dbname = 'taxogeno')

    ##Proteomes list
    proteomeListDfReactive<-reactive({
        jobIdInput<-input$jobid
        if(is.na(jobIdInput) || is.null(jobIdInput) || jobIdInput==""){
            jobIdInput<-NA
        }

        if(is.na(jobIdInput))
            proteomeList(conn)
        else
            userProteomeList(conn,jobIdInput)
    })
    
    
    output$proteomeList<-DT::renderDataTable({
        DT::datatable (
                # proteomeListDF[, input$proteomeListColumnsShown, drop = FALSE],
                proteomeListDfReactive(),
                selection="single",
                rownames=proteomeListDfReactive()$proteomeid
            )
    })

    ## reactive para tener el proteoma seleccionado de la lista que viene de la base de datos
    dbProteomeIdSelectedReactive<-reactive({
        proteomeIdSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        proteomeIdSelected
    })

    ##proteomeId selected reactive
    ## puede venir o bien de un proteoma de usuario o bien de un proteoma de la base de datos
    userProteomeIdSelectedReactive<-reactive({
        proteomeIdSelected<-proteomeIdFromJobId(input$jobId)
        proteomeIdSelected
    })
    
    ##Proteome detail
    proteomeInfoReactive<-eventReactive(input$proteomeSubsection,{
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        proteomeInfo(conn,proteomeidSelected)
    })

    output[["proteomeInfo.proteomeid"]]<-renderText(proteomeInfoReactive()[["proteomeid"]])
    output[["proteomeInfo.basename"]]<-renderText(proteomeInfoReactive()[["basename"]])
    output[["proteomeInfo.creationtimestamp"]]<-renderText(proteomeInfoReactive()[["creationtimestamp"]])
    output[["proteomeInfo.is_userproteome"]]<-renderText(proteomeInfoReactive()[["is_userproteome"]])
    output[["proteomeInfo.genecount"]]<-renderText(proteomeInfoReactive()[["genecount"]])
    output[["proteomeInfo.ncbitaxid"]]<-renderText(proteomeInfoReactive()[["ncbitaxid"]])
    output[["proteomeInfo.scientific_name"]]<-renderText(proteomeInfoReactive()[["scientific_name"]])
    output[["proteomeInfo.gcaid"]]<-renderText(proteomeInfoReactive()[["gcaid"]])
    output[["proteomeInfo.dbname"]]<-renderText(proteomeInfoReactive()[["dbname"]])
    output[["proteomeInfo.sourceurl"]]<-renderText(proteomeInfoReactive()[["sourceurl"]])

    ## ############
    ## Gene
    ## ############
    geneListReactive<-eventReactive(input$genesInfoTabsetPanel,{
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        geneList(conn,proteomeidSelected)
    })
    output$geneList<-DT::renderDataTable({
        DT::datatable(
                geneListReactive(),
                selection="single",
                rownames=geneListReactive()$geneid
            ) #end datatable
    })

    ## ############
    ## Annotations
    ## ############
    geneGoAnnotationListReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,"geneid"]
        geneGoAnnotationList(conn,geneidSelected) ## !!
    })
    output$geneGoAnnotationList<-DT::renderDataTable({
        DT::datatable(
                geneGoAnnotationListReactive(),
                selection="single",
                rownames=geneGoAnnotationListReactive()$geneid
            )
    })
    ##
    geneKeywordAnnotationListReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,"geneid"]
        geneKeywordAnnotationList(conn,geneidSelected) ## !!
    })
    output$geneKeywordAnnotationList<-DT::renderDataTable({
        DT::datatable(
                geneKeywordAnnotationListReactive(),
                selection="single",
                rownames=geneKeywordAnnotationListReactive()$geneid
            )
    })
    ##
    geneEnzymeAnnotationListReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,"geneid"]
        geneEnzymeAnnotationList(conn,geneidSelected)
    })
    output$geneEnzymeAnnotationList<-DT::renderDataTable({
        DT::datatable(
                geneEnzymeAnnotationListReactive(),
                selection="single",
                rownames=geneEnzymeAnnotationListReactive()$geneid
            )
    })
    ##gene pathway annotation list
    genePathwayAnnotationListReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,"geneid"]
        genePathwayAnnotationList(conn,geneidSelected)
    })
    output$genePathwayAnnotationList<-DT::renderDataTable({
        DT::datatable(
                genePathwayAnnotationListReactive(),
                selection="single",
                rownames=genePathwayAnnotationListReactive()$geneid
            ) #end datatable
    })
    
    ##Gene detail
    geneInfoReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,][["geneid"]]
        geneInfo(conn,geneidSelected)
    })

    output[["geneInfo.geneid"]]<-renderText(geneInfoReactive()[["geneid"]])
    output[["geneInfo.genename"]]<-renderText(geneInfoReactive()[["genename"]])
    output[["geneInfo.genedescription"]]<-renderText(geneInfoReactive()[["genedescription"]])
    output[["geneInfo.fastaheader"]]<-renderText(geneInfoReactive()[["fastaheader"]])
    output[["geneInfo.aasequence"]]<-renderText(str_replace_all(geneInfoReactive()[["aasequence"]],"(.{60})", "\\1\n"))
    
    ## Tree
    output$tree <- renderTree({
        getTaxonomyNestedStructure(conn,1)
    })
    
    referenceProteomeIdsForTaxonomySelectedNodesReactive<-reactive({
        getReferenceProteomeIdsForTaxonomySelectedNodes( conn, input[["tree"]] )
    })

    ## Similarity
    similarityListReactive<-reactive({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,"proteomeid"]
        referenceProteomeIds<-referenceProteomeIdsForTaxonomySelectedNodesReactive()
        getEuclideanInterproteomeDistances(conn,'generated_goslim', proteomeidSelected, referenceProteomeIds) ## !!
    })
    output$similarityList<-DT::renderDataTable({
        DT::datatable(
                similarityListReactive() ##,
                ## selection="single" ##,
                ## rownames=genePathwayAnnotationListReactive()$geneid
            ) #end datatable
    })
    ## Clust
##    output$euclideanDistancesDendogramPlot<-renderForceNetwork({
    output$euclideanDistancesDendogramPlot<-renderPlot({
        referenceProteomeIds<-referenceProteomeIdsForTaxonomySelectedNodesReactive()
        distancesDf<-dbGetQuery(conn, "SELECT * FROM taxogeno.euclidean_distances_keyword WHERE (proteomeid_greatest  IN ($1) OR proteomeid_least IN ($1))",params=list(referenceProteomeIds))
        distancesDf<-merge(distancesDf,dbGetQuery(conn,"select genecount as genecount_greatest, proteomeid as proteomeid_greatest, basename as basename_greatest from taxogeno.proteome where proteomeid in ($1)",list(referenceProteomeIds)), by="proteomeid_greatest")
        distancesDf<-merge(distancesDf,dbGetQuery(conn,"select genecount as genecount_least, proteomeid as proteomeid_least, basename as basename_least from taxogeno.proteome where proteomeid in ($1)",list(referenceProteomeIds)), by="proteomeid_least")
        distancesDf[,"basename_greatest"]<-gsub("^(.*)_GCA.*$","\\1",distancesDf[,"basename_greatest"])
        distancesDf[,"basename_least"]<-gsub("^(.*)_GCA.*$","\\1",distancesDf[,"basename_least"])

        distancesDf<-distancesDf[,c("basename_greatest","basename_least","distance")]
        
        uniqueBasenamesVec<-unique(c(distancesDf$basename_least,distancesDf$basename_greatest))
        reflexiveDistancesVec<-rep(0,length(uniqueBasenamesVec))
        reflexiveDistancesDf<-data.frame(basename_least=uniqueBasenamesVec,
                                         basename_greatest=uniqueBasenamesVec,
                                         distance=reflexiveDistancesVec)
        distancesDf<-rbind(distancesDf,reflexiveDistancesDf)
        distancesDfSpecular<-data.frame(basename_least=distancesDf$basename_greatest,
                                basename_greatest=distancesDf$basename_least,
                                distance=distancesDf$distance)

        distancesDf<-rbind(distancesDf,distancesDfSpecular)
        rm(distancesDfSpecular)

        scalar1 <- function(x) { x / sqrt(sum(x^2))}

        distancesDf$similarity<-rep(1,length(distancesDf$distance))-scalar1(distancesDf$distance)

        distMatrix<-xtabs(distance~basename_greatest+basename_least,distancesDf)
        
        heatmap(distMatrix)
        
        ## distancesdf[,"similarity"]<-1/distancesdf[,"distance"]
        ## distancesdf<-distancesdf[order(distancesdf$distance),]
        ## distmatrix<-xtabs(similarity~basename_greatest+basename_least,distancesdf,sparse=FALSE)
        ## #distobj<-as.dist(xtabs(similarity~basename_greatest+basename_least,distancesdf,sparse=TRUE))

        ## ##euclideanDistancesClust<-hclust(distobj,method="complete")
        
        
        ## uniqueProteomeIdVec<-unique(c(distancesdf[,"proteomeid_least"],distancesdf[,"proteomeid_greatest"]))
        ## seqIdProteomeIdMap<-data.frame(id=seq_along(uniqueProteomeIdVec)-1,proteomeid=uniqueProteomeIdVec)
        ## rm(uniqueProteomeIdVec)

        ## sourceNodesDf<-merge(seqIdProteomeIdMap, distancesdf, by.x="proteomeid",by.y="proteomeid_least")
        ## colnames(sourceNodesDf)[colnames(sourceNodesDf) == 'basename_least'] <- 'name'
        ## targetNodesDf<-merge(seqIdProteomeIdMap, distancesdf, by.x="proteomeid",by.y="proteomeid_greatest")
        ## colnames(targetNodesDf)[colnames(targetNodesDf) == 'basename_greatest'] <- 'name'
        
        ## d3ForceGraphLinksDf<-data.frame(source=c(sourceNodesDf[,"id"], targetNodesDf[,"id"]),
        ##                                 target=c(targetNodesDf[,"id"], sourceNodesDf[,"id"]),
        ##                                 value=c(sourceNodesDf[,"distance"],targetNodesDf[,"distance"])*50)
        ## d3ForceGraphLinksDf<-unique(d3ForceGraphLinksDf[order(d3ForceGraphLinksDf[,"source"]),])
        
        ## d3ForceGraphNodesDf<-unique(data.frame(id=c(sourceNodesDf[,"id"],targetNodesDf[,"id"]),
        ##                                        name=c(sourceNodesDf[,"name"],targetNodesDf[,"name"]),
        ##                                        group=rep(1,nrow(sourceNodesDf))))
        ## d3ForceGraphNodesDf<-d3ForceGraphNodesDf[order(d3ForceGraphNodesDf[,"id"]),]

        ## print(d3ForceGraphNodesDf)
        ## print(d3ForceGraphLinksDf)
        
        ## rm(sourceNodesDf)
        ## rm(targetNodesDf)

        ## forceNetwork(Links = d3ForceGraphLinksDf, Nodes = d3ForceGraphNodesDf, Source = "source", Target = "target", Value = "value", NodeID = "name", Group="group", linkDistance=JS("function(d){return d.value * 30}"), opacity = 0.8, zoom=TRUE)


        ## distancesdf[,"similarity"]<-distancesdf[,"distance"]
        ## distobj<-as.dist(xtabs(similarity~basename_greatest+basename_least,distancesdf,sparse=TRUE))
        ## euclideanDistancesClust<-hclust(distobj,method="complete")

        ## dendroNetwork(euclideanDistancesClust, zoom=TRUE)
    })
    
    ## Global annotations
    molecularFunctionGoslimAnnotationPointReactive <- reactive({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,"proteomeid"]
        molecularFunctionGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, proteomeidSelected, "molecular_function" ) ## !!
        molecularFunctionGoslimAnnotationDF$golabel<-str_wrap(molecularFunctionGoslimAnnotationDF$golabel,width=15)
        molecularFunctionGoslimAnnotationDF
    })
    molecularFunctionGoslimAnnotationListReactive <- reactive({
        referenceProteomeIds<- referenceProteomeIdsForTaxonomySelectedNodesReactive()
        molecularFunctionGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, referenceProteomeIds, "molecular_function" ) ## !!
        molecularFunctionGoslimAnnotationDF$golabel<-str_wrap(molecularFunctionGoslimAnnotationDF$golabel,width=15)
        molecularFunctionGoslimAnnotationDF
    })
    molecularFunctionGoslimAnnotationBoxplotDataReactive <- reactive({
        referenceProteomeIds<- referenceProteomeIdsForTaxonomySelectedNodesReactive()
        molecularFunctionGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, referenceProteomeIds, "molecular_function" ) ## !!
        molecularFunctionGoslimAnnotationDF <- molecularFunctionGoslimAnnotationDF[order(molecularFunctionGoslimAnnotationDF$goid),]
        goDescriptionDf<-unique(molecularFunctionGoslimAnnotationDF[,c("goid","goaspect","golabel")])
        data.frame(
            goid=goDescriptionDf$goid,
            goaspect=goDescriptionDf$goaspect,
            golabel=goDescriptionDf$golabel,
            Q1=aggregate(relval~goid,molecularFunctionGoslimAnnotationDF,quantile,0.25)$relval,
            Q2=aggregate(relval~goid,molecularFunctionGoslimAnnotationDF,quantile,0.50)$relval,
            Q3=aggregate(relval~goid,molecularFunctionGoslimAnnotationDF,quantile,0.75)$relval,
            IQR=aggregate(relval~goid,molecularFunctionGoslimAnnotationDF,IQR)$relval
        )
    })
    output$molecularFunctionGoslimAnnotationList <- DT::renderDataTable({
        DT::datatable(
                molecularFunctionGoslimAnnotationBoxplotDataReactive(),
                selection="single"
            ) #end datatable
    })
    output$molecularFunctionGoslimAnnotationBoxplot <- renderGirafe({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,"proteomeid"]
        girafe(
            ggobj = ggplot(aes(y=relval,x=golabel),
                           data=molecularFunctionGoslimAnnotationListReactive()) +
                geom_boxplot() +
                geom_point(data=molecularFunctionGoslimAnnotationPointReactive(),aes(y=relval,x=golabel),color="red", size=5) +
                coord_flip(),
            options = list(opts_sizing(rescale = FALSE)),
            width_svg = 10, height_svg = 30
        )
    })

    ##biologicalProcessGoslim
    biologicalProcessGoslimAnnotationPointReactive <- reactive({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        biologicalProcessGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, proteomeidSelected, c("biological_process") )
        biologicalProcessGoslimAnnotationDF$golabel<-sapply(biologicalProcessGoslimAnnotationDF$golabel,str_wrap,width=20)
        biologicalProcessGoslimAnnotationDF
    })

    biologicalProcessGoslimAnnotationListReactive <- reactive({
        referenceProteomeIds<- referenceProteomeIdsForTaxonomySelectedNodesReactive()
        biologicalProcessGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, referenceProteomeIds, c("biological_process") )
        ## biologicalProcessGoslimAnnotationDF$iqrByGoid<-sapply(
        ##         biologicalProcessGoslimAnnotationDF$goid,
        ##         function(goidParam){IQR(subset(biologicalProcessGoslimAnnotationDF, goid==goidParam)$relval)}
        ## )
        biologicalProcessGoslimAnnotationDF$golabel<-sapply(biologicalProcessGoslimAnnotationDF$golabel,str_wrap,width=20)
        biologicalProcessGoslimAnnotationDF
    })
    output$biologicalProcessGoslimAnnotationList <- DT::renderDataTable({
        DT::datatable(
                biologicalProcessGoslimAnnotationListReactive(),
                selection="single"
            ) #end datatable
    })
    output$biologicalProcessGoslimAnnotationBoxplot <- renderGirafe({
        girafe(
            ggobj = ggplot(aes(y=relval,x=golabel),
                           data=biologicalProcessGoslimAnnotationListReactive()) +
                geom_boxplot() +
                geom_point(data=biologicalProcessGoslimAnnotationPointReactive(),aes(y=relval,x=golabel),color="red", size=5) +
                coord_flip(),
            options = list(opts_sizing(rescale = FALSE)),
            width_svg = 10, height_svg = 30
        )
    })
    ##cellularComponentGoslim
    cellularComponentGoslimAnnotationPointReactive <- reactive({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        cellularComponentGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, proteomeidSelected, c("cellular_component") )
        cellularComponentGoslimAnnotationDF$golabel<-sapply(cellularComponentGoslimAnnotationDF$golabel,str_wrap,width=20)
        cellularComponentGoslimAnnotationDF
    })    

    cellularComponentGoslimAnnotationListReactive <- reactive({
        referenceProteomeIds<- referenceProteomeIdsForTaxonomySelectedNodesReactive()
        cellularComponentGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, referenceProteomeIds, c("cellular_component") )
        ## cellularComponentGoslimAnnotationDF$iqrByGoid<-sapply(
        ##         cellularComponentGoslimAnnotationDF$goid,
        ##         function(goidParam){IQR(subset(cellularComponentGoslimAnnotationDF, goid==goidParam)$relval)}
        ## )
        cellularComponentGoslimAnnotationDF$golabel<-sapply(cellularComponentGoslimAnnotationDF$golabel,str_wrap,width=20)
        cellularComponentGoslimAnnotationDF
    })
    output$cellularComponentGoslimAnnotationList <- DT::renderDataTable({
        DT::datatable(
                cellularComponentGoslimAnnotationListReactive(),
                selection="single"
            ) #end datatable
    })
    output$cellularComponentGoslimAnnotationBoxplot <- renderGirafe({
        girafe(
            ggobj = ggplot(aes(y=relval,x=golabel),
                           data=cellularComponentGoslimAnnotationListReactive()) +
                geom_boxplot() +
                geom_point(data=cellularComponentGoslimAnnotationPointReactive(),aes(y=relval,x=golabel),color="red", size=5) +
                coord_flip(),
            options = list(opts_sizing(rescale = FALSE)),
            width_svg = 10, height_svg = 30
        )
    })

    ## user proteome upload
    tsvFilePath<-reactive({
        fileTsv <- input$tsv
        path <- fileTsv$datapath
        path
    })
    tsvFileName<-reactive({
        fileTsv <- input$tsv
        path <- fileTsv$name
        path
    })
    multifastaFilePath<-reactive({
        fileMultifasta <- input$multifasta
        path <- fileMultifasta$datapath
		path
	})
    tags<-reactive({
    	unlist(strsplit(input$tags,","))
    })
	
    jobIdOutputReactive<-eventReactive(input$sendButton,{        
        if(!is.null(input$tsv) && !is.null(input$multifasta)){
            jobId<-saveSma3sFileSet(conn, tsvFilePath(), multifastaFilePath(),tags(),isUserProteome=TRUE,webFile=TRUE,webSma3sAnnotationFileName=tsvFileName())
            jobId
        }else{
            "no proteome"
        }
    })
    output$jobId <- renderText(jobIdOutputReactive())

    ## Annotation info download
    
    output$downloadData <- downloadHandler(
       filename = function() {
         paste('data-', Sys.Date(), '.csv', sep='')
       },
       content = function(con) {
         write.csv(similarityListReactive(), con)
       })
    
}
