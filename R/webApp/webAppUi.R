## Página con tres pestañas
## Cada pestaña contiene un sidebarLayout
ui <- navbarPage(
    title = "Taxogeno",
    tabPanel(
        title= "Step 0: Upload proteome",
        sidebarLayout(
            sidebarPanel(
                helpText("Upload here your proteome annotated with Sma3s"),
                width=2
            ),
            mainPanel(
                column(
                    width=4,
                    h3("Upload form"),
                    textInput("tags", "Insert comma separated tags for this proteome"),
                    fileInput("tsv", "Choose Sma3s output TSV file", accept = ".tsv"),
                    fileInput("multifasta", "Choose multifasta file", accept = ".fasta"),
                    actionButton("sendButton","Insert proteome")
                ),
                column(
                    width=4,
                    h3("Warnings and errors")
                    ## Errores y advertencias podrían ser que los ficheros tuvieran nombres distintos, que los ficheros no se ajustaran al formato (por ejemplo que se subiera el summary en lugar del tsv con todos los datos) o que el multifasta no fuera un multifasta, o el resumen de si cuadra todo o no (esto no está implementado, esto depende de modificar taxogeno.R)
                ),
                column(
                    width=4,
                    h3("Upload result"),
                    h4("Job ID for uploaded proteome"),
                    verbatimTextOutput("jobId",placeholder=TRUE)
                )
            )
        )
    ),
    tabPanel(
        title = "Step 1: Select a proteome to compare",
        sidebarLayout(    
            sidebarPanel(
                conditionalPanel(
                    condition = 'input.proteomeSubsection === "Proteome selection"',
                    helpText("Proteome selection.")
                ),
                conditionalPanel(
                    condition = 'input.proteomeSubsection === "General info"',
                    helpText("General information about proteome.")
                ),
                conditionalPanel(
                    condition = 'input.proteomeSubsection === "Annotation info"',
                    helpText("Information about functional annotation of proteome.")
                ),
                conditionalPanel(
                    condition =
                        'input.proteomeSubsection === "Genes info"',
                    helpText("Genes info conditional panel."),
                    conditionalPanel(
                        condition = 'input.genesInfoTabsetPanel === "Genes table"'
                    )
                ),
                conditionalPanel(
                    condition = 'input.proteomeSubsection === "Similarity info"',
                    helpText("Display 5 records by default.")
                )
               ,width=2),
            
            mainPanel(
                tabsetPanel(
                    id = 'proteomeSubsection',
                    tabPanel(
                        title = "Proteome selection",
                        h3("Proteome list"),
                        p("If you have uploaded a proteome, insert your jobid below. Your proteome is not visible for anyone. A valid jobid is a condition for it to appear."),
                        textInput("jobid", "Insert your jobid"),
                        DT::dataTableOutput("proteomeList")
                    ),
                    tabPanel(
                        title = "General info",
                        h3("General info"),
                        column(
                            6,
                            h4("Internal proteome Identificator in this database"),
                            verbatimTextOutput("proteomeInfo.proteomeid")
                        ),
                        column(
                            6,
                            h4("Common name of uploaded files"),
                            verbatimTextOutput("proteomeInfo.basename")
                        ),
                        column(
                            6,
                            h4("Date and time of insertion in database"),
                            verbatimTextOutput("proteomeInfo.creationtimestamp")
                        ),
                        column(
                            6,
                            h4("Has been this proteome uploaded by an user?"),
                            verbatimTextOutput("proteomeInfo.is_userproteome")
                        ),
                        column(
                            6,
                            h4("Number of genes"),
                            verbatimTextOutput("proteomeInfo.genecount")
                        ),
                        column(
                            6,
                            h4("Assembly"),
                            verbatimTextOutput("proteomeInfo.gcaid")
                        ),                        
                        column(
                            6,
                            h4("NCBI Taxonomy Identificator"),
                            verbatimTextOutput("proteomeInfo.ncbitaxid")
                        ),
                        column(
                            6,
                            h4("Scientific name"),
                            verbatimTextOutput("proteomeInfo.scientific_name")
                        ),
                        column(
                            6,
                            h4("Database of origin"),
                            verbatimTextOutput("proteomeInfo.dbname")
                        ),
                        column(
                            6,
                            h4("Proteome URL resource"),
                            verbatimTextOutput("proteomeInfo.sourceurl")
                        ) 
                    ), #end tabPanel
                    ## ###################################
                    ## GENES INFO
                    ## ###################################
                    tabPanel(
                        title = "Genes info",
                        ##Los elementos de la pestaña sep. por comas
                        h3("Genes info"),
                        column(
                            12,
                            tabsetPanel(
                                id='genesInfoTabsetPanel',
                                tabPanel(
                                    title = "Genes table",
                                    ##Los elementos de la pestaña sep. por comas
                                    h4("Gene table"),
                                    DT::dataTableOutput("geneList")
                                ), #end tabPanel
                                tabPanel(
                                    title="Gene info",
                                    h4("Gene info"),
                                    h5("geneid"),
                                    verbatimTextOutput("geneInfo.geneid"),
                                    h5("genename"),
                                    verbatimTextOutput("geneInfo.genename"),
                                    h5("genedescription"),
                                    verbatimTextOutput("geneInfo.genedescription"),
                                    h5("fastaheader"),
                                    verbatimTextOutput("geneInfo.fastaheader"),
                                    h5("aasequence"),
                                    verbatimTextOutput("geneInfo.aasequence")
                                ), #end tabPanel
                                tabPanel(
                                    title="Gene annotation info",
                                    h4("Gene GO annotations"),
                                    DT::dataTableOutput("geneGoAnnotationList"),
                                    h4("Gene Keyword annotations"),
                                    DT::dataTableOutput("geneKeywordAnnotationList"),
                                    h4("Gene EC annotations"),
                                    DT::dataTableOutput("geneEnzymeAnnotationList"),
                                    h4("Gene Pathway annotations"),
                                    DT::dataTableOutput("genePathwayAnnotationList")
                                )
                            ) #end tabsetPanel
                        ) #end column
                    ) #end tabPanel
                ) #end tabsetPanel
               ,width=10) #end mainPanel
        ) #end sidebarLayout
    ), #end tabPanel
    tabPanel(
        title = "Step 2: Select a set of reference proteomes",
        sidebarLayout(
            sidebarPanel(
                helpText("Step 1: Select proteomes to compare from taxonomic tree")
            ),
            mainPanel(
                shinyTree("tree", checkbox = TRUE, search = TRUE)
            )
        )
    ),
    tabPanel(
        title = "Step 3: Check annotation info",
        sidebarLayout(
            sidebarPanel(
                helpText("Here there is only information about generated GO Slim")
            ),
            mainPanel(
                h3("Annotation info"),
                column(
                    12,
                    tabsetPanel(
                        id='annotationInfoTabsetPanel',
                        tabPanel(
                            title="GOslim nnotation",
                            h3("GOslim annotation"),
                            column(
                                12,
                                tabsetPanel(
                                    tabPanel(
                                        title="GOslim Annotation, boxplot",
                                        ##Los elementos de la pestaña sep. por comas
                                        h4("Molecular function GOslim Annotation, boxplot"),
                                        girafeOutput("molecularFunctionGoslimAnnotationBoxplot"),
                                        h4("Biological process GOslim Annotation, boxplot"),
                                        girafeOutput("biologicalProcessGoslimAnnotationBoxplot"),
                                        h4("Cellular component GOslim Annotation, boxplot"),
                                        girafeOutput("cellularComponentGoslimAnnotationBoxplot")
                                    ),
                                    tabPanel(
                                        title="GOslim Annotation, table",
                                        h3("Molecular Function GOslim Annotation, table"),
                                        DT::dataTableOutput("molecularFunctionGoslimAnnotationList"),
                                        h3("Biological process GOslim Annotation, table"),
                                        DT::dataTableOutput("biologicalProcessGoslimAnnotationList"),
                                        h3("Cellular component GOslim Annotation, table"),
                                        DT::dataTableOutput("cellularComponentGoslimAnnotationList")
                                    )
                                ) ##end tabset panel goslim boxplot/list
                            ) ## end column
                        ) ## end tabpanel goslim
                    ) ##end tabset panel goslim/kw/...
                ) ## end column
            )
        )
    ), #end tabPanel
    
    tabPanel(
        title = "Step 4: Check similarity info",
        sidebarLayout(
            sidebarPanel(
                helpText("For this funcionality to work, you must select a proteome in step 1 and a set of proteomes from taxonomic tree in step 2"),
                helpText("You can review generated data in tabulated format and download it. Also there is a graphical representation of similarities between selected proteomes."),
                ## Para hacer la similitud debe seleccionarse un proteoma en el paso 1 y un conjunto de proteomas de referencia. Sólo cuando ambos estén seleccionados deben salir los elementos gráficos
                width=2
            ),
            mainPanel(
                h3("Similarity info between selected proteome and reference proteomes"),
                tabsetPanel(
                    tabPanel(
                        title="Similarity Info, tabulated",
                        DT::dataTableOutput("similarityList"),
                        downloadButton("downloadData", label = "Download table")
                    ),
                    tabPanel(
                        title="Similarity Info, dendrogram",
                        plotOutput("euclideanDistancesDendogramPlot",height="2000px")
                    )
                )
            )
        )
    ) #end tabPanel
) #end fluidPage
