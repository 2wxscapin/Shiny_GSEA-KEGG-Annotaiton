# shiny tool for GSEA-KEGG analysis based on clusterProfiler
## ui




library(shiny)
library(shinythemes)
library(DT)


ui <- fluidPage(theme = shinytheme("paper"),
                titlePanel("GSEA-KEGG"),
                hr(),
                h3("Get GSEA Info Table"),
                sidebarLayout(
                  
                  sidebarPanel(
                    fileInput("file1","Upload a RANKED GENE LIST (decreasing order) file",
                              multiple = FALSE,
                              accept = c('text/csv',
                                         'text/comma-separated-values',
                                         'text/tab-separated-values',
                                         'text/plain',
                                         '.csv')),
                    helpText("Note1. Your genelist should contain gene ID with corresponding numeric variable (e.g., FC)."),
                    helpText("Note2. The preferred format of gene ID is ENTREZID."),
                    radioButtons("categorycd", "Select a Category Code from ", 
                                 c("H"="hallmark","C1"="c1","C2"="c2","C3"="c3",
                                   "C4"="c4","C5"="c5","C6"="c6","C7"="c7"),
                                 selected = 'hallmark'),
                    br(),
                    textInput("annodb","Type your Annotation Package","org.Hs.eg.db"),
                    br(),
                    actionButton("action", "RUN"),
                    br(),br(),
                    actionButton("reset", "Rest the Category")
                    # hr(),
                    # downloadButton("downID","Gene Symbol Table")
                    # verticalLayout(downloadButton("downID","Gene Symbol Table"),
                    #                br(),
                    #                downloadButton("downMSig","MSig Table"),
                    #                br(),
                    #                downloadButton("downKEGG","KEGG Table"))
                    # 
                  ),
                  
                  mainPanel(
                    tabsetPanel(type = "tabs", 
                                tabPanel("Your Genes", DT::dataTableOutput("genetable")),
                                tabPanel("MSigDb Table", DT::dataTableOutput("msigtable")),
                                tabPanel("KEGG pathways", DT::dataTableOutput("keggoutput")))
                    
                    
                    
                  )
                  
                ),
                  
                hr(),
                h3("Pathview Download"),
                flowLayout(textInput("selectedpath","Type one KEGG ID","hsa00270"),
                br(),br(),
                actionButton("downloadimage", "Pathview Image",
                             icon("project-diagram")),
                hr()
                  
                )
)
                



# server <- function(input, output) {}
# shinyApp(ui = ui, server = server)
