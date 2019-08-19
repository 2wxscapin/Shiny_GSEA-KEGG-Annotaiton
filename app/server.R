# shiny tool for GSEA-KEGG analysis based on clusterProfiler
## server

 
if(!require(clusterProfiler)){
  BiocManager::install("clusterProfiler")
  library(clusterProfiler)
}
if(!require(dplyr)){
  BiocManager::install("dplyr")
  library(dplyr)
}
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}
if(!require(msigdf)){
  devtools::install_github("ToledoEM/msigdf")
  library(msigdf)
}

library(pathview)
        


server <- function(input,output) {
  
  annodbpkg <- reactive({
    input$annodb
  })
  
  # symbtb <- reactive({
    output$genetable <- DT::renderDataTable({
      if (is.null(input$file1))
        return(NULL)
      yrgene <- read.csv(input$file1$datapath,header = F,sep = "\t")
      yrsymbol <- bitr(yrgene[,1], fromType = "ENTREZID", 
                       toType = "SYMBOL",
                       OrgDb = annodbpkg(),
                       drop = F) 
      genedf <- cbind(yrsymbol,yrgene[,2])
      colnames(genedf) <- c("ENTREZID","SYMBOL","Value")
      return(genedf)
    })
  # })
  
  
  genelist <- reactive({
    yrgene <- read.csv(input$file1$datapath,header = F,sep = "\t")
    entrezwithFC <- yrgene[,2]
    names(entrezwithFC) <- yrgene[,1]
    entrezwithFC
  })
  
  symblist <- reactive({
    yrgene <- read.csv(input$file1$datapath,header = F,sep = "\t")
    yrsymbol <- bitr(yrgene[,1], fromType = "ENTREZID", 
                     toType = "SYMBOL",
                     OrgDb = annodbpkg()) 
    filtered <- intersect(row.names(yrgene),row.names(yrsymbol))
    symbls <- yrgene[filtered,2]
    names(symbls) <- yrsymbol[,2]
    symbls
  })
  
  selectcd <- reactive({
    if(input$categorycd == "hallmark"){
      selectcd <- msigdf.human %>% 
        filter(category_code == "hallmark") %>% 
        dplyr::select(geneset, symbol) %>% as.data.frame}
    else if(input$categorycd == "c1"){
      selectcd <- msigdf.human %>% 
        filter(category_code == "c1") %>% 
        dplyr::select(geneset, symbol) %>% as.data.frame}
    else if(input$categorycd == "c1"){
      selectcd <- msigdf.human %>% 
        filter(category_code == "c1") %>% 
        dplyr::select(geneset, symbol) %>% as.data.frame}
    else if(input$categorycd == "c1"){
      selectcd <- msigdf.human %>% 
        filter(category_code == "c1") %>% 
        dplyr::select(geneset, symbol) %>% as.data.frame}
    else if(input$categorycd == "c4"){
      selectcd <- msigdf.human %>% 
        filter(category_code == "c4") %>% 
        dplyr::select(geneset, symbol) %>% as.data.frame}
    else if(input$categorycd == "c5"){
      selectcd <- msigdf.human %>% 
        filter(category_code == "c5") %>% 
        dplyr::select(geneset, symbol) %>% as.data.frame}
    else if(input$categorycd == "c6"){
      selectcd <- msigdf.human %>% 
        filter(category_code == "c6") %>% 
        dplyr::select(geneset, symbol) %>% as.data.frame}
    else {
      selectcd <- msigdf.human %>% 
        filter(category_code == "c7") %>% 
        dplyr::select(geneset, symbol) %>% as.data.frame}
    })
  
  # msigtb <- reactive({
    observeEvent(input$action,{
      output$msigtable <- DT::renderDataTable({
        res <- GSEA(symblist(), TERM2GENE = selectcd())
        res@result
      })
    })
    
    observeEvent(input$rest,{
      output$msigtable <- NULL
    })
  # })
  
  # keggtb <- reactive({
    observeEvent(input$action, {
      output$keggoutput <- DT::renderDataTable({
        if (is.null(input$file1))
          return(NULL)
        gseaKEGGres <- gseKEGG(geneList     = genelist(),
                               organism     = 'hsa',
                               nPerm        = 1000,
                               minGSSize    = 10,
                               pvalueCutoff = 1,
                               verbose      = FALSE)
        
        gseaKEGGres@result
        
      })
    })
  # })
  
  
  
  observeEvent(input$downloadimage, {
    pathview(gene.data  = genelist(),
             pathway.id = input$selectedpath,
             species    = "hsa",
             limit      = list(gene=max(abs(genelist())), cpd=1),
             kegg.dir   = "../www") 
  })
  
  # output$downID <- downloadHandler(
  #   filename = function() {
  #     paste("GeneSymb_", Sys.Date(), ".csv", sep="")
  #   },
  #   content = function(file) {
  #     write.csv(symbtb(),file, row.names = F)
  #   }
  # )
  # 
  # output$downMSig <- downloadHandler(
  #   filename = function() {
  #     paste("MSig_", Sys.Date(), ".csv", sep="")
  #   },
  #   content = function(file) {
  #     write.csv(msigtb(), file, row.names = F)
  #   }
  # )
  # 
  # output$downKEGG <- downloadHandler(
  #   filename = function() {
  #     paste("KEGG_", Sys.Date(), ".csv", sep="")
  #   },
  #   content = function(file) {
  #     write.csv(keggtb(), file, row.names = F)
  #   }
  # )
}
