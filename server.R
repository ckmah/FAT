library(shiny)
library(DT)
library(shinyjs)
library(ggplot2)
library(cummeRbund)
library(plotly)
library(gage)

if (!exists("keggAnalysis", mode = "function")) {
  source("functional-analysis.R")
}

# server dev options
options(shiny.autoreload = TRUE,
        shiny.maxRequestSize = 50 * 1000 * 1024 * 2)

shinyServer(function(input, output) {
  # ---------------------------------------- INITIALIZATION ----------------------------------------
  ##### Initialize Disabled Elements
  js$disableTabs()
  disable('viewFile')
  disable('submitButton')
  # set paths
  path = getwd()
  outDir = 'output'
  
  # ---------------------------------------- LOAD DATA ----------------------------------------
  
  observeEvent(input$file, {
    enable('submitButton')
  })
  
  # load Cuffdiff data
  cufflinksData <- reactive({
    progress <- shiny::Progress$new()
    progress$set(message = "Processing data...", value = 0.8)
    on.exit(progress$close())
    
    # WORKING UPLOAD
    zipFile <- input$file
    unzip(zipFile$datapath, exdir = outDir)
    readCufflinks(dir = file.path(path, outDir))
    # TEST DATA FOR NOW
    # readCufflinks(dir = system.file("extdata", package = "cummeRbund"))
    # readCufflinks(dir = 'data/cuffdiff_out_pLKO1_vs_shLPA2/', rebuild = T)
  })
  
  # FPKM W/OUT REPLICATES
  fpkmData <- reactive({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    allfpkm <- fpkmMatrix(genes(cufflinksData()))
    return(cleanUpData(allfpkm))
  })
  
  # FPKM NORMALIZED W/OUT REPLICATES
  fpkmNormData <- reactive({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    allfpkm <- countMatrix(genes(cufflinksData()))
    return(cleanUpData(allfpkm))
  })
  
  # FPKM WITH REPLICATES
  fpkmRepData <- reactive({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    allfpkm <- repFpkmMatrix(genes(cufflinksData()))
    return(cleanUpData(allfpkm))
  })
  
  # FPKM NORMALIZED WITH REPLICATES
  fpkmNormRepData <- reactive({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    allfpkm <- repCountMatrix(genes(cufflinksData()))
    return(cleanUpData(allfpkm))
  })
  
  # Remove genes with 0 counts across all samples
  # Add pseudocount of 1 (default)
  # return log10 of counts
  cleanUpData <- function(dat, pseudocount = 1) {
    dat = dat[rowSums(dat) != 0,]
    dat = dat + 1
    return(log10(dat))
  }
  
  # ---------------------------------------- START ----------------------------------------
  output$dataSetup <- renderUI({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    samp <- cummeRbund::samples(genes(cufflinksData()))
    radioButtons("reference", "Indicate reference sample:",
                 choices = samp)
  })
  
  output$fileList <- renderUI({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    fileNames <- list.files(outDir)
    selectInput("fileDropdown", "View selected file:", fileNames)
  })
  
  # raw data preview
  output$rawDataTable <- DT::renderDataTable({
    # dependency on submitButton
    if (input$submitButton == 0 |
        is.null(input$fileDropdown) | is.null(input$file))
      return()
    
    # enable functions
    js$enableTabs()
    enable('viewFile')
    
    # TODO: hide irrelevant files
    tab <-
      read.table(file = paste(outDir, input$fileDropdown, sep = "/"),
                 header = TRUE)
    
    DT::datatable(
      # style = 'bootstrap',
      tab,
      extensions = list('Scroller'),
      selection = 'none',
      options = list(
        scrollY = 300,
        scroller = TRUE,
        scrollX = TRUE,
        dom = 'frtip',
        pageLength = 100,
        fixedColumns = list(leftColumns = 1)
      )
    )
  })
  
  # ---------------------------------------- QUALITY CONTROL ----------------------------------------
  ##### Disperson Tab #####
  output$dispersionPlot <- renderPlotly({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    # bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') # QUOTES
    
    aes = aes(
      x = ind,
      y = values,
      fill = ind,
      colour = ind
    )
    scale = scale_fill_hue(c = 45, l = 80)
    colors = geom_boxplot(outlier.colour = NULL)
    
    if (input$dispNorm == TRUE) {
      if (input$dispRep == TRUE) {
        gg <-
          ggplot(stack(fpkmNormRepData()), aes) + colors + scale + labs(x = "Samples with Replicates", y = "log10FPKM") + ggtitle("")
      }
      else {
        gg <-
          ggplot(stack(fpkmNormData()), aes) + colors + scale + labs(x = "Samples", y = 'log10FPKM')
      }
    }
    else {
      if (input$dispRep == TRUE) {
        gg <-
          ggplot(stack(fpkmRepData()), aes) + colors + scale + labs(x = "Samples with Replicates", y = "log10FPKM")
      } else {
        gg <-
          ggplot(stack(fpkmData()), aes) + colors + scale + labs(x = "Samples", y = 'log10FPKM')
      }
    }
    ggplotly(gg)
  })
  
  output$dispersionTable <- DT::renderDataTable({
    DT::datatable(
      # style = 'bootstrap',
      repFpkm(genes(cufflinksData())),
      extensions = list('Scroller', 'Buttons'),
      options = list(
        scrollY = 300,
        scroller = TRUE,
        scrollX = TRUE,
        pageLength = 100,
        # pageLength = length(repFpkm(genes(
        #   cufflinksData()
        # ))[[1]]),
        fixedColumns = list(leftColumns = 1),
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
    )
  })
  
  ##### Variability #####
  output$variabilityPlot <- renderPlotly({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    plot <-
      fpkmSCVPlot(genes(cufflinksData()), showPool = input$varRep) + labs(x =
                                                                            'log10FPKM', y = 'CV2')
    
    ggplotly(plot)
    
  })
  
  output$variabilityTable <- DT::renderDataTable({
    DT::datatable(
      # style = 'bootstrap',
      repFpkm(genes(cufflinksData())),
      extensions = list('Scroller'),
      options = list(
        scrollY = 300,
        scroller = TRUE,
        scrollX = TRUE,
        pageLength = 100,
        # pageLength = length(repFpkm(genes(
        #   cufflinksData()
        # ))[[1]]),
        dom = 'frtip',
        fixedColumns = list(leftColumns = 1)
      )
    )
  })
  
  ##### Clustering #####
  output$pcaPlot <- renderPlotly({
    # dependency on submitButton
    if (input$submitButton == 0)
      return()
    
    xAxis <- paste("PC", sep = "", input$clustX)
    yAxis <- paste("PC", sep = "", input$clustY)
    plot <-
      PCAplot(
        genes(cufflinksData()),
        x = xAxis,
        y = yAxis,
        replicates = input$clustRep,
        showPoints = FALSE
        # logMode = input$clustLog
      )
    ggplotly(plot)
  })
  
  output$clustTable <- DT::renderDataTable({
    DT::datatable(
      # style = 'bootstrap',
      repFpkm(genes(cufflinksData())),
      extensions = list('Scroller'),
      options = list(
        scrollY = 300,
        scroller = TRUE,
        scrollX = TRUE,
        pageLength = 100,
        # pageLength = length(repFpkm(genes(
        #   cufflinksData()
        # ))[[1]]),
        dom = 'frtip',
        fixedColumns = list(leftColumns = 1)
      )
    )
  })
  
  # ---------------------------------------- FUNCTIONAL ANALYSIS ----------------------------------------
  ##### Kegg Analysis #####
  kegg.geneset <- reactive({
    if (is.null(input$organism))
      return()
    # get mouse annotations
    if (input$organism == "mmu")  {
      kg.mouse <- kegg.gsets("mouse")
      return(kg.mouse$kg.sets[kg.mouse$sigmet.idx])
    } else {
      # default to human
      data(kegg.gs)
      return(kegg.gs)
    }
  })
  
  ##### KEGG Analysis with GAGE #####
  # Returns list with mapped matrix and gage analysis pathways
  keggResults <- reactive({
    progress <- shiny::Progress$new()
    progress$set(message = "Performing KEGG Analysis...", value = 0.8)
    on.exit(progress$close())
    result <-
      keggAnalysis(fpkmData(),
                   cufflinksData(),
                   kegg.geneset(),
                   input$reference)
    return(result)
  })
  
  mapped <- reactive({
    keggResults()$mapped
  })
  
  gageAnalysis <- reactive({
    keggResults()$kegg.gage
  })
  
  # Process selected geneset
  geneset <- reactive({
    s <- input$kgTable_rows_selected
    if (!length(s))
      return()
    
    return(s[length(s)])
  })
  
  output$test <- renderPrint({
    return(geneset())
  })
  
  genesetDown <- reactive({
    s <- input$keggLessTable_rows_selected
    if (!length(s))
      return()
    
    return(s[length(s)])
  })
  
  output$test <- renderPrint({
    return(geneset())
  })
  ##### STATISTICS #####
  output$keggBubblePlot <- renderPlotly({
    statData <-
      as.data.frame(gageAnalysis()$stats[!is.na(gageAnalysis()$stats[, 1]),])
    statData$set.size <-
      gageAnalysis()$greater[rownames(statData), 'set.size']
    
    g <-
      ggplot(statData, aes(
        x = rownames(statData),
        y = stat.mean,
        size = set.size
      )) + geom_point(
        fill = "#56B4E9",
        colour = "#FFFFFF",
        alpha = 0.5,
        pch = 21
      ) + scale_size_area(max_size = 10) + theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none"
      ) + labs(x = "KEGG Pathway", y = "Mean Statistic theme")
    gp <- plotly_build(g)
    gp$data[[1]]$text <-
      paste(
        rownames(statData),
        paste("stat.mean:", statData$stat.mean),
        paste("set.size", statData$set.size),
        sep = '<br>'
      )
    gp
  })
  
  output$keggStatsTable <- DT::renderDataTable({
    dat <- gageAnalysis()$stats[!is.na(gageAnalysis()$stats[, 1]),]
    # s = input$x1_rows_selected
    DT::datatable(
      dat,
      # style = 'bootstrap',
      extensions = list('Scroller'),
      selection = 'single',
      options = list(
        scrollY = 300,
        scroller = TRUE,
        scrollX = TRUE,
        pageLength = nrow(dat),
        dom = 'frti',
        fixedColumns = list(leftColumns = 1)
      )
    )
  })
  
  ##### UPREGULATED #####
  # HEATMAP
  output$keggHeatmap <- renderPlotly({
    s <- geneset()
    heatmap.g <- mapped()[kegg.geneset()[[s]],]
    heatmap.g <- heatmap.g[!is.na(heatmap.g[1]), ]
    rownames(heatmap.g) <-
      eg2id(rownames(heatmap.g),
            org = input$organism,
            category = "symbol")[, 2]
    
    # A <- heatmap.g[[input$reference]]
    # B <- heatmap.g[[names(heatmap.g)[names(heatmap.g) != input$reference]]]
    # heatmap.fold <- (B - A)/A
    
    # convert to z-scores and melt for plotting
     
    scaled <- scale(as.numeric(as.matrix(melt(heatmap.g))[, 2]))
    heatmap.z <- cbind(as.matrix(melt(heatmap.g))[, 1], scaled)
    rownames(heatmap.z) <- melt(as.matrix(heatmap.g))[['Var1']]
    g <- ggplot(as.data.frame(heatmap.z),
                aes(
                  x = Var2,
                  y = as.factor(Var1),
                  fill = value
                )) + geom_raster() + scale_fill_gradient(low = "red", high = "blue")+ labs(x = "Samples", y = "Genes")
    ggplotly(g)
  })
  
  # SCATTERPLOT
  output$keggScatterPlot <- renderPlotly({
    s <- geneset()
    heatmap.g <- mapped()[kegg.geneset()[[s]],]
    heatmap.g <- heatmap.g[!is.na(heatmap.g[1]), ]
    rownames(heatmap.g) <-
      eg2id(rownames(heatmap.g),
            org = input$organism,
            category = "symbol")[, 2]
    
    g <-
      ggplot(heatmap.g, aes(x = heatmap.g[[input$keggUpS1]],
                            y = heatmap.g[[input$keggUpS2]])) + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("") + labs(x = input$keggUpS1, y = input$keggUpS2)
    gp <- plotly_build(g)
    gp$data[[1]]$text <- rownames(heatmap.g)
    gp
  })
  
  output$keggUpOptions <- renderUI({
    inputPanel(
      selectInput('keggUpS1', 'Sample 1', choices = names(mapped())),
      selectInput('keggUpS2', 'Sample 2', choices = names(mapped()))
    )
  })
  
  output$kgTable <- DT::renderDataTable({
    dat <-
      gageAnalysis()$greater[!is.na(gageAnalysis()$greater[, 1]),]
    DT::datatable(
      dat,
      # style = 'bootstrap',
      extensions = list('Scroller'),
      selection = 'single',
      options = list(
        scrollY = 300,
        scroller = TRUE,
        scrollX = TRUE,
        pageLength = nrow(dat),
        dom = 'frti',
        fixedColumns = list(leftColumns = 1)
      )
    )
  })
  
  ##### DOWNREGULATED #####
  
  output$keggDownHeatmap <- renderPlotly({
    s <- genesetDown()
    heatmap.g <- mapped()[kegg.geneset()[[s]],]
    heatmap.g <- heatmap.g[!is.na(heatmap.g[1]), ]
    rownames(heatmap.g) <-
      eg2id(rownames(heatmap.g),
            org = input$organism,
            category = "symbol")[, 2]
    
    g <- ggplot(melt(as.matrix(heatmap.g)), aes(
      x = Var2,
      y = as.factor(Var1),
      fill = value
    )) + geom_raster() + labs(x = "Samples", y = "Genes")
    ggplotly(g)
  })
  
  # SCATTERPLOT
  output$keggDownScatterPlot <- renderPlotly({
    s <- genesetDown()
    heatmap.g <- mapped()[kegg.geneset()[[s]],]
    heatmap.g <- heatmap.g[!is.na(heatmap.g[1]), ]
    rownames(heatmap.g) <-
      eg2id(rownames(heatmap.g),
            org = input$organism,
            category = "symbol")[, 2]
    
    g <-
      ggplot(heatmap.g, aes(x = heatmap.g[[input$keggUpS1]],
                            y = heatmap.g[[input$keggUpS2]])) + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("") + labs(x = input$keggDownS1, y = input$keggDownS2)
    gp <- plotly_build(g)
    gp$data[[1]]$text <- rownames(heatmap.g)
    gp
  })
  
  output$keggDownOptions <- renderUI({
    inputPanel(
      selectInput('keggDownS1', 'Sample 1', choices = names(mapped())),
      selectInput('keggDownS2', 'Sample 2', choices = names(mapped()))
    )
  })
  
  output$keggLessTable <- DT::renderDataTable({
    # s = input$x1_rows_selected
    dat <- gageAnalysis()$less[!is.na(gageAnalysis()$less[, 1]),]
    DT::datatable(
      dat,
      # style = 'bootstrap',
      extensions = list('Scroller'),
      selection = 'single',
      options = list(
        scrollY = 300,
        scroller = TRUE,
        scrollX = TRUE,
        pageLength = nrow(dat),
        dom = 'frti',
        fixedColumns = list(leftColumns = 1)
      )
    )
  })
})
