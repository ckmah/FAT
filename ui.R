library(shiny)
library(DT)
library(shinyjs)
library(cummeRbund)
library(plotly)

##### Extending shinjs #####
jscode <- "
shinyjs.disableTabs = function(name) {
var tab = $('.nav li a[data-toggle=tab]');
tab.bind('click.tab', function(e) {
e.preventDefault();
return false;
});
tab.addClass('disabled');
}

shinyjs.enableTabs = function(name) {
var tab = $('.nav li a[data-toggle=tab]');
tab.unbind('click.tab');
tab.removeClass('disabled');
}

"
css <- "
.nav li a.disabled {
background-color: #aaa !important;
color: #333 !important;
cursor: not-allowed !important;
border-color: #aaa !important;
}"


# ---------------------------------------- START ----------------------------------------
startPage <- tabPanel(
  "Start",
  useShinyjs(),
  extendShinyjs(text = jscode),
  inlineCSS(css),
  fixedPage(
    fileInput('file',
              'Upload Cuffdiff Analysis Output', accept = c('.zip')),
    selectInput('organism', 'Select organism', c("Human" = "hsa", "Mouse" = "mmu")),
    actionButton('submitButton', 'Submit'),
    uiOutput('dataSetup'),
    uiOutput('fileList'),
    DT::dataTableOutput('rawDataTable')
  )
)

# ---------------------------------------- QUALITY CONTROL ----------------------------------------
qcPage <- tabPanel(
  "Quality Control",
  useShinyjs(),
  extendShinyjs(text = jscode),
  inlineCSS(css),
  fixedPage(
    tabsetPanel(
      type = 'tabs',
      position = 'above',
      tabPanel(
        'Dispersion',
        fluidRow(column(
          4,
          inputPanel(
            checkboxInput('dispNorm', 'VIEW NORMALIZED', value = FALSE),
            checkboxInput('dispRep', 'SHOW REPLICATES', value = TRUE)
          )
        ),
        column(8, plotlyOutput('dispersionPlot'))),
        DT::dataTableOutput('dispersionTable')
      ),
      tabPanel(
        'Variability',
        fluidRow(column(4,
                        inputPanel(
                          checkboxInput('varRep', 'COMBINE SAMPLES', value = FALSE)
                        )),
                 column(8, plotlyOutput('variabilityPlot'))),
        DT::dataTableOutput('variabilityTable')
      ),
      tabPanel(
        'Clustering',
        fluidRow(column(
          4,
          inputPanel(
            checkboxInput('clustRep', 'SHOW REPLICATES', value = TRUE),
            # checkboxInput('clustLog', 'LOG10 TRANSFORM FPKM', value = FALSE),
            numericInput(
              'clustX',
              'X-axis PC',
              value = 1,
              min = 1,
              max = 3
            ),
            numericInput(
              'clustY',
              'Y-axis PC',
              value = 2,
              min = 1,
              max = 3
            )
            
            
          )
        ),
        column(8, fluidRow(
          plotlyOutput('pcaPlot')
        ))),
        DT::dataTableOutput('clustTable')
      )
    )
  )
)

# ---------------------------------------- FUNCTIONAL ANALYSIS ----------------------------------------
funcPage <- tabPanel(
  'Functional Analysis',
  useShinyjs(),
  extendShinyjs(text = jscode),
  inlineCSS(css),
  fixedPage(
    tabsetPanel(
      type = 'tabs',
      position = 'above',
      tabPanel(
        'Statistics',
        plotlyOutput('keggBubblePlot'),
        DT::dataTableOutput('keggStatsTable')
      ),
      tabPanel(
        'Upregulated',
        fluidRow(
          tabsetPanel(
            type = 'tabs',
            position = 'above',
            tabPanel('Heatmap', plotlyOutput('keggHeatmap')),
            tabPanel('Scatter Plot',
                     column(7, plotlyOutput('keggScatterPlot')),
                     column(5, uiOutput('keggUpOptions')))
          )
        ),
        DT::dataTableOutput('kgTable')
      ),
      tabPanel(
        'Downregulated',
        fluidRow(
          tabsetPanel(
            type = 'tabs',
            position = 'above',
            tabPanel('Heatmap', plotlyOutput('keggDownHeatmap')),
            tabPanel('Scatter Plot',
                     column(7, plotlyOutput('keggDownScatterPlot')),
                     column(5, uiOutput('keggDownOptions')))
          )
        ),
        DT::dataTableOutput('keggLessTable')
      )
      
    )
  )
)

shinyUI(fluidPage(
  ##### Material Bootstrap theme and custom CSS #####
  includeCSS("www/bootstrap.min.css"),
  includeCSS("www/styles.css"),
  
  navbarPage('FAT: RNA-seq Functional Analysis Tool',
             startPage,
             qcPage,
             funcPage)
))