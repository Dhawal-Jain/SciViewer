tabBox(title = "",width = NULL,
       tabPanel(title="1. Cluster Details",
                icon = icon("object-group"),
                shiny::HTML("<i>Each point on the plot represents single cells grouped based on their expression
                patterns. t-SNE or UMAP algoriths are used for embedding the cells into low dimensional space. 
                Clusters are defined in 2D space using graph-based methods and labelled using the expression 
                patterns of marker genes.</i><br><br>"),
                fluidRow(column(width=6, 
                                div(style="display: inline-block;vertical-align:top; width: 500px;",
                                    selectInput("sc_celltype1a",label="select celltype",choices = NULL,multiple = T)),
                                tags$br(),
                                tags$br(),
                                plotlyOutput("sc_pcaplot1a") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                tags$br()
                                ),
                         column(width=6, 
                                fluidRow(column(width=6,selectInput("sc_sel_1",label="feature",choices = NULL,multiple = F)),
                                         column(width=6,selectInput("sc_sel_1a",label="sub-feature",choices = NULL,multiple = T)),
                                         style = "height:75px; background-color: white;"
                                         ),
                                tags$br(),
                                plotlyOutput("sc_selplot1a") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                tags$br()
                                ),
                         style = "height:550px; background-color: white;"
                ),
                fluidRow(column(width=6, 
                                tags$br(),
                                htmlOutput("scpca1a_hover"),
                                tags$br()
                                ),
                         column(width=6,
                                tags$br(),
                                htmlOutput("scsel1a_hover"),
                                tags$br()
                                ), 
                         style = "height:50px; background-color: white;"
                ),  
                tags$br(),
                tags$hr(),
                tabsetPanel(type = 'tabs',
                  tabPanel(title="1.1 Cell Counts",
                           icon=icon("signal", lib = "glyphicon"),
                           shiny::HTML("<i>Following barplot summarizes fraction of total cells per identified cell type in the experiment. </i><br><br>"),
                           highchartOutput("sc_barplot",width = "1000px",height = "500px") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                           tags$br()
                  ),
                  tabPanel(title="1.2 Summary",
                           icon=icon("th", lib = "glyphicon"),
                           shiny::HTML("<i>The spreadsheet summarizes cell type proportion in the selected experiment. </i><br><br>"),
                           DT::dataTableOutput("sc_summaryTable"),
                           tags$br()
                  ),
                  tabPanel(title="1.3 QC plots",
                           icon=icon("eye-open", lib = "glyphicon"),
                           shiny::HTML("<i>Quality checks for individual cells are displayed below. </i><br><br>")
                  )
                ),
                tags$br(),
                tags$br()
       ),
       tabPanel(title="2.Single Gene",
                icon = icon("object-group"),
                shiny::HTML("<i>Gene expression for the <i><u>first selected gene</i></u> is summarized on the 
                             cluster map. Select clusters from the dropdown menu or zoom on the left-hand plot for interactive
                             visualization.</i><br><br>"), 
                fluidRow(column(width=6, 
                                div(style="display: inline-block;vertical-align:top; width: 500px;",
                                    selectInput("sc_celltype2a",label="select celltype",choices = NULL,multiple = T)),
                                plotlyOutput("sc_pcaplot2a") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                tags$br()
                                ),
                         column(width=6,
                                div(style="display: inline-block;vertical-align:top; width: 500px;",
                                    selectInput("c1t2expnscale",label="select color-scale",choices = NULL,multiple = F)),
                                plotlyOutput("sc_pcaplot2b") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                tags$br()
                                ),
                         style = "height:550px; background-color: white;"
                ),
                fluidRow(column(width=6, 
                                tags$br(),
                                htmlOutput("scpca2a_hover"),
                                tags$br()
                                ),
                         column(width=6),
                         style = "height:50px; background-color: white;"
                ),
                tags$br(),
                tags$hr(),
                tabsetPanel(type = 'tabs',
                            tabPanel(title="2.1 Expression facets",
                                     icon=icon("signal", lib = "glyphicon"),
                                     shiny::HTML('<i>Following barplot summarizes gene expression across the study features. 
                       Please select one or more attributes for faceting the log-normalized average expression values across selected groups.
                       The error bars represent 95% confidence intervals.
                       By default, values that lie beyond 1.5IQR are labelled as outliers and consequently are ommited. 
                       Thus, if only handful of cells express a gene, the boxplot will likely show median expression
                      at 0 and the few cells will be ommitted from display likely because of their expression falls beyond the 1.5IQR.</i><br><br>'),
                                     radioButtons("t2radio", label = h5("Select graph type"),
                                                  choices = list("Barplot" = 1, "Boxplot (outlier omitted)" = 2), 
                                                  selected = 1,inline = T),
                                     div(style="display: inline-block;vertical-align:top; width: 300px;",
                                         selectInput("sc_sel_2a",label="feature",choices = NULL,multiple = T)),
                                     highchartOutput("sc_barplot_2a",width = "1000px",height = "500px") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                     tags$br()                
                            )
                ),
                tags$br()
       ),
       tabPanel(title="3.Multi Gene",
                icon = icon("object-group"),
                shiny::HTML("<i>Gene expression is summarized on the cluster map. Select clusters from the dropdown menu 
                         or zoom on the top plot for interactive linked-visualization.</i><br><br>"),
                fluidRow(column(width=6,
                                div(style="display: inline-block;vertical-align:top; width: 500px;",
                                    selectInput("sc_celltype3a",label="select celltype",choices = NULL,multiple = T)),
                                plotlyOutput("sc_pcaplot3a") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                tags$br()
                         ),
                         column(width = 6,
                                tags$br(),
                                tags$br(),
                                tags$br(),
                                htmlOutput("scpca3a_hover"),
                                tags$br()
                         ),
                         style = "height:500px; background-color: white;"
                ),
                tags$br(),
                tags$br(),
                fluidRow(style = "height:10px; background-color: white;"),
                tags$br(),
                tags$hr(),
                tabsetPanel(type = 'tabs',
                            tabPanel(title="3.1 Feature plots",
                                     icon=icon("th", lib = "glyphicon"),
                                     shiny::HTML("<b>NOTE:</b> The expression values on the UMAP plot below are not directly comparable. <br>
                  The Application uses aggregation method for summarizing the points within a pixel.Please use these plots to 
                                                 <u>qualitatively</u> compare expression across cell types.<br><br>"),
                                     div(style='width:1200px;overflow-x: scroll;height:800px;overflow-y: scroll;',
                                         uiOutput("sc_pcaplot3b")),
                                     tags$br()
                            ),
                            tabPanel(title="3.2 Dot plot",
                                     icon=icon("th", lib = "glyphicon"),
                                     shiny::HTML("<i>Interactive dotplot summarize average gene expression for all selected genes.
                              Size of the points represent % of cells in the cell type expressing the genes. While
                              color represents expression scale.The views on this tab are linked to the selection of clusters/zoom 
                              on above 2D plot.</i><br><br>"), 
                                     div(style='width:600px;overflow-x: scroll;height:500px;overflow-y: scroll;',
                                         uiOutput("sc_dotplot3a")%>% shinycssloaders::withSpinner(color="#0dc5c1")),
                                     tags$br()
                            ),
                            tabPanel(title="3.3 Correlaton plot",
                                     icon=icon("th", lib = "glyphicon"),
                                     shiny::HTML("<i>This section summarizes correlation between gene expression of the selected genes. 
                <b>Pearson correlation coefficients</b> are displayed in below heatmap. You can select one or 
                multiple cell types from the drop-down menu to further assess the gene expression correlations for selected cell types.<br> 
                <b>Note:</b> At least 2 genes are required for producing the correlation plot.</i><br><br>"),
                                     tags$br(),
                                     fluidRow(column(width = 6,
                                                     selectInput("sc_sel_3a",label="select cell types",choices = NULL,multiple = T)),
                                              column(width = 6,
                                                     radioButtons(inputId = 'dfgacvf',label = 'data selection:',choices = c('zero-ommitted','all'),selected = 'zero-ommitted',inline = T)),
                                              style = "height:70px; background-color: white;"
                                     ),
                                     tags$br(),
                                     uiOutput("sc_corr_3a") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                     tags$br()
                            )
                ),
                tags$br()
       ),
       tabPanel(title="4. Celltype markers",
                icon = icon("object-group"),
                shiny::HTML("<i>Cell type markers are computed using FindMarkers function in Seurat package. 
                Select one or more genes from the table to visualize their expression. Genes with FDR <0.05 
                are displayed.</i><br><br>"),
                fluidRow(
                  column(width = 4,
                         tags$br(),
                         selectInput("sc_celltype4a",label="select celltype",choices = NULL,multiple = F)
                  ),
                  column(width = 4,
                         tags$br(),
                         selectInput("sc_sel_4a",label="select DE test",choices = NULL,multiple = F)
                  ),
                  column(width = 4,
                         sliderInput("slider_4a", label = h5("FDR Cutoff"), min=0,max=0.2, value = 0.05,step = 0.05)
                  ),
                  style = "height:75px; background-color: white;"
                ),
                tags$br(),
                tags$hr(),
                fluidRow(column(width = 6,
                                plotlyOutput("sc_pcaplot4a") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                tags$br()
                        ),
                        column(width=6,
                               #plotlyOutput("sc_volcano") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                               div(style='width:500px;overflow-x: scroll;height:500px;overflow-y: scroll;',
                                   uiOutput("sc_markerdotplot")),
                               tags$br()
                        ),
                        style = "height:500px; background-color: white;"
                ),
                tags$br(),
                fluidRow(
                  column(width = 6,
                         htmlOutput("scpca4a_hover"),
                         tags$br()
                  ),
                  style = "height:50px; background-color: white;"
                ),
                tags$br(),
                DT::dataTableOutput("marker_tab"),
                tags$br(),
                tags$hr(),
                tabsetPanel(type = 'tabs',
                            tabPanel(title="4.1 Feature plots",
                                     icon=icon("th", lib = "glyphicon"),
                                     shiny::HTML("<i><b>NOTE:</b> The expression values on the UMAP plot below are not directly comparable.<br>
                                     The Application uses aggregation method for summarizing the points within a pixel.Please use these plots
                                     to <u>qualitatively</u> compare expression across cell types.</i><br><br>"),
                                     div(style='width:1200px;overflow-x: scroll;height:800px;overflow-y: scroll;',
                                         uiOutput("sc_pcaplot4marker")),
                                     tags$br()
                            ),
                            tabPanel(title="4.2 Enrichments",
                                     icon=icon("th", lib = "glyphicon"),
                                     shiny::HTML("<i>The genes at FDR cutoff as selected above are used for downstream enrichment analyses. The app uses 
                                                 <a href=\"https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html\",target=\"_blank\">
                                                  gProfiler2 </a>for calculating the enrichments and displays.</i><br><br>"),
                                     radioButtons("t4aradio", label = h5("Direction of change"),
                                                  choices = list("Overexpressed" = 1, "Underexpressed" = 2), 
                                                  selected = 1,inline = T),
                                     uiOutput("enrich_marker4a")%>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                     tags$br(),
                                     DT::dataTableOutput("enrichmarker_tab"),
                                     tags$br()
                            )
                ),
                tags$br()
       ),
       tabPanel(title="5. DE",
                icon = icon("object-group"),
                shiny::HTML("<i>Differential Expression (DE) for the genes is computed using pseudobulk method and ajusting for 
                            covariates/surrogate variables -whenever possible. Please contact Data science team, 
                            if you need more information or help in interpreting the results.This section will populated only when
                            the results are computed elsewhere and provided in the input data object.</i><br><br>"),
                fluidRow(
                  column(width = 4,
                         tags$br(),
                         selectInput("sc_celltype5a",label="select celltype/group",choices = NULL,multiple = F)
                  ),
                  column(width = 4,
                         tags$br(),
                         selectInput("sc_sel_5a",label="select DE test",choices = NULL,multiple = F)
                  ),
                  column(width = 4,
                         sliderInput("slider_5a", label = h5("FDR Cutoff"), min=0,max=0.2, value = 0.05,step = 0.05)
                  ),
                  style = "height:75px; background-color: white;"
                ),
                tags$br(),
                tags$hr(),
                fluidRow(column(width = 6,
                                plotlyOutput("sc_pcaplot5a") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                tags$br()
                ),
                column(width=6,
                       div(style='width:500px;overflow-x: scroll;height:500px;overflow-y: scroll;',
                           uiOutput("sc_degdotplot")),
                       tags$br()
                ),
                style = "height:500px; background-color: white;"
                ),
                tags$br(),
                fluidRow(
                  column(width = 6,
                         htmlOutput("scpca5a_hover"),
                         tags$br()
                  ),
                  style = "height:50px; background-color: white;"
                ),
                tags$br(),
                DT::dataTableOutput("deg_tab"),
                tags$br(),
                tags$hr(),
                tabsetPanel(type = 'tabs',
                            tabPanel(title="5.1 Feature plots",
                                     icon=icon("th", lib = "glyphicon"),
                                     shiny::HTML("<i><b>NOTE:</b> The expression values on the UMAP plot below are not directly comparable.<br>
                                     The Application uses aggregation method for summarizing the points within a pixel.Please use these plots
                                     to <u>qualitatively</u> compare expression across cell types.</i><br><br>"),
                                     div(style='width:1200px;overflow-x: scroll;height:800px;overflow-y: scroll;',
                                         uiOutput("sc_pcaplot5deg")),
                                     tags$br()
                            ),
                            tabPanel(title="5.2 Enrichments",
                                     icon=icon("th", lib = "glyphicon"),
                                     shiny::HTML("<i>The genes at FDR cutoff as selected above are used for downstream enrichment analyses. The app uses 
                                                 <a href=\"https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html\",target=\"_blank\">
                                                  gProfiler2 </a>for calculating the enrichments and displays.</i><br><br>"),
                                     radioButtons("t5aradio", label = h5("Direction of change"),
                                                  choices = list("Overexpressed" = 1, "Underexpressed" = 2), 
                                                  selected = 1,inline = T),
                                     uiOutput("enrich_deg5a")%>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                     tags$br(),
                                     DT::dataTableOutput("enrichdeg_tab"),
                                     tags$br()
                            )
                ),
                tags$br()
       )
) 
  
