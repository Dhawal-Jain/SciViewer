

tabBox(title = "",width = NULL,
       tabPanel(title="Study info",
                icon = icon("info"),
                tags$br(), 
                fluidRow(
                  htmlOutput("sc_study_attribs1"),
                  tags$br(),
                  tags$br(),
                  htmlOutput("comp_study_attribs"),
                  tags$br()
                ),
                tags$br()
       ),
       tabPanel(title="Comparison",
                icon = icon("object-group"),
                tags$br(), 
                shiny::HTML("<i>Two dimensional representation of the cluster plots are displayed side-by-side.</i><br><br>"),
                
                fluidRow(column(width=6,
                                selectInput("sc_celltype21",label="select celltype",choices = NULL,multiple = T)
                         ),
                         column(width=6,
                                selectInput("sc_celltype22",label="select celltype",choices = NULL,multiple = T)
                         ),
                         style = "height:70px; background-color: white;"
                ),
                tags$br(),
                fluidRow(column(width=6,
                                plotlyOutput("sc_pcaplot21") %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                         column(width=6,
                                plotlyOutput("sc_pcaplot22") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
                ),
                tags$br(),
                fluidRow(column(width=6,
                                tags$br(),
                                htmlOutput("scpca21_hover"),
                                tags$br()
                          ),
                         column(width=6,
                                tags$br(),
                                htmlOutput("scpca22_hover"),
                                tags$br()
                         ),
                         style = "height:100px; background-color: white;"
                ),
                tags$br(),
                tags$br(),
                tags$hr(),
                tabsetPanel(type = 'tabs',
                            tabPanel(title="Expression",
                                     icon=icon("th", lib = "glyphicon"),
                                     tags$br(),
                                     fluidRow(column(width=6, 
                                                     plotlyOutput("sc_pcaplot23") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                              ),
                                              column(width=6, 
                                                     plotlyOutput("sc_pcaplot24") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                              )
                                     ),
                                     tags$br()
                              )
                 )
       )
)