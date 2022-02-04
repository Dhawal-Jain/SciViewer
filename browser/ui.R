source("vars.R",local = TRUE,encoding = "UTF-8")
options(shiny.maxRequestSize = 100*1024^2)

# javascript code to collapse box
jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

## New theme with shiny dashboard
shinyUI(
  tagList(
    tags$style(HTML(".selectize-input { font-size: 13px; }")),
    tags$style(HTML(".select-dropdown { font-size: 13px; }")),
    tags$style(HTML(".main-sidebar { font-size: 13px; }")),
    tags$head(HTML('
                 <!-- Global site tag (gtag.js) - Google Analytics -->
                 <script async src="https://www.googletagmanager.com/gtag/js?id=UA-18253988-3"></script>
                 <script>
                 window.dataLayer = window.dataLayer || [];
                 function gtag(){dataLayer.push(arguments);}
                 gtag("js", new Date());
                 
                 gtag("config", "UA-18253988-3");
                 </script>')),
    
    dashboardPage(
      
      dashboardHeader(
        titleWidth = "96%",
        title = shiny::HTML(
          "<div style = 'background-color:#ffffff;color: #6a51a3;font-size:30px;font-weight:bold; vertical-align:middle'>
           <img src = 'logo.png' align = 'left'  height = '55px' width = '200px'>
            Single cell Interactive Viewer 
           </div>")
      ),
      
      dashboardSidebar(collapsed = F,
        sidebarMenu(id = 'AppTab',
                    tags$br(),
                    menuItem(text = shiny::HTML("<b>1)</b>   Select Study"),tabName = "scsel"),
                    menuItem(text = shiny::HTML("<b>2)</b>   scStudy"),tabName = "scstudy"),
                    menuItem(text = "...Celltypes",tabName = "",
                             div(style='width:225px;overflow-x: scroll;overflow-y: scroll;',
                                 uiOutput("scclust_legend"))
                             ),
                    menuItem(text = shiny::HTML("<b>3)</b>   Compare"),tabName = "compare"),
                    menuItem(text = shiny::HTML("<b>4)</b>   Help"),tabName = "help"),
                    tags$br(),
                    tags$br(),
                    tags$br(),
                    tags$br(),
                    tags$br()
        )
      ),
      dashboardBody(
        #tags$hr(),
        tags$style(HTML("
           .box.box-solid.box-primary>.box-header {
           color:#fff;
           background:#e5f5e0
          }
          .box.box-solid.box-primary{
           border-bottom-color:#a1d99b;
           border-left-color:#a1d99b;
           border-right-color:#a1d99b;
           border-top-color:#a1d99b;
           background:#e5f5e0
          }")),
        
        shinyDashboardThemes(theme = "onenote"),
        tabItems(
          tabItem(tabName = "scsel",
                  h4("Welcome to SciView!"),
                  tags$br(),
                  shiny::HTML("This visualization platform allows for browsing
                            multiple single cell RNA-seq, CITE-seq as well as TCR-seq studies in an 
                            interactive manner.<br><br>
                            
                            Use the drop-don menu below to select the study of interest. If the study of
                            interest is not listed in the dropdown menu, browse the folder where .db files
                            are located. For additional details on the how to convert single cell studies 
                            (from .h5ad or Seurat Objects) to .db file visit <u><i>SciViewIn</u></i> app.<br> "),
                  tags$br(),
                  fluidRow(column(width=8,textInput("datadir", "browse directory path", NULL)),
                           column(width=4,tags$br(),actionButton("ab573dfg3","Go!",icon = icon("play-circle")))
                  ),
                  verbatimTextOutput('direxists'),
                  shiny::HTML('<i>Please provide the directory path where .db files are available. this path is
                added to the <u> Select scRNA-seq study </u> section below. If left blank, following directories are 
                            available for browsing (/home; .; C:/  ).</i><br>'),
                  tags$br(),
                  shinyFilesButton(id = 'infile', label = '+ Add more studies!',title =  'Please select a dataset', multiple = T),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$hr(),
                  h4("Select scRNA-seq study"),
                  selectInput("sc_study",label="",choices = NULL,multiple = F),
                  tags$br(),
                  tags$hr(),
                  h4("Details of selected study"),
                  shiny::HTML("<i> After the study is selected, the genes section on sidebar menu will populate.
                            Please visit the <u>scStudy</u> menu for exploring the study. Also, use <u>Compare</u> 
                            menu to compare the selected study with other available studies. </i>"),
                  tags$br(),
                  tags$br(),
                  htmlOutput("sc_study_attribs"),
                  tags$br()
          ),## End of scsel
          tabItem(tabName = "scstudy",
                  fluidRow(column(width=4,
                                  selectizeInput("geneA",label="select gene/s of interest",choices = NULL,multiple = T),
                           ),
                           column(width=4,
                                  tags$br(),
                                  actionButton("geneAgo","Make selection!",icon = icon("play-circle"))
                           ),
                           column(width=4)),
                  tags$br(),
                  source(file = 'ui_sc.R',local = T,encoding = "UTF-8")$value,
                  tags$br()
          ),## End of scstudy
          tabItem(tabName = "compare",
                  fluidRow(
                    column(width=4,
                           tags$br(),
                           htmlOutput('firstsc')),
                    column(width=4,
                           selectInput("comp_study",label="select study to compare",choices = NULL,multiple = F)),
                    column(width=1,
                           tags$br(),
                           actionButton("t2action1","Hit run!",icon = icon("play-circle")))
                  ),
                  source(file = 'ui_compare.R',local = T,encoding = "UTF-8")$value,
                  tags$br()  
          ), ## End of compare
          tabItem(tabName = "help",
                  source(file = 'ui_help.R',local = T,encoding = "UTF-8")$value
          ) ## End of help
        ) # end of tabItems
      ) # end of dashboard body
    ),# end of dashboard
    
    tags$footer(align="center",
                tags$hr(),
                shiny::HTML(
                  "<div vertical-align:middle'>
           <img src = 'PDD_logo.png' align = 'left'  height = '100px' width = '130px'>
           </div>"),
                tags$p("Copyright (c) 2020",color="black"), 
                tags$a(" Bayer US LLC.; ", href = "http://www.bayer.us/",target="_blank"), 
                tags$a(" PDD labs, Strategic collaboration - Speciality lung diseases", href = "https://www.brighamhealth.org/home", target="_blank"), 
                tags$p("Version 1.0"),
                #tags$p("App logo created using www.designevo.com"),
                style = "bottom:0;width:100%;color: #6a51a3;padding: 10px;
                background-color: #ffffff;z-index: 1000;"
    )
  )  #End of taglist
) #End




