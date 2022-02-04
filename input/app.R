source("input_functions.R",local = TRUE,encoding = "UTF-8")
source("vars.R",local = TRUE,encoding = "UTF-8")


##----------------- extra code
# javascript code to collapse box
jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"


##----------------- application
server <- function(input, output,session) {
  VARS <- reactiveValues(
    roots = ROOTS,
    outdir=NULL
  )
  observeEvent({input$ab573dfg3},{
    if(isTruthy(input$ab573dfg3)){
      if(dir.exists(input$datadir)){
        VARS$roots = ROOTS
        inputdir=ifelse(isTRUE(grep("/$",input$datadir)),
                        as.character(input$datadir),
                        as.character(paste0(input$datadir,"/")))
        names(inputdir) = 'inputdir'
        cat("inputdir:",inputdir,"\n")
        VARS$roots = append(inputdir,VARS$roots)
        output$direxists <- NULL
        rm(inputdir)
      }else{
        output$direxists <- renderPrint({ cat('Directory does not exits!\n') })
      }
    }
  })
  observeEvent({VARS$roots},{
    if(isTruthy(VARS$roots)){
      shinyFileChoose(input,'infile', session=session,roots=VARS$roots, 
                      filetypes=c('', 'rds','h5ad'))
    }
  })
  
  variables = reactiveValues(
    so = NULL,
    redmap = NULL,
    si_markerfile = NULL,
    si_defile = NULL,
    simarker = NULL,
    mat = NULL,
    meta = NULL,
    db = NULL,
    studyname = NULL,
    descr = NULL,
    tissue = NULL,
    pmid = NULL,
    geo = NULL,
    status = NULL,
    rating = NULL,
    run = NULL,
    cell_type = NULL,
    barcode = NULL,
    clusterx = NULL,
    clustery = NULL,
    donor = NULL,
    disease = NULL,
    de=NULL,
    covariates = NULL,
    contvars = NULL,
    catvars = NULL
  )
  
  observeEvent({input$infile},{
    inFile <- parseFilePaths(roots=VARS$roots, input$infile)
    if(isTruthy(inFile$datapath)){
      variables$so = NULL
      variables$redmap = NULL
      variables$si_markerfile = NULL
      variables$si_defile = NULL
      variables$simarker = NULL
      variables$mat = NULL
      variables$meta = NULL 
      variables$db = NULL
      variables$studyname = NULL
      variables$descr = NULL
      variables$tissue = NULL
      variables$pmid = NULL
      variables$geo = NULL
      variables$status = NULL
      variables$rating = NULL
      variables$run = NULL
      variables$cell_type = NULL
      variables$barcode = NULL
      variables$clusterx = NULL
      variables$clustery = NULL
      variables$donor = NULL
      variables$de = NULL
      variables$disease = NULL
      variables$covariates = NULL
      variables$contvars = NULL
      variables$catvars = NULL
      
      ext <- tools::file_ext(unname(inFile$datapath))
      VARS$outdir <- dirname(unname(inFile$datapath))
      message("observed file extention: ",ext,"\n")
      validate(need(ext %in% c("h5ad","rds"), "Please upload a correct file!"))
      output$direxists <- renderPrint({ cat('outdir: ',VARS$outdir) })
      so = NULL
      if(ext=='rds'){
        scProg <- shiny::Progress$new()
        on.exit(scProg$close())
        scProg$set(message = "reading the file..", value = 0)
        so = readRDS(unname(inFile$datapath))
        scProg$set(message = "done!", value = 1)
      }else if(ext =='h5ad'){
        scProg <- shiny::Progress$new()
        on.exit(scProg$close())
        scProg$set(message = "reading the file..", value = 0)
        so <- local({
          ad <- anndata::read_h5ad(filename = unname(inFile$datapath))
          meta <- ad$obs
          meta$V1 <- meta$V2 <- NULL
          mat <- ad$X
          qq <- dimnames(mat)
          mat <- as(mat,'matrix.csr')
          mat <- as(mat,'dgCMatrix')
          dimnames(mat) <- qq
          rm(qq)
          mat <- t(mat)
          cat(dim(mat)[1],"\n")
          so = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0,meta.data = meta)
          scvis <- ad$obsm
          if(length(scvis)>1) message('Multiple data reductions found. Will populate them.\n')
          for(i in 1:length(scvis)){
            message(" ..adding coordinates: ",gsub("^X_","",names(scvis)[i]),"\n")
            so[[gsub("^X_","",names(scvis)[i])]] = CreateDimReducObject(embeddings = scvis[[i]], key = gsub("^X_","",names(scvis)[i]))
          }
          so
        })
        scProg$set(message = "done!", value = 1)
      }else{
        stop('Input data is not the correct format.\n')
      }
      if(!is.null(so)){
        variables$meta <- local({
          coords <- data.frame(SAMPID=rownames(so@meta.data))
          for(f in names(so@reductions)){
            qq <- as.data.frame(so@reductions[[f]]@cell.embeddings)
            qq <- qq[,1:2]
            names(qq) <- c('V1','V2')
            names(qq) <- paste0(names(qq),'_',f)
            coords <- cbind(coords, qq)
            rm(qq)
          }
          my.meta <- so@meta.data
          my.meta <- cbind(my.meta,coords)
          my.meta
        })
        output$redmap <- renderUI({
          selectInput("redmap",label="Default reduction map to display",choices = names(so@reductions),multiple = F)
        })
        variables$so = so
      }
    }
  })
  
  observeEvent({input$infile_marker},{
    inFile <- parseFilePaths(roots=VARS$roots, input$infile_marker)
    ext <- tools::file_ext(unname(inFile$datapath))
    validate(need(ext %in% c("txt","csv"), "Please upload txt/csv file!"))
    cf <- NULL
    if(ext=='txt'){
      cf = read.delim(unname(inFile$datapath),as.is = T,colClasses = "#",sep="\t")
    }else if(ext == 'csv'){
      cf = read.csv2(unname(inFile$datapath),as.is = T,colClasses = "#")
    }
    variables$si_markerfile = cf
  })
  
  observeEvent({input$infile_de},{
    inFile <- parseFilePaths(roots=VARS$roots, input$infile_de)
    ext <- tools::file_ext(unname(inFile$datapath))
    validate(need(ext %in% c("txt","csv"), "Please upload txt/csv file!"))
    cf <- NULL
    if(ext=='txt'){
      cf = read.delim(unname(inFile$datapath),as.is = T,colClasses = "#",sep="\t")
    }else if(ext == 'csv'){
      cf = read.csv2(unname(inFile$datapath),as.is = T,colClasses = "#")
    }
    variables$si_defile = cf
  })
  
  
  output$vars <- renderUI({
    orderInput('variables', 'Available variables in the meta-data file', items = colnames(variables$meta),
               item_class = 'info',
               connect = c('cell_type','donor','disease','covariates','clusterx','clustery')) #'barcode',
  })
  output$metatable = DT::renderDataTable({
    DT::datatable(variables$meta,rownames = F,
                  options = list(searching = TRUE,pageLength = 5,
                                 scrollX = TRUE,scrollCollapse = TRUE))
  })
  output$metasummary = renderPrint({
    summary(variables$meta)
  })
  output$cell_type <- renderUI({
    orderInput('cell_type', 'Cell cluster annotation column', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables')
  })
  output$clusterx <- renderUI({
    orderInput('clusterx', 'X-coordinate of reduction method', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables')
  })
  output$clustery <- renderUI({
    orderInput('clustery', 'Y-coordinate of reduction method', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables')
  })

  output$vars1 <- renderUI({
    qq <- colnames(variables$meta)
    if(isTruthy(input$cell_type_order[[2]])){
      qq <- qq[qq!=input$cell_type_order[[2]]]
    }
    if(isTruthy(input$clusterx_order[[2]])){
      qq <- qq[qq!=input$clusterx_order[[2]]]
    }
    if(isTruthy(input$clustery_order[[2]])){
      qq <- qq[qq!=input$clustery_order[[2]]]
    }
    orderInput('variables1', 'Available variables in the meta-data file', items = qq,
               item_class = 'warning',
               connect = c('contvars','catvars'))
  })
  output$contvars <- renderUI({
    orderInput('contvars', 'Continuous variables in the data, that you want to use in the display', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables1')
  })
  output$catvars <- renderUI({
    orderInput('catvars', 'Categorical variables in the data, that you want to use in the display', items = input$cell_type_order[[2]], placeholder = 'Drag item/s here...', connect = 'variables1')
  })
  
  output$vars2 <- renderUI({
    orderInput('variables2', 'Available variables in the meta-data file', items = colnames(variables$meta),
               item_class = 'info',
               connect = c('disease','donor','covariates')) #'barcode',
  })
  output$donor <- renderUI({
    orderInput('donor', 'Donor annotation column', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables2')
  })
  output$disease <- renderUI({
    orderInput('disease', 'Disease/Test/Experiment annotation column', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables2')
  })
  output$covariates <- renderUI({
    orderInput('covariates', 'Variables for computing marker genes', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables2')
  })
  
  
  ##--------- calculations
  observeEvent({input$run},{
    if(isTruthy(input$run)){
      
      if(isTruthy(input$db)) { variables$db = input$db}
      if(isTruthy(input$studyname)) { variables$studyname = input$studyname}
      if(isTruthy(input$redmap)) { variables$redmap = input$redmap}
      if(isTruthy(input$descr)) { variables$descr = input$descr}
      if(isTruthy(input$tissue)) { variables$tissue = input$tissue}
      if(isTruthy(input$pmid)) { variables$pmid = input$pmid}
      if(isTruthy(input$geo)) { variables$geo = input$geo}
      if(isTruthy(input$status)) { variables$status = input$status}
      if(isTruthy(input$rating)) { variables$rating = input$rating}
      if(isTruthy(input$de)) { variables$de = input$de}
      if(isTruthy(input$simarker)) { variables$simarker = input$simarker}
      if(isTruthy(input$cell_type_order[[2]])) { variables$cell_type = input$cell_type_order[[2]][1]}
      if(isTruthy(input$clusterx_order[[2]])) { variables$clusterx = input$clusterx_order[[2]][1]}
      if(isTruthy(input$clustery_order[[2]])) { variables$clustery = input$clustery_order[[2]][1]}
      if(isTruthy(input$donor_order[[2]])) { variables$donor = input$donor_order[[2]][1]}
      if(isTruthy(input$disease_order[[2]])) { variables$disease = input$disease_order[[2]][1]}
      if(isTruthy(input$covariates_order[[2]])) { variables$covariates = paste0(input$covariates_order[[2]],collapse = ",")}
      if(isTruthy(input$contvars_order[[2]])) { variables$contvars = paste0(input$contvars_order[[2]],collapse = ",")}
      if(isTruthy(input$catvars_order[[2]])) { variables$catvars = paste0(input$catvars_order[[2]],collapse = ",")}
      

      output$finalSel <- renderText({
        paste0("<b> Output database: </b> ",variables$db," </font><br>",
               "<b> Study identifier: </b> ",variables$studyname," </font><br>",
               "<b> Description: </b> ",variables$descr," </font><br>",
               "<b> Tissue source: </b> ",variables$tissue," </font><br>",
               "<b> PMID: </b> ",variables$pmid," </font><br>",
               "<b> GEO: </b> ",variables$geo," </font><br>",
               "<b> Availability status: </b> ",variables$status," </font><br>",
               "<b> User rating: </b> ",variables$rating," </font><br>",
               "<b> Whether to log normalize the data: </b> ",variables$lognorm," </font><br>",
               "<b> ColumnID for cell types: </b> ",variables$cell_type," </font><br>",
               "<b> ColumnID for cell barcodes: </b> ",variables$barcode," </font><br>",
               "<b> ColumnID for donor: </b> ",variables$donor," </font><br>",
               "<b> ColumnID for disease: </b> ",variables$disease," </font><br>",
               "<b> ColumnID for X/Y cluster coordinates: </b> ",variables$clusterx,",",variables$clustery," </font><br>",
               "<b> ColumnID for covariates in DE modeling: </b> ",variables$covariates," </font><br>",
               "<b> ColumnID for continuous variables: </b> ",variables$contvars," </font><br>",
               "<b> ColumnID for categorical variables: </b> ",variables$catvars," </font><br>"
        )
      })
      
      ##---------------- rename some columns
      #if(isTruthy(variables$cell_type)){
      #  names(variables$meta) <- ifelse(names(variables$meta)==variables$cell_type,as.character("cell_type"),as.character(names(variables$meta)))
      #}
      #if(isTruthy(variables$clusterx)){
      #  names(variables$meta) <- ifelse(names(variables$meta)==variables$clusterx,as.character("V1"),as.character(names(variables$meta)))
      #}
      #if(isTruthy(variables$clustery)){
      #  names(variables$meta) <- ifelse(names(variables$meta)==variables$clustery,as.character("V2"),as.character(names(variables$meta)))
      #}
      #if(isTruthy(variables$catvars)){
      #  qq <- unlist(strsplit(variables$catvars,","))
      #  qq <- ifelse(qq==variables$cell_type,as.character("cell_type"),as.character(qq))
      #  variables$catvars <- paste0(qq,collapse = ",")
      #  rm(qq)
      #}
      
      ##------------------------
      write.scstudy2.sqlitedb(so = variables$so,
                              db_address=variables$db,
                              StudyName=variables$studyname,
                              Celltype = variables$cell_type,
                              Reduction_map = variables$redmap,
                              Donors_VariableName = variables$donor,
                              DE_Calc = variables$de,
                              DE_Precomputed=variables$si_defile,
                              Disease_VariableName=variables$disease,
                              Marker_Calc = variables$simarker,
                              Marker_Covariates = variables$covariates,
                              Marker_Precomputed=variables$si_markerfile,
                              StudyDescr=variables$descr,
                              Tissue=variables$tissue,
                              PMID=variables$pmid,
                              GEO=variables$geo,
                              StudyStatus=variables$status,
                              StudyRating=variables$rating,
                              Continuous_Vars=variables$contvars,
                              Categorical_Vars=variables$catvars,
                              OUTDIR = VARS$outdir)
      tryCatch({
        sendSweetAlert(session = session,
                       title = "DONE",
                       text = "Generated SQLite database file for provided study!",
                       type = "success")
      },error = function(e){
        sendSweetAlert(session = session,
                       title = "Failed",
                       text = "Unable to generate database file at this point!",
                       type = "error")
      })
      
      
    }
  })
}

ui <- fluidPage(
  shiny::HTML(
    "<div style = 'background-color:#ffffff;color: #6a51a3;font-size:30px;font-weight:bold; vertical-align:middle'>
     <img src = 'logo.png' align = 'left'  height = '55px' width = '200px'>
        Single-cell Interactive Viewer: Data Input 
     </div>"
  ),
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
          }")
  ),
  tags$hr(),
  tags$br(),
  h2('1.Input file'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  shiny::HTML('<i>The application accepts scRNA-seq files in .h5ad or .rds formats.</i><br> 
              <b><i>.h5ad file: </b></i> should contain at least the RNA expression layer with metadata and cluster coordinate information.<br>
              <b><i>.rds file: </b></i> this is a seurat object with count data, metadata and cluster coordinates<br>'),
  tags$br(),
  
  shinyFilesButton('infile', 'File Select', 'Please select a dataset', FALSE),
  actionButton("filenotfound", "Can't locate the file",icon =icon("question-circle"),width = '300px',class = "btn-warning"),
  
  conditionalPanel(
    condition = "input.filenotfound>'0'",
    box(id="sdf",title = "", 
        solidHeader = T,collapsed = F,width = 12,status = 'primary',collapsible = T,
        shiny::HTML('<i>please provide the directory path where .h5ad/.rds files are available. this path is added to the file browser option above.<br>
              if left blank, following directories are available for browsing as default (<b>/home </b>; </b> .</b>; <b>C:/</b>  ).</i><br>'),
        fluidRow(column(width=8,textInput("datadir", "directory path", NULL)),
                 column(width=2,tags$br(),actionButton("ab573dfg3","Go!",icon = icon("play-circle")))
        ),
        verbatimTextOutput('direxists')
    )
  ),
    
  tags$br(),
  
  tags$br(),
  tags$br(),
  
  h2('2. Sneakpeak into the metadata'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  conditionalPanel(
    condition = "input.infile>'0'",
    box(id="m1",title = "metadata table", 
        solidHeader = T,collapsed = T,width = 12,status = 'primary',collapsible = T,
        DT::dataTableOutput("metatable")
    ),
    box(id="m2",title = " metadata summary", 
        solidHeader = T,collapsed = T,width = 12,status = 'primary',collapsible = T,
        verbatimTextOutput('metasummary')
    ),
    tags$br(),
    tags$br(),
  ),
  tags$br(),
  
  h2('3.Defining variables'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  uiOutput("vars"), 
  tags$br(),
  fluidRow(column(width = 3,uiOutput("cell_type"),tags$br()),
           column(width = 3,uiOutput("redmap"),tags$br()),
           column(width = 3),
           column(width = 3)
           #column(width = 3,uiOutput("clusterx"),tags$br()),
           #column(width = 3,  uiOutput("clustery"),tags$br())
  ),
  tags$br(),
  
  h4('3.1.Continuous/Categorical variables'),
  tags$hr(style = "border-top: 1px solid #d9d9d9;"),
  shiny::HTML('<i> continuous variables are grouped into 5 distinct groups. while, for each categorical variable a separate color is assigned</i>'),
  uiOutput("vars1"),
  tags$br(),
  fluidRow(column(width = 6,uiOutput("contvars"),tags$br()),
           column(width = 6,uiOutput("catvars"),tags$br())
  ),
  tags$br(),
  
  
  h2('4.Calculating markers/ differential expression'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  shiny::HTML('<i>please note calculating the markers or differentially expressed genes can take a significant amount of time.
              </i><br><br>'),
  tags$br(),
  uiOutput("vars2"), 
  tags$br(),
  shiny::HTML("<i>Markers are computed using <u>FindMarkers</u> function from Seurat package.
                           The default non-parametric Wilcox test is used. The cells can be grouped by cell type, disease or any other covariate/s.
                           The provided function will group the cells based on a selected covariate and compute differences against the rest of the cells.</i><br>"),
  fluidRow(column(width = 4,  
                  selectInput("simarker",label="Compute Marker genes",choices = c('TRUE','FALSE'),multiple = F)
           ),
           column(width = 4,
             conditionalPanel(
               condition = "input.simarker == 'TRUE'",
               uiOutput("covariates"),
               shiny::HTML("<br><br>")
             )
           ),
           column(width = 4)
  ),
  conditionalPanel(
    condition = "input.simarker == 'FALSE'",
    shiny::HTML("<i> Do you consider uploading the user defined markers? If yes, please use the file link below to upload the marker data. 
                 You can also <b> skip</b> this step.</i><br>"),
    shinyFilesButton('infile_marker', 'File Select', 'Please select pre-computed result file', TRUE),
    tags$br()
  ),
  shiny::HTML("<br><br>"),

  shiny::HTML("<i>Differentially expressed genes are computed using <u>FindMarkers</u> function from Seurat package.
                           The default non-parametric Wilcox test is used. The cells can be grouped by cell type, disease or any other covariate/s.
                           The provided function will group the cells based on a selected cell types and compute differences using the user-provided contrast.</i><br>"),
  fluidRow(column(width = 4,  
                  selectInput("de",label="Compute differential expression",choices = c('TRUE','FALSE'),multiple = F)
           ),
           column(width = 4,
                  conditionalPanel(
                    condition = "input.de == 'TRUE'",
                    uiOutput("disease"),
                  )
           ),
           column(width = 4)
           ),
  conditionalPanel(
    condition = "input.de == 'FALSE'",
    shiny::HTML("<i> Do you consider uploading the user defined differential genes? If yes, please use the file link below to upload the data. 
                 You can also <b> skip</b> this step.</i><br>"),
    shinyFilesButton('infile_de', 'File Select', 'Please select pre-computed result file', TRUE),
    tags$br()
    #fluidRow(column(width = 4,uiOutput("donor"),tags$br()),
    #         column(width = 4,uiOutput("disease"),tags$br()),
    #         column(width = 4,uiOutput("covariates"),tags$br())
    #)
  ),
  #verbatimTextOutput('covariatesout'),
  tags$br(),
  
  
  h2('5.Study details'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  fluidRow(column(width = 6,textInput("db", "output database name", NULL)),
           column(width = 6,textInput("studyname", "Study name", NULL))
  ),
  tags$br(),
  
  h4('5.1 Optional study details'),
  tags$hr(style = "border-top: 1px solid #d9d9d9;"),
  shiny::HTML('<i> these details are listed next to the study in the visualizations.</i>'),
  fluidRow(column(width = 6,textInput("descr", "Brief description of the study/Abstract", NULL)),
           column(width = 6,textInput("tissue", "Tissue used for the scRNA-seq", NULL))
  ),
  tags$br(),
  fluidRow(column(width = 6,textInput("pmid", "If the study is published, then the PMID #", NULL)),
           column(width = 6,textInput("geo", "If the study is published, link to data (e.g. GEO)", NULL))
  ),
  tags$br(),
  fluidRow(column(width = 6,
                  selectInput("status",label="Select status of the study",choices = c('Internal','Publicly Available'),multiple = F)),
           column(width = 6,
                  selectInput("rating",label="Select your rating for the study",choices = c('High','Medium','Poor'),multiple = F))
  ),
  tags$br(),
  tags$br(),
  
  tags$hr(style = "border-top: 1px solid #000000;"),
  fluidRow(column(width = 4),
           column(width = 4,  actionButton("run", "Go!",icon =icon("play"),width = '300px',class = "btn-warning")),
           column(width = 4)
           ),
  tags$hr(style = "border-top: 1px solid #000000;"),
  
  tags$br(),
  h2('6.Run details'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  htmlOutput("finalSel"),
  
  
  tags$br()
  
)

shinyApp(ui, server)