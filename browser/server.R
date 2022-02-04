source("vars.R",local = TRUE,encoding = "UTF-8")
source("plot_functions.R",local = TRUE,encoding = "UTF-8")
source("MMdbLib.R",local = TRUE,encoding = "UTF-8")

##----------------------------------------------
warningMessage <- paste0("<font color=\"#e6550d\">
                           <br><br><br>
                           Requested plot/table is unavailable. <br>
                           Reasons:<br> 
                           1) Gene not present in the underlying database <br>
                           2) Gene locus is unavaialbe in Ensembl GRCh37 transcript database. <br>
                           3) A technical glitch with plot rendering? <br>
                           Try refreshing the page or contact the Data science team/ developers for details.
                           </font><br><br><br>")

##----------------------------------------------
js <- "
    function(el, x, inputName){
      var id = el.getAttribute('id');
      var gd = document.getElementById(id);
      var xScrollPos = el.scrollLeft || document.documentElement.scrollLeft;
      var yScrollPos = el.scrollTop || document.documentElement.scrollTop;
      var d3 = Plotly.d3;
      Plotly.plot(id).then(attach);
        function attach() {
          var xaxis = gd._fullLayout.xaxis;
          var yaxis = gd._fullLayout.yaxis;
          var l = gd._fullLayout.margin.l;
          var t = gd._fullLayout.margin.t;
          var coordinates = [null, null]
          
          gd.addEventListener('mousemove', 
          function(evt) {
              var coordinates = [xaxis.p2c(evt.pointerX-l), 
                                 yaxis.p2c(evt.pointerY-t)];
              Shiny.setInputValue(inputName, coordinates);
          });

        };
  }
  "

##----------------------------------------------

server <- function(input, output, session) {
  
  ROOTS=c(workdir='.',
          datadir='C:/',
          home='/home/',
          shinydata='/srv/shiny-server/data/')
  
  ###---- 0) meta
  VARS <- reactiveValues(
    connGenes = list(),
    connList = list(),
    dbb = NULL,
    roots = ROOTS,
    sc_studies = c(),
    availGenes = NULL
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
                      filetypes=c('', 'db'))
    }
  })
  
  
  ###---  1) genes database
  VARS$connGenes <- ({
    RSQLite::dbConnect(RSQLite::SQLite(),GENEFILE)
  })
  
  ###---- 2) set the sqlite databases with access 
  VARS$dbb <- ({
    tryCatch({
      read.delim(file = 'DBs.txt',header = F,comment.char = "#",stringsAsFactors = F)$V1
    },error= function(e){
    })
  })
  observeEvent({input$infile},{
    VARS$dbb = NULL
    VARS$dbb = tryCatch({
      read.delim(file = SCDBFILE,header = F,comment.char = "#",stringsAsFactors = F)$V1
    },error= function(e){
    })
    inFile <- parseFilePaths(roots=VARS$roots, input$infile)
    if(isTruthy(inFile$datapath)){
      ext <- tools::file_ext(unname(inFile$datapath))
      validate(need(ext %in% c("db"), "Please select files with .db extension!"))
      VARS$dbb <- append(VARS$dbb,unname(inFile$datapath))
      VARS$dbb <- unique(VARS$dbb)
      rm(ext)
    }
  })
  
  ##--- if remote databases are allowed
  observe({
    if(USE_REMOTE_DB==TRUE){
      VARS$dbb = NULL
      cf <- GetStudyName(AppName=REPO_NAME)
      cf <- as.character(unlist(cf[,1]))
      cf <- cf[grep("^scRNA_",cf)]
      if(isTruthy(cf)){
        VARS$dbb <- unique(cf)
      }
      rm(cf)
    }else{
      print("Currently working with local databases\n")
    }
  })
  
  
  observeEvent({VARS$dbb},{
    if(isTruthy(VARS$dbb)){
      VARS$connList <- ({
        connList <- list()
        if(length(VARS$dbb)>0){
          for(f in VARS$dbb){
            id = pddExpn::randomStringGenerator(n = 1,lenght = 8)
            if(USE_REMOTE_DB==TRUE){
              connList[[id]] <- f #gsub("^scRNA_","",f)
            }else{
              connList[[id]] <- RSQLite::dbConnect(RSQLite::SQLite(),f)
            }
            rm(id,f)
          }
        }
        connList
      })
      VARS$sc_studies <- ({
        sc_studies <- c()
        for(i in 1:length(connList)){
          if(USE_REMOTE_DB==F){
            tablist <- RSQLite::dbListTables(connList[[i]])
            tablist <- tablist[grep("_study",tablist)]
          }else{
            tablist <- paste0(connList[[i]],"_study")
            #tablist <- gsub("^scRNA_","",tablist)
          }
          for(j in tablist){
            cat(i,"\t",j,"\n")
            query <- paste0("SELECT * FROM ",j)
            sc_studies <- rbind(sc_studies,
                                cbind(
                                  queryDB(HANDLER=connList[[i]], QUERY=query,
                                          REPO_NAME=REPO_NAME,
                                          USE_REMOTE_DB=USE_REMOTE_DB),
                                  ObjID=names(connList)[i]
                                  )
                                )
          }
          rm(i,j)
        }
        sc_studies
      })
    }
  })
  
  

  ##-----------------------------------------------------------------------
  ###------------ variables
  variables = reactiveValues(
    ## general
    gene = NULL,
    geneA = NULL,
    sc_study = NULL,
    comp_study = NULL,
    ## sc_study
    connID_1 = NULL,
    feature = NULL,
    add_features = NULL,
    cellTypes = NULL,
    add_feature_list = list(),
    sc_attribute = NULL,
    sc_louvain = NULL,
    sc_df = NULL,
    pldf = list(),
    cordf = NULL,
    corrdf_sel=NULL,
    marker_volc = NULL,
    marker_dot = c(),
    deg_volc = NULL,
    deg_dot = c(),
    pldf_marker = list(),
    pldf_deg = list(),
    ## enrichments
    markenrich_up = list(),
    markenrich_dn = list(),
    degenrich_up = list(),
    degenrich_dn = list(),
    markerTests = NULL,
    degTests = NULL,
    ## comparing study
    comp_feature = NULL,
    comp_add_features = NULL,
    comp_cellTypes = NULL,
    comp_add_feature_list = list(),
    comp_attribute = NULL,
    comp_louvain = NULL,
    comp_df = NULL
  )
  gc()
  Prog <- shiny::Progress$new()
  on.exit(Prog$close())
  Prog$set(message = "initializing the app..", value = 0)
  
  observe({ cat("tab item: ",input$AppTab,"\n") })
  Prog$set(message = "populating the dropdowns..", value = 0.3)
  
  ##---------- sc studies
  observeEvent({VARS$sc_studies},{
    if(isTruthy(VARS$sc_studies)){
      sc_studies <- VARS$sc_studies
      updateSelectizeInput(session, 'sc_study', choices = VARS$sc_studies$Database,
                           server = TRUE,
                           selected =VARS$sc_studies$Database[1])
    }
  })
  ##----- comp study
  observe({
    sc_studies <- VARS$sc_studies
    comp_studies <- sc_studies[sc_studies$Database!=input$sc_study,]
    updateSelectizeInput(session, 'comp_study', choices = comp_studies$Database, 
                         server = TRUE,selected =comp_studies$Database[1] )
  })
  
  
  ##----- List of available genes
  gGenes <- reactive({
    query <- paste0("SELECT * FROM gAliasProtein")
    gGenes <- RSQLite::dbGetQuery(VARS$connGenes, query)
  })
  #observe({
  #  availGenes <- as.character(unique(gGenes()$geneSymbol))
  #  #availGenes <- c('WNT5A','ACF','A1CF','POSTN','NPPA','NPPB','TTN','SFTPD')
  #  updateSelectizeInput(session, 'geneA', choices = availGenes, 
  #                       server = TRUE,selected ="NPPB")
  #})
  observeEvent({input$sc_study},{
    if(isTruthy(input$sc_study)){
      sc_studies <- VARS$sc_studies
      qq = sc_studies[sc_studies$Database==input$sc_study,]$ObjID %>% 
        as.character
      query <- paste0("SELECT geneSymbol FROM ",input$sc_study,"_data")
      cat(" getting gene lists for: ", input$sc_study, " objID: ", qq,"\n")
      availGenes <- queryDB(HANDLER=VARS$connList[[qq]], 
                            QUERY=query,REPO_NAME=REPO_NAME,
                            USE_REMOTE_DB=USE_REMOTE_DB)
      updateSelectizeInput(session, 'geneA', choices = availGenes[,1], 
                           server = TRUE,selected ="WNT5A")
    }
  })

  
  ###-------- alias genes
  gAlias <- reactive({
    gAlias <- list()
    for(f in input$geneA){
      hgnc <- gGenes()[gGenes()$geneSymbol==f,]$HGNC %>% as.character
      gg <- as.character(gGenes()[gGenes()$HGNC%in%hgnc,]$geneSymbol)
      #gg <- gg[gg!=f]
      gAlias[[f]] <- unique(c(f,gg))
    }
    gAlias
  })
  observe({
    cat("Alias gene names: ", paste0(unlist(gAlias()),collapse = ";"),"\n")
  })
  

  
  Prog$set(message = "done!", value = 1)
  
  source(file = "server_df.R",local = TRUE,encoding = "UTF-8")$value

  observeEvent({input$geneA
    input$sc_study
    input$comp_study},{
      cat("object size\t", object.size(variables),"\n")
    })
  
  source(file = "server_sc.R",local = TRUE,encoding = "UTF-8")$value
  
  source(file = "server_compare.R",local = TRUE,encoding = "UTF-8")$value

}