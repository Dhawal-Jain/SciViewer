## server_df.R

###----------------------------------------------------------
###---------------------- Study 1
###----------------------------------------------------------
  ## sc_louvain, assign variables based on study name
observeEvent({input$sc_study},{
  if(isTruthy(input$sc_study)  &
     (is.null(variables$sc_study) | isTruthy(input$sc_study != variables$sc_study))
     ){
    
    sc_studies = VARS$sc_studies 
    variables$sc_study = input$sc_study
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "resetting/initializing DF calculations..", value = 0)
    
    cat("###----------------------------------\n")
    cat(" Inside bracket-1\n setting the study values to null\n\n")
    
    variables$connID_1 = NULL
    variables$feature = NULL
    variables$add_features = NULL
    variables$cellTypes = NULL
    variables$add_feature_list = list()
    variables$sc_attribute = NULL
    variables$sc_louvain = NULL
    variables$sc_df = NULL
    variables$pldf = list()
    variables$pldf_marker = list()
    variables$pldf_deg = list()
    variables$cordf = NULL
    variables$corrdf_sel=NULL
    variables$marker_dot = c()
    variables$marker_volc = NULL
    variables$gene = NULL
    variables$geneA = NULL
    variables$markenrich_up = list()
    variables$markenrich_dn = list()
    variables$degenrich_up = list()
    variables$degenrich_dn = list()
    variables$markerTests = NULL
    variables$degTests = NULL
    gc()
    
    variables$connID_1 <- ({
      sc_studies[sc_studies$Database==input$sc_study,]$ObjID %>% 
        as.character
    })
    cat(" study1: ", input$sc_study,"\n")
    cat(" study1 connID: ", variables$connID_1,"\n")
    
    variables$sc_louvain <- ({
      query <- paste0("SELECT * FROM ",input$sc_study,"_metaFeatures")
      louvain <- queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
              QUERY=query,REPO_NAME=REPO_NAME,
              USE_REMOTE_DB=USE_REMOTE_DB)
      louvain$value <- 0
      louvain
    })
    scProg$set(message = "Calculated louvain DF..", value = 0.5)

    ####################################################
    variables$sc_attribute <- 
      paste("<font color=\"#FF0000\"><b> STUDY STATUS: </b> ",sc_studies[sc_studies$Database==input$sc_study,]$STATUS,"! </font><br>",
            "<b> # of cells: </b>",format(nrow(variables$sc_louvain) ,big.mark = ",",scientific = F),"<br>",
            "<b> TISSUE: </b>", sc_studies[sc_studies$Database==input$sc_study,]$TISSUES," ",
            "<b> SAMPLE SIZE: </b>", sc_studies[sc_studies$Database==input$sc_study,]$SampleSize,"<br>",
            "<b> PUBMED: </b>", a(sc_studies[sc_studies$Database==input$sc_study,]$PMID,href=paste0("https://pubmed.ncbi.nlm.nih.gov/",sc_studies[sc_studies$Database==input$sc_study,]$PMID), target="_blank"),
            "<b> ACCSSION: </b>", a("COVID19",href=paste0(sc_studies[sc_studies$Database==input$sc_study,]$GEO), target="_blank"),"<br>",
            "<b> STUDY ABSTRACT: </b> <font color=\"#bdbdbd\">", sc_studies[sc_studies$Database==input$sc_study,]$Description,"</font><br>"
      )
    
    x <- unlist(strsplit(as.character(sc_studies[sc_studies$Database==input$sc_study,]$CATVAR),","))
    if(!is.null(x)){
      variables$feature=x[1]
    }else{
      variables$feature='cell_type'
    }
    if(length(x)>1){
      variables$add_feature=x[2:length(x)]
    }
    
    ## messages
    cat(" Input geneA ",input$geneA,"\n")
    cat(" Variables geneA ",input$geneA,"\n")
    cat(" variables gene ",input$gene,"\n")
    cat(" Input feature ",variables$feature,"\n")
    cat(" Input additional features ",variables$add_feature,"\n")
    cat(" Input scStudy_2 (for comparison) ",input$comp_study,"\n")
    cat("\n")
    scProg$set(message = "done..", value = 1)
  }
})

  ## generate data frame for genes
observeEvent({input$geneA
  variables$sc_louvain},{

  if(isTruthy(input$sc_study) & isTruthy(input$geneA) & 
     isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
    
    if(is.null(variables$gene) | isTruthy(variables$gene!=input$geneA[1])
    ){
      variables$sc_df = NULL
      variables$gene=input$geneA[1]
      variables$marker_dot = c()
      variables$deg_dot = c()
      
      cat(" Inside bracket-2\nvariables$gene=", variables$gene,
          "\nconnID: ", variables$connID_1,"\n")
      
      z <- paste0(input$sc_study,'_',variables$gene)
      ##-- use all gene Aliases to search the database 
      cntr = 1
      cat(" xxx: ",length(gAlias()[[variables$gene]])," ", isTRUE(nrow(variables$sc_df))," ",gAlias()[[variables$gene]][cntr],"\n")
      
      while(cntr <= length(gAlias()[[variables$gene]]) &&
            !isTRUE(nrow(variables$sc_df)>0)){
        cat("  gene/alias name: ", gAlias()[[variables$gene]][cntr],"\n")
        variables$sc_df <- variables$pldf[[z]] <- 
                   get_plot_df_sc04(connSc = VARS$connList[[variables$connID_1]],
                                    study = input$sc_study,
                                    genename = gAlias()[[variables$gene]][cntr],
                                    louvain = variables$sc_louvain,
                                    REPO_NAME=REPO_NAME,
                                    USE_REMOTE_DB=USE_REMOTE_DB)
        cntr = cntr+1
      }
      #if(is.null(variables$sc_df)){
      #  variables$sc_df <- NULL
      #  variables$pldf[[z]] <- NULL
      #}
      cat(" names(variables$pldf):",names(variables$pldf),"\n")
      
      ## populate cell type selection bar
      variables$cellTypes <- unique(variables$sc_df[,variables$feature])
      
      if(!is.null(variables$add_feature)){
        for(f in variables$add_feature){
          if(f %in% names(variables$sc_df)){
            variables$add_feature_list[[f]] <- unique(variables$sc_df[,f])
          }else{
            variables$add_feature <- variables$add_feature[variables$add_feature!=f]
          }
        }
      }
      #####
      cat("calculated first-gene DF for plot\n\n")
      #scProg$set(message = "done..", value = 1)
    }
    
    if(is.null(variables$geneA) | 
        sum(variables$geneA!=input$geneA)>0
    ){
      
      variables$geneA=input$geneA  
      variables$marker_dot = c()
      variables$deg_dot = c()
      
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = "calculating DF for multiple genes..", value = 0)
      
      cat(" Inside bracket-3\ncalculating for multiple genes..\n")
      ##------------ generate and update pldf
      for(g in input$geneA){
        y <- paste0(input$sc_study,'_',g)
        if(!y %in% names(variables$pldf)){
          pl <- NULL
          cntr =1 
          while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
            cat("  gene/alias name: ", gAlias()[[g]][cntr],"\n")
            pl <-  get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                                            input$sc_study,
                                            gAlias()[[g]][cntr],
                                            variables$sc_louvain,
                                            REPO_NAME=REPO_NAME,
                                            USE_REMOTE_DB=USE_REMOTE_DB)
            cntr = cntr +1 
          }
          if(isTRUE(nrow(pl)>0)){
            cat('  getting plot df for gene ',g, ":",y,":",names(variables$pldf),"\n")
            variables$pldf[[y]] <- pl
          }
        }else{
          cat('  plot df for gene ',g, " is already calculated\n")
        }
      }
      rm(y,cntr,g)
      z <- paste0(input$sc_study,'_',input$geneA)
      for(y in names(variables$pldf)){
        if(!y%in%z){
          cat('  removing plot df for gene ',y, "\n")
          if(y %in%names(variables$pldf)){
            variables$pldf[[y]] <- NULL
          }
        }
      }
      
      ###------- multigene corr 
      variables$corrdf <- NULL
      variables$corrdf_sel <- NULL
      if(length(variables$pldf)>1){
        scProg$set(message = "Calculating multigene correlation matrix..", value = 0.8)
        cat('  Calculating multigene correlation matrix. Available cell group: ', variables$feature,"\n")
        cat(" ", names(variables$pldf),"\n")
        qw <- data.frame(cell_type=variables$sc_df[,'cell_type'])
        if(nrow(qw)>0){
          for(i in 1:length(variables$pldf)){
            if(nrow(variables$pldf[[i]])>0){
              cat(" corr matrix: i=",i," gene=",names(variables$pldf)[i],"\n")
              qw <- cbind(qw,variables$pldf[[i]][,'value'])
            }
          }
          gnames <- gsub(paste0(input$sc_study,"_"),"",names(variables$pldf))
          names(qw)[2:ncol(qw)] <- gnames
          variables$corrdf <- qw
          variables$corrdf_sel <- qw[,2:ncol(qw)]
        }
      }
      cat("\n")
      
      ###---------- get marker results
      scProg$set(message = "Retrieving Markers..", value = 0.85)
      for(g in input$geneA){
        pl <- NULL
        cntr <-1
        while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
          cat("  Marker gene/alias name: ", gAlias()[[g]][cntr],"\n")
          query <- paste0("SELECT * FROM ",input$sc_study,"_Marker WHERE geneSymbol = '",gAlias()[[g]][cntr],"'")
          pl <- tryCatch({
            queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                    QUERY=query,REPO_NAME=REPO_NAME,
                    USE_REMOTE_DB=USE_REMOTE_DB)
            
          },error= function(e){
            return(NULL)
          })
          cntr = cntr +1 
        }
        if(isTRUE(nrow(pl)>0)){
          cat( " Marker rows for dotplot: ",nrow(pl),"\n")
          variables$marker_dot <- rbind(variables$marker_dot,pl)
        }
      }
      rm(cntr,g,pl)
          
      ###---------- get deg results
      scProg$set(message = "Retrieving DEG", value = 0.9)
      for(g in input$geneA){
        pl <- NULL
        cntr <-1
        while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
          cat("  Marker gene/alias name: ", gAlias()[[g]][cntr],"\n")
          query <- paste0("SELECT * FROM ",input$sc_study,"_DEG WHERE geneSymbol = '",gAlias()[[g]][cntr],"'")
          pl <- tryCatch({
            queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                    QUERY=query,REPO_NAME=REPO_NAME,
                    USE_REMOTE_DB=USE_REMOTE_DB)
            
          },error= function(e){
            return(NULL)
          })
          cntr = cntr +1 
        }
        if(isTRUE(nrow(pl)>0)){
          cat( " DEG rows for dotplot: ",nrow(pl),"\n")
          variables$deg_dot <- rbind(variables$deg_dot,pl)
        }
      }
      rm(cntr,g,pl)
      
      ###---------- end
      scProg$set(message = "done..", value = 1)
      
    }
  }
})

 ## select cell types for correlation matrices
observeEvent({input$sc_sel_3a},{
  if(isTruthy(input$sc_sel_3a) & isTruthy(variables$corrdf)){
    cat("Inside bracket-4\n\n")
    qw <- variables$corrdf
    qw <- qw[qw$cell_type%in%input$sc_sel_3a,]
    cat("sc_sel_3a ", input$sc_sel_3a,"  nrow:",nrow(qw),"\n")
    cat(paste0(qw[1,],collapse = ";"),"\n")
    variables$corrdf_sel <- qw[,2:ncol(qw)]
  }else if(isTruthy(variables$corrdf)){
    qw <- variables$corrdf
    cat("NULLVal: sc_sel_3a ", input$sc_sel_3a,"  nrow:",nrow(qw),"\n")
    variables$corrdf_sel <- qw[,2:ncol(qw)]
  }
},ignoreNULL = F)


###----------------------------------------------------------
###---------------------- Enrichment analyses for Markers/DEG
###----------------------------------------------------------
  ## select cell types for displaying Marker volcano plots
observeEvent({input$sc_celltype4a},{
  if(isTruthy(input$sc_celltype4a)){
    variables$pldf_marker = list()
    variables$markerTests = NULL
    query <- paste0("SELECT * FROM ",input$sc_study,"_Marker WHERE Tag = '",input$sc_celltype4a,"'")
    pl <- tryCatch({
      queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
              QUERY=query,REPO_NAME=REPO_NAME,
              USE_REMOTE_DB=USE_REMOTE_DB)
    },error=function(e){
      return(NULL)
    })
    if(isTRUE(nrow(pl)>0)){
      cat("  getting the Marker data frame for volcano plot\t", input$sc_celltype4a,"\t",nrow(pl),"\n")
      variables$markerTests = unique(pl$Test)
      variables$marker_volc <- pl
    }
  }
})

 ## select cell types for displaying DEG volcano plots
observeEvent({input$sc_celltype5a},{
  if(isTruthy(input$sc_celltype5a)){
    variables$pldf_deg = list()
    variables$degTests = NULL
    query <- paste0("SELECT * FROM ",input$sc_study,"_DEG WHERE Tag = '",input$sc_celltype5a,"'")
    pl <- tryCatch({
      queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
              QUERY=query,REPO_NAME=REPO_NAME,
              USE_REMOTE_DB=USE_REMOTE_DB)
    },error=function(e){
      return(NULL)
    })
    if(isTRUE(nrow(pl)>0)){
      cat("  getting the DEG data frame for volcano plot\t", input$sc_celltype5a,"\t",nrow(pl),"\n")
      variables$degTests = unique(pl$Test)
      variables$deg_volc <- pl
    }
  }
})


###----------------------------------------------------------
###---------------------- Option selection in the app
###----------------------------------------------------------
  ## Update input choices across the app
observe({
  sc_studies = VARS$sc_studies
  cat(" Inside bracket-5\n")
  cat(" SelectInputbar population_1..\n")
  cat("\n")
  updateSelectizeInput(session, 'sc_celltype1a', 
                       choices = variables$cellTypes, 
                       server = TRUE,selected =NULL)
  updateSelectizeInput(session, 'sc_celltype2a', 
                       choices = variables$cellTypes, 
                       server = TRUE,selected =NULL)
  updateSelectizeInput(session, 'sc_celltype3a', 
                       choices = variables$cellTypes, 
                       server = TRUE,selected =NULL)
  updateSelectizeInput(session, 'sc_celltype4a', 
                       choices = variables$cellTypes, 
                       server = TRUE,selected =NULL)
  updateSelectizeInput(session, 'sc_celltype5a', 
                       choices = variables$cellTypes, 
                       server = TRUE,selected =NULL)
  updateSelectizeInput(session, 'sc_celltype6a', 
                       choices = variables$cellTypes, 
                       server = TRUE,selected =NULL)
  updateSelectizeInput(session, 'sc_sel_3a', 
                       choices = variables$cellTypes, 
                       server = TRUE,selected =NULL)
  ##
  feat <- unlist(strsplit(as.character(sc_studies[sc_studies$Database==input$sc_study,]$CATVAR),","))
  updateSelectizeInput(session, 'sc_sel_2a',  
                       choices = feat, 
                       server = TRUE,selected =feat[1])
  #updateSelectizeInput(session, 'sc_sel_1n',  
  #                     choices = feat, 
  #                     server = TRUE,selected =feat[1])

  ## additional pca display
  if(!is.null(variables$add_feature)){
    updateSelectizeInput(session, 'sc_sel_1',  
                         choices = variables$add_feature, 
                         server = TRUE,selected =variables$add_feature[1])
  }
})
  
## select 1a
observe({
  cat(" Inside bracket-6\n")
  cat(" SelectInputbar population_2..\n")
  cat("\n")
  updateSelectizeInput(session, 'sc_sel_1a', 
                       choices = variables$add_feature_list[[input$sc_sel_1]], 
                       server = TRUE,selected =NULL)
})
### expression color scale choices
observe({
  cat("Inside bracket-7\n\n")
  updateSelectizeInput(session, 'c1t2expnscale', 
                       choices = c("GrayBlue","Clay","magma","Viridis",
                                   "GreenYellow","Purple","Reds",
                                   "Oranges"), 
                       server = TRUE,selected =NULL)
})
observe({
  cat("Inside bracket-7a\n\n")
  updateSelectizeInput(session, 'sc_sel_4a', 
                       choices = variables$markerTests, 
                       server = TRUE,selected =NULL)
})
observe({
  cat("Inside bracket-7b\n\n")
  updateSelectizeInput(session, 'sc_sel_5a', 
                       choices = variables$degTests, 
                       server = TRUE,selected =NULL)
})

###----------------------------------------------------------
###---------------------- Multigene display for celltype marker genes
###----------------------------------------------------------
observeEvent({input$marker_tab_rows_selected
  variables$sc_louvain},{
    
    if(isTruthy(input$marker_tab_rows_selected) & 
       isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
      
      ID = variables$marker_volc[input$marker_tab_rows_selected,]$geneSymbol %>%
        as.character()
      
      cat(" Inside bracket-6Marker\ncalculating for multiple genes..\n")
      ##------------ generate and update pldf
      for(g in ID){
        y <- paste0(input$sc_study,'_',g)
        if(!y %in% names(variables$pldf_marker)){
          pl <- NULL
          pl <- get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                                          input$sc_study,
                                          g,
                                          variables$sc_louvain,
                                          REPO_NAME=REPO_NAME,
                                          USE_REMOTE_DB=USE_REMOTE_DB)
          if(isTRUE(nrow(pl)>0)){
            cat('  getting plot df for gene ',g, ":",y,":",names(variables$pldf_marker),"\n")
            variables$pldf_marker[[y]] <- pl
          }
        }else{
          cat('  plot df for gene ',g, " is already calculated\n")
        }
      }
      rm(y,g)
      z <- paste0(input$sc_study,'_',ID)
      for(y in names(variables$pldf_marker)){
        if(!y%in%z){
          cat('  removing plot df for gene ',y, "\n")
          if(y %in%names(variables$pldf_marker)){
            variables$pldf_marker[[y]] <- NULL
          }
        }
      }
      
    }
})

###----------------------------------------------------------
###---------------------- Multigene display for celltype DEG genes
###----------------------------------------------------------
observeEvent({input$deg_tab_rows_selected
  variables$sc_louvain},{
    
    if(isTruthy(input$deg_tab_rows_selected) & 
       isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
      
      ID = variables$deg_volc[input$deg_tab_rows_selected,]$geneSymbol %>%
        as.character()
      
      cat(" Inside bracket-6DEG\ncalculating for multiple genes..\n")
      ##------------ generate and update pldf
      for(g in ID){
        y <- paste0(input$sc_study,'_',g)
        if(!y %in% names(variables$pldf_deg)){
          pl <- NULL
          pl <- get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                                 input$sc_study,
                                 g,
                                 variables$sc_louvain,
                                 REPO_NAME=REPO_NAME,
                                 USE_REMOTE_DB=USE_REMOTE_DB)
          if(isTRUE(nrow(pl)>0)){
            cat('  getting plot df for gene ',g, ":",y,":",names(variables$pldf_deg),"\n")
            variables$pldf_deg[[y]] <- pl
          }
        }else{
          cat('  plot df for gene ',g, " is already calculated\n")
        }
      }
      rm(y,g)
      z <- paste0(input$sc_study,'_',ID)
      for(y in names(variables$pldf_deg)){
        if(!y%in%z){
          cat('  removing plot df for gene ',y, "\n")
          if(y %in%names(variables$pldf_deg)){
            variables$pldf_deg[[y]] <- NULL
          }
        }
      }
      
    }
  })

###----------------------------------------------------------
###---------------------- Study 2
###----------------------------------------------------------

observeEvent({input$comp_study
  input$t2action1},{
    if(isTruthy(input$comp_study) & isTruthy(variables$gene) &
       input$AppTab=='compare'
    ){
      
      sc_studies = VARS$sc_studies
      cat("setting the comp-study values to null\n")
      variables$comp_study = input$comp_study
      variables$comp_feature = NULL
      variables$comp_add_features = NULL
      variables$comp_cellTypes = NULL
      variables$comp_add_feature_list = list()
      variables$comp_attribute = NULL
      variables$comp_louvain = NULL
      variables$comp_df = NULL
      
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = "initializing DF generation for comparing study..", value = 0)

      cat("Inside bracket-8\ncalculating df for comparison study..\n")

      connID2 <- reactive({
        sc_studies[sc_studies$Database==input$comp_study,]$ObjID %>%
          as.character
      })
      
      variables$comp_louvain <- local({
        query <- paste0("SELECT * FROM ",input$comp_study,"_metaFeatures")
        louvain <- queryDB(HANDLER=VARS$connList[[connID2()]], 
                QUERY=query,REPO_NAME=REPO_NAME,
                USE_REMOTE_DB=USE_REMOTE_DB)
        louvain$value <- 0
        louvain
      })
      scProg$set(message = "Calculated louvain DF..", value = 0.5)
      
      variables$comp_attribute <- 
        paste("<font color=\"#FF0000\"><b> STUDY STATUS: </b> ",sc_studies[sc_studies$Database==input$comp_study,]$STATUS,"! </font><br>",
              "<b> # of cells: </b>",format(nrow(variables$comp_louvain) ,big.mark = ",",scientific = F),"<br>",
              "<b> TISSUE: </b>", sc_studies[sc_studies$Database==input$comp_study,]$TISSUES," ",
              "<b> SAMPLE SIZE: </b>", sc_studies[sc_studies$Database==input$comp_study,]$SampleSize,"<br>",
              "<b> PUBMED: </b>", a(sc_studies[sc_studies$Database==input$comp_study,]$PMID,href=paste0("https://pubmed.ncbi.nlm.nih.gov/",sc_studies[sc_studies$Database==input$comp_study,]$PMID), target="_blank"),
              "<b> ACCSSION: </b>", a("COVID19",href=paste0(sc_studies[sc_studies$Database==input$comp_study,]$GEO), target="_blank"),"<br>",
              "<b> STUDY ABSTRACT: </b> <font color=\"#bdbdbd\">", sc_studies[sc_studies$Database==input$comp_study,]$Description,"</font><br>"
        )

            ## feature column, this part should go! update the database column name  
      x <- unlist(strsplit(as.character(sc_studies[sc_studies$Database==input$comp_study,]$CATVAR),","))
      if(!is.null(x)){
        variables$comp_feature=x[1]
      }else{
        variables$comp_feature='cell_type'
      }
      
      #####-----------------------------------------------
      ##-- use all gene Aliases to search the database 
      cntr = 1
      while(cntr <= length(gAlias()[[variables$gene]]) &&
            !isTRUE(nrow(variables$comp_df)>0)){
        cat("  compStudy gene/alias name: ", gAlias()[[variables$gene]][cntr],"\n")
        variables$comp_df <- 
          get_plot_df_sc04(connSc = VARS$connList[[connID2()]],
                           study = input$comp_study,
                           genename = gAlias()[[variables$gene]][cntr],
                           louvain = variables$comp_louvain,
                           REPO_NAME=REPO_NAME,
                           USE_REMOTE_DB=USE_REMOTE_DB)
        cntr = cntr+1
      }
      rm(cnt)
      scProg$set(message = "Calculated plot DF..", value = 0.75)
      
      ### cell type in study 1
      updateSelectizeInput(session, 'sc_celltype21', 
                           choices = variables$cellTypes, 
                           server = TRUE,selected =NULL)
      ### cell type in study 2
      cat( 'study b feature: ', variables$comp_feature,'\n')
      b_celltypes <- unique(variables$comp_df[,variables$comp_feature])
      updateSelectizeInput(session, 'sc_celltype22', 
                           choices = b_celltypes, 
                           server = TRUE,selected =NULL)
      cat("\n")
      scProg$set(message = "done..", value = 1)
      
    }
})



