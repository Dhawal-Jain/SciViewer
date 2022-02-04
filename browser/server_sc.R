# server_sc.R

observe({
  cat( "ssc: printing study-1 attribute\n")
  output$sc_study_attribs1 <- output$sc_study_attribs <- renderText({ 
    variables$sc_attribute
  })
})

##----------------------------------------------------------------------
## tab page 1
##----------------------------------------------------------------------

observeEvent({variables$sc_df},{
    if(isTruthy(variables$sc_df)){
      
      cat("ssc: tab1 displays\n")

      ### legend
      output$scclust_legend <- renderUI({ 
        suppressWarnings({
          scclust_legend <- plot_umap_raster(pl=variables$sc_df,
                             feature = variables$feature,
                             genename = "Cell type clusters",
                             legendPlot = T,
                             selectedFeature = input$sc_celltype1a)
          
          output$dfgertawef <- renderPlotly(scclust_legend %>%
                                              layout(height = 400, width = 300)
                                           )
          plotlyOutput("dfgertawef")
        })
      })
  
      ## tab1 PCA plot 
      sc_pcaplot1a <- reactive({
        plot_umap_raster(pl=variables$sc_df,
                         feature = variables$feature,
                         genename = "Cell type clusters",
                         legendPlot = F,
                         selectedFeature = input$sc_celltype1a,
                         source.val='sc_pca1',js=js)
      })
      output$scpca1a_hover <- renderText({
        getHoverText(variables$sc_df, 
                     colids=c('V1','V2',variables$feature), 
                     hovercoord = input$sc_pca1 )
      })
      output$sc_pcaplot1a <- renderPlotly({ 
        suppressWarnings({
          sc_pcaplot1a() %>% layout(height = 450, width = 450)
        })
      })
      
      ## tab1 barplot 
      output$sc_barplot <- renderHighchart({
        dhc_columnPlot(pl=variables$sc_df,
                       features = variables$feature,
                       ycol = 'value',plotType = 'count',
                       log.transform = F,main = "Cell type proportions",
                       xlab = "",ylab = "% of total",
                       sourceref=SOURCEREF)
      })
      
      ## summary table
      output$sc_summaryTable = DT::renderDataTable({
        cf <-  cellCountSummary(pl = variables$sc_df,
                                features = variables$feature,
                                ycol = 'value')
        
        datatable(cf,rownames = F,extensions = 'Buttons',
                  options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                                 buttons = list(list(extend = 'csv',  filename = paste0("CellCountSummary_",input$sc_study,"_",input$sc_study,"")),
                                                list(extend ='excel', filename = paste0("CellCountSummary_",input$sc_study,"_",input$sc_study,""))), 
                                 scrollX = TRUE,scrollCollapse = TRUE
                  )) %>% formatStyle('Number of cells',
                                     background = styleColorBar(range(cf$`Number of cells`,na.rm=T), 'lightblue'),
                                     backgroundSize = '98% 88%',
                                     backgroundRepeat = 'no-repeat',
                                     backgroundPosition = 'center')
      })  
      #####  
    }
})

## tab1 selected pca
observeEvent({input$sc_sel_1},{
    if(input$AppTab=='scstudy' &
       isTruthy(input$sc_sel_1) & 
       isTruthy(variables$sc_df)){
      
      cat('ssc: tab1 hover..\n')
      cat('ssc: input$sc_sel_1 ',input$sc_sel_1,"\n")
      cat('ssc: input$sc_sel_1a ',input$sc_sel_1a,"\n")
      sc_selplot1a <- reactive({
        plot_umap_raster(pl=variables$sc_df,
                         x='V1',y='V2',feature = input$sc_sel_1,
                         genename = paste0(input$sc_sel_1, " clusters"),
                         legendPlot = F,
                         selectedFeature = input$sc_sel_1a,
                         source.val='scsel1',js=js)
      })
      output$scsel1a_hover <- renderText({
        getHoverText(variables$sc_df, 
                     colids=c('V1','V2',input$sc_sel_1), 
                     hovercoord = input$scsel1 )
        
      })
      output$sc_selplot1a <- renderPlotly({ 
        suppressWarnings({
          sc_selplot1a() %>% layout(height = 450, width = 450)
        })
      })
    }
})


##----------------------------------------------------------------------
## tab page 2
##----------------------------------------------------------------------

observeEvent({variables$sc_df
  input$c1t2expnscale},{
    if(input$AppTab=='scstudy' &
       isTruthy(variables$sc_df)){
      
      cat("ssc: variables$gene: ", variables$gene,"\n")
      ## pca plot2a    
      sc_pcaplot2a <- reactive({
        plot_umap_raster(pl=variables$sc_df,
                         feature = variables$feature,
                         genename = "Cell type clusters",
                         legendPlot = F,
                         selectedFeature = input$sc_celltype2a,
                         source.val='sc_pca2a',js=js)
      })
      range2a <- reactive({
        zoom <- event_data("plotly_relayout", "sc_pca2a")
        getZoomrange(zoom)
      })
      cat("range: ",range2a()[[1]],"\n")
      output$scpca2a_hover <- renderText({
        getHoverText(variables$sc_df, 
                     colids=c('V1','V2',variables$feature), 
                     hovercoord = input$sc_pca2a )
        
      })
      output$sc_pcaplot2a <- renderPlotly({ 
        suppressWarnings({
          sc_pcaplot2a() %>%
            layout(height = 450, width = 450)
        })
      })
      
      ## pca2b
      output$sc_pcaplot2b <- renderPlotly({ 
        ql <- select_colScale(input$c1t2expnscale)
        
        suppressWarnings({
          plot_multi_umap_raster(list(variables$sc_df),
                                 genenames = variables$gene,
                                 legendPlot = F,xrange =range2a()[[1]],
                                 yrange=range2a()[[2]],
                                 colScale=ql[[1]],
                                 minValCol=ql[[2]]) %>%
            layout(height = 450, width = 450)
        })
      })
      #########    
    }
  })

observeEvent({input$t2radio
  input$sc_sel_1},{
    if(input$AppTab=='scstudy' & 
       isTruthy(input$sc_sel_2a) & isTruthy(variables$sc_df)){
      
      cat("ssc: input$sc_sel_2a: ", input$sc_sel_2a, "\n")
      output$sc_barplot_2a <- renderHighchart({
        if(input$t2radio==1){
          dhc_columnPlot(pl=variables$sc_df,
                         features = input$sc_sel_2a,
                         plotType = 'expn',min.cell.number = 10,
                         ycol = 'value',
                         log.transform = F,main = paste0(input$geneA[1]," expression"),
                         xlab = "",ylab = "Normalized expression",
                         sourceref=SOURCEREF)
        }else if(input$t2radio==2){
          dhc_boxPlot(pl=variables$sc_df,
                      features = input$sc_sel_2a,
                      ycol = 'value',
                      log.transform = F,
                      main = paste0(input$geneA[1]," expression"),
                      xlab = "",ylab = "Normalized expression",
                      sourceref=SOURCEREF)
        }
      })
    }
    
  })


##----------------------------------------------------------------------
## tab page 3
##----------------------------------------------------------------------

observeEvent({variables$pldf
  variables$sc_df},{
    if(input$AppTab=='scstudy' &
       isTruthy(input$sc_study)& 
       isTruthy(variables$pldf) & 
       isTruthy(variables$sc_df) ){
      
      cat("ssc: tab3/4 \n")
      sc_pcaplot3a <- reactive({
        plot_umap_raster(pl=variables$sc_df,
                         feature = variables$feature,
                         genename = "Cell type clusters",
                         legendPlot = F,
                         selectedFeature = input$sc_celltype3a,
                         source.val='scpca2',js=js)
      })
      range1 <- reactive({
        zoom <- event_data("plotly_relayout", "scpca2")
        getZoomrange(zoom)
      })
      output$sc_pcaplot3a <- renderPlotly({ 
        suppressWarnings({
          sc_pcaplot3a() %>% layout(height = 400, width = 400)
        })
      })
      output$scpca3a_hover <- renderText({
        getHoverText(variables$sc_df, 
                     colids=c('V1','V2',variables$feature), 
                     hovercoord = input$scpca2 )
      })
      
      ## multigene expression plot 
      sc_expnlenged <- reactive({
        plot_multi_umap_raster(variables$pldf,
                               genenames = input$geneA,
                               legendPlot = T)
      })
      output$sc_pcaplot3b <- renderUI({ 
        gnames <- gsub(paste0(input$sc_study,"_"),"",names(variables$pldf))
        sc_pcaplot3b <- plot_multi_umap_raster(variables$pldf,
                               genenames = gnames,
                               legendPlot = F,xrange =range1()[[1]],
                               yrange=range1()[[2]])
        output$tmpa <- renderPlotly(sc_pcaplot3b)
        plotlyOutput("tmpa")
      })
      
      ## interactive dotplot/heatmap 
      sc_heatplot1 <- reactive({
        gnames <- gsub(paste0(input$sc_study,"_"),"",names(variables$pldf))
        plot_multigene_grouped_heatmap(variables$pldf,
                                       genenames = gnames,
                                       feature = variables$feature,
                                       x = 'V1',y = 'V2',ycol = 'value',
                                       log.transform = F,
                                       xrange = range1()[[1]],
                                       yrange = range1()[[2]])
      })
      output$sc_dotplot3a <- renderUI({ 
        rows <- 800
        columns <- 600 
        rows <- paste0(rows,"px")
        columns <- paste0(columns,'px')
        output$tempc <- renderPlotly(sc_heatplot1()[[2]])
        plotlyOutput("tempc",width = columns,height = rows)
      })

      
    }
})

  ####### correlation plot

observeEvent({variables$corrdf_sel
  input$dfgacvf
  variables$sc_df},{
    if(input$AppTab=='scstudy' &
       isTruthy(input$sc_study) & 
       isTruthy(variables$sc_df)){
      
      if(isTruthy(variables$corrdf_sel)){
        z <- variables$corrdf_sel
        z <- as.data.frame(z)
        z <- as.matrix(z)
        z[!is.finite(z)] <- 0
        z <- as.data.frame(z)
        if(input$dfgacvf=='all'){
        }else{
          flt <- apply(z, 1, function(x){
            ifelse(sum(x==0)==length(x),F,T)
          })
          z <- z[flt,]
          rm(flt)
        }
        cat("ssc: corr matrix dataframe dim: ", nrow(cor(z)),"\n")
        cat("  -->", cor(z,use = 'pairwise.complete.obs'),"\n")
        hc <- hchart_cor(cor(z,use = 'pairwise.complete.obs'),sourceref = SOURCEREF)
        rows <- max(c(300,ncol(z) *100))
        rows <- paste0(rows,"px")
        output$sc_corr_3a <- renderUI({ 
          output$tempd <- renderHighchart(hc)
          highchartOutput("tempd",width = rows,height = rows)
        })
      }else{
        output$sc_corr_3a <- renderUI({})
      }
    }
  })


##----------------------------------------------------------------------
## tab page 4
##----------------------------------------------------------------------
observeEvent({input$AppTab
  variables$sc_df},{
    if(input$AppTab=='scstudy' &
       isTruthy(input$sc_study) &
       isTruthy(variables$sc_df)){
      
      cat("ssc: tab6\n")
      sc_pcaplot4a <- reactive({
        plot_umap_raster(pl=variables$sc_df,
                         feature = variables$feature,
                         genename = "Cell type clusters",
                         legendPlot = F,
                         selectedFeature = input$sc_celltype4a,
                         source.val='scpca4',js=js)
      })
      output$sc_pcaplot4a <- renderPlotly({ 
        suppressWarnings({
          sc_pcaplot4a() %>% layout(height = 450, width = 450)
        })
      })
      output$scpca4a_hover <- renderText({
        getHoverText(variables$sc_df, 
                     colids=c('V1','V2',variables$feature), 
                     hovercoord = input$scpca4 )
      })
    }
})

  ##-- dotplot
observeEvent({input$AppTab
  variables$marker_dot
  input$sc_sel_4a},{
    if(input$AppTab=='scstudy' & isTruthy(input$sc_sel_4a) &
       isTRUE(nrow(variables$marker_dot)>0)){
      
      cat(" ssc: preparing the marker dotplot\t", input$sc_sel_4a,"\t",nrow(variables$marker_dot),"\n")
      output$sc_markerdotplot <- renderUI({ 
        output$tempc <- renderPlotly({
          qq <- variables$marker_dot
          qq <- qq[qq$Test==input$sc_sel_4a,]
          sc_dge_dotplot(qq,
                         x='geneSymbol',y='Tag',
                         colCol='logFC',sizeCol='adj.P.Val',
                         testCol='Test',configlayout=T)
        })
        plotlyOutput("tempc")
      })
    }else{
      output$sc_markerdotplot <- renderUI({})
    }
})

  ##--- spreadsheet and enrichment plots
observeEvent({input$AppTab
  input$slider_4a
  input$sc_sel_4a
  variables$marker_volc},{
    
    if(input$AppTab=='scstudy' & isTRUE(nrow(variables$marker_volc)>0) &
       isTruthy(input$sc_sel_4a) & isTruthy(input$slider_4a)){
      
      variables$markenrich_up = list()
      variables$markenrich_dn = list()
      cat(" ssc: printing marker data table\n")
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = "ssc: printing marker data table..", value = 0)
      ## spreadsheet
      output$marker_tab = DT::renderDataTable({
        cf <-  variables$marker_volc
        cf <- cf[cf$Test==input$sc_sel_4a,]
        cf <- cf[as.numeric(cf$adj.P.Val) <= as.numeric(input$slider_4a),]
        cf$logFC = round(cf$logFC,2)
        cf$AveExpr = round(cf$AveExpr,2)
        cf$t = round(cf$t,2)
        cf$B = round(cf$B,2)
        datatable(cf,rownames = F,extensions = 'Buttons',
                  options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                                 buttons = list(list(extend = 'csv',  filename = paste0("Markers_FDR_",input$slider_4a,"_",input$sc_study,"_",input$sc_celltype4a,"")),
                                                list(extend ='excel', filename = paste0("Markers_FDR_",input$slider_4a,"_",input$sc_study,"_",input$sc_celltype4a,""))), 
                                 scrollX = TRUE,scrollCollapse = TRUE
                  )) %>% formatSignif(columns = c('adj.P.Val','P.Value'),digits = 2)
      })
      
      ## enrichment
      scProg$set(message = "calculating enrichment for upregulated markers..", value = 0.2)
      variables$markenrich_up <- ({
        qf <-  variables$marker_volc
        qf <- qf[qf$Test==input$sc_sel_4a,]
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$slider_4a) & qf$logFC>0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Upregulated (wrt rest cells)'
          gprofiler2::gost(query = qf[1],organism = 'hsapiens',
                                 ordered_query = F,exclude_iea = T,
                                 significant = F,
                                 user_threshold = 0.05,correction_method = 'g_SCS')
        }else{ NULL }
      })
      scProg$set(message = "calculating enrichment for downregulated markers..", value = 0.6)
      variables$markenrich_dn <- ({
        qf <-  variables$marker_volc
        qf <- qf[qf$Test==input$sc_sel_4a,]
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$slider_4a) & qf$logFC < 0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Downregulated (wrt rest cells)'
          gprofiler2::gost(query = qf[1],organism = 'hsapiens',
                           ordered_query = F,exclude_iea = T,
                           significant = F,
                           user_threshold = 0.05,correction_method = 'g_SCS')
        }else{ NULL }
      })
      scProg$set(message = "done..", value = 1)
      
    }else{
      output$marker_tab = DT::renderDataTable({})
      output$enrichmarker_tab = DT::renderDataTable({})
      output$enrich_marker4a  = renderPlotly({})
    }
  })

observeEvent({input$AppTab
  variables$markenrich_up
  variables$markenrich_dn
  input$t4aradio},{
    if(input$AppTab=='scstudy' & isTruthy(input$t4aradio)){
      if(input$t4aradio==1){
        rr1 <- reactive({ variables$markenrich_up  })
      }else{
        rr1 <- reactive({ variables$markenrich_dn  })
      }
      
      output$enrich_marker4a <- renderUI({ 
        if(isTruthy(rr1())){
          output$ewascjy1 <- renderPlotly({
            gostplot(gostres = rr1(),capped = T,interactive = T)
          })
          plotlyOutput("ewascjy1")
        }else{
          renderPlotly({})
        }  
      })
      output$enrichmarker_tab <- DT::renderDataTable({ 
        if(isTruthy(rr1())){
          rr1()[[1]] %>%
            select(c("source","term_id" ,"term_name","term_size", "intersection_size","p_value")) %>%
            filter(intersection_size>1 & p_value<0.05) %>%
            arrange(p_value) %>%
            datatable(rownames = F,extensions = 'Buttons',
                      options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                                     buttons = list(list(extend = 'csv',  filename = paste0("Enrichments_FDR_",input$slider_4a,"_",input$sc_study,"_",input$sc_celltype4a,"")),
                                                    list(extend ='excel', filename = paste0("Enrichments_FDR_",input$slider_4a,"_",input$sc_study,"_",input$sc_celltype4a,""))), 
                                     scrollX = TRUE,scrollCollapse = TRUE)) %>% 
            formatSignif(columns = c('p_value'),digits = 2)
        }else{ DT::renderDataTable({}) }  
      })
    }
})

 ## multigene feature plot
observeEvent({variables$pldf_marker},{
    if(input$AppTab=='scstudy' &
       isTruthy(variables$pldf_marker)  ){
      
      output$sc_pcaplot4marker <- renderUI({ 
        gnames <- gsub(paste0(input$sc_study,"_"),"",names(variables$pldf_marker))
        qq <- plot_multi_umap_raster(variables$pldf_marker,
                               genenames = gnames,
                               legendPlot = F)
        output$tmpa_marker <- renderPlotly(qq)
        plotlyOutput("tmpa_marker")
      })
      
    }
})

##----------------------------------------------------------------------
## tab page 5
##----------------------------------------------------------------------
observeEvent({input$AppTab
  variables$sc_df},{
    if(input$AppTab=='scstudy' &
       isTruthy(input$sc_study) &
       isTruthy(variables$sc_df)){
      cat("ssc: tab7\n")
      sc_pcaplot5a <- reactive({
        plot_umap_raster(pl=variables$sc_df,
                         feature = variables$feature,
                         genename = "Cell type clusters",
                         legendPlot = F,
                         selectedFeature = input$sc_celltype5a,
                         source.val='scpca5',js=js)
      })
      output$sc_pcaplot5a <- renderPlotly({ 
        suppressWarnings({
          sc_pcaplot5a() %>% layout(height = 450, width = 450)
        })
      })
      output$scpca5a_hover <- renderText({
        getHoverText(variables$sc_df, 
                     colids=c('V1','V2',variables$feature), 
                     hovercoord = input$scpca5 )
      })
    }
})

##-- dotplot
observeEvent({input$AppTab
  variables$deg_dot
  input$sc_sel_5a},{
    if(input$AppTab=='scstudy' & isTruthy(input$sc_sel_5a) &
       isTRUE(nrow(variables$deg_dot)>0)){
      
      cat(" ssc: preparing the DEG dotplot\t", input$sc_sel_5a,"\t",nrow(variables$deg_dot),"\n")
      output$sc_degdotplot <- renderUI({ 
        output$asodksc <- renderPlotly({
          qq <- variables$deg_dot
          qq <- qq[qq$Test==input$sc_sel_5a,]
          sc_dge_dotplot(qq,
                         x='geneSymbol',y='Tag',
                         colCol='logFC',sizeCol='adj.P.Val',
                         testCol='Test',configlayout=T)
        })
        plotlyOutput("asodksc")
      })
    }else{
      output$sc_degdotplot <- renderUI({})
    }
  })

##--- spreadsheet and enrichment plots
observeEvent({input$AppTab
  input$slider_5a
  input$sc_sel_5a
  variables$deg_volc},{
    
    if(input$AppTab=='scstudy' & isTRUE(nrow(variables$deg_volc)>0) &
       isTruthy(input$sc_sel_5a) & isTruthy(input$slider_5a)){
      
      variables$degenrich_up = list()
      variables$degenrich_dn = list()
      
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = "ssc: printing deg data table..", value = 0)
      cat(" ssc: printing deg data table\n")
      ## spreadsheet
      output$deg_tab = DT::renderDataTable({
        cf <-  variables$deg_volc
        cf <- cf[cf$Test==input$sc_sel_5a,]
        cf <- cf[as.numeric(cf$adj.P.Val) <= as.numeric(input$slider_5a),]
        cf$logFC = round(cf$logFC,2)
        cf$AveExpr = round(cf$AveExpr,2)
        cf$t = round(cf$t,2)
        cf$B = round(cf$B,2)
        datatable(cf,rownames = F,extensions = 'Buttons',
                  options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                                 buttons = list(list(extend = 'csv',  filename = paste0("DEG_FDR_",input$slider_5a,"_",input$sc_study,"_",input$sc_celltype5a,"")),
                                                list(extend ='excel', filename = paste0("DEG_FDR_",input$slider_5a,"_",input$sc_study,"_",input$sc_celltype5a,""))), 
                                 scrollX = TRUE,scrollCollapse = TRUE
                  )) %>% formatSignif(columns = c('adj.P.Val','P.Value'),digits = 2)
      })
      
      ## enrichment
      scProg$set(message = "calculating enrichment of upregulated deg..", value = 0.2)
      variables$degenrich_up <- ({
        qf <-  variables$deg_volc
        qf <- qf[qf$Test==input$sc_sel_5a,]
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$slider_5a) & qf$logFC>0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- paste0('Upregulated in ',input$sc_celltype5a)
          gprofiler2::gost(query = qf[1],organism = 'hsapiens',
                           ordered_query = F,exclude_iea = T,
                           significant = F,
                           user_threshold = 0.05,correction_method = 'g_SCS')
        }else{ NULL }
      })
      scProg$set(message = "calculating enrichment of downregulated deg..", value = 0.6)
      variables$degenrich_dn <- ({
        qf <-  variables$deg_volc
        qf <- qf[qf$Test==input$sc_sel_5a,]
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$slider_5a) & qf$logFC < 0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- paste0('Downregulated in ',input$sc_celltype5a)
          gprofiler2::gost(query = qf[1],organism = 'hsapiens',
                           ordered_query = F,exclude_iea = T,
                           significant = F,
                           user_threshold = 0.05,correction_method = 'g_SCS')
        }else{ NULL }
      })
      scProg$set(message = "done..", value = 1)
    }else{
      output$deg_tab = DT::renderDataTable({})
      output$enrichdeg_tab = DT::renderDataTable({})
      output$enrich_deg5a  = renderPlotly({})
    }
  })

observeEvent({input$AppTab
  variables$degenrich_up
  variables$degenrich_dn
  input$t5aradio},{
    if(input$AppTab=='scstudy' & isTruthy(input$t5aradio)){
      if(input$t5aradio==1){
        rr1 <- reactive({ variables$degenrich_up  })
      }else{
        rr1 <- reactive({ variables$degenrich_dn  })
      }
      
      output$enrich_deg5a <- renderUI({ 
        if(isTruthy(rr1())){
          output$ewascjyx <- renderPlotly({
            gostplot(gostres = rr1(),capped = T,interactive = T)
          })
          plotlyOutput("ewascjyx")
        }else{
          renderPlotly({})
        }  
      })
      output$enrichdeg_tab <- DT::renderDataTable({ 
        if(isTruthy(rr1())){
          rr1()[[1]] %>%
            select(c("source","term_id" ,"term_name","term_size", "intersection_size","p_value")) %>%
            filter(intersection_size>1 & p_value<0.05) %>%
            arrange(p_value) %>%
            datatable(rownames = F,extensions = 'Buttons',
                      options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                                     buttons = list(list(extend = 'csv',  filename = paste0("Enrichments_FDR_",input$slider_5a,"_",input$sc_study,"_",input$sc_celltype5a,"")),
                                                    list(extend ='excel', filename = paste0("Enrichments_FDR_",input$slider_5a,"_",input$sc_study,"_",input$sc_celltype5a,""))), 
                                     scrollX = TRUE,scrollCollapse = TRUE)) %>% 
            formatSignif(columns = c('p_value'),digits = 2)
        }else{ DT::renderDataTable({}) }  
      })
    }
})

## multigene feature plot
observeEvent({variables$pldf_deg},{
  if(input$AppTab=='scstudy' &
     isTruthy(variables$pldf_deg)  ){
    
    output$sc_pcaplot5deg <- renderUI({ 
      gnames <- gsub(paste0(input$sc_study,"_"),"",names(variables$pldf_deg))
      qq <- plot_multi_umap_raster(variables$pldf_deg,
                                   genenames = gnames,
                                   legendPlot = F)
      output$tmpa_deg <- renderPlotly(qq)
      plotlyOutput("tmpa_deg")
    })
    
  }
})


