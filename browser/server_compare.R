# server_compare.R

observeEvent({variables$sc_df},{
    if(isTruthy(variables$sc_df)){
      
      cat( " comparestudy: bracket-1\n\n")
      
      output$firstsc <- renderText({ 
        paste("<b> Primary study </b>", input$sc_study,"<br>")
      })
      ###----------------------------------------------------
      ## study 1
      sc_pcaplot21 <- reactive({
        plot_umap_raster(pl=variables$sc_df,
                         feature = variables$feature,
                         genename = 'Cluster map',
                         legendPlot = F,
                         selectedFeature = input$sc_celltype21,
                         source.val='sc_pca21',js=js)
      })
      output$scpca21_hover <- renderText({
        getHoverText(variables$sc_df, 
                     colids=c('V1','V2',variables$feature), 
                     hovercoord = input$sc_pca21 )
      })
      range21 <- reactive({ 
        zoom <- event_data("plotly_relayout", "sc_pca21")
        getZoomrange(zoom)
      })
      cat("range 2_1: ",range21()[[1]],"\n")
      ### expression plots
      sc_pcaplot23 <- reactive({
        plot_multi_umap_raster(list(variables$sc_df),
                               genenames = variables$gene,
                               legendPlot = F,xrange =range21()[[1]],
                               yrange=range21()[[2]],configlayout=F)
      })
      cat(" got plot-1 \n\n")
      output$sc_pcaplot21 <- renderPlotly({ 
        if(isTruthy(sc_pcaplot21() )){
          suppressWarnings({
            sc_pcaplot21() %>%layout(height = 400, width = 400)
          })
        }
      })
      output$sc_pcaplot23 <- renderPlotly({ 
        if(isTruthy(sc_pcaplot23())){
          suppressWarnings({
            sc_pcaplot23() %>% layout(height = 400, width = 400)
          })
        }
      })
      
    }
})


# &isTruthy(input$t2action1)
observeEvent({input$AppTab
  variables$comp_df},{
  if(input$AppTab=='compare' & 
     isTruthy(variables$comp_df)){
    
    cat( " comparestudy: bracket-2\n\n")
    output$comp_study_attribs <- renderText({ 
      variables$comp_attribute
    })
    
    ##----------------- study 2
    sc_pcaplot22 <- reactive({
      plot_umap_raster(pl=variables$comp_df,
                       feature = variables$comp_feature,
                       genename = 'Cluster map',
                       legendPlot = F,
                       selectedFeature = input$sc_celltype22,
                       source.val='sc_pca22',js=js)
    })
    output$scpca22_hover <- renderText({
      getHoverText(variables$comp_df, 
                   colids=c('V1','V2',variables$comp_feature), 
                   hovercoord = input$sc_pca22 )
    })
    cat("variables$comp_feature ", variables$comp_feature,"\n")
    range22 <- reactive({ 
      zoom <- event_data("plotly_relayout", "sc_pca22")
      getZoomrange(zoom)
    })
    cat("range 2_2: ",range22()[[1]],"\n")
    
    ### expression plots
    sc_pcaplot24 <- reactive({
      plot_multi_umap_raster(list(variables$comp_df),
                             genenames = variables$gene,
                             legendPlot = F,xrange =range22()[[1]],
                             yrange=range22()[[2]],configlayout=F)
    })
    cat(" got plot-2 \n\n")
    output$sc_pcaplot22 <- renderPlotly({ 
      if(isTruthy(sc_pcaplot22() )){
        suppressWarnings({
          sc_pcaplot22() %>% layout(height = 400, width = 400)
        })
      }
    })
    output$sc_pcaplot24 <- renderPlotly({ 
      if(isTruthy(sc_pcaplot24())){
        suppressWarnings({
          sc_pcaplot24() %>% layout(height = 400, width = 400)
        })
      }
    })
    

  }
})

