## following functions are highly customized, they should become part of pddExpn eventually

queryDB <- function(HANDLER=NULL, QUERY=NULL,
                    REPO_NAME=NULL,USE_REMOTE_DB=FALSE){
  
  if(is.null(HANDLER) | is.null(QUERY)){
    return(NULL)
  }
  if(USE_REMOTE_DB==FALSE){
    return(RSQLite::dbGetQuery(HANDLER, QUERY)) 
  }
  else if(USE_REMOTE_DB==TRUE){
    #QUERY = paste0("SELECT * FROM    dfg WHRE x = 'from'")
    #HANDLER= 'qqq'
    #gsub("FROM\\s+",paste0("FROM ",HANDLER,"."),QUERY)
    QUERY = gsub("FROM\\s+",paste0("FROM ",HANDLER,"."),QUERY)
    return(GetData(AppName=REPO_NAME, sSQL=QUERY))
  }
  else{
    return(NULL)
  }
}

get_marker_df <- function(connSc,study){
  
  markers <- matrix(NA,nrow = 0,ncol=11)
  markers <- as.data.frame(markers)
  names(markers) <- c("cluster", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", 
                      "gene", "cell_type", "test", "id", "top10")
  
  query <- paste0('SELECT * FROM ',study,"_markers")
  cf <- tryCatch({
    queryDB(HANDLER=connSc, QUERY=query,REPO_NAME=REPO_NAME,USE_REMOTE_DB=USE_REMOTE_DB)
    #RSQLite::dbGetQuery(connSc,query)
  },error = function(e){
    cf <- markers
    return(cf)
  })
  return(cf)
  
}

get_plot_df_sc04 <- function(connSc,study="Madisson_LungTissue",
         genename="WNT4",louvain,REPO_NAME=NULL,USE_REMOTE_DB=FALSE){
  #louvain must have a value column
  if(!'value'%in%names(louvain)){
    louvain$value <- 0
  }
  
  query <- paste0("SELECT * FROM ",study,"_data WHERE geneSymbol = '", genename,"'")
  gdf <- queryDB(HANDLER=connSc, 
          QUERY=query,REPO_NAME=REPO_NAME,
          USE_REMOTE_DB=USE_REMOTE_DB)
  
  if(nrow(gdf)>0){
    cnta <- matrix(0,nrow = nrow(louvain),ncol=nrow(gdf)) %>% as.data.frame
    for(i in 1:nrow(gdf)){
      gdfa <- data.frame(col_index= as.numeric(unlist(strsplit(gdf[i,]$col_index,","))),
                         value = as.numeric(unlist(strsplit(gdf[i,]$value,","))))
      gdfa <- gdfa[order(gdfa$col_index,decreasing = F),]
      cnta[gdfa$col_index,i] <- gdfa$value
      rm(gdfa,i)
    }
    if(ncol(cnta)>1){
      cnta <- rowMeans(cnta)
    }else{
      cnta <- cnta$V1
    }
    louvain$value <- cnta
  }else{
    louvain <- louvain[0,]
  }
  return(louvain)
}

hchart_cor <- function(mcor,sourceref='pDDSingleCellApp') {
  
  mcor <- round(mcor,2)
  mcor[is.na(mcor)] = 0
  
  #fntltp <- JS("function(){
  #                return this.series.xAxis.categories[this.point.x] + ' ~ ' +
  #                       this.series.yAxis.categories[this.point.y] + ': <b>' +
  #                       Highcharts.numberFormat(this.point.value, 2)+'</b>';
  #             ; }")
  cor_colr <- list( list(0, '#FF5733'),
                    list(0.5, '#F8F5F5'),
                    list(1, '#2E86C1')
  )
  
  hs <- hchart(mcor) %>%
    hc_plotOptions(
      series = list(
        boderWidth = 0,
        dataLabels = list(enabled = TRUE)
      ))   %>%
    hc_legend(align = "right", layout = "vertical") %>%
    hc_colorAxis(  stops= cor_colr,min=-1,max=1)
  hs <- hs %>%
    hc_credits(
      enabled = TRUE,
      text = paste0("Source:",sourceref),
      href = "",
      style = list(fontSize = "11px")
    )%>%
    hc_exporting(enabled=T
    ) %>%
    hc_scrollbar(
      barBackgroundColor = "gray",
      barBorderRadius = 7,
      barBorderWidth = 0,
      buttonBackgroundColor = "gray",
      buttonBorderWidth = 0,
      buttonArrowColor = "yellow",
      buttonBorderRadius = 7,
      rifleColor = "yellow",
      trackBackgroundColor = "white",
      trackBorderWidth = 1,
      trackBorderColor = "silver",
      trackBorderRadius = 7
    ) %>%
    hc_chart(zoomType = "xy",
             borderColor = "#EBBA95",
             borderRadius = 10,
             borderWidth = 3,
             plotShadow = F,
             allowForce = T,
             turboThreshold =1,
             allowForce=T,
             animation=T,
             boostThreshold = 1,
             usePreallocated = T,
             useGPUTranslations =T,
             seriesThreshold = 2
    )
  
  return(hs)
}



