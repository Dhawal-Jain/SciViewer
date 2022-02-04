
##------------------------------------------------------------------------
## Input functions
##------------------------------------------------------------------------

matrix.reformat <- function(mat,meta){
  
  require(data.table)|| stop("Can't find data.table package. Please install and restart! \n")
  require(Matrix)|| stop("Can't find Matrix package. Please install and restart! \n")
  
  if(!'SAMPID'%in%names(meta)){
    stop('SAMPID column is not found in the metadata file. Please make sure to label the cell id column as SAMPID\n')
  }
  if(length(which(is(mat)=="sparseMatrix"))==0){
    warning(' provided input matrix is not in sparsematrix format!\n Converting data to sparse matrix format.')
    mat <- as.matrix(mat)
    mat <- as(mat,'dgCMatrix')
  }
  genes <- data.frame(gene_id=rownames(mat),geneSymbol=rownames(mat))
  
  
  meta <- meta[meta$SAMPID%in%colnames(mat),]
  mat <- mat[,colnames(mat)%in%meta$SAMPID]
  if(ncol(mat)!=sum(colnames(mat)%in%meta$SAMPID)){
    stop('Provided number of cells in the count matrix do not match to those in the metadata file')
  }
  
  message('Matrix dimensions: ',paste0(dim(mat),collapse = ","),"\n")
  message('Generating indexes..\n')
  rnames <- rownames(mat)
  rnames <- data.frame(geneID=rnames,id=1:length(rnames))
  cnames <- colnames(mat)
  cnames <- data.frame(cellID=cnames,id=1:length(cnames))
  mat <- as.data.frame(summary(mat))
  
  message("Regrouping the data..\n")
  message('  -->this step takes relatively long time to finish, please wait..')
  mat <- data.table(mat)
  sumfun <- function(x){ paste0(x,collapse = ",") }
  mat <- mat[, lapply(.SD, sumfun), by=list(i) ]
  mat <- as.data.frame(mat)
  names(mat)[1:3] <- c("row_index","col_index",'value')
  
  ## match the gene-dataframe with matrix index
  rnames <- rnames[match(mat$row_index,rnames$id),]
  genes <- genes[match(rnames$geneID,genes$gene_id),]
  if(nrow(mat)!=sum(genes$gene_id==rnames$geneID)){
    stop('error, the rows do not match! \n')
  }
  
  ## binding the matrix information
  mat <- cbind(genes,mat)
  rownames(mat) <- NULL
  
  ## re-arrange meta dataframe rows to match original matrix index
  meta <- meta[match(cnames$cellID,meta$SAMPID),]
  meta$col_index <- cnames$id
  
  return(list(A=mat,B=meta))
  
}


seurat2sqlite <- function(so=NULL,si_study=NULL,
                          si_reduction=NULL,
                          si_compute_cellmarker = F,si_cell_markers=NULL,si_cell_col = NULL,
                          si_compute_diseasemarker = F,si_disease_markers=NULL,si_disease_col=NULL,si_celltype=NULL,
                          db_address=NULL,OUTDIR=NULL){
  
  require(Seurat) || stop("Can't find Seurat package. Please install and restart! \n")
  require(RSQLite) || stop("Can't find RSQLite package. Please install and restart! \n")
  
  #-- dependencies
  if(is.null(so)) stop(" Please provide Seurat object for recording the study\n")
  if(is.null(si_study)) stop(" Please provide study details for recording the study\n")
  if(is.null(db_address)) stop(" Please provide DB name for recording the study\n")
  
  #-- collect variables
  si_assays<- Seurat::Assays(so)
  si_coords <- names(so@reductions)
  if(!is.null(si_reduction) & si_reduction %in% names(so@reductions)){
    si_default_coord <- si_reduction
  }else{
    si_default_coord <- ifelse(grep('^umap$|^UMAP$',si_coords),
                               si_coords[grep('^umap$|^UMAP$',si_coords)],
                               si_coords[1])
  }
  si_study$COORDS <- paste0(si_coords,collapse = ",")
  si_study$DEFAULT_COORDS <- si_default_coord
  
  #-- data normalization etc.
  if(length(si_assays)>1){
    warning('The provided Seurat object contains multiple assay objects.Please make sure that for non-RNA assays you have appropriately normalized the data.RNA assay object will be tested by default and log-normalized in case non-normalized data is found.\n')
  }
  if(range(so[[si_assays[1]]]@data)[2]>25){
    message('Normalizing the RNA data\n')
    so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)  
  }
  
  #-- create metadata file
  if(!'percent.mt'%in%names(so@meta.data)) so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  if(!'percent.RPS'%in%names(so@meta.data)) so[["percent.RPS"]] <- PercentageFeatureSet(so, pattern = "^RPS")
  if(!'percent.RPL'%in%names(so@meta.data)) so[["percent.RPL"]] <- PercentageFeatureSet(so, pattern = "^RPL")
  my.meta <- local({
    coords <- data.frame(SAMPID=rownames(so@meta.data))
    qq <- as.data.frame(so@reductions[[si_default_coord]]@cell.embeddings)
    qq <- qq[,1:2]
    names(qq) <- c('V1','V2')
    coords <- cbind(coords, qq)
    rm(qq)
    for(f in si_coords){
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
  
  #-- loop over assays 
  for(my_assay in si_assays){
    message("working on assay:", my_assay,"\n")
    DefaultAssay(so) <- my_assay
    si_cell_markers_tmp <- c()
    si_disease_markers_tmp <- c()
    
    #-- compute cell type markers
    if(si_compute_cellmarker==T & !is.null(si_cell_col)){
      if(is.null(si_cell_col)) stop(" Column name/s of the variable against which markers are to be computed is not provided!\n")
      if(sum(si_cell_col%in%names(so@meta.data))!=length(si_cell_col)){
        stop(' one or more marker column names are not present in the Seurat object\n')
      }
      message(" Computing markers for the variable/s:", si_cell_col,"\n")
      if(length(si_cell_col)>1){
        message(' -> Note: multiple varialbes are provided for computing the markers. Will proceed!\n')
      }
      for(mycolm in si_cell_col){
        message("  ..computing markers for the variable:", mycolm,"\n")
        Idents(so) <- mycolm
        my.markers <- tryCatch({FindAllMarkers(so,  
                                               min.pct = 0.25, 
                                               logfc.threshold = 0.1,
                                               return.thresh = 0.1,
                                               test.use = 'wilcox')
        },error=function(e){
          return(c())
        })
        if(!is.null(my.markers)){
          names(my.markers) <- c("P.Value", "logFC", "pct.1", "pct.2", 
                                 "adj.P.Val", "Tag", "geneSymbol")
          my.markers$Test <- paste0(my.markers$Tag,"-Rest")
          my.markers$AveExpr <- NA
          my.markers$t <- NA
          my.markers$B <- NA
          my.markers$Assay <- my_assay
          my.markers$Method <- 'Wilcox'
          my.markers <- my.markers[,c("logFC", "AveExpr", "t", 
                                      "P.Value", "adj.P.Val", "B", "Test","Tag",
                                      "geneSymbol",'Method','Assay',"pct.1","pct.2")]
          si_cell_markers_tmp <- rbind(si_cell_markers_tmp,my.markers)
        }
        rm(my.markers)
      }
      rm(mycolm)
    }else if(si_compute_cellmarker ==F & isTruthy(si_cell_markers)){
      mycols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Test","Tag","geneSymbol")
      if(sum(mycols%in%names(si_cell_markers))!=length(mycols)){
        warning("One or more columns are missing in the provided cell marker file. Please ensure that correct file is provided. Skipping this file and moving on!")
      }else{
        mycols <- append(mycols,
                         names(si_cell_markers)[!names(si_cell_markers)%in%mycols]
        )
        si_cell_markers_tmp <- si_cell_markers[,mycols]
      }
    }
    
    #-- compute disease markers
    if(si_compute_diseasemarker==T &  !is.null(si_disease_col) & !is.null(si_celltype)){
      my.disease <- levels(so@meta.data[,si_disease_col])
      my.celltypes <- levels(so@meta.data[,si_celltype])
      if(length(my.disease)<2) stop(" the disease/exp group doesn't have contrast\n")
      if(length(si_celltype)>1) stop(" please define unique cell type annotation column\n")
      if(length(si_disease_col)>1) stop(" please define unique disease/experimental condition annotation column\n")
      if(!si_disease_col%in%names(so@meta.data)) stop(' the disease/experimental condition column is not found in the metadata\n')
      my.tests <- expand.grid(my.disease,my.disease,stringsAsFactors = F)
      my.tests <- my.tests[my.tests$Var1>my.tests$Var2,]
      message(" Performing the differential expression tests using provided disease/experiment contrast.\n")
      message("  ..total ",nrow(my.tests)," contrasts are found. Starting the execution!\n")
      message("  ..this can take a while.\n")
      
      so$temp_celltype <- so@meta.data[,si_celltype]
      Idents(so) <- si_disease_col
      for(i in 1:nrow(my.tests)){
        sa <- subset(x = so, idents = c(my.tests$Var1[i],my.tests$Var2[i]))
        for(f in my.celltypes){
          sb <- tryCatch({
            subset(x = sa, subset = temp_celltype == f)
          },error = function(e){
            return(c())
          })
          if(is.null(sb)) next
          message('  ..working on celltype: ',f ,"  (",my.tests$Var1[i],"-", my.tests$Var2[i],")\n")
          q <- tryCatch({
            FindMarkers(sb, ident.1 = my.tests$Var1[i], 
                        ident.2 = my.tests$Var2[i],print.bar = T,
                        min.pct = 0.25, 
                        logfc.threshold = 0.1,
                        return.thresh = 0.1,min.cells.group = 3,
                        test.use = 'wilcox')
          },error=function(e){
            return(c())
          })
          if(!is.null(q)){
            names(q) <- c("P.Value", "logFC", "pct.1", "pct.2", 
                          "adj.P.Val")
            q$Test <- paste0(my.tests$Var1[i],"-", my.tests$Var2[i])
            q$Tag <- f
            q$geneSymbol <- rownames(q)
            q$Method <- 'Wilcox'
            q$Assay <- my_assay
            q$AveExpr <- NA
            q$t <- NA
            q$B <- NA
            q <- q[,c("logFC", "AveExpr", "t","P.Value", "adj.P.Val",
                      "B", "Test","Tag","geneSymbol",'Method','Assay',
                      "pct.1","pct.2")]
            si_disease_markers_tmp <- rbind(si_disease_markers_tmp,q)
          }
          rm(q,sb)
        }
      }
      so$temp_celltype <- NULL
      
    }else if(si_compute_diseasemarker==F & !is.null(si_disease_col) & !is.null(si_celltype)){
      mycols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Test","Tag","geneSymbol")
      if(sum(mycols%in%names(si_disease_markers))!=length(mycols)){
        warning("One or more columns are missing in the provided disease marker file. Please ensure that correct file is provided. Skipping this file and moving on!")
      }else{
        mycols <- append(mycols,
                         names(si_disease_markers)[!names(si_disease_markers)%in%mycols]
        )
        si_disease_markers_tmp <- si_disease_markers[,mycols]
      }
    }
    
    #-- matrix reformatting
    message('Reformatting the matrix for writing down\n')
    l<- matrix.reformat(mat=so@assays[[my_assay]]@data,
                        meta=my.meta)
    
    #-- write down
    message('Writing down the matrix and study details into SQLite database\n')
    if(isTRUE(nrow(si_study)>0) & isTRUE(nrow(l[[1]])>0) & isTRUE(nrow(l[[2]])>0) ){
      StudyName = paste0(my_assay,"_",si_study$Database)
      si_study$Database = StudyName
      if(length(grep("^scRNA_",db_address))==0){
        db_address <- paste0("scRNA_",my_assay,"_",db_address)
      }
      if(!is.null(OUTDIR)){
        connSc <- RSQLite::dbConnect(RSQLite::SQLite(),paste0(OUTDIR,"/",db_address))
      }else{
        connSc <- RSQLite::dbConnect(RSQLite::SQLite(),paste0(db_address))
      }
      RSQLite::dbWriteTable(connSc,paste0(StudyName,"_study"),si_study,overwrite=T)
      RSQLite::dbWriteTable(connSc,paste0(StudyName,"_metaFeatures"),l[[2]],overwrite=T)
      RSQLite::dbWriteTable(connSc,paste0(StudyName,"_data"),l[[1]],overwrite=T)
      if(isTRUE(nrow(si_cell_markers_tmp)>0)){
        RSQLite::dbWriteTable(connSc,paste0(StudyName,"_Marker"),si_cell_markers_tmp,overwrite=T)
      }
      if(isTRUE(nrow(si_disease_markers_tmp)>0)){
        RSQLite::dbWriteTable(connSc,paste0(StudyName,"_DEG"),si_disease_markers_tmp,overwrite=T)
      }
      RSQLite::dbListTables(connSc)
      RSQLite::dbDisconnect(connSc)
    }
  }
  return(1)
}


write.scstudy2.sqlitedb <- function(so,
                                    db_address=NULL,
                                    StudyName=NULL,
                                    Celltype = NULL,
                                    Reduction_map = NULL,
                                    Donors_VariableName=NULL,
                                    Disease_VariableName=NULL,
                                    Marker_Calc = F,
                                    Marker_Covariates = NULL,
                                    Marker_Precomputed=NULL,
                                    DE_Calc=F,
                                    DE_Precomputed = NULL,
                                    StudyDescr=NULL,
                                    Tissue=NULL,
                                    PMID=NULL,
                                    GEO=NULL,
                                    StudyStatus='Internal',
                                    StudyRating='High',
                                    Continuous_Vars=NULL,
                                    Categorical_Vars=NULL,
                                    OUTDIR=NULL){
  
  ##---- dependencies
  require(dplyr) || stop("Can't find dplyr package. Please install and restart! \n")
  require(Matrix) || stop("Can't find Matrix package. Please install and restart! \n")
  require(data.table) || stop("Can't find data.table package. Please install and restart! \n")
  require(RSQLite) || stop("Can't find RSQLite package. Please install and restart! \n")
  require(Seurat) || stop("Can't find Seurat package. Please install and restart! \n")
  
  scProg <- shiny::Progress$new()
  on.exit(scProg$close())
  scProg$set(message = "assessing dependancies..", value = 0)
  
  ##-- prechecks
  if(is.null(db_address)){
    stop('Please provide valid database address where the data is to be written')
  }
  if(length(grep("[.]db",db_address))==0){
    #db_address <- gsub("_sc$","",db_address)
    #db_address <- paste0(db_address,"_sc.db")
    db_address <- paste0(db_address,".db")
  }else{
    #db_address <- gsub(".db","_sc.db",db_address)
  }
  message("database address: ",db_address,"\n")
  if(is.null(StudyName)){
    stop('Please provide StudyName using which you would like to access this particualr data in the SingleCellApp')
  }
  if(is.null(StudyDescr)){
    warning('No study description is provided. You will see this section as blank in the Application\n')  
  }else{
    message(StudyDescr,"\n\n")
  }
  if(is.null(Donors_VariableName)){
    warning('Donors variable name is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Identifier for donors: ',Donors_VariableName,"\n\n")
  }
  if(is.null(Tissue)){
    warning('Tissue used for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Identifier for tissue: ',Tissue,"\n\n")
  }
  if(is.null(Disease_VariableName)){
    warning('Disease variable (e.g. experimental condition that is not control/healthy etc.) for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Identifier for disease: ',Disease_VariableName,"\n\n")
  }
  if(is.null(PMID)){
    message('(If public), the PMID for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('PMID: ',PMID,"\n\n")
  }
  if(is.null(GEO)){
    message('(If public), the data link for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Data link (if public): ',GEO,"\n\n")
  }
  if(is.null(StudyStatus)){
    warning('If public/Internal, the status for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Study status: ',StudyStatus,"\n\n")
  }
  if(is.null(StudyRating)){
    warning('If public/Internal, the status rating (e.g. high confidence/ poor study etc.) for the experiment is not provided. This information in the summary section will remain blank\n')  
  }else{
    message('Study rating: ',StudyRating,"\n\n")
  }
  if(is.null(Continuous_Vars)){
    warning('No continuous varibles in the metadata file are identified. This information in the summary section will remain blank\n')  
  }else{
    qq <- unlist(strsplit(as.character(Continuous_Vars),","))
    if(sum(qq%in%names(so@meta.data))!=length(qq) ){
      stop(" could not find provided continuous variables in the metadata file\n") 
    }
    message('Continuous variables in the metadata: ',Continuous_Vars,"\n\n")
  }
  if(is.null(Celltype)){
    stop("cell type column is not declared in the input. Please declare the column and rerun the function\n")
  }else{
    message(' cell type variable in the metadata: ', Celltype, "\n\n")
    so@meta.data <- local({
      meta <- so@meta.data
      names(meta) <- ifelse(names(meta)==Celltype,as.character("cell_type"),as.character(names(meta)))
      meta
    })
  }
  if(is.null(Categorical_Vars)){
    stop('No categorical varibles in the metadata file are identified. This information in the summary section will remain blank\n')  
  }else{
    qq <- unlist(strsplit(as.character(Categorical_Vars),","))
    qq <- ifelse(qq==Celltype,as.character("cell_type"),as.character(qq))
    if(sum(qq%in%names(so@meta.data))!=length(qq) ){
      #cat(paste0(qq,collapse = ","),"\n",
      #    paste0(names(so@meta.data),","),"\n")
      stop(" could not find provided categorical variables in the metadata file\n") 
    }
    if(!'cell_type'%in%qq){
      stop(" could not find column with name 'cell_type' in the metadata file\n") 
    }
    Categorical_Vars = qq
    message('Categorical variables in the metadata: ',Categorical_Vars,"\n\n")
  }
  if(!is.null(Marker_Covariates)){
    qq <- unlist(strsplit(as.character(Marker_Covariates),","))
    qq <- ifelse(qq==Celltype,as.character("cell_type"),as.character(qq))
    Marker_Covariates = qq
    message(" Covariates for which to calculate marker genes: ", Marker_Covariates,"\n")
    rm(qq)
  }
  
  study <- data.frame(Database=ifelse(!is.null(StudyName),as.character(StudyName),NA),
                      Description=ifelse(!is.null(StudyDescr),as.character(StudyDescr),NA),
                      SampleSize=ifelse(!is.null(Donors_VariableName),length(unique(meta[,Donors_VariableName])),NA),
                      TISSUES=ifelse(!is.null(Tissue),as.character(Tissue),NA),
                      DISEASE_Variable=ifelse(!is.null(Disease_VariableName),as.character(Disease_VariableName),NA), 
                      PMID=ifelse(!is.null(PMID),as.character(PMID),NA),
                      GEO=ifelse(!is.null(GEO),as.character(GEO),NA),
                      STATUS=ifelse(!is.null(StudyStatus),as.character(StudyStatus),NA),
                      RATING=ifelse(!is.null(StudyRating),as.character(StudyRating),NA),
                      CONTVAR=ifelse(!is.null(Continuous_Vars),as.character(Continuous_Vars),NA),
                      CATVAR=ifelse(!is.null(Categorical_Vars),as.character(Categorical_Vars),NA))
  
  
  
  seurat2sqlite(so=so,
                si_study=study,
                si_reduction=Reduction_map,
                si_compute_cellmarker = Marker_Calc,
                si_cell_markers=Marker_Precomputed,
                si_cell_col = Marker_Covariates,
                si_compute_diseasemarker = DE_Calc,
                si_disease_markers=DE_Precomputed,
                si_disease_col=Disease_VariableName,
                si_celltype='cell_type',
                db_address=db_address,
                OUTDIR=OUTDIR)
  
  scProg$set(message = 'done!', value = 1)
  
  return(1)
}

