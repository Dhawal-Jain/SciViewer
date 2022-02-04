if(F){
  source("input_functions.R",local = TRUE,encoding = "UTF-8")
  source("vars.R",local = TRUE,encoding = "UTF-8")
  
  so <- readRDS("/home/dhawal_jain/s3shared/data/PDD_Pilot1/PDD_SingleCellPilot/SeuratObjects/PDD_pilot1_qc00_wclusters_s2_Lymphoid.rds")
  so@active.ident <- factor(so$IntAnno)
  sampled.cells <- sample(x = colnames(so), 
                          size = 10000, replace = F)
  so <- subset(so, cells = sampled.cells)
  rm(sampled.cells)
  
  so <- local({
    ad <- anndata::read_h5ad(filename = "C:/Dhawal/scRNASeq_data/PCL/tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad")
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
  
  
  meta <- local({
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
  
  
  study <- data.frame(Database="abc",
             Description=NA,
             SampleSize=NA,
             TISSUES=NA,
             DISEASE_Variable=NA, 
             PMID=NA,
             GEO=NA,
             STATUS='Public',
             RATING='High',
             CONTVAR=NA,
             CATVAR='cell_ontology_class')
  
    seurat2sqlite(so=so,
                  si_study=study,
                  si_reduction='umap',
                  si_compute_cellmarker = F,
                  si_cell_markers=NULL,
                  si_cell_col = NULL,
                  si_compute_diseasemarker = F,
                  si_disease_markers=NULL,
                  si_disease_col=NULL,
                  si_celltype='cell_ontology_class',
                  db_address='asd',
                  OUTDIR=NULL)
  
    isTRUE(ncol(si_cell_markers)>0)
  
}