#!/usr/bin/env Rscript
require("optparse")
require("getopt")
source("input_functions.R",local = TRUE,encoding = "UTF-8")
source("vars.R",local = TRUE,encoding = "UTF-8")

option_list = list(
  
  ##------- Input data files
  make_option(c("-Z", "--h5ad"), type="character", default=NULL, 
              help="If the data file is in the H5AD format, then the file path [default= NULL]", metavar="character"),
  make_option(c("-W", "--seurat"), type="character", default=NULL, 
              help="If the data file is in the Seurat format, then the file path [default= NULL]", metavar="character"),

  ##------- Mandatory inputs
  make_option(c("-a", "--db_address"), type="character", default=NULL, 
              help="(*) full address where the database file is to be written [default= NULL]", metavar="character"),
  make_option(c("-b", "--StudyName"), type="character", default=NULL, 
              help="(*) name using which you would like to identify the study in the application [default= NULL]", metavar="character"),
  make_option(c("-d", "--celltypeCol"), type="character", default=NULL, 
              help="(*) name of cell type annotation column in the metadata file [default= NULL]", metavar="character"),
  make_option(c("-e", "--UMAPs"), type="character", default=NULL, 
              help="(*) comma-separated names of column where the X/Y coodinate information is stored in the metadata file  [default= NULL]", metavar="character"),
  
  ##------- Inputs for DE calculation
  make_option(c("-f", "--lognormalize"), type="logical", default=F, 
              help="whether to log-normalize the data [default= %default]", metavar="logical"),
  make_option(c("-g", "--donorCol"), type="character", default=NULL, 
              help="Name of donor annotation column in the metadata file [default= NULL]", metavar="character"),
  make_option(c("-i", "--diseaseVarCol"), type="character", default=NULL, 
              help="Name of disease variable column in the metadata file [default= NULL]", metavar="character"),
  make_option(c("-j", "--covariates"), type="character", default=NULL, 
              help="Comma-separated list of covariates from the metadata file to be used for modelling the counts and calculating of differential expression [default= NULL]", metavar="character"),
  
  ##------- Study details
  make_option(c("-k", "--StudyDescr"), type="character", default=NULL, 
              help="detailed description of the study, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("-l", "--Tissue"), type="character", default=NULL, 
              help="tissue used in the study, this info will be reflected in the summary section", metavar="character"),
  make_option(c("-m", "--PMID"), type="character", default=NULL, 
              help="publication id if available, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("-n", "--GEO"), type="character", default=NULL, 
              help="data repository link if available, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("-o", "--StudyStatus"), type="character", default='Internal', 
              help="status of the study (e.g. public/internal), this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("-p", "--StudyRating"), type="character", default='High', 
              help="your rating for the study (e.g. high/poor/average), this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("-q", "--Continuous_Vars"), type="character", default=NULL, 
              help="continuous variables identified from the metadata file, this info will be used for data display [default= %default]", metavar="character"),
  make_option(c("-r", "--Categorical_Vars"), type="character", default=NULL, 
              help="categorical variables identified from the metadata file, this info will be used for data display. column name 'cell_type' is required [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


###--------------------- Input tests
if(is.null(opt$db_address)){
  stop('path to write database is not provided!')
}
if(is.null(opt$UMAPs)){
  stop('comma-separated x/y column names for clusters not provided!')
}
if(is.null(opt$barcodeCol)){
  stop('column name for cell barcodes not provided!')
}
if(is.null(opt$celltypeCol)){
  stop('column name for cell type annotations not provided!')
}
if(is.null(opt$Categorical_Vars)){
  stop('categorical variables in the metadata (at least cell_type column) not specified!')
}

###--------------------- H5AD to data objects
if(!is.null(opt$h5ad)){
  require(anndata)
  require(dplyr)
  require(Matrix)
  h5adpath=opt$h5ad
  #h5adpath="/gpfs01/home/glxaw/data/scRNASeq_datasets/Tabula_Muris/tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad"
  #h5adpath='/gpfs01/home/glxaw/data/scRNASeq_datasets/Barbry_HealthyAirways/Barbry_HealthyAirways_Lung_fromDhawal_mahmoudConvert_seuratDisk.h5ad'
  ad <- anndata::read_h5ad(filename = h5adpath)
  meta <- ad$obs
  meta$V1 <- meta$V2 <- NULL
  mat <- ad$X
  mat <- t(as.matrix(mat))
  mat <- as(mat,'dgCMatrix')
  cat(dim(mat)[2],"\n")
  scvis <- ad$obsm
  if(length(grep('umap',names(scvis)))>0){
    cat('.. using UMAP coordinates\n')
    scvis <- scvis[[grep('umap',names(scvis))]] %>% as.data.frame()
    names(scvis) <- c('V1','V2')
  }else if(length(grep('tsne',names(scvis)))>0){
    cat('.. using t-sne coordinates\n')
    scvis <- scvis[[grep('tsne',names(scvis))]] %>% as.data.frame()
    names(scvis) <- c('V1','V2')
  }else{
    stop(' x/y coordinates for the cluster do not exist!')
  }
  meta <- cbind(meta,scvis)
  meta$SAMPID = rownames(meta)
  rownames(meta) = NULL
  rm(scvis,ad)
  if(sum(colnames(mat)==meta$SAMPID)!=nrow(meta)){
    stop(' sample ids in the metadata file and count matrix do not match! ')
  }
}else if(!is.null(opt$seurat)){
  so <- readRDS(file = opt$seurat)
  l <- seurat2objects(so)
  mat <- l[[1]]
  meta <- l[[2]]
  rm(l,so)
}else{
  load(file = opt$mat)
  load(file = opt$metadata)
  if(!'meta'%in%ls() | !'mat'%in%ls()){
    stop(' provided input does not have mattrix (mat) and metadata (meta) files!')
  }
  qq <- unlist(strsplit(opt$UMAPs,","))
  names(meta) <- ifelse(names(meta)==qq[1],as.character("V1"),as.character(names(meta)))
  names(meta) <- ifelse(names(meta)==qq[2],as.character("V2"),as.character(names(meta)))
  rm(qq)
}
names(meta) <- ifelse(names(meta)==opt$celltypeCol,as.character("cell_type"),as.character(names(meta)))
names(meta) <- ifelse(names(meta)==opt$barcodeCol,as.character("SAMPID"),as.character(names(meta)))
qq <- unlist(strsplit(opt$Categorical_Vars,","))
qq <- ifelse(qq==opt$celltypeCol,as.character("cell_type"),as.character(qq))
opt$Categorical_Vars <- paste0(qq,collapse = ",")
rm(qq)


###--------------------- generate database
write.scstudy2.sqlitedb(mat = mat,
                        meta = meta,
                        lognormalize=opt$lognormalize,
                        db_address=opt$db_address,
                        StudyName=opt$StudyName,
                        Donors_VariableName = opt$donorCol,
                        DE_Covariates = opt$covariates,
                        Disease_VariableName=opt$Disease_VariableName,
                        StudyDescr=opt$StudyDescr,
                        Tissue=opt$Tissue,
                        PMID=opt$PMID,
                        GEO=opt$GEO,
                        StudyStatus=opt$StudyStatus,
                        StudyRating=opt$StudyRating,
                        Continuous_Vars=opt$Continuous_Vars,
                        Categorical_Vars=opt$Categorical_Vars)


if(F){
  
  setwd("~/data/scRNASeq_datasets/LungAtlas_Nasal")
  load("temp.RData")
  load("anno.RData")
  lognormalize=F
  db_address='qq.db'
  StudyName='a'
  StudyDescr=NULL
  Tissue=NULL
  Disease_VariableName=NULL
  PMID=NULL
  GEO=NULL
  StudyStatus='Internal'
  StudyRating='High'
  Continuous_Vars=NULL
  Categorical_Vars='cell_type'
  Donors_VariableName='mouse.id'
  DE_Covariates='age,sex'
  
  write.scstudy2.sqlitedb(mat,meta,lognormalize=F,db_address=db_address,
                          StudyName=StudyName,StudyDescr=StudyDescr,
                          Donors_VariableName=Donors_VariableName,
                          Tissue=Tissue,
                          Disease_VariableName=Disease_VariableName,
                          DE_Covariates = DE_Covariates,
                          PMID=PMID,GEO=GEO,
                          StudyStatus=StudyStatus,
                          StudyRating=StudyRating,
                          Continuous_Vars=Continuous_Vars,
                          Categorical_Vars=Categorical_Vars)
  
  
  #Rscript --vanilla write_sc2sqlite.R -Z /gpfs01/home/glxaw/data/scRNASeq_datasets/Barbry_HealthyAirways/Barbry_HealthyAirways_Lung_fromDhawal_mahmoudConvert_seuratDisk.h5ad -a qq.db -b a -c SAMPID -d cell_type -e 'V1,V2' -f F -g Donor -r 'cell_type'

  lognormalize=F
  db_address='qq.db'
  StudyName='a'
  StudyDescr=NULL
  Tissue=NULL
  Disease_VariableName='Disease'
  PMID=NULL
  GEO=NULL
  StudyStatus='Internal'
  StudyRating='High'
  Continuous_Vars=NULL
  Categorical_Vars='cell_type'
  Donors_VariableName='Donor'
  DE_Covariates=NULL
  
  use_python("C:/Users/glxaw/AppData/Local/r-miniconda/python.exe")
  library(anndata)
  h5adpath="C:/Dhawal/tmp/tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad"
  ad <- anndata::read_h5ad(filename = h5adpath)
  meta <- ad$obs
  meta$V1 <- meta$V2 <- NULL
  mat <- ad$X
  
  
}






