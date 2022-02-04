library(shiny)
library(shinyjqui)
library(shinydashboard)
library(shinyFiles)
library(shinyWidgets)

#library(anndata)
library(Seurat)
library(reticulate)

library(data.table)
library(Matrix)
library(SparseM)
library(dplyr)

library(limma)
library(sva)
library(edgeR)

library(RSQLite)


ROOTS=c(workdir='.',
        datadir='C:/',
        home='/home/',
        shinydata='/srv/shiny-server/data/')

