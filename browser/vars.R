library(shiny)
library(shinyjs)
library(dashboardthemes)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(shinyBS)
library(shinyFiles)
#library(shinyjqui) ## can add draggable/movable/resizable plots 

library(highcharter)
library(ggvis)
library(plotly)
library(rasterly)
library(png)
library(data.table)
library(DT)
library(pddExpn)
library(optparse)
library(gprofiler2)

library(RSQLite)
library(RMySQL)
library(DBI)

#DATADIR="C:/Dhawal/SHINYAPP_DATA/"
GENEFILE="C:/Dhawal/SHINYAPP_DATA/hs_genes.db"
SCDBFILE="DB.txt"
SOURCEREF='SciView'

##-- if using remote database
USE_REMOTE_DB = FALSE
REPO_NAME="Minuteman"


