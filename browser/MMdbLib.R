library(RMySQL)
library(DBI)

ConnectDB<-function(AppName)
{
  db_user <- 'admin'
  db_password <- 'HpgcDev=2021msql'
  db_name <- 'config'
  db_host <- 'hpgc-mysql-dev.crhjphok688d.us-east-1.rds.amazonaws.com'
  db_port <- 3306
  mydb <- dbConnect(MySQL(), user = db_user, password = db_password,
                    dbname = db_name, host = db_host, port = db_port)
}


DisConnectDB<-function(connObject)
{
  dbDisconnect(connObject)
}


GetStudyName<-function(AppName)
{
  connObject <- ConnectDB(AppName)
  query <- ("select schema_name as StudyName from information_schema.schemata where schema_name not in ('sys', 'mysql', 'config') order by schema_name;")
  rs <- dbSendQuery(connObject, query)
  
  df <- fetch(rs)
  DisConnectDB(connObject)
  return(df)
}


GetData<-function(AppName, sSQL)
{
  connObject <- ConnectDB(AppName)
  rs <- dbSendQuery(connObject, sSQL)
  df <- fetch(rs)
  DisConnectDB(connObject)
  return(df)
}
