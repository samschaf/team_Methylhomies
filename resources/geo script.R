library(GEOmetadb) 
getSQLiteFile() 
con <- dbConnect(SQLite(), "GEOmetadb.sqlite") 

dbListFields(con, "gsm")

x<-dbGetQuery(con, "select title,description,series_id,gsm,source_name_ch1,characteristics_ch1 from gsm where gpl='GPL13534'")

sampl<-x