# Heat shock proteins and molecular chaperones from DESeq2 Results table.
# Patrick Campbell
# März 2021

rm(list=ls())

################################################
##~ User Data Entry ~###########################
################################################

## INSERT YOUR PATH DIRECTORY WHERE DEOBJECT IS STORED HERE ##
# DEObject_name is the name of your saved LFC DEObject,
#
#   i.e:   lfc_Tumordata <- lfcShrink(dds_Tumordata, coef=2, res=res_Tumordata, type = 'apeglm')
#          save(lfc_Tumordata, file=paste(fig_save_path,'lfc_Tumordata',sep=''))
#          DEObject_name would then be: 'lfc_Tumordata'

path_dir <- '/path/to/LFC_Object/'   # enter dir here
DEObject_name <- 'lfc_Tumordata'     # enter name here

################################################
##~ Dependencies and Functions ~################
################################################
load_libraries <- function() {
  library(hash)
  library(DESeq2)
  library(conflicted)
  conflict_prefer("keys", "hash")
}

get_heatshocks <- function(DEObject_name){
  # Gets Heatshock Proteins,
  # DEObject_name is the name of your saved LFC DEObject,
  #   i.e:   lfc_Tumordata <- lfcShrink(dds_Tumordata, coef=2, res=res_Tumordata, type = 'apeglm')
  #          save(lfc_Tumordata, file=paste(fig_save_path,'lfc_Tumordata',sep=''))
  #       DEObject_name would then be: 'lfc_Tumordata'
  
  loadRData <- function(DEObject_name){
    #loads an RData file, and returns it
    # DEObject_name: name of DEObject that was saved in your directory, ie: 'lfc'
    file_name <- paste(path_dir, DEObject_name, sep='')
    load(file_name)
    DEObject <- get(ls()[!(ls() %in% c("file_name","DEObject_name"))])
    return(DEObject)
  }
  
  X <- loadRData(DEObject_name)
  print(X)
  # https://www.genenames.org
  dict <- hash()
  dict[['90kDA']] <- c('HSP90', 'TRAP1')
  dict[['70kDA']] <- c('HSPA', 'HSPH1','HYOU1','HSP70','GRP78','GRP170')  # GRP are aliasis names
  dict[['40kDA']] <- c('DNAJ', 'HSCB', 'SEC63', 'GAK', 'SACS')
  dict[['Small']] <- c('HSB','CRYAA','CRYAB','ODF1')
  dict[['Chaperones']] <- c('BBS10', 'BBS12', 'TCP1', 'CCT', 'CLPB', 'HSPD1', 'HSPE1', 'MKKS')
  
  X$type <- NA
  df_channels <- X[FALSE,]
  grepper <- function(grepee){
    # grepee is the search term
    grepee <- paste('^',grepee,sep='')
    return(X[grep(grepee, rownames(X)),])
  }
  
  for (k in keys(dict)){
    for (v in dict[[k]]){
      grepped <- grepper(v)
      if (nrow(grepped) != 0){
        grepped$type <- k
      }
      df_channels <- data.frame(rbind(as.matrix(df_channels), as.matrix(grepped)))
    }
  }
  
  remove_bad_greps <- df_channels[grep('^GLRX', rownames(df_channels)),] #
  df_channels <- df_channels[!rownames(df_channels) %in% rownames(remove_bad_greps),]#
  # To remove future bad greps, just double the last 2 lines but modified
  
  df_channels <- df_channels[order(df_channels$log2FoldChange),]
  write.table(df_channels, file=paste(path_dir, DEObject_name, 'HEATSHOCKS.csv', sep=''), col.names = NA)
  
  return(df_channels)
}

################################################
##~ Call Functions ~############################
################################################
load_libraries()
Y <- get_heatshocks(DEObject_name)


