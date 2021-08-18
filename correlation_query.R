
##################################functions to query and format correlation pairs###############################

target <- function(table1, corr = .6, drugname = NULL, gene = NULL){
  library(dplyr)
  table1 <- subset(table1, select=which(!duplicated(names(table1))))
  high_correlation <- table1 %>% select_if(~any(. > corr | . < -corr))
  if(is.null(drugname) & is.null(gene)){
    return(targethelp(table1, high_correlation, corr))
  }else if(!is.null(drugname) & is.null(gene)){
    return(targethelp2(table1, high_correlation, corr, drugname))
  }else if(is.null(drugname) & !is.null(gene)){
    return(targethelp3(table1, high_correlation, corr, gene))
  }else{
    #given all parameters simply generate a list of correlations(given repeated names)
    rows = which(high_correlation[,1]==gene)
    cols = which(colnames(high_correlation)==drugname)
    corlist <- list()
    for(i in rows){
      for(j in cols){
        correlation <- high_correlation[i, ..j]
        colnum <- which(colnames(table1)==drugname)
        corlist[[paste(i, "," , colnum, sep="")]] <- correlation
      }
    }
    return(corlist)
  }
  
}


#given 0 or one parameters
targethelp <- function(table1, table2, corr = .6){
  corlist <- data.table(matrix(ncol = 5, nrow = 0))
  names(corlist) = c("column", "drug", "row", "gene", "correlation")
  for(drg in 2:ncol(table2)){
    high <- which(abs(table2[, ..drg]) > corr)
    for(g in high){
      gene <- as.character(table2[g, 1])
      colnum <- which(colnames(table1)==names(table2[, ..drg]))
      entry <- data.frame(column = colnum, drug = names(table2[, ..drg]),
                          row = g, gene = gene, correlation = as.numeric(table2[g, ..drg]))
      corlist <- rbind(corlist, entry, fill = TRUE)
    }
  }
  return(corlist)
}

#given the drug name
targethelp2 <- function(table1, table2, corr, drugname = NULL){
  corlist <- data.table(matrix(ncol = 5, nrow = 0))
  names(corlist) = c("column", "drug", "row", "gene", "correlation")
  cols <- which(colnames(table2)==drugname)
  for(i in cols){
    high <- which(abs(table2[, ..i]) > corr)
    for(g in high){
      gene <- as.character(table2[g, 1])
      colnum <- which(colnames(table1)==drugname)
      entry <- data.frame(column = colnum, drug = drugname, 
                          row = g, gene = gene, correlation = as.numeric(table2[g, ..i]))
      corlist <- rbind(corlist, entry, fill = TRUE)
    }
  }
  return(corlist)
}

#given the gene name
targethelp3 <- function(table1, table2, corr, gene = NULL){
  corlist <- data.table(matrix(ncol = 5, nrow = 0))
  names(corlist) = c("column", "drug", "row", "gene", "correlation")
  rows <- which(table2[, 1]==gene)
  for(i in rows){
    row <- table2[i]
    for(val in 2:ncol(row)){
      if(abs(table2[i, ..val]) > corr){
        colnum <- which(colnames(table1)==names(table2[, ..val]))
        drugname <- as.character(names(table2[, ..val]))
        #drug <- paste("(", colnum,") ", names(table2[, ..val]), sep = "")
        #gene <- paste("(", i, ") ", gene, "--->", sep = "")
        #pair <- as.character(c(drug, gene, table2[i, ..val]))
        entry <- data.frame(column = colnum, drug = drugname, 
                          row = i, gene = gene, correlation = as.numeric(table2[i, ..val]))
        corlist <- rbind(corlist, entry, fill = TRUE)
      }
    }
  }
  return(corlist)
}


#query a table based on target genes for each drug in the treatment table
targetgene <- function(cortable, treatment){
  corlist <- data.table(matrix(ncol = 5, nrow = 0))
  names(corlist) = c("column", "drug", "row", "gene", "correlation")
  for(row in 1:nrow(treatment)){
    drugname <- treatment[row,2] 
    tgt <- treatment[row,6]
    if(!is.null(tgt)){
      tgt <- str_split(treatment[row, 6], ", ")
      tgt <- tgt[[1]]
      for(genename in tgt){
        rownum <- which(cortable[, 1]==genename)
        colnum <- which(colnames(cortable)==drugname)
        for(rn in rownum){
          for(cn in colnum){
            corr <- cortable[rn, ..cn]
            entry <- data.frame(column = cn, drug = drugname, 
                            row = rn, gene = genename, correlation = as.numeric(corr))
            corlist <- rbind(corlist, entry, fill = TRUE)
          }
        }
      }
    }
  }
  return(corlist)
}


#version given to different correlation cutoff values
target2 <- function(table1, corr1 = -0.6, corr2 = 0.6, drugname = NULL, gene = NULL){
  library(dplyr)
  table1 <- subset(table1, select=which(!duplicated(names(table1))))
  high_correlation <- table1 %>% select_if(~any(. > corr2 | . < corr1))
  if(is.null(drugname) & is.null(gene)){
    return(target2help(table1, high_correlation, corr1, corr2))
  }else if(!is.null(drugname) & is.null(gene)){
    return(target2help2(table1, high_correlation, corr1, corr2, drugname))
  }else if(is.null(drugname) & !is.null(gene)){
    return(target2help3(table1, high_correlation, corr, gene))
  }else{
    #given all parameters simply generate a list of correlations(given repeated names)
    rows = which(high_correlation[,1]==gene)
    cols = which(colnames(high_correlation)==drugname)
    corlist <- list()
    for(i in rows){
      for(j in cols){
        correlation <- high_correlation[i, ..j]
        colnum <- which(colnames(table1)==drugname)
        corlist[[paste(i, "," , colnum, sep="")]] <- correlation
      }
    }
    return(corlist)
  }
  
}


#given 0 or one parameters
target2help <- function(table1, table2, corr1, corr2){
  corlist <- data.table(matrix(ncol = 5, nrow = 0))
  names(corlist) = c("column", "drug", "row", "gene", "correlation")
  for(drg in 2:ncol(table2)){
    high <- which(table2[, ..drg] > corr2 | table2[, ..drg] < corr1)
    for(g in high){
      gene <- as.character(table2[g, 1])
      colnum <- which(colnames(table1)==names(table2[, ..drg]))
      entry <- data.frame(column = colnum, drug = names(table2[, ..drg]),
                          row = g, gene = gene, correlation = as.numeric(table2[g, ..drg]))
      corlist <- rbind(corlist, entry, fill = TRUE)
    }
  }
  return(corlist)
}

#given the drug name
target2help2 <- function(table1, table2, corr1, corr2, drugname = NULL){
  corlist <- data.table(matrix(ncol = 5, nrow = 0))
  names(corlist) = c("column", "drug", "row", "gene", "correlation")
  cols <- which(colnames(table2)==drugname)
  for(i in cols){
    high <- which(table2[, ..i] > corr2) | which(table2[, ..drg] < corr1)
    for(g in high){
      gene <- as.character(table2[g, 1])
      colnum <- which(colnames(table1)==drugname)
      entry <- data.frame(column = colnum, drug = drugname, 
                          row = g, gene = gene, correlation = as.numeric(table2[g, ..i]))
      corlist <- rbind(corlist, entry, fill = TRUE)
    }
  }
  return(corlist)
}

#given the gene name
target2help3 <- function(table1, table2, corr1, corr2, gene = NULL){
  corlist <- data.table(matrix(ncol = 5, nrow = 0))
  names(corlist) = c("column", "drug", "row", "gene", "correlation")
  rows <- which(table2[, 1]==gene)
  for(i in rows){
    row <- table2[i]
    for(val in 2:ncol(row)){
      if(table2[i, ..val] > corr2 | table2[i, ..val] < corr1){
        colnum <- which(colnames(table1)==names(table2[, ..val]))
        drugname <- as.character(names(table2[, ..val]))
        #drug <- paste("(", colnum,") ", names(table2[, ..val]), sep = "")
        #gene <- paste("(", i, ") ", gene, "--->", sep = "")
        #pair <- as.character(c(drug, gene, table2[i, ..val]))
        entry <- data.frame(column = colnum, drug = drugname, 
                            row = i, gene = gene, correlation = as.numeric(table2[i, ..val]))
        corlist <- rbind(corlist, entry, fill = TRUE)
      }
    }
  }
  return(corlist)
}


#removes duplicate gene and drug pairs that appear in the table
remove_dupes <- function(targets){
  t1 <- unique( targets[,c("drug", "gene")] ) 
  c1 <- character()
  c2 <- character()
  for(i in 1:nrow(targets)){
    pair <- targets[i,]
    drug <- as.character(pair[,2])
    gene <- as.character(pair[,4])
    druggene <- paste(drug, gene, sep = "")
    c1 <- c(c1, druggene)
  }
  
  for(i in 1:nrow(t1)){
    pair <- t1[i,]
    drug <- as.character(pair[,1])
    gene <- as.character(pair[,2])
    druggene <- paste(drug, gene, sep = "")
    c2 <- c(c2, druggene)
  }
  
  m1 <- match(c2, c1)
  column <- targets$column[m1]
  row <- targets$row[m1]
  correlation <- targets$correlation[m1]
  t1 <- cbind( column, "drug" = t1$drug, row, "gene" = t1$gene, correlation)
  t1 <- data.table(t1)
  return(t1)
}


