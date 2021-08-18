
#########################functions to format the initial raw data##############################


#alters the cell numbers to the actual name of the cell line
cell_name <- function(data, names, name_col = 2){
  update <- data
  for(i in 1:length(row.names(data))){
    row <- which(row.names(names)==row.names(data)[i])
    if(length(row) > 0){
      line_name <- names[row, name_col]
    }else{
      line_name <- row.names(data)[i]
    }
    row.names(update)[i] <- line_name
  }
  return(update)
}



#simplify the name of the gene
gene_name <- function(data){
  for(i in 1:ncol(data)){
    full <- colnames(data)[i]
    index <- 0
    for(j in 0:nchar(full)){
      if(substr(full, j, j)=="."){
        break
      }
      index <- j
    }
    colnames(data)[i] <- substr(full, 0, index)
  }
  return(data)
}


#remove non alphanumeric characters from name
splice <- function(names){
  library(stringr)
  for(i in 1:length(names)){
    names[i] <- str_replace_all(names[i], "[^[:alnum:]]", "")
  }
  return(names)
}


#simplify the names of the drugs
drug_name <- function(data, names, name_col = 2){
  splice1 <- splice(colnames(data))
  splice2 <- splice(row.names(names))
  for(i in 1:length(splice1)){
    full <- splice1[i]
    row <- which(splice2==full)
    if(length(row) > 0){
      col_name <- names[row, name_col]
    }else{
      col_name <- colnames(data)[i]
    }
    colnames(data)[i] <- col_name
  }
  return(data)
}