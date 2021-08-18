

#################functions to generate plots based on cor tables, not really used in the server############


#function to plot a barplot of the high correlation drugs
drug_corplot <- function(table, corval = "0.5"){
  library(ggplot2)
  library(ggpubr)
  library(grid)
  cl <- setDT(table)[,freq := .N, by = c("drug")]
  cl <- cl[order(freq, decreasing = T),] 
  theme_set(theme_pubr())
  count <- unique(table$drug)
  grob <- grobTree(textGrob(paste("count:", length(count)), x=0.7,  y=0.1, hjust=0,
                    gp=gpar(col="black", fontsize=13, fontface="italic")))
  plot1 <- ggplot(cl, aes(reorder(drug, freq))) +
    geom_bar(fill = "#0073C2FF") +
    theme_pubclean()+coord_flip() + 
    labs(title = paste("drugs that have correlation greater than ", corval) ,y="number of pairs", x = "drugs") +
    annotation_custom(grob) +
    geom_text(stat='count', aes(label=..count..), vjust = 0.3, hjust = -1) +
    theme(plot.title = element_text(hjust = 0.5)) 
  return(plot1)
}


#function to plot a barplot of the high correlation drugs
gene_corplot <- function(table, corval = "0.5"){
  library(ggplot2)
  library(ggpubr)
  library(grid)
  cl <- setDT(table)[,freq := .N, by = c("gene")]
  cl <- cl[order(freq, decreasing = T),]
  theme_set(theme_pubr())
  count <- unique(table$gene)
  grob <- grobTree(textGrob(paste("count:", length(count)), x=0.7,  y=0.1, hjust=0,
                            gp=gpar(col="black", fontsize=13, fontface="italic")))
  plot1 <- ggplot(cl, aes(reorder(gene, freq))) +
    geom_bar(fill = "#2c9e21") +
    theme_pubclean()+coord_flip() + 
    labs(title = paste("genes that have correlation greater than ", corval) ,y="number of pairs", x = "genes") +
    annotation_custom(grob) +
    geom_text(stat='count', aes(label=..count..), vjust = 0.3, hjust = -1) +
    theme(plot.title = element_text(hjust = 0.5)) 
  return(plot1)
}

#generate a histogram of a random sample of specified length. 
cor_sample <- function(table, targets){
  library(ggplot2)
  val1 <- c(1:nrow(table))
  # val1 <- c(val1, targets$row)
  # var.names <- names(table(val1))[table(val1) == 1]
  # val1 <- match(var.names, val1)
  # 
  val2 <- c(1:ncol(table))
  # val2 <- c(val2, targets$column)
  # var.names2 <- names(table(val2))[table(val2) == 1]
  # val2 <- match(var.names2, val2)
  
  sample <- numeric()
  for(i in 1:nrow(targets)){
    rows <- sample(val1, 1, replace = TRUE)
    cols <- sample(val2, 1, replace = TRUE)
    r <- rows[1]
    c <- cols[1]
    val <- as.numeric(table[r, ..c])
    #since there are nas in the prepared table, when we find NA that means it is a target pair, which we do not want to add to the sample.
    while (is.na(val)) {
      rows <- sample(val1, 1, replace = TRUE)
      cols <- sample(val2, 1, replace = TRUE)
      r <- rows[1]
      c <- cols[1]
      val <- as.numeric(table[r, ..c])
    }
    sample <- c(sample, val)
  }
  plot <- ggplot()+aes(sample)+geom_histogram() + xlim(-0.7,0.7)
  return(plot)
}


#plots the distributions of both sample and target
combine_plots <- function(table1, plot1){
  dat <- data.frame(correlation = c(as.numeric(table1$correlation), 
                           as.numeric(plot1$plot_env$sample)),Data = rep(c("Target", "Random"),each = nrow(table1)))
  sum1 <- summary(table1$correlation)
  sum2 <- summary(plot1$plot_env$sample)
  grob <- grobTree(textGrob(paste("mean of target cors:", round(sum1[4], 10)), x=0.6,  y=0.9, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))
  grob2 <- grobTree(textGrob(paste("mean of random cors:", round(sum2[4], 10)), x=0.6,  y=0.8, hjust=0,
                            gp=gpar(col="green", fontsize=13, fontface="italic")))
  plot <- ggplot(dat,aes(x=xx)) + 
    geom_histogram(data=subset(dat,Data == 'Target'), bins = 40, fill = "red", color = "red", alpha = 0.2, size = 2) +
    geom_histogram(data=subset(dat,Data == 'Random'), bins = 40, fill = "green", color = "green",  alpha = 0.2, size = 1) +
    theme_gray()+
    labs(title = "Target correlation dist(red) vs. Random correlation dist(green)", x = "correlation", y = "frequency") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotation_custom(grob) +
    annotation_custom(grob2)
  return(plot)
}

#density plot
combine_plots2 <- function(table, plot){
  dat <- data.frame(correlation = c(as.numeric(table$correlation), 
                                    as.numeric(plot$plot_env$sample)),Data = rep(c("Target", "Random"),each = nrow(table)))
  plot <- ggplot(dat, aes(x=correlation, color=Data, fill = Data)) + 
    geom_density(alpha = 0.4, size=2) +
    theme_gray() + scale_fill_brewer(palette = "Pastel2") + 
    scale_alpha(range = c(0.4, 0.8)) +
    labs(title = "Target correlation dist vs. Random correlation dist", x = "Correlation", y = "Density") +
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}

#absolute value density plot
combine_plots3 <- function(table, plot){
  dat <- data.frame(correlation = c(as.numeric(abs(table$correlation)), 
                                    abs(as.numeric(plot$plot_env$sample))),Data = rep(c("Target", "Random"),each = nrow(table)))
  plot <- ggplot(dat, aes(x=correlation, color=Data, fill = Data)) + 
    geom_density(alpha = 0.4, size=2) +
    theme_gray() + scale_fill_brewer(palette = "Pastel2") + 
    scale_alpha(range = c(0.4, 0.8)) +
    labs(title = "Target correlation dist vs. Random correlation dist (absolute value)", x = "Correlation", y = "Density") +
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}

#function that sets all values in the correlation matrix that have a corresponding value in the target table to
#nas to be used by the sampling function
difference <- function(table, targets){
  columns <- as.numeric(targets$column)
  rows <- as.numeric(targets$row)
  for(i in 1:length(columns)){
    column <- columns[i]
    row <- rows[i]
    table[row, column] <- NA
  }
  return(table)
}


#create a adjacency matrix
adjacency <- function(pairs){
  drugs <- unique(pairs$drug)
  genes <- unique(pairs$gene)
  factors <- c(drugs, genes)
  adj_matrix <- matrix(nrow = length(factors), ncol = length(factors))
  #sadj_matrix <- data.table(adj_matrix)
  colnames(adj_matrix) <- factors
  adj_matrix[is.na(adj_matrix)] <- 0
  rownames(adj_matrix) <- factors
  for(i in 1:nrow(pairs)){
    drug <- as.character(pairs[i, 2])
    gene <- as.character(pairs[i, 4])
    column <- as.numeric(which(factors==drug))
    row <- as.numeric(which(factors==gene))
    adj_matrix[row, column] <- adj_matrix[row, column] + 1
    adj_matrix[column, row] <- adj_matrix[column, row] + 1
  }
  return(adj_matrix)
}

#graph attempt
dgnet <- function(matrix, valuable_targets){
  library(statnet)
  set.seed(12345)
  num_nodes <- 10
  diag(adj_matrix) <- 0
  
  net <- as.network(x = matrix, # the network object
                    directed = TRUE, # specify whether the network is directed
                    loops = FALSE, # do we allow self ties (should not allow them)
                    matrix.type = "adjacency" # the type of input
  )
  
  drugs <- unique(valuable_targets$drug)
  genes <- unique(valuable_targets$gene)
  factors <- (c(drugs, genes))
  network.vertex.names(net) <- factors
  # Create the variable
  type <- c(rep("Drug",length(drugs)),rep("Gene",length(genes)))
  # Take a look at our variable
  # Add it to the network object
  set.vertex.attribute(net, # the name of the network object
                       "Type", # the name we want to reference the variable by in that object
                       type # the value we are giving that variable
  ) 
  
  # age <- round(rnorm(num_nodes,20,3))
  # set.vertex.attribute(net,"Age",age)
  # summary.network(net, # the network we want to look at
  #                 print.adj = FALSE # if TRUE then this will print out the whole adjacency matrix.
  # )
  
  node_colors <- rep("",length(factors))
  for(i in 1:length(factors)){
    if(get.node.attr(net,"Type")[i] == "Drug"){
      node_colors[i] <- "lightblue"
    }else{
      node_colors[i] <- "maroon"
    }
  }
  plot.network(net, # our network object
               vertex.col = node_colors, # color nodes by gender
               vertex.cex = 1, # size nodes by their age
               displaylabels = T, # show the node names
               label.pos = 5 # display the names directly over nodes
  )
  return(net)
}
