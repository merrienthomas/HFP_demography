setwd(getwd())

#### LIBRARY ####

library(stringi)
library(tidyr)
library(phytools)
library(ape)
library(phangorn)
library(geiger)
library(dplyr)
library(car)
library(FactoMineR)
library(mgcv)
library(MCMCglmm)
library(geosphere)
library(DHARMa)
library(ade4)
library(geometry)
library(lme4)
library(nlme)
library(alphahull)
library(alphashape3d)
library(ggplot2)
library(ggthemes)
library(ggfortify)
library(data.table)
library(Rcompadre)


#### FUNCTIONS ####
at.node <-
  function(phylogeny, location.node, tip.label)
  {
    phylo <- reorder(phylogeny)
    if (!is.numeric(location.node))
      location.node <- which(phylo$node.label == location.node) + length(phylo$tip.label)
    a <- location.node - length(phylo$tip.label)
    EL <- branching.times(phylo)[a]
    a0 <- a + length(phylo$tip.label)
    a1 <- which(phylo$edge[,1] == a0)[1]
    aa <- length(which(phylo$edge[1 : (a1 - 1), 2] <= length(phylo$tip.label)))
    eG0 <- matrix(c(a0 + 1, aa + 1), nrow = 1)
    eG <- matrix(phylo$edge[a1 : dim(phylo$edge)[1], ], ncol = 2)
    eG[, 1] <- eG[, 1] + 1
    s <- which(eG[, 2] > aa)
    eG[, 2][s] <- eG[, 2][s] + 1
    eG <- rbind(eG0, eG)    
    eL <- c(EL, phylo$edge.length[a1 : length(phylo$edge.length)])
    tL <- c(tip.label, phylo$tip.label[(aa + 1) : length(phylo$tip.label)])
    if (a1 > 1) 
    { 
      eGn <- matrix(phylo$edge[1 : (a1 - 1), ], ncol = 2)
      eGn[, 1] <- eGn[, 1] + 1
      s <- which(eGn[, 2] > aa)
      eGn[, 2][s] <- eGn[, 2][s] + 1
      eG <- rbind(eGn, eG) 
      eL <- c(phylo$edge.length[1 : (a1 - 1)], eL)
    }
    if (aa > 0)
    {
      tL <- c(phylo$tip.label[1 : aa], tL)
    }
    phylo$edge <- eG
    phylo$tip.label <- tL
    phylo$edge.length <- eL
    return(phylo)
  }

correcting_data <- function(col_trait, list_trait, value_checking, df, size_correction, spatial_correction, log_transform, inverse_transform, tree){
  
  index <- grep(col_trait,list_trait)
  
  #Cleaning the data by removing the outliers
  
  if(value_checking[index]==1){
    newlist <- Removing_outliers_and_false_values(df,col_trait)
  } else {
    newlist <- Removing_outliers(df,col_trait)
  }
  dfbis <- newlist[[1]]
  storage <- newlist[[2]]
  
  treebis<-matching_tree_data(storage,tree)
  
  #Changing node labels
  treebis$node.label<-unique(treebis$edge[,1])
  
  #Transformation of the data
  
  if (length(grep(TRUE, log_transform == col_trait)) > 0){
    dfbis[, col_trait] <- log(dfbis[, col_trait] + 0.0000001)
  } else {
    if (length(grep(TRUE, inverse_transform == col_trait)) > 0){
      dfbis[,col_trait] <- (1/(dfbis[, col_trait]) + 0.0000001)
    } else {
      dfbis[,col_trait] <- dfbis[,col_trait]
    }
  }
  
  
  #Correction for size and spatial autocorrelation
  
  if (size_correction[index] == 1) {
    
    if (spatial_correction[index] == 1) {
      
      dfbis[,col_trait] <- gam(dfbis[,col_trait] ~ dfbis$size + te(dfbis$Lat,dfbis$Lon))$residual
      
    } else {
      
      dfbis[,col_trait] <- lm(dfbis[,col_trait] ~ dfbis$size)$residual
      
    }
    
  } else {
    
    if (spatial_correction[index] == 1) {
      
      dfbis[,col_trait] <- gam(dfbis[,col_trait] ~ 1 + te(dfbis$Lat,dfbis$Lon))$residual
      
    } else {
      
      dfbis[,col_trait] <- dfbis[,col_trait]
      
    }
  }
  
  return(list(dfbis,treebis))
  
}

#function from V.phylomaker package

ext.node <-
  function(phylogeny, location.tip, tip.label, node.label = NULL, position = 0.5)
  {
    phylo <- reorder(phylogeny)
    if (!is.numeric(location.tip))
    {
      location.tip <- which(phylo$tip.label == location.tip)
    }
    a <- location.tip
    a1 <- which(phylo$edge[, 2] == a)
    h <- phylo$edge[a1, 1]         
    if (is.null(phylo$node.label))
    {
      if (!is.null(node.label))
      {
        nL <- rep(NA, phylo$Nnode + 1)
        n <- h - length(phylo$tip.label)
        nL[n + 1] <- node.label
      }
      if (is.null(node.label))
      {
        nL <- NULL
      }
    }
    if (!is.null(phylo$node.label))
    {
      if (!is.null(node.label)) 
      {
        n <- h - length(phylo$tip.label)
        nL <- c(phylo$node.label[1 : n], node.label)
        if (n < phylo$Nnode) 
        {
          nL <- c(nL, phylo$node.label[(n + 1) : phylo$Nnode]) 
        }
      }
      if (is.null(node.label)) 
      {
        n <- h - length(phylo$tip.label)
        nL <- c(phylo$node.label[1 : n], NA)
        if (n < phylo$Nnode) 
        {
          nL <- c(nL,phylo$node.label[(n + 1) : phylo$Nnode]) 
        }
      }
    }
    eG0 <- matrix(c(h + 1, h + 2, h + 2, a, h + 2, a + 1), nrow = 3, byrow = T)
    eG <- matrix(phylo$edge[1 : (a1 - 1), ], ncol = 2)
    eG[, 1] <- eG[, 1] + 1
    s <- which(eG[, 1] > (h + 1))
    eG[, 1][s] <- eG[, 1][s] + 1
    s <- which(eG[, 2] > a)
    eG[, 2][s] <- eG[, 2][s] + 1
    s <- which(eG[, 2] > (h + 1))
    eG[, 2][s] <- eG[, 2][s] + 1
    eG <- rbind(eG, eG0)
    tL <- c(phylo$tip.label[1 : a], tip.label)
    eL <- c(phylo$edge.length[1 : (a1 - 1)],  phylo$edge.length[a1] * (1 - position), phylo$edge.length[a1] * position, phylo$edge.length[a1] * position)
    if (a < length(phylo$tip.label)) 
    { 
      eGn <- matrix(phylo$edge[(a1 + 1) : (dim(phylo$edge)[1]), ], ncol = 2)
      eGn[,1] <- eGn[,1] + 1
      s <- which(eGn[,1] > (h + 1))
      eGn[, 1][s] <- eGn[, 1][s] + 1
      s <- which(eGn[, 2] > a)
      eGn[, 2][s] <- eGn[, 2][s] + 1
      s <- which(eGn[, 2] > (h + 1))
      eGn[, 2][s] <- eGn[, 2][s] + 1
      eG <- rbind(eG, eGn) 
      tL <- c(tL, phylo$tip.label[(a + 1) : length(phylo$tip.label)])
      eL <- c(eL,phylo$edge.length[(a1 + 1) : length(phylo$edge.length)])
    }
    phylo$edge <- eG
    phylo$tip.label <- tL
    phylo$edge.length <- eL
    phylo$Nnode <- phylo$Nnode + 1
    phylo$node.label <- nL
    return(phylo)
  }


#Code provided by Liam Revell

force.ultrametric<-function(tree,method=c("nnls","extend")){ 
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

matching_tree_data<-function(storage,tree){
  if ((length(storage)>0)==TRUE){
    for (i in 1:length(storage)){
      tree<-drop.tip(tree, tip = tree$tip.label[which(tree$tip.label == storage[i])])
    }
  }
  tree$node.label<-unique(tree$edge[,1])
  return(tree)
}

Moran.Index<-function(df,trait){
  
  ## Create a matrix of distance between point
  
  Distance <- matrix(nrow = length(df$X), ncol = length(df$X))
  Long <- df$Lon
  Lati <- df$Lat
  
  coords <- as.matrix(cbind(Long, Lati))
  
  # Generate the distance matrix
  coords.dists <- as.matrix(distm(coords))
  
  # Take the inverse of the distance matrix
  coords.dists.inv <- 1/coords.dists
  
  #Replace infinite value by 0
  coords.dists.inv[coords.dists.inv=="Inf"] <- 0
  
  return(Moran.I(trait,coords.dists.inv,na.rm = TRUE))
  
}

#function from V.phylomaker package

phylo.maker <-
  
  function (sp.list, tree = GBOTB.extended, nodes = nodes.info.1, output.sp.list = TRUE, output.tree = FALSE, scenarios = "S3", r = 1)
    
  {
    
    # options(scipen = 999)
    
    treeX <- tree
    
    if (is.null(tree$node.labels))
      
      tree$node.label <- rep("", tree$Nnode)
    
    dimnames(sp.list)[[2]][1:3] <- c("species", "genus", "family")
    
    sp.list[sapply(sp.list, is.factor)] <- lapply(sp.list[sapply(sp.list,
                                                                 
                                                                 is.factor)], as.character)
    
    if (any(duplicated(sp.list$species))) {
      
      print("Duplicated species detected and removed.")
      
      print(sp.list$species[duplicated(sp.list$species)])
      
    }
    
    sp.list <- sp.list[!duplicated(sp.list$species), ]
    
    sp.list.original <- sp.list
    
    sp.list$species <- gsub(" ", "_", sp.list$species)
    
    sp.list$species <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$species,
                            
                            perl = TRUE)
    
    sp.list$genus <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$genus,
                          
                          perl = TRUE)
    
    sp.list$family <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$family,
                           
                           perl = TRUE)
    
    rnN <- data.frame(node.label = paste("N", 1:length(tree$node.label),
                                         
                                         sep = ""), oriN = tree$node.label, stringsAsFactors = FALSE)
    
    nodes[, c("level", "family", "genus", "rn", "bn", "taxa")] <- lapply(nodes[,
                                                                               
                                                                               c("level", "family", "genus", "rn", "bn", "taxa")], as.character)
    
    tree$node.label <- paste("N", 1 : length(tree$node.label),
                             
                             sep = "")
    
    kk <- c()
    
    for (i in 1 : length(tree$tip.label)) {
      
      kk <- c(kk, substring(tree$tip.label[i], 1, gregexpr("_",
                                                           
                                                           tree$tip.label[i])[[1]][1] - 1))
      
    }
    
    m <- data.frame(num = 1 : length(kk), genus = kk, species = tree$tip.label)
    
    m <- merge(m, nodes[, c("genus", "family")])
    
    mX <- m
    
    m <- m[, c("genus", "family")]
    
    m <- m[!duplicated(m$genus), ]
    
    dimnames(m)[[2]][2] <- "family_in_tree"
    
    m<-m[, c("genus", "family_in_tree")]
    
    m0 <- sp.list[!duplicated(sp.list$genus), c("genus", "family")]
    
    dimnames(m0)[[2]][2] <- "family_in_sp.list"
    
    mm<-merge(m0, m)
    
    g<-mm[which(is.na(match(paste(mm$genus, mm$family_in_sp.list, sep = "_"),
                            
                            paste(mm$genus, mm$family_in_tree, sep = "_")))), ]
    
    if (dim(g)[1] > 0)
      
    {
      
      print("Taxonomic classification not consistent between sp.list and tree.")
      
      print(g)
      
    }
    
    add.tip <- sp.list[which(is.na(match(sp.list$species, tree$tip.label))), ]
    
    status <- rep("prune", dim(sp.list)[1])
    
    status[which(is.na(match(sp.list$species, tree$tip.label)))] <- "bind"
    
    if (dim(add.tip)[1] == 0 & length(na.omit(match(sp.list$species, tree$tip.label))) == 0)
      
      stop("Incorrect format of species list.")
    
    if (length(setdiff(sp.list$species, treeX$tip.label)) == 0 & length(na.omit(match(sp.list$species, tree$tip.label))) > 0)
      
    {
      
      print("All species in sp.list are present on tree.")
      
      splis <- sp.list.original
      
      treeX <- drop.tip(treeX, setdiff(treeX$tip.label, sp.list$species))
      
      splis$status <- "prune"
      
      phyloX <- list(scenario.1 = NULL, scenario.2 = NULL, scenario.3 = NULL, species.list = splis)
      
      if ("S1" %in% scenarios) {phyloX$scenario.1 <- treeX}
      
      if ("S2" %in% scenarios) {phyloX$scenario.2 <- treeX}
      
      if ("S3" %in% scenarios) {phyloX$scenario.3 <- treeX}
      
      phyloX[sapply(phyloX, is.null)] <- NULL
      
      return(phyloX)
      
      stop()
      
    }
    
    add.tip$sort <- ""
    
    add.tip$sort[which(!is.na(match(add.tip$genus, nodes[nodes$level ==
                                                           
                                                           "G", ]$genus)))] <- "G1"
    
    add.tip$sort[which(is.na(match(add.tip$genus, nodes[nodes$level ==
                                                          
                                                          "G", ]$genus)) & !is.na(match(add.tip$family, nodes[nodes$level ==
                                                                                                                
                                                                                                                "F", ]$family)))] <- "F1"
    
    add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort ==
                                                            
                                                            "F1", ]$genus)] <- "F2"
    
    a <- which(add.tip$sort == "")
    
    if (length(a) > 0) {
      
      print(paste("Note:", length(a), "taxa fail to be binded to the tree,",
                  
                  sep = " "))
      
      print(add.tip$species[a])
      
      status[match(add.tip$species[a], sp.list$species)] <- "fail to bind"
      
    }
    
    sp.list.original$status <- status
    
    if ("S1" %in% scenarios) {
      
      t1 <- tree
      
      rnN1 <- rnN
      
      nG <- nodes[nodes$level == "G", ]
      
      nF <- nodes[nodes$level == "F", ]
      
      data <- add.tip[add.tip$sort == "F1" | add.tip$sort ==
                        
                        "F2", ]
      
      if (dim(data)[1] > 0) {
        
        for (i in 1:dim(data)[1]) {
          
          num <- nF$bn[match(data$family[i], nF$family)]
          
          t1 <- at.node(t1, location.node = num, tip.label = data$species[i])
          
        }
        
      }
      
      data <- add.tip[add.tip$sort == "G1", ]
      
      if (dim(data)[1] > 0) {
        
        for (i in 1:dim(data)[1]) {
          
          num <- nG$bn[match(data$genus[i], nG$genus)]
          
          t1 <- at.node(t1, location.node = num, tip.label = data$species[i])
          
        }
        
      }
      
      t1$edge.length <- as.numeric(t1$edge.length)
      
      tree1 <- t1
      
      tree1$node.label <- rnN1$oriN[match(tree1$node.label, rnN1$node.label)]
      
      toDrop <- setdiff(1:length(t1$tip.label), which(!is.na(match(t1$tip.label,
                                                                   
                                                                   sp.list$species))))
      
      t1 <- drop.tip(t1, tip = toDrop)
      
      Re <- which(!is.na(match(t1$node.label, rnN1$node.label)))
      
      noRe <- which(is.na(match(t1$node.label, rnN1$node.label)))
      
      t1$node.label[Re] <- rnN1$oriN[match(t1$node.label, rnN1$node.label)[Re]]
      
      t1$node.label[noRe] <- ""
      
    }
    
    else {
      
      t1 <- NULL
      
      tree1 <- NULL
      
    }
    
    if ("S2" %in% scenarios) {
      
      t2r <- replicate(r, list())
      
      names(t2r) <- paste("run", 1 : r, sep = ".")
      
      tree2r <- replicate(r, list())
      
      names(tree2r) <- paste("run", 1 : r, sep = ".")
      
      for (o in 1 : r)
        
      {
        
        t2 <- tree
        
        rnN2 <- rnN
        
        nG <- nodes[nodes$level == "G", ]
        
        nF <- nodes[nodes$level == "F", ]
        
        data <- add.tip[add.tip$sort == "F1", ]
        
        if (dim(data)[1] > 0) {
          
          for (i in 1:dim(data)[1]) {
            
            n <- match(data$family[i], nF$family)
            
            g <- nF$gen.n[n]
            
            s <- nF$sp.n[n]
            
            if (g == 1 & s == 1) {
              
              num <- match(nF$taxa[n], t2$tip.label)
              
              nlabel <- paste("N", t2$Nnode + 1, sep = "")
              
              t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                             
                             node.label = nlabel, position = 2/3)
              
              nF$gen.n[n] <- g + 1
              
              nF$sp.n[n] <- s + 1
              
              x <- which(t2$node.label == nlabel)
              
              xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel,
                                                        
                                                        oriN = ""))
              
              xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
              
              rnN2 <- xx
              
              nF$bn[n] <- nlabel
              
            }
            
            else {
              
              num <- sample(nG$bn[which(nG$family %in% data$family[i])],
                            
                            1)
              
              t2 <- at.node(t2, location.node = num, tip.label = data$species[i])
              
              nF$gen.n[n] <- g + 1
              
              nF$sp.n[n] <- s + 1
              
            }
            
          }
          
        }
        
        data <- add.tip[add.tip$sort == "F2", ]
        
        if (dim(data)[1] > 0) {
          
          for (i in 1:dim(data)[1]) {
            
            n <- grep(paste(data$genus[i], "_", sep = ""),
                      
                      t2$tip.label)
            
            nlabel <- paste("N", t2$Nnode + 1, sep = "")
            
            if (length(n) == 1) {
              
              num <- n
              
              t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                             
                             node.label = nlabel, position = 1/2)
              
              x <- which(t2$node.label == nlabel)
              
              xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel,
                                                        
                                                        oriN = ""))
              
              xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
              
              rnN2 <- xx
              
            }
            
            if (length(n) > 1) {
              
              num <- t2$edge[which(t2$edge[, 2] == n[1]),
                             
                             1]
              
              t2 <- at.node(t2, location.node = num, tip.label = data$species[i])
              
            }
            
          }
          
        }
        
        data <- add.tip[add.tip$sort == "G1", ]
        
        if (dim(data)[1] > 0) {
          
          for (i in 1:dim(data)[1]) {
            
            n0 <- match(data$genus[i], nG$genus)
            
            n <- nG$sp.n[n0]
            
            nlabel <- paste("N", t2$Nnode + 1, sep = "")
            
            if (n == 1) {
              
              num <- t2$tip.label[match(nG$taxa[n0], t2$tip.label)]
              
              t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                             
                             node.label = nlabel, position = 1/2)
              
              x <- which(t2$node.label == nlabel)
              
              xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel,
                                                        
                                                        oriN = ""))
              
              xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
              
              rnN2 <- xx
              
              nG$sp.n[n0] <- nG$sp.n[n0] + 1
              
            }
            
            if (n > 1) {
              
              num <- which(t2$node.label == nG$bn[n0]) +
                
                length(t2$tip.label)
              
              num1 <- which(t2$edge[, 1] %in% num)
              
              part1 <- t2$edge[1:min(num1), ]
              
              n1 <- max(which(part1[, 1] < num), 0) + 1
              
              part2 <- t2$edge[max(num1):dim(t2$edge)[1],
                               
                               ]
              
              n2 <- min(which(part2[, 1] < num), dim(part2)[1] +
                          
                          1) + max(num1) - 2
              
              sect <- t2$edge[n1:n2, ]
              
              sect <- sort(unique(c(sect[, 1], sect[, 2])))
              
              sect <- sect[which(sect > length(t2$tip.label))]
              
              num2 <- sect[sample(1:length(sect), 1)]
              
              t2 <- at.node(t2, location.node = num2, tip.label = data$species[i])
              
            }
            
          }
          
        }
        
        t2$edge.length <- as.numeric(t2$edge.length)
        
        tree2 <- t2
        
        tree2$node.label <- rnN2$oriN[match(tree2$node.label, rnN2$node.label)]
        
        toDrop <- setdiff(1:length(t2$tip.label), which(!is.na(match(t2$tip.label,
                                                                     
                                                                     sp.list$species))))
        
        t2 <- drop.tip(t2, tip = toDrop)
        
        Re <- which(!is.na(match(t2$node.label, rnN2$node.label)))
        
        noRe <- which(is.na(match(t2$node.label, rnN2$node.label)))
        
        t2$node.label[Re] <- rnN2$oriN[match(t2$node.label, rnN2$node.label)[Re]]
        
        t2$node.label[noRe] <- ""
        
        t2r[[o]] <- t2
        
        tree2r[[o]] <- tree2
        
      }
      
    }
    
    else {
      
      t2r <- NULL
      
      tree2r <- NULL
      
    }
    
    if ("S3" %in% scenarios) {
      
      t3 <- tree
      
      rnN3 <- rnN
      
      nG <- nodes[nodes$level == "G", ]
      
      nF <- nodes[nodes$level == "F", ]
      
      data <- add.tip[add.tip$sort == "F1", ]
      
      if (dim(data)[1] > 0) {
        
        for (i in 1:dim(data)[1]) {
          
          n <- match(data$family[i], nF$family)
          
          g <- nF$gen.n[n]
          
          s <- nF$sp.n[n]
          
          if (g == 1 & s == 1) {
            
            num <- match(nF$taxa[n], t3$tip.label)
            
            nlabel <- paste("N", t3$Nnode + 1, sep = "")
            
            t3 <- ext.node(t3, location.tip = num, tip.label = data$species[i],
                           
                           node.label = nlabel, position = 2/3)
            
            nF$gen.n[n] <- g + 1
            
            nF$sp.n[n] <- s + 1
            
            x <- which(t3$node.label == nlabel)
            
            xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel,
                                                      
                                                      oriN = ""))
            
            xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
            
            rnN3 <- xx
            
            nF$bn[n] <- nlabel
            
          }
          
          if (g == 1 & s > 1) {
            
            nlabel <- paste("N", t3$Nnode + 1, sep = "")
            
            if ((2/3) * nF$rn.bl[n] <= nF$bn.bl[n]) {
              
              len <- (nF$rn.bl[n] - nF$bn.bl[n])/2
              
            }
            
            if ((2/3) * nF$rn.bl[n] > nF$bn.bl[n]) {
              
              len <- nF$rn.bl[n] * 2/3 - nF$bn.bl[n]
              
            }
            
            port <- len/(nF$rn.bl[n] - nF$bn.bl[n])
            
            t3 <- int.node(t3, location.node = nF$bn[n],
                           
                           tip.label = data$species[i], node.label = nlabel,
                           
                           position = port)
            
            nF$gen.n[n] <- g + 1
            
            nF$sp.n[n] <- s + 1
            
            x <- which(t3$node.label == nlabel)
            
            xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel,
                                                      
                                                      oriN = ""))
            
            xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
            
            rnN3 <- xx
            
            nF$bn[n] <- nlabel
            
          }
          
          if (g > 1) {
            
            t3 <- at.node(t3, location.node = nF$bn[n],
                          
                          tip.label = data$species[i])
            
          }
          
        }
        
      }
      
      data <- add.tip[add.tip$sort == "F2", ]
      
      if (dim(data)[1] > 0) {
        
        for (i in 1:dim(data)[1]) {
          
          n <- grep(paste(data$genus[i], "_", sep = ""),
                    
                    t3$tip.label)
          
          nlabel <- paste("N", t3$Nnode + 1, sep = "")
          
          if (length(n) == 1) {
            
            t3 <- ext.node(t3, location.tip = n, tip.label = data$species[i],
                           
                           node.label = nlabel, position = 1/2)
            
            nG$sp.n[match(data$genus[i], nG$genus)] <- length(n) +
              
              1
            
            x <- which(t3$node.label == nlabel)
            
            xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel,
                                                      
                                                      oriN = ""))
            
            xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
            
            rnN3 <- xx
            
            nG$bn[match(data$genus[i], nG$genus)] <- nlabel
            
          }
          
          if (length(n) > 1) {
            
            num <- ancestor(t3, min(n), max(n))
            
            t3 <- at.node(t3, location.node = num, tip.label = data$species[i])
            
          }
          
        }
        
      }
      
      data <- add.tip[add.tip$sort == "G1", ]
      
      if (dim(data)[1] > 0) {
        
        for (i in 1:dim(data)[1]) {
          
          n0 <- match(data$genus[i], nG$genus)
          
          s <- nG$sp.n[n0]
          
          nlabel <- paste("N", t3$Nnode + 1, sep = "")
          
          if (s == 1) {
            
            num <- t3$tip.label[match(nG$taxa[n0], t3$tip.label)]
            
            t3 <- ext.node(t3, location.tip = num, tip.label = data$species[i],
                           
                           node.label = nlabel, position = 1/2)
            
            nG$sp.n[n0] <- nG$sp.n[n0] + 1
            
            x <- which(t3$node.label == nlabel)
            
            xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel,
                                                      
                                                      oriN = ""))
            
            xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
            
            rnN3 <- xx
            
            nG$bn[n0] <- nlabel
            
          }
          
          if (s > 1) {
            
            t3 <- at.node(t3, location.node = nG$bn[n0],
                          
                          tip.label = data$species[i])
            
          }
          
        }
        
      }
      
      t3$edge.length <- as.numeric(t3$edge.length)
      
      tree3 <- t3
      
      tree3$node.label <- rnN3$oriN[match(tree3$node.label, rnN3$node.label)]
      
      toDrop <- setdiff(1:length(t3$tip.label), which(!is.na(match(t3$tip.label,
                                                                   
                                                                   sp.list$species))))
      
      t3 <- drop.tip(t3, tip = toDrop)
      
      Re <- which(!is.na(match(t3$node.label, rnN3$node.label)))
      
      noRe <- which(is.na(match(t3$node.label, rnN3$node.label)))
      
      t3$node.label[Re] <- rnN3$oriN[match(t3$node.label, rnN3$node.label)[Re]]
      
      t3$node.label[noRe] <- ""
      
    }
    
    else {
      
      t3 <- NULL
      
      tree3 <- t3
      
    }
    
    if (output.sp.list == FALSE) {
      
      sp.list.original <- NULL
      
    }
    
    if (output.tree == FALSE) {
      
      tree1 <- tree2r <- tree3 <- NULL
      
    }
    
    phylo <- list(scenario.1 = t1, scenario.2 = t2r, scenario.3 = t3,
                  
                  species.list = sp.list.original, tree.scenario.1 = tree1,
                  
                  tree.scenario.2 = tree2r, tree.scenario.3 = tree3)
    
    phylo[sapply(phylo, is.null)] <- NULL
    
    return(phylo)
    
  }

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.025, .975), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

Removing_outliers_and_false_values <- function (df,trait_column){
  
  dfbis<-df #Creation of a copy in order not to change the original one
  dfbis[,trait_column]<-remove_outliers(df[,trait_column]) # Removing outliers (Replace the outlier value with a NA)
  storage<-grep("TRUE",is.na(dfbis[,trait_column])) #Saving the location of the outliers and NA in order to remove them in the variance-covariance matrice
  if ((length(grep("TRUE",dfbis[,trait_column]<0))>0)==TRUE){ #Removing value that are not supposed to exist
    storage<-c(storage,grep("TRUE",dfbis[,trait_column]<0)) 
  } 
  dfbis <- dfbis[!(dfbis[,trait_column]<0),]
  dfbis <- drop_na(dfbis,colnames(dfbis)[trait_column]) #Removing all the NA value
  storage <- sort(storage)
  storage <- rownames(df)[storage]
  newlist <- list(dfbis,storage)
  return(newlist)
}

Removing_outliers <- function(df,trait_column){
  
  dfbis <- df #Creation of a copy in order not to change the original one
  dfbis[,trait_column]<-remove_outliers(df[,trait_column]) # Removing outliers (Replace the outlier value with a NA)
  storage <- grep("TRUE",is.na(dfbis[,trait_column])) #Saving the location of the outliers and NA in order to remove them in the variance-covariance matrice
  dfbis <- drop_na(dfbis,colnames(dfbis)[trait_column]) #Removing all the NA value
  storage <- sort(storage)
  storage <- rownames(df)[storage]
  newlist <- list(dfbis,storage)
  return(newlist)
}

Removing_outliers_and_zero <- function (df,trait_column){
  
  dfbis<-df #Creation of a copy in order not to change the original one
  dfbis[,trait_column]<-remove_outliers(df[,trait_column]) # Removing outliers (Replace the outlier value with a NA)
  storage<-grep("TRUE",is.na(dfbis[,trait_column])) #Saving the location of the outliers and NA in order to remove them in the variance-covariance matrice
  if ((length(grep("TRUE",dfbis[,trait_column] <= 0))>0)==TRUE){ #Removing value that are not supposed to exist
    storage<-c(storage,grep("TRUE",dfbis[,trait_column]<0)) 
  } 
  dfbis <- dfbis[!(dfbis[,trait_column] <= 0),]
  dfbis <- drop_na(dfbis,colnames(dfbis)[trait_column]) #Removing all the NA value
  storage <- sort(storage)
  storage <- rownames(df)[storage]
  newlist <- list(dfbis,storage)
  return(newlist)
}

Running_MCMCglmm <- function(df, trait_list, value_checking, list_transform, tree){
  
  results_1 <- data.frame(matrix(0,nrow = length(trait_list), ncol = 7))
  colnames(results_1) <- c("coord1", "p-value", "coord2", "p-value", "size", "p-value", "DIC")
  row.names(results_1) <- colnames(df[,trait_list])
  
  results_2 <- results_1
  
  spatial_size_correction <- data.frame(matrix(0,nrow = length(trait_list), ncol = 2))
  colnames(spatial_size_correction) <- c("size correction", "spatial correction")
  row.names(spatial_size_correction) <- colnames(df[,trait_list])
  
  
  for (k in 1:length(trait_list)){
    
    #How should the data be cleaned
    if (value_checking[k] == 1){
      
      list <- Removing_outliers_and_false_values(df, trait_list[k])
      
      dfbis <- list[[1]]
      storage <- list[[2]]
      
    } else {
      
      if (value_checking[k] == 2){
        
        list <- Removing_outliers(df, trait_list[k])
        
        dfbis <- list[[1]]
        storage <- list[[2]]
        
      } else {
        
        #if (value_checking[k] == 3){
        
        list <- Removing_outliers_and_zero(df, trait_list[k])
        
        dfbis <- list[[1]]
        storage <- list[[2]]
        
        #}
        
      }
    }
    
    #How to transform the data
    if (list_transform[k] == 1){
      
      dfbis$transform <- dfbis[, trait_list[k]]
      
    } else {
      
      if (list_transform[k] == 2) {
        
        dfbis$transform <- log(dfbis[, trait_list[k]] + 0.000001)
        
      } else {
        
        if (list_transform[k] == 3) {
          
          dfbis$transform <- sqrt(dfbis[, trait_list[k]] + 0.000001)
          
        } else {
          
          if (list_transform[k] == 4) {
            
            dfbis$transform <- (1/(dfbis[, trait_list[k]] + 0.000001))
            
          } else {
            
            #if (list_transform[k] == 5) {
            
            dfbis$transform <- (sqrt(dfbis[, trait_list[k]] + 0.0001))
            
            #}
            
          }
          
        }
        
      }
      
    }
    
    
    #Do we correct for size?
    size_lm <- summary(lm(dfbis$transform ~ dfbis$logsize))
    
    if (size_lm$coefficients[8] <= 0.05){
      
      dfbis$res <- size_lm$residuals
      
      spatial_size_correction[k,1] <- "YES"
      
    } else {
      
      spatial_size_correction[k,1] <- "NO"
      
      dfbis$res <- dfbis$transform
      
    }
    
    
    #Do we correct for spatial autocorrelation?
    
    MoranI <- Moran.Index(dfbis, dfbis$res)
    
    if (abs(MoranI$observed) >= 0.1){
      
      if (spatial_size_correction[k,1]=="YES"){
        
        spatial_size_lm <- gam((dfbis$transform) ~ dfbis$logsize + te(dfbis$Lon, dfbis$Lat))
        
        dfbis$res <- spatial_size_lm$residuals
        
      } else {
        
        spatial_size_lm <- gam((dfbis$transform) ~ 1 + te(dfbis$Lon, dfbis$Lat))
        
        dfbis$res <- spatial_size_lm$residuals
        
      }
      
      spatial_size_correction[k,2] <- "YES"
      
    } else {
      
      spatial_size_correction[k,2] <- "NO"
      
    }
    
    
    #MCMCglmm
    treebis <- matching_tree_data(storage, tree)
    dfbis$phylo <- rownames(dfbis)
    
    prior2 <- list(G=list(G1=list(V=3,nu=20), G2=list(V=1,nu=2)),R=list(V=1,nu=0.02))
    
    pglm2 <- MCMCglmm(res ~ coord1dim + coord2dim, random=~phylo + SpeciesAccepted,
                      ginverse = list(phylo = inverseA(treebis,nodes = "TIPS")$Ainv),
                      data=dfbis,
                      prior = prior2,
                      nitt = 25000,
                      burnin = 1000,
                      thin = 10)
    
    if (spatial_size_correction[k,1]=="YES"){
      
      pglm1 <- MCMCglmm(transform ~ coord1dim + coord2dim + logsize, random=~phylo+SpeciesAccepted,
                        ginverse = list(phylo = inverseA(treebis,nodes = "TIPS")$Ainv),
                        data=dfbis,
                        prior = prior2,
                        nitt = 25000,
                        burnin = 1000,
                        thin = 10)
    } else {
      pglm1 <- MCMCglmm(transform ~ coord1dim + coord2dim, random=~phylo+SpeciesAccepted,
                        ginverse = list(phylo = inverseA(treebis,nodes = "TIPS")$Ainv),
                        data=dfbis,
                        prior = prior2,
                        nitt = 25000,
                        burnin = 1000,
                        thin = 10)
    }
    
    
    
    results_1[k, 7] <- summary(pglm1)$DIC
    if (spatial_size_correction[k,1]=="YES"){
      results_1[k, 1:6] <- summary(pglm1)$solution[c(2, 18 ,3, 19,4, 20)]
    } else {
      results_1[k, 1:4] <- summary(pglm1)$solution[c(2, 14, 3, 15)]
    }
    
    
    results_2[k, 7] <- summary(pglm2)$DIC
    results_2[k, 1:4] <- summary(pglm2)$solution[c(2, 14, 3, 15)]
    
  }
  
  output <- list(results_1, spatial_size_correction, results_2)
  return(output)
}


tricking_tree<-function(df,tree){
  epsilon=0.0000001 #set the distance that will be used between populations of the same species
  
  for (i in 1:length(unique(df$SpeciesAccepted))){ #For each unique species
    
    nb_pop <- length(grep(unique(df$SpeciesAccepted)[i], df$SpeciesAccepted)) #find the number of replicates
    
    if ((nb_pop > 1) == TRUE){ #if there is more than one population
      
      if ((nb_pop > 2) == TRUE){
        
        ordre <- sample(2:(nb_pop))#create a random vector with the order of the population
        
      } else {
        
        ordre <- 2
        
      }
      
      pop_names <- paste(unique(df$SpeciesAccepted)[i], ordre, sep = "_") #create a list of names
      
      pop_names <- c(unique(df$SpeciesAccepted)[i], pop_names)
      
      for (j in 1:length(ordre)){
        
        tree<-bind.tip(tree, pop_names[j+1], edge.length = ordre[j] * epsilon, 
                       where=which(tree$tip.label==pop_names[j])) #Add the new population in the tree
        
        index <- mrca(tree)[pop_names[j+1], pop_names[j]] #last nodes between two populations
        
        tree$edge.length[which(tree$edge[,1] == index)] <- (nb_pop - j) * epsilon #Set the right length for the new branch
        
      }
    }
  }
  tree<-force.ultrametric(tree) #To ensure that the tree remain ultrametric after all the additions
  
  return(tree)
}

tricking_tree<-function(df,tree){
  epsilon=0.0000001 #set the distance that will be used between populations of the same species
  
  for (i in 1:length(unique(df$SpeciesAccepted))){ #For each unique species
    
    nb_pop <- length(grep(unique(df$SpeciesAccepted)[i], df$SpeciesAccepted)) #find the number of replicates
    
    if ((nb_pop > 1) == TRUE){ #if there is more than one population
      
      if ((nb_pop > 2) == TRUE){
        
        ordre <- sample(2:(nb_pop))#create a random vector with the order of the population
        
      } else {
        
        ordre <- 2
        
      }
      
      pop_names <- paste(unique(df$SpeciesAccepted)[i], ordre, sep = "_") #create a list of names
      
      pop_names <- c(unique(df$SpeciesAccepted)[i], pop_names)
      
      for (j in 1:length(ordre)){
        
        tree<-bind.tip(tree, pop_names[j+1], edge.length = ordre[j] * epsilon, 
                       where=which(tree$tip.label==pop_names[j])) #Add the new population in the tree
        
        index <- mrca(tree)[pop_names[j+1], pop_names[j]] #last nodes between two populations
        
        tree$edge.length[which(tree$edge[,1] == index)] <- (nb_pop - j) * epsilon #Set the right length for the new branch
        
      }
    }
  }
  tree<-force.ultrametric(tree) #To ensure that the tree remain ultrametric after all the additions
  
  return(tree)
}

#### 1- LOADING RAW DATA ####

load("COMADRE_v.3.0.0.Rdata")
load("COMPADRE_v.5.0.0.Rdata")

cleancompadre <- read.delim(file = "cleanCOM_P_ADRE_Output 4 Oct 2019.csv", sep = ",", header = TRUE) #data compiling the different traits value
species_and_HFP_complete <- read.delim(file = "COMPADRE_HFP.csv", sep=",", header = TRUE) #data with the value of HFP for the different species in COM(P)ADRE

animal_trait <-  read.csv("Animal_traits.csv", header = TRUE, sep = ",")

amniote_data <- read.csv("Amniote_Life_History_Database_10_feb_2020.csv", header = TRUE, sep = ",")
phylacine_data <- read.csv("Phylacine_trait_list.csv", header = TRUE, sep = ",")

plant_trait <-  read.csv("Plant_traits.csv", header = TRUE, sep = ",")

plant_size1 <- read.csv("mismatch_TRY_mydata_plant_height.csv", header = TRUE, sep = ",")
plant_size2 <- fread("8609.txt",header = T, sep = "\t", dec = ".", quote = "", data.table = T, showProgress = FALSE)

animal_tree<-read.newick (file = "animal_timetree.nwk")

load(file = "GBOTB.extended.rda")
load(file = "nodes.info.1.rda")

HFP_layer_10_new <- read.csv("HFP_10kmAVG_coords.csv")
animal_migration <- read.csv("Migration_distance.csv")


#### 2- CLEANING COMPADRE_HFP ####

#Transforming some colums from numeric to factor
cleancompadre[,'StartYear'] <- as.numeric(as.character(cleancompadre[,'StartYear']))
cleancompadre[,'EndYear'] <- as.numeric(as.character(cleancompadre[,'EndYear']))

#Here we remove all the data that do not have a full set of HFP values
species_and_HFP_complete <- na.omit(species_and_HFP_complete)

#### 3- CLEANING CLEANCOMPADRE ####

#Removing marine species
cleancompadre <- cleancompadre[!(cleancompadre$Habitat=="Marine"),]

#Removing non-mean matrices
cleancompadre <- cleancompadre[(cleancompadre$MatrixComposite=="Mean"),]

#Removing growth rate above 3
cleancompadre <- cleancompadre[!(cleancompadre$Lambda>3),]

#Removing data that are not unmanipulated
cleancompadre <- cleancompadre[(cleancompadre$Treatment=="Unmanipulated"),]

#Removing non ergodic, non irreducible and non primitive matrices
cleancompadre <- cleancompadre[(cleancompadre$Ergodic==TRUE),]
cleancompadre <- cleancompadre[(cleancompadre$Primitive==TRUE),]
cleancompadre <- cleancompadre[(cleancompadre$Irreducible==TRUE),]

#Removing studies that start before 1985
cleancompadre <- cleancompadre[!(cleancompadre$StartYear < 1985),]

#### 4- CREATION OF A NEW DATAFRAME THAT GATHER THE DATA WE WANT ####
df <-data.frame(matrix(ncol = 89,nrow = 1005))
colnames(df)<-c("SpeciesAccepted",
                "Kingdom",
                "Class",
                "Family",
                "Matrix_name",
                "Treatment",
                "Matrixtype",
                "Ergodic",
                "Primitive",
                "Irreducible",
                "StudyStart",
                "StudyEnd",
                "Lat",
                "Lon",
                "Indexcleancompadre",
                "Built1994",
                "Built2009",
                "Built",
                "Croplands1992",
                "Croplands2005",
                "Croplands",
                "HFP1993_int",
                "HFP2009_int",
                "HFP_int",
                "HFP1993",
                "HFP2009",
                "HFP",
                "Lights1994",
                "Lights2009",
                "Lights",
                "Navwater1994",
                "Navwater2009",
                "Navwater",
                "Pasture1993",
                "Pasture2009",
                "Pasture",
                "Popdensity1990",
                "Popdensity2010",
                "Popdensity",
                "Railways",
                "Roads",
                "HFP_date",
                "Lambda",
                "MatrixDimension",
                "R0",
                "S",
                "GenT",
                "Lalpha",
                "H",
                "pRep",
                "Rho",
                "Pi",
                "Reactivity",
                "FirstStepAtt",
                "WignerDisk.mean",
                "SurvSSD",
                "GrowSSD",
                "ShriSSD",
                "RepSSD",
                "CloSSD",
                "ShriSSDPreRep",
                "SurvSSDRep",
                "GrowSSDRep",
                "Esurv",
                "Eshri",
                "Eclo",
                "Egrow",
                "Erep",
                "EsurvPreRep",
                "EshriPreRep",
                "Ssurv",
                "Sgrow",
                "Sshri",
                "Srep",
                "SsurvRep",
                "habitat",
                "size",
                "diet",
                "locomotion",
                "coord1dim",
                "coord2dim",
                "coord3dim",
                "MaxAttenuation",
                "InertiaLow",
                "resilience",
                "resistance_FSA",
                "resistance_MA",
                "resistance_IL",
                "MaxAmp")


#Here we will find the corresponding rows between species_and_HFP_complete and cleancompadre
dfrow=1 
for (i in 1:length(species_and_HFP_complete$cat)){
  select <- cleancompadre[cleancompadre$SpeciesAuthor==stri_c(species_and_HFP_complete$SpeciesAuthor[i]),]
  select <- select[grep(TRUE,select$Population==stri_c(species_and_HFP_complete$MatrixPopulation[i])),]
  select <- select[grep(TRUE,select$Lon==species_and_HFP_complete$Lon[i]),]
  select <- select[grep(TRUE,select$Lat==species_and_HFP_complete$Lat[i]),]
  if ((length(select$X) > 0)==TRUE){
    for (selectrow in 1:length(select$X)){
      df$SpeciesAccepted[dfrow] <- stri_c(select$SpeciesAccepted[selectrow])
      df$Kingdom[dfrow] <- stri_c(select$Kingdom[selectrow])
      df$Class[dfrow] <- stri_c(select$Class[selectrow])
      df$Family[dfrow] <- stri_c(select$Family[selectrow])
      df$Matrix_name[dfrow] <- stri_c(select$Population[selectrow])
      df$Treatment[dfrow] <- stri_c(select$Treatment[selectrow])
      df$Matrixtype[dfrow] <- stri_c(select$MatrixComposite[selectrow])
      df$Ergodic[dfrow] <- stri_c(select$Ergodic[selectrow])
      df$Primitive[dfrow] <- stri_c(select$Primitive[selectrow])
      df$Irreducible[dfrow] <- stri_c(select$Irreducible[selectrow])
      df$StudyStart[dfrow] <- select$StartYear[selectrow]
      df$StudyEnd[dfrow] <- select$EndYear[selectrow]
      df$Lat[dfrow] <- select$Lat[selectrow]
      df$Lon[dfrow] <- select$Lon[selectrow]
      df$Indexcleancompadre[dfrow] <- select$X[selectrow]
      df$Built1994[dfrow] <- species_and_HFP_complete$Built1994[i]
      df$Built2009[dfrow] <- species_and_HFP_complete$Built2009[i]
      df$Croplands1992[dfrow] <- species_and_HFP_complete$croplands1992[i]
      df$Croplands2005[dfrow] <- species_and_HFP_complete$croplands2005[i]
      df$HFP1993_int[dfrow] <- species_and_HFP_complete$HFP1993_int[i]
      df$HFP2009_int[dfrow] <- species_and_HFP_complete$HFP2009_int[i]
      #df$HFP_int[dfrow]
      df$HFP1993[dfrow] <- species_and_HFP_complete$HFP1993[i]
      df$HFP2009[dfrow] <- species_and_HFP_complete$HFP2009[i]
      df$Lights1994[dfrow] <- species_and_HFP_complete$Lights1994[i]
      df$Lights2009[dfrow] <- species_and_HFP_complete$Lights2009[i]
      df$Navwater1994[dfrow] <- species_and_HFP_complete$NavWater1994[i]
      df$Navwater2009[dfrow] <- species_and_HFP_complete$NavWater2009[i]
      df$Pasture1993[dfrow] <- species_and_HFP_complete$Pasture1993[i]
      df$Pasture2009[dfrow] <- species_and_HFP_complete$Pasture2009[i]
      df$Popdensity1990[dfrow] <- species_and_HFP_complete$Popdensity1990[i]
      df$Popdensity2010[dfrow] <- species_and_HFP_complete$Popdensity2010[i]
      df$Railways[dfrow] <- species_and_HFP_complete$Railways[i]
      df$Roads[dfrow] <- species_and_HFP_complete$Roads[i]
      #df$HFP_date[dfrow]
      df$Lambda[dfrow] <- select$Lambda[selectrow]
      df$MatrixDimension[dfrow] <- select$MatrixDimension[selectrow]
      df$R0[dfrow] <- select$R0[selectrow]
      df$S[dfrow] <- select$S[selectrow]
      df$GenT[dfrow] <- select$GenTfun[selectrow]
      df$Lalpha[dfrow] <- select$La[selectrow]
      df$H[dfrow] <- select$H[selectrow]
      df$pRep[dfrow] <- select$pRep[selectrow]
      df$Rho[dfrow] <- select$Rho[selectrow]
      df$Pi[dfrow] <- select$Pi[selectrow]
      df$Reactivity[dfrow] <- select$Reactivity[selectrow]
      df$FirstStepAtt[dfrow] <- select$FirstStepAtt[selectrow]
      df$WignerDisk.mean[dfrow] <- select$WignerDisk.mean[selectrow]
      df$SurvSSD[dfrow] <- select$SurvSSD[selectrow]
      df$GrowSSD[dfrow] <- select$GrowSSD[selectrow]
      df$ShriSSD[dfrow] <- select$ShriSSD[selectrow]
      df$RepSSD[dfrow] <- select$RepSSD[selectrow]
      df$CloSSD[dfrow] <- select$CloSSD[selectrow]
      df$ShriSSDPreRep[dfrow] <- select$ShriSSDPreRep[selectrow]
      df$SurvSSDRep[dfrow] <- select$SurvSSDRep[selectrow]
      df$GrowSSDRep[dfrow] <- select$GrowSSDRep[selectrow]
      df$Esurv[dfrow] <- select$Esurv[selectrow]
      df$Eshri[dfrow] <- select$Eshri[selectrow]
      df$Eclo[dfrow] <- select$Eclo[selectrow]
      df$Egrow[dfrow] <- select$Egrow[selectrow]
      df$Erep[dfrow] <- select$Erep[selectrow]
      df$EsurvPreRep[dfrow] <- select$EsurvPreRep[selectrow]
      df$EshriPreRep[dfrow] <- select$EshriPreRep[selectrow]
      df$Ssurv[dfrow] <- select$Ssurv[selectrow]
      df$Sgrow[dfrow] <- select$Sgrow[selectrow]
      df$Sshri[dfrow] <- select$Sshri[selectrow]
      df$Srep[dfrow] <- select$Srep[selectrow]
      df$SsurvRep[dfrow] <- select$SsurvRep[selectrow]
      df$habitat[dfrow] <- stri_c(select$Habitat[selectrow])
      df$MaxAttenuation[dfrow] <- select$MaxAtt[selectrow]
      df$InertiaLow[dfrow] <- select$InertiaLow[selectrow]
      df$MaxAmp[dfrow] <- select$MaxAmp[selectrow]
      dfrow=dfrow+1
    }
  }
}

# Adding the HFP layer relative to the study ----------------------------------------------------------
#To add the HFP date relative to the study, we compare the dates of the study with the years of the two sets of
#HFP layers (1993 and 2009). 2001 is the median between the two dates.

for (k in 1:length(df$StudyStart)){
  if (is.na(df$StudyEnd[k]) == FALSE & is.na(df$StudyStart[k]) == FALSE){
    
    if (((abs(df$StudyStart[k] - 2001) >= abs(df$StudyEnd[k] - 2001)) & (df$StudyStart[k] < 2001)) == TRUE){
      df$HFP_date[k] <- 1993
    } else {
      df$HFP_date[k] <- 2009
    }
    
  } else {
    
    if (is.na(df$StudyEnd[k]) == TRUE & is.na(df$StudyStart[k]) == TRUE){
      df$HFP_date[k] <- NA
      
    } else {
      
      compare <- sum(df$StudyEnd[k],df$StudyStart[k], na.rm = TRUE)
      if (compare < 2001){
        df$HFP_date[k] <- 1993
      } else {
        
        df$HFP_date[k] <- 2009
      }
      
    }
  }
}

df <- drop_na(df,colnames(df)[42])

for (i in 1:length(df$SpeciesAccepted)){
  if (df$HFP_date[i] == 1993){
    df$HFP[i] <- df$HFP1993[i]
    df$Built[i] <- df$Built1994[i]
    df$Croplands[i] <- df$Croplands1992[i]
    df$Pasture[i] <- df$Pasture1993[i]
    df$Lights[i] <- df$Lights1994[i]
    df$Navwater[i] <- df$Navwater1994[i]
    df$Popdensity[i] <- df$Popdensity1990[i]
  } else {
    df$HFP[i] <- df$HFP2009[i]
    df$Built[i] <- df$Built2009[i]
    df$Croplands[i] <- df$Croplands2005[i]
    df$Pasture[i] <- df$Pasture2009[i]
    df$Lights[i] <- df$Lights2009[i]
    df$Navwater[i] <- df$Navwater2009[i]
    df$Popdensity[i] <- df$Popdensity2010[i]
  }
}

#### 5- HARMONIZATION OF THE NAMES ####
df$SpeciesAccepted<-gsub(" ", "_", df$SpeciesAccepted)

#Removing subspecies and false species like XXX sp
list<-grep("subsp", df$SpeciesAccepted) #to find subspecies
list<-c(list, grep("_sp.", df$SpeciesAccepted)) # to find false species (like XXXX sp.)
list<-sort(list, decreasing = FALSE)

index=0
for (i in 1:length(list)){
  df <- df[-(list[i]-index),]
  index <- index + 1
}

#### 6- POST-PROCESS CORRECTION ####
df <- df[!(df$SpeciesAccepted == "Cyprinodon_diabolis"),]
df <- df[!(df$SpeciesAccepted == "Ammocrypta_pellucida"),]
df <- df[!(df$SpeciesAccepted == "Arctodiaptomus_salinus"),]

df <- df[!(df$SpeciesAccepted == "Echinospartum_ibericum_algibicum"),]
df$SpeciesAccepted <- gsub("Vella_pseudocytisus_paui","Vella_pseudocytisus", df$SpeciesAccepted)

df <- drop_na(df, 1)

#### 7- GIVING ROW NAMES TO OUR DATA ####
#We are gonna create population names with a number to identify them clearly
pop_name <- unique(df$SpeciesAccepted)
row_name <- df$SpeciesAccepted

for (i in 1:length(pop_name)){
  location <- grep(pop_name[i], df$SpeciesAccepted)
  if (length(location) > 1){
    index <- 2:length(location) #number associated with the populations
    new_names <- paste(pop_name[i], index, sep = "_")
    new_names <- c(pop_name[i], new_names)
    for(j in 2:length(location)){
      row_name[location[j]] <- new_names[j]
    }
  }
}

row.names(df) <- row_name

#### 8- PCA WITH THE HFP ####

#PCA using the mean HFP value ------------------------------------------------------------------------
#Creation of a PCA using random sample allover the HFP map
samplesize <- 60000
set.seed(59)
rd_select <- sample_n(HFP_layer_10_new, samplesize)

PCA1 <- data.frame(matrix(nrow = samplesize, ncol = 8))
colnames(PCA1) <- c("Built", "Croplands", "Lights", "Navwater", "Pasture", "Popdensity", "Railways", "Roads")
Binom <- sample(c(0,1), samplesize, replace = TRUE)
for (i in 1:length(Binom)){
  if (Binom[i]==1){
    PCA1[i,] <- rd_select[i,c(3,5,11,13,15,17,19,20)]
  } else {
    PCA1[i,] <- rd_select[i,c(4,6,12,14,16,18,19,20)] 
  }
}

acp.res <- PCA(PCA1,scale.unit = TRUE, ncp=8, graph = FALSE)

#Interpolating the coordinates of our population on the PCA axis
species_HFP <- df[,c(18,21,30,33,36,39,40,41)]
colnames(species_HFP) <- c("Built", "Croplands", "Lights", "Navwater", "Pasture", "Popdensity", "Railways", "Roads")

PCA_prediction <- predict.PCA(acp.res, species_HFP)

df$coord1dim <- PCA_prediction$coord[,1]
df$coord2dim <- PCA_prediction$coord[,2]
df$coord3dim <- PCA_prediction$coord[,3]

#### 9- CREATION OF THE ANIMAL AND PLANT SUB-SELECTION ####

#First you need to run the PCA with all the species together from the Analysis section (3_Analysis)

dfanimal <- df[(df$Kingdom == "Animalia"),]
dfplant <- df[(df$Kingdom == "Plantae"),]

#### 10- FINDING ANIMAL SPECIES BODY MASS ####

animal_trait$body_mass_g <- NA
animal_trait$bodymass_source <- NA

# Looking for species body mass in amniote life history database ------------------------------------------

amniote_data$SpeciesAccepted <- paste(amniote_data$genus, amniote_data$species, sep = "_")

for (i in 1:length(animal_trait$SpeciesAccepted)){
  index <- grep(animal_trait$SpeciesAccepted[i], amniote_data$SpeciesAccepted)
  if (length(index)>0){
    animal_trait$body_mass_g[i] <- amniote_data$adult_body_mass_g[index[1]]
    animal_trait$bodymass_source[i] <- "Amniote"
  }
}

# Looking for species body mass in Phylacine -------------------------------------------------------------

for (i in 1:length(animal_trait$SpeciesAccepted)){
  index <- grep(animal_trait$SpeciesAccepted[i], phylacine_data$Binomial.1.2)
  if (length(index)>0 & is.na(animal_trait$body_mass_g[i]) == TRUE){
    animal_trait$body_mass_g[i] <- phylacine_data$Mass.g[index[1]]
    animal_trait$bodymass_source[i] <- "Phylacine"
  }
}

# Adding the body mass of Antrhopoides paradiseus --------------------------------------------------------
#According to the Shaw et al., 2010 the bodymass of A. paradiseus vary from 4600 to 5090 g. We decided to take the
#mean value: 4845
#Jessica M Shaw , Andrew R Jenkins , Peter G Ryan & Jon J Smallie (2010) A preliminary survey of avian mortality on power lines in the Overberg, South Africa, OSTRICH, 81:2, 109-113, DOI:10.2989/00306525.2010.488421

animal_trait$body_mass_g[which(animal_trait$SpeciesAccepted == "Anthropoides_paradiseus")] <- 4845
animal_trait$bodymass_source[which(animal_trait$SpeciesAccepted == "Anthropoides_paradiseus")] <- "Shaw et al., 2010"

#### 11- ADDING PLANT SIZE ####

#part of script from Pol Capdevila to get plant height

compadreUse <- cdb_fetch("compadre")

compadreUse <- cdb_flatten(compadreUse)

plant_size2 <- plant_size2 %>%
  filter(AccSpeciesName%in%compadreUse$SpeciesAccepted,
         TraitName=="Plant height vegetative") %>% 
  group_by(AccSpeciesName) %>%
  summarise(mean_height=mean((StdValue)),
            max_height=max((StdValue))) %>% 
  left_join(plant_size2, ., by = c('AccSpeciesName')) %>%
  group_by(AccSpeciesName) %>%
  filter(StdValue==max_height) %>%
  dplyr::select(AccSpeciesName, mean_height, max_height, Reference)%>%
  distinct(AccSpeciesName, .keep_all = T)

plant_trait$Growth_form <- compadre$metadata$OrganismType[grep(plant_trait$SpeciesAccepted,compadre$metadata$SpeciesAccepted)[1]]

for (i in 1:length(plant_trait$SpeciesAccepted)){
  select <- compadre$metadata[gsub(" ", "_", compadre$metadata$SpeciesAccepted) == stri_c(plant_trait$SpeciesAccepted[i]), ]
  plant_trait$Growth_form[i] <- unique(drop_na(select, OrganismType)$OrganismType)
  select <- grep(plant_trait$SpeciesAccepted[i],plant_size1$x)
  if (length(select)>0) {
    plant_trait$size.m[i] <- plant_size1$height..cm.[select]*0.01
  } else {
    select <- grep(plant_trait$SpeciesAccepted[i],gsub(" ", "_", plant_size2$AccSpeciesName))
    if (length(select)>0){
      plant_trait$size.m[i] <- plant_size2$mean_height[select]
    }
  }
}

#### 12- ADDING ADDITIONAL TRAITS ####

#First you need to run the script "1,5_Getting_morphological_and_functionnal_trait"
for (i in 1:length(animal_trait[,1])){
  
  index <- grep(animal_trait$SpeciesAccepted[i],dfanimal$SpeciesAccepted)
  dfanimal$size[index] <- animal_trait$body_mass_g[i]
  
}

plant_trait$SpeciesAccepted <- gsub(" ", "_", plant_trait$SpeciesAccepted)
for (i in 1:length(plant_trait[,1])){
  index <- grep(plant_trait$SpeciesAccepted[i],dfplant$SpeciesAccepted)
  dfplant$size[index] <- plant_trait$size.m[i]
}


#### 13- CLASSIC PHYLOGENY FOR ANIMALS AND PLANTS ####

animal_pop_tree <- read.tree(file = "animal_pop_phylo.tre")

#is.binary(animal_pop_tree)
#is.ultrametric(animal_pop_tree)
#is.rooted(animal_pop_tree)

plant_pop_tree <- read.tree(file = "plant_pop_phylo.tre")

#is.binary(plant_pop_tree)
#is.ultrametric(plant_pop_tree)
#is.rooted(plant_pop_tree)


#### 14- DATA TRANSFORMATION AND PREREQUISITE FOR THE MAIN LHT ANALYSES ####
#We are gonna stock the variables that need to be transform later analysis

list_traits <- c(44:75)#column index of the traits we are interested in

#For plants

list_traits_plant <- list_traits[-c(1,12,21,22,23,24,25,26,27,28,29,30,31,32)]

value_checking_plant <- c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)

#For animals
list_traits_animal <- list_traits[-c(1, 15, 17, 18, 22, 23, 27, 30)]
list_traits_animal <- list_traits_animal[-c(11,17,18,19,20,21,22,23,24)]
value_checking_animal <- c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)


#Log transformation of the body mass
dfanimal$logsize <- log(dfanimal$size)
dfplant$logsize <- log(dfplant$size)


#Assessing the transformation that need to be done
#A value of 1 means no transformation, a value of 2 means log transformation, 3 means square root transformation and 4 means inverse transformation
list_transform_animal <- c(2,1,1,1,1,3,2,4,2,1,1,1,2,1,1)
list_transform_plant <- c(2,1,2,3,1,3,2,2,2,1,1,1,2,2,1,3,2,3)


#### 15- EFFECT OF HFP ON DEMOGRAPHIC TRAITS WITH CORRECTIONS ####

selection_with_interaction <- c(2, 14, 3, 15)
selection_without_interaction <- c(2, 14, 3, 15)

# Model for animals with mean HFP value -------------------------------------------------------------
set.seed(15)
output_animal <- Running_MCMCglmm(dfanimal, list_traits_animal, value_checking_animal, list_transform_animal, animal_pop_tree)

# Model for plants with mean HFP value --------------------------------------------------------------
set.seed(16)
output_plant <- Running_MCMCglmm(dfplant, list_traits_plant, value_checking_plant, list_transform_plant, plant_pop_tree)

#### 16-1 - OUTPUTS OF THE SIMULATION ####

write.csv(output_animal[[1]], file = "Animal_phylo_results.csv")
write.csv(output_animal[[2]], file = "Animal_results_correction.csv")
write.csv(output_animal[[3]], file = "Animal_residual_mixed_model.csv")

write.csv(output_plant[[1]], file = "Plant_phylo_results.csv")
write.csv(output_plant[[2]], file = "Plant_results_correction.csv")
write.csv(output_plant[[3]], file = "Plant_residual_mixed_model.csv")


#### 16-2- SIMULATION FOR ANIMALS DEPENDING ON THE MIGRATION DISTANCE ####

for (i in 1:length(dfanimal$SpeciesAccepted)){
  dfanimal$migration[i] <- animal_migration$Distance[animal_migration$SpeciesAccepted==dfanimal$SpeciesAccepted[i]]
}

set.seed(15)
output_animal_long_dist <- Running_MCMCglmm(dfanimal[dfanimal$migration=="Long",], list_traits_animal, value_checking_animal, list_transform_animal, animal_pop_tree)
output_animal_short_dist <- Running_MCMCglmm(dfanimal[dfanimal$migration=="Short",], list_traits_animal, value_checking_animal, list_transform_animal, animal_pop_tree)

write.csv(output_animal_long_dist[[1]], file = "Animal_long_phylo_results.csv")
write.csv(output_animal_long_dist[[2]], file = "Animal_long_results_correction.csv")
write.csv(output_animal_long_dist[[3]], file = "Animal_long_residual_mixed_model.csv")

write.csv(output_animal_short_dist[[1]], file = "Animal_short_phylo_results.csv")
write.csv(output_animal_short_dist[[2]], file = "Animal_short_results_correction.csv")
write.csv(output_animal_short_dist[[3]], file = "Animal_short_residual_mixed_model.csv")

#### 16-3- TRAIT ANALYSIS WITHOUT CORRECTION ####

#Besoin de stocker les résultats des différentes régression dans un seul objet : glmm_stock. 
glmm <- as.data.frame()

for (i in 1:length(list_traits_animal)){

  if (value_checking_animal[k] == 1){
    
    list <- Removing_outliers_and_false_values(dfanimal, list_traits_animal[k])
    
    dfbis <- list[[1]]
    storage <- list[[2]]
    
  } else {
    
    if (value_checking[k] == 2){
      
      list <- Removing_outliers(dfanimal, list_traits_animal[k])
      
      dfbis <- list[[1]]
      storage <- list[[2]]
      
    } else {
      
      #if (value_checking[k] == 3){
      
      list <- Removing_outliers_and_zero(dfanimal, list_traits_animal[k])
      
      dfbis <- list[[1]]
      storage <- list[[2]]
      
      #}
      
    }
  }
  
  #How to transform the data
  if (list_transform_animal[k] == 1){
    
    dfbis$transform <- dfbis[, list_traits_animal[k]]
    
  } else {
    
    if (list_transform_animal[k] == 2) {
      
      dfbis$transform <- log(dfbis[, list_traits_animal[k]] + 0.000001)
      
    } else {
      
      if (list_transform_animal[k] == 3) {
        
        dfbis$transform <- sqrt(dfbis[, list_traits_animal[k]] + 0.000001)
        
      } else {
        
        if (list_transform_animal[k] == 4) {
          
          dfbis$transform <- (1/(dfbis[, list_traits_animal[k]] + 0.000001))
          
        } else {
          
          #if (list_transform[k] == 5) {
          
          dfbis$transform <- (sqrt(dfbis[, list_traits_animal[k]] + 0.0001))
          
          #}
          
        }
        
      }
      
    }
    
  }
  
  glmm[k] <- glmer(dfbis[,list_traits_animal[k]] ~ dfbis$coord1dim * dfbis$coord2dim + (1|dfbis$SpeciesAccepted), family = gaussian)
}

#### 17- RESILIENCE - RESISTANCE FRAMEWORK ANALYSIS ####

# Setting the variables -----------------------------------------------------------------------------
#Creation of the resilience and resistance variables

for (i in 1:length(df$SpeciesAccepted)){
  df$resilience[i] <- log(df$Rho[i] + 0.000001)
  if (is.na(df$FirstStepAtt[i]) == FALSE){
    if (df$FirstStepAtt[i] == 1){
      df$resistance_FSA[i] <- NA
    } else {
      df$resistance_FSA[i] <- log(1-df$FirstStepAtt[i])
    }
  } else {
    df$resistance_FSA[i] <- NA
  }
  if (is.na(df$MaxAttenuation[i]) == FALSE){
    if (df$MaxAttenuation[i] == 1){
      df$resistance_MA[i] <- NA
    } else {
      df$resistance_MA[i] <- log(1-df$MaxAttenuation[i])
    }
  } else {
    df$resistance_MA[i] <- NA
  }
  if (is.na(df$InertiaLow[i]) == FALSE){
    if (df$InertiaLow[i] == 1){
      df$resistance_IL[i] <- NA
    } else {
      df$resistance_IL[i] <- log(1-df$InertiaLow[i])
    }
  } else {
    df$resistance_IL[i] <- NA
  }
}

for (i in 1:length(df$SpeciesAccepted)){
  df$resonance[i] <- log(df$MaxAmp[i] + 0.000001)
}

#For animals

for (i in 1:length(dfanimal$SpeciesAccepted)){
  dfanimal$resonance[i] <- log(dfanimal$MaxAmp[i] + 0.000001)
}

for (i in 1:length(dfanimal$SpeciesAccepted)){
  dfanimal$resilience[i] <- log(dfanimal$Rho[i] + 0.000001)
  if (is.na(dfanimal$FirstStepAtt[i]) == FALSE){
    if (dfanimal$FirstStepAtt[i] == 1){
      dfanimal$resistance_FSA[i] <- NA
    } else {
      dfanimal$resistance_FSA[i] <- log(1-dfanimal$FirstStepAtt[i])
    }
  } else {
    dfanimal$resistance_FSA[i] <- NA
  }
  if (is.na(dfanimal$MaxAttenuation[i]) == FALSE){
    if (dfanimal$MaxAttenuation[i] == 1){
      dfanimal$resistance_MA[i] <- NA
    } else {
      dfanimal$resistance_MA[i] <- log(1-dfanimal$MaxAttenuation[i])
    }
  } else {
    dfanimal$resistance_MA[i] <- NA
  }
  if (is.na(dfanimal$InertiaLow[i]) == FALSE){
    if (dfanimal$InertiaLow[i] == 1){
      dfanimal$resistance_IL[i] <- NA
    } else {
      dfanimal$resistance_IL[i] <- log(1-dfanimal$InertiaLow[i])
    }
  } else {
    dfanimal$resistance_IL[i] <- NA
  }
}

for (i in 1:length(dfanimal$SpeciesAccepted)){
  dfanimal$resonance[i] <- log(dfanimal$MaxAmp[i] + 0.000001)
}

#for plants

for (i in 1:length(dfplant$SpeciesAccepted)){
  dfplant$resonance[i] <- log(dfplant$MaxAmp[i] + 0.000001)
}

for (i in 1:length(dfplant$SpeciesAccepted)){
  dfplant$resilience[i] <- log(dfplant$Rho[i] + 0.000001)
  if (is.na(dfplant$FirstStepAtt[i]) == FALSE){
    if (dfplant$FirstStepAtt[i] == 1){
      dfplant$resistance_FSA[i] <- NA
    } else {
      dfplant$resistance_FSA[i] <- log(1-dfplant$FirstStepAtt[i])
    }
  } else {
    dfplant$resistance_FSA[i] <- NA
  }
  if (is.na(dfplant$MaxAttenuation[i]) == FALSE){
    if (dfplant$MaxAttenuation[i] == 1){
      dfplant$resistance_MA[i] <- NA
    } else {
      dfplant$resistance_MA[i] <- log(1-dfplant$MaxAttenuation[i])
    }
  } else {
    dfplant$resistance_MA[i] <- NA
  }
  if (is.na(dfplant$InertiaLow[i]) == FALSE){
    if (dfplant$InertiaLow[i] == 1){
      dfplant$resistance_IL[i] <- NA
    } else {
      dfplant$resistance_IL[i] <- log(1-dfplant$InertiaLow[i])
    }
  } else {
    dfplant$resistance_IL[i] <- NA
  }
}

for (i in 1:length(dfplant$SpeciesAccepted)){
  dfplant$resonance[i] <- log(dfplant$MaxAmp[i] + 0.000001)
}


# Creation of the convexhull of the resistance-resilience framework ---------------------------------

framework <- data.frame(df$coord1dim ,df$resistance_FSA, df$resilience, df$resonance)
framework2 <- data.frame(df$coord2dim ,df$resistance_FSA, df$resilience, df$resonance)
framework <- drop_na(framework)
framework2 <- drop_na(framework2)
pc <- (as.matrix(framework))
pc2 <- (as.matrix(framework2))

framework3 <- data.frame(df$coord1dim , df$coord2dim, df$resistance_FSA, df$resilience, df$resonance)
framework3 <- drop_na(framework3)
colnames(framework3) <- c("Human_presence", "Agricultural_landuse", "Resistance", "Speed_of_recovery", "Resonance")
pc3 <- as.matrix(framework3)

frames <- framework
pcs <- pc

hull <- convhulln(pcs, "FA")
obs.vol <- hull$vol

# Creation of the null hypothesis of a uniform distribution of the traits ---------------------------
# Code from : Díaz, S., Kattge, J., Cornelissen, J.H.C., Wright, I.J., Lavorel, S., Dray, S., Reu, B., Kleyer, M., Wirth, C., Colin Prentice, I., Garnier, E., Bönisch, G., Westoby, M., Poorter, H., Reich, P.B., Moles, A.T., Dickie, J., Gillison, A.N., Zanne, A.E., Chave, J., Joseph Wright, S., Sheremet’ev, S.N., Jactel, H., Baraloto, C., Cerabolini, B., Pierce, S., Shipley, B., Kirkup, D., Casanoves, F., Joswig, J.S., Günther, A., Falczuk, V., Rüger, N., Mahecha, M.D., Gorné, L.D., 2016. The global spectrum of plant form and function. Nature 529, 167–171. https://doi.org/10.1038/nature16489



sim.uniform.vol <- rep(0,1000)
sim.uniform.overlap <- rep(0,1000)
set.seed(446)
for(i in 1:1000){
  newpc_unif1 <- matrix(pcs[,1])
  newpc_unif <- (apply(pcs[,-1], 2, function(x) runif(nrow(pcs), min = min(x), max = max(x))))
  newpc_unif <- cbind(newpc_unif1, newpc_unif)
  hull_unif <- convhulln(newpc_unif,"FA")
  #uni_overlap <- geometry::intersectn(hull$hull, hull_unif$hull)
  #sim.uniform.overlap[i] <- uni_overlap$ch$vol/(uni_overlap$ch1$vol)
  sim.uniform.vol[i] <- hull_unif$vol
}

# Creation of the null hypothesis of a normal distribution of the traits ----------------------------
# NB: The distribution of the PC1 location is still uniform
sim.norm.vol <- rep(0,1000)
sim.norm.overlap <- rep(0,1000)
set.seed(447)
for(i in 1:1000){
  newpc1 <- matrix(pcs[,1])
  #newpc2 <- matrix(runif(nrow(pcs), min = min(pcs[,2]), max = max(pcs[,2])))
  newpc3 <- matrix(rnorm(nrow(pcs), mean = mean(frames$df.resistance_FSA), sd = sd(frames$df.resistance_FSA)), nrow(pcs), 1)
  newpc4 <- matrix(rnorm(nrow(pcs), mean = mean(frames$df.resilience), sd = sd(frames$df.resilience)), nrow(pcs), 1)
  newpc5 <- matrix(rnorm(nrow(pcs), mean = mean(frames$df.resonance), sd = sd(frames$df.resonance)), nrow(pcs), 1)
  new_pc_norm <- cbind(newpc1, newpc3, newpc4, newpc5)
  hull_norm <- convhulln(new_pc_norm,"FA")
  #norm_overlap <- geometry::intersectn(hull$hull, hull_norm$hull, fp = )
  #sim.norm.overlap[i] <- norm_overlap$ch$vol/(norm_overlap$ch1$vol)
  sim.norm.vol[i] <- hull_norm$vol
}

#Creation of the null hypothesis of a gamma distribution of the traits ------------------------------

sim.gamma.vol <- rep(0,1000)
sim.gamma.overlap <- rep(0,1000)

set.seed(447)
for(i in 1:1000){
  newpc1 <- matrix(pcs[,1])
  newpc2 <- matrix(rgamma(nrow(pcs), shape = ((mean(frames$df.resistance_FSA)^2)/sd(frames$df.resistance_FSA)), scale = sd(frames$df.resistance_FSA)/((mean(frames$df.resistance_FSA))^2)), nrow(pcs), 1)
  newpc3 <- matrix(rgamma(nrow(pcs), shape = ((mean(frames$df.resilience)^2)/sd(frames$df.resilience)), scale = sd(frames$df.resilience)/((mean(frames$df.resilience))^2)), nrow(pcs), 1)
  newpc4 <- matrix(rgamma(nrow(pcs), shape = ((mean(frames$df.resonance)^2)/sd(frames$df.resonance)), scale = sd(frames$df.resonance)/((mean(frames$df.resonance))^2)), nrow(pcs), 1)
  new_pc_gamma <- cbind(newpc1, newpc2, newpc3, newpc4)
  hull_gamma <- convhulln(new_pc_gamma,"FA")
  #gamma_overlap <- geometry::intersectn(hull$hull, hull_gamma$hull)
  #sim.gamma.overlap[i]<- gamma_overlap$ch$vol/(gamma_overlap$ch1$vol)
  sim.gamma.vol[i] <- hull_gamma$vol
}

# Comparison of the different model -----------------------------------------------------------------

test.unif <- as.randtest(obs=obs.vol,sim=sim.uniform.vol, alter="greater")
test.norm <- as.randtest(obs=obs.vol,sim=sim.norm.vol, alter="greater")
test.gamma <- as.randtest(obs=obs.vol,sim=sim.gamma.vol, alter="less")

100 - mean(obs.vol/test.unif$expvar[2]) * 100
100 - mean(obs.vol/sim.uniform.vol) * 100
mean(sim.uniform.vol)/obs.vol

100 - mean(obs.vol/test.norm$expvar[2])*100
100 - mean(obs.vol/sim.norm.vol)*100
mean(sim.norm.vol)/obs.vol

100 - mean(obs.vol/test.gamma$expvar[2])*100
100 - mean(obs.vol/sim.gamma.vol)*100
mean(sim.gamma.vol)/obs.vol

# Overlap between volumes ----------------------------------------------------------------------------

uni_overlap <- geometry::intersectn(hull$hull, hull_unif$hull)
uni_overlap$ch$vol/(uni_overlap$ch1$vol)

norm_overlap <- geometry::intersectn(hull$hull, hull_norm$hull)
norm_overlap$ch$vol/(norm_overlap$ch1$vol)

gamma_overlap <- geometry::intersectn(hull$hull, hull_gamma$hull)
gamma_overlap$ch$vol/(gamma_overlap$ch1$vol)


