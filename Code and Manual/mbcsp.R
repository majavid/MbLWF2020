############################# An LWF CG from the LCD package ###################################
toy.graph <- matrix(c(0,1,1,0,0,0,0,0,0,0,0,
                      1,0,0,1,0,0,0,0,0,0,0,
                      1,0,0,0,1,0,0,0,0,0,0,
                      0,1,0,0,0,1,1,0,0,0,0,
                      0,0,0,0,0,1,0,0,1,0,0,
                      0,0,0,0,1,0,0,0,0,0,1,
                      0,0,0,0,0,0,0,1,0,0,1,
                      0,0,0,0,0,0,1,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,1,0,
                      0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0),
                    nrow = 11, byrow = TRUE)
rownames(toy.graph) <- colnames(toy.graph) <- c("a","b","c","d","e","f","g","h","i","j","k")
#############################################################################################################
#############################################################################################################
############################# Plot LWF CGs ##################################################################
#############################################################################################################
#############################################################################################################
# Input
# -amat: adjacency matrix of the LWF CG

## ----------------------------------------------------------------------
## Author: ***************** August, 2019 (adapted from the "lcd")
################################################################################
################################################################################
################################################################################
# library(igraph)
`plotCG` <- function(amat, plain = TRUE){
  `skeletonCG` <- function(amat)
  {
    0 + (amat + t(amat) > 0)
  }
    e.total <- length(which(skeletonCG(amat) == 1))/2
    if(e.total == 0)
      return("The graph has no edge!")
    ## el <- matrix("", nrow = e.total, ncol = 2)
    el <- matrix(NA, nrow = e.total, ncol = 2)
    vset <- sort(rownames(amat))
    amat <- amat[vset,vset]
    p <- nrow(amat)
    arrow <- c()
    k <- 1
    for(i in 1:(p-1))
      for(j in (i+1):p){
        if(amat[i,j] == 1){
          ## el[k,1] <- i-1
          el[k,1] <- i
          ## el[k,2] <- j-1
          el[k,2] <- j
          k <- k + 1
          if(amat[j,i] == 0){
            arrow <- c(arrow, ">")
          } else {
            arrow <- c(arrow, "-")
          }
        } else {
          if(amat[j,i] == 1){
            ## el[k,1] <- j-1
            el[k,1] <- j
            ## el[k,2] <- i-1
            el[k,2] <- i
            k <- k + 1
            arrow <- c(arrow, ">")
          }
        }
      }
    g <- graph.empty()
    g <- add.vertices(g, p)
    V(g)$name <- vset
    g <- add.edges(g, t(el))
    if(plain){
      plot.igraph(g, #layout = layout.reingold.tilford,
                  edge.arrow.mode = arrow,
                  vertex.size = min(15, 500/p),
                  edge.arrow.size = min(1, 30/p))
    } else{
      tkplot(g, layout = layout.reingold.tilford,
             edge.arrow.mode = arrow,
             vertex.size = min(15, 500/p),
             edge.arrow.size = min(1, 30/p))
    }
    V(g)
  }






#############################################################################################################
#############################################################################################################
############################# MBC-CSP Algorithm for LWF CGs ################################################
#############################################################################################################
#############################################################################################################
# -dataset: a data frame containing the variables in the model.

# -cluster: an optional cluster object from package parallel.

# -test: a character string, the label of the conditional independence test to be used in the
# algorithm. If none is specified, the default test statistic is the mutual information
# for categorical variables, the Jonckheere-Terpstra test for ordered factors and
# the linear correlation for continuous variables.

# -alpha: a numeric value, the target nominal type I error rate.

# -max.sx: a positive integer, the maximum allowed size of the conditioning sets used in
# conditional independence tests. The default is that there is no limit on size 

## ----------------------------------------------------------------------
## Author: ***************** August, 2019 (adapted from the "bnlearn" and "lcd")
################################################################################
################################################################################
################################################################################

`mbcsp` <- function(dataset,alpha = 0.05 , test = "zf", cluster = NULL,
                    max.sx = ncol(dataset), debug = FALSE){
  data.info = bnlearn:::check.data(dataset, allow.missing = TRUE)
  complete=data.info$complete.nodes
  corr = cor(dataset)
  #####################################################
  ##################### Computing Adjs & sepset  ######
  #####################################################
  
  ##  Input: T: the target variable; alpha: significance level
  ##  Output: adj(T): the set of variables adjacent to T; sepset: separation set
  
  ## how to build a sepset with names
  ## note that dataset should be in the data frame format
  
  `adjecencies` <- function(dataset, alpha = alpha){
    sepset <- lapply(colnames(dataset), function(.) {vec<-vector("list",ncol(dataset))
    names(vec)<-colnames(dataset)
    return(vec)})
    names(sepset)<-colnames(dataset)
    Adjs <-vector("list",ncol(dataset))
    names(Adjs)<-colnames(dataset)
    
    
    
    for (T in colnames(dataset)){
      #message("T = ", T)
      S = character(0) # empty set
      adjT = character(0) # empty set
      for (v in setdiff(colnames(dataset),T)) {
        pval <- bnlearn:::indep.test(v, T, sx = S, test = test, data = dataset,
                                     B = 0L, alpha = alpha, complete = complete)
        # #message("v = ",v)
        # #message("pval = ",pval)
        if (abs(pval[1]) >= alpha){
          sepset[[T]][[v]] <- S
        }else{
          adjT <- append(adjT, v)
        }
      }
      # Sort adj(T ) in increasing corr(Vi, T ) value
      sorted_adj <- sort(corr[,T][adjT])
      # #print(sorted_adj)
      # #print(names(sorted_adj))
      sorted_adj <- names(sorted_adj)
      k = 1
      nNum = min(length(adjT),max.sx)
      if(nNum < 1){
        Adjs[[T]] <- character(0)
      }else{
        while(k <= nNum){
          # all subsets S \subset adj(T) needs to be tested whose
          # sizes do not exceed k
          for (v in sorted_adj) {
            #message("v = ", v)
            condset=setdiff(sorted_adj,v) 
            len <- length(condset)
            if(k > len){
              k <- nNum
              break
            }
            #for (j in 1:k) {
            a <- bnlearn:::allsubs.test(v, T, sx=condset, fixed = character(0), data=dataset, test=test, B = 0L,
                                        alpha = alpha, min = k, max = k, complete=complete, debug = FALSE)
            #print(a)
            if(a["p.value"] >= alpha){
              sorted_adj <- setdiff(sorted_adj,v)
              sepset[[T]][[v]] <- attr(a,"dsep.set")
            }
            #} 
          }
          k <- k+1
          nNum <- length(sorted_adj)
        }
        Adjs[[T]] <- sorted_adj
      }
    }# End of Alg 2
    # message("before symmetry correction")
    # print(Adjs)
    ## symmetry correction for removing false positives from adjV(T)
    ## i.e., if X \in adjV(Y) but Y \not\in adjV(X) then 
    ## adjV(Y) = adjV(Y)\{X} and sepset[[Y]][[X]] <- sepset[[X]][[Y]]
    for (var in colnames(dataset)) {
      for (u in Adjs[[var]]) {
        if(!(var %in% Adjs[[u]])){
          Adjs[[var]] <- setdiff(Adjs[[var]],u)
          sepset[[var]][[u]] <- sepset[[u]][[var]]
        }
      }
    }
    # message("after symmetry correction")
    # print(Adjs)
    return(list(adjs=Adjs,sepset = sepset))
  }#End adjecencies
  #############################################################################################
  ##################### MMBC-CSP Algorithm for Mb recovery ####################################
  #############################################################################################
  ## Input: T : the target variable; alpha: significant level
  ## Output: Mb(T )
  
  `findMb` <- function(T, dataset, adjs, sepset,alpha = alpha){
    adjT <- adjs[[T]]
    # bd(T) U ch(T) \subsetq Mb(T)
    mmbT <- adjT
    # complex-spouses recovery phase
    # candidates of complex-spouses
    candids <- setdiff(colnames(dataset),union(T,adjT))
    for (ch in adjT) {
      for (csp in candids) {
        pval1 <- bnlearn:::indep.test(csp, T, sx = sepset[[T]][[csp]], test = test, data = dataset,
                                      B = 0L, alpha = alpha, complete = complete)
        pval2 <- bnlearn:::indep.test(csp, T, sx = union(sepset[[T]][[csp]],ch), test = test, data = dataset,
                                      B = 0L, alpha = alpha, complete = complete)
        if(abs(pval1[1]) > alpha & abs(pval2[1]) <= alpha){
          mmbT <- union(mmbT,csp)
        }
      }
    }
    ## false positives phase (shrinking phase)
    continue <- TRUE
    delcand <- mmbT    # only those added later is allowed to be deleted
    if (length(delcand) == 0)
      continue <- FALSE
    while (continue) {       # shrink the Markov blanket
      p.val <- sapply(delcand, function(x)
        bnlearn:::indep.test(T, x, sx = setdiff(mmbT,x), test = "zf", data = dataset,
                             B = 0L, alpha = alpha, complete = complete))
      # multinom.ci.test(freq.tb, var, x, setdiff(mmbT,x))$p.value)
      ## this step could be speeded up significantly!!!
      p.val.max <- max(p.val)
      idx <- which(p.val == p.val.max)[1]
      if(p.val.max > alpha) {
        mmbT <- setdiff(mmbT, delcand[idx])
        delcand <- delcand[-idx]
      } else {
        continue <- FALSE
      }
    }
    return(mmbT)
  }#End findMb
  
  # 1. [Compute Markov Blankets]
  result <- adjecencies(dataset = dataset,alpha = alpha)
  nam <- names(dataset)
  mb = bnlearn:::smartSapply(cluster, as.list(nam), findMb, dataset = dataset, 
                             adjs = result$adjs, sepset = result$sepset, alpha=alpha)
  names(mb) = nam
  p <- ncol(dataset)
  G <- matrix(0, p, p)
  colnames(G)<-rownames(G)<-nam
  # a list similar to sepset in pcalg package
  # initialize all elements = 0
  # sepset<-apply(G, 1, as.list)
  
  # 2. [Compute Graph Structure]
  ###############################################################
  # # check symmetry in the output of the algorithm
  # message("before symmetry correction")
  # print(mb)
  ## symmetry correction for removing false positives from Mb(T)
  ## i.e., if X \in Mb(Y) but Y \not\in Mb(X) then 
  ## Mb(Y) = Mb(Y)\{X} and sepset[[Y]][[X]] <- sepset[[X]][[Y]]
  for (var in colnames(dataset)) {
    for (u in mb[[var]]) {
      if(!(var %in% mb[[u]])){
        mb[[var]] <- setdiff(mb[[var]],u)
        #sepset[[var]][[u]] <- sepset[[u]][[var]]
      }
    }
  }#End check symmetry

  return(mb)
}