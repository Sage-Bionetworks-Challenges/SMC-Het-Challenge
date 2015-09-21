# Simulation from the V-measure paper to test the desirable
# properties of the different metrics empirically.

#### PREAMBLE #################################################
library(BoutrosLab.plotting.general)

#### SIMULATION PARAMETERS ####################################
# flag for whether to print information about the simulation status
output = F
# number of iterations to perform for each set of parameters
n.iter = 1
# number of elements in each class
n.C = 200
# number of useful classes
C = 5
# number of useful clusters
K.u = 5
# number of noise clusters
K.n = 4
# percentage of elements that are useful but assigned to an 
# incorrect useful cluster
epsilon1 = 0.1
# percentage of elements that are useful but assigned to 
# noise clusters
epsilon2 = 0.3
# percentage of elements that are incorrectly assigned
epsilon = epsilon1 + epsilon2

# names of parameters that can be specified for in the simulation
# used in the Data Output/Analysis section
sim.params = c("K.u", "K.n", "epsilon1", "epsilon2")

#### HELPER FUNCTIONS ########################################

#### Data Creation Helper Functions ##########################

#### assign.useful ###########################################
# Assign elements from classes to useful clusters
#
# INPUT
# class.pairs - vector of length K.u where the i'th entry is 1
#       if cluster i is paired with the class in question, and 
#       is 0 otherwise
# OUTPUT
# class.elem - vector of length K.u where the i'th entry is
#       the number of elements of the class in question that are
#       in the i'th useful cluster
assign.useful <- function(class.pairs){
  # number of elements in each class that are incorrectly assigned 
  # to the wrong useful cluster
  num.class.incorr <- floor(n.C * (epsilon1))
  num.clust.unpair <- sum(1 - class.pairs) # number of clusters paired with this class
  
  num.class.corr <- n.C * (1 - epsilon) # number of elements in each class that are correctly assigned
  num.clust.pair <- sum(class.pairs) # number of clusters paired with this class
  
  class.elem <- (class.pairs * floor(num.class.corr / num.clust.pair)) + # correctly assigned elements
    ((1 - class.pairs) * floor(num.class.incorr / num.clust.unpair)) # incorrectly assigned elements
  
  toassign.corr = num.class.corr %% num.clust.pair
  toassign.incorr = num.class.incorr %% num.clust.unpair
  # start looking at a random index so that the extra elements are
  # assigned randomly and not all to the same classes
  index <- sample(1:K.u, 1)
  # check to ensure that elements are not "lost"
  while(toassign.corr > 0 || toassign.incorr > 0){ # assign "extra" elements
    #print(paste("Correct", toassign.corr, "Incorrect", toassign.incorr))
    if(class.pairs[index] && toassign.corr){
      class.elem[index] <- class.elem[index] + 1
      toassign.corr <- toassign.corr - 1
    } else if(!class.pairs[index] && toassign.incorr){
      class.elem[index] <- class.elem[index] + 1
      toassign.incorr <- toassign.incorr - 1
    }
    index <- (index %% K.u) + 1
  }
  return(class.elem)
}

#### assign.noise ##########################################
# Assign elements from classes to noise clusters
#
# OUTPUT
# class.elem - matrix of size (K.n x C) where the (i,j)'th 
#       entry is the number of elements of the j'th class 
#       that are in the i'th noise cluster
assign.noise <- function(){
  # number of elements in each class that are incorrectly assigned 
  # to a noise cluster
  num.class.incorr <- n.C * (epsilon2) 
  
  class.elem <- matrix(rep(floor(n.C * epsilon2 / K.n), K.n * C), 
                       nrow = K.n,
                       ncol = C)
  
  toassign = (n.C * epsilon2) %% K.n
  class.elem <- apply(class.elem, 2, assign.noise.help, toassign=toassign)
  return(class.elem)
}

#### assign.noise.help ######################################
# Helper function for assigning elements to noise clusters
#
# INPUT:
# class.elem - vector of length K.n where the i'th entry is
#       the number of elements of the class in question that are
#       currently in the i'th noise cluster
# toassign - number of extra elements that need to be assigned
# OUTPUT:
# class.elem - vector of length K.n where the i'th entry 
#       is the updated number of elements of the class in question 
#       that are in the i'th noise cluster
assign.noise.help <- function(class.elem, toassign){
  indices <- sample(1:K.n, toassign)
  cur <- class.elem[1]
  class.elem[indices] <- cur+1
  return(class.elem)
}

#### Metric Helper Functions ####################################

#### get.ccm ####################################################
# Calculate the true and predicted Co-clustering matrices given 
# data as a matrix of cluster-class pairs with the number of 
# elements in each such pair
#
# INPUT:
#     data - data from the "create.data" function for a given
#         parameter setting
# OUTPUT (list of values):
#     ccm.t - the true (classes) co-clustering matrix
#     ccm.p - the predicted (clusters) co-clustering matrix
get.ccm <- function(data){
  k <- dim(data)[[1]]
  c <- dim(data)[[2]]
  n.c <- sum(data[,1])
  ccm.t <- ccm.p <- data.frame(matrix(0,nrow=n.c*c,ncol=n.c*c))
  
  # Assign values to ccm.t
  for(i in 1:c){
    ccm.t[(((i-1)*n.c)+1):(i*n.c),(((i-1)*n.c)+1):(i*n.c)] <- 1
  }
  
  # Assign values to ccm.p
  # Need to figure out which elements belong to which clusters
  # By default, elements in K1 from a given class are the first
  # elements in the matrix from that class, K2 are next, etc.
  for(i in 1:k){
    # vector of indices for elements in the cluster
    k.elem <- c()
    for(j in 1:c){ # class index
      if(data[i,j] > 0){
        if(i > 1){
          start <- sum(data[1:(i-1),j]) + 1
        } else{
          start <- 1
        }
        k.elem <- c(k.elem, 
                    ((j-1) * n.c + start):
                      ((j-1) * n.c + start + data[i,j] - 1))
      }
    }
    ccm.p[k.elem,k.elem] <- 1
  }
  
  return(list(ccm.t=ccm.t, ccm.p=ccm.p))
}

#### tp #################################################
# Returns the number of True Positives for the given clustering 
# assignment
#
# INPUT:
# data - matrix of size (K.u+K.n x C) that gives the number of
#     elements of each class that are assigned to each cluster
# OUTPUT:
# tp - # of True Positives
tp <- function(data){
  n.elem <- sum(data)
  # number of pairs of elements that are correctly assigned
  # to the same cluster
  tp <- (sum(apply(data, c(1,2), function(x){x*(x-1)})) / 2) + n.elem
  
  return(tp)
}

#### fp #################################################
# Returns the number of False Positives for the given clustering 
# assignment
#
# INPUT:
# data - matrix of size (K.u+K.n x C) that gives the number of
#     elements of each class that are assigned to each cluster
# OUTPUT:
# fp - # of False Positives
fp <- function(data){
  # find the total number of pairs of elements that are clustered together
  total.pos <- apply(data, 1, sum)
  #total.pos <- total.pos * (total.pos-1)
  total.pos <- sum(total.pos^2)
  
  fp <- total.pos - tp(data)
  return(fp)
}

#### tn #################################################
# Returns the number of True Negatives for the given clustering 
# assignment
#
# INPUT:
# data - matrix of size (K.u+K.n x C) that gives the number of
#     elements of each class that are assigned to each cluster
# OUTPUT:
# tn - # of True Negatives
tn <- function(data){
  n.elem <- sum(data) # total number of elements
  # find the total number of pairs of elements that are not clustered together
  total.neg <- apply(data, 1, sum)
  total.neg <- total.neg * (n.elem - total.neg)
  total.neg <- sum(total.neg)
  
  tn <- total.neg - fn(data)
  return(tn)
}

#### fn  ##################################################
# Returns the number of False Negatives for the given clustering 
# assignment
#
# INPUT:
# data - matrix of size (K.u+K.n x C) that gives the number of
#     elements of each class that are assigned to each cluster
# OUTPUT:
# fn - number of True Negatives
fn <- function(data){
  # number of pairs of elements that are incorrectly assigned
  # to the same clusters
  fn <- sum(apply(data, 1, fn.help)) / 2
  return(fn)
}

# Takes in a row vector associated with a cluster and returns the
# number of False Negatives associated with that cluster
fn.help <- function(row){
  #number of elements in the given cluster
  n.cluster <- sum(row)
  # returns (elements both in the class and cluster) * [(total # elements) 
  #         - (# elements in the cluster) - (# elements in the class outside the cluster)]
  fn <- sum(sapply(row,function(x){x * (n.cluster-x)}))
  return(fn)
}

#### tpr  ###############################################
# Returns the True Positive Rate for the given clustering 
# assignment
#
# INPUT:
# data - matrix of size (K.u+K.n x C) that gives the number of
#     elements of each class that are assigned to each cluster
# OUTPUT:
# tpr - True Positive Rate

tpr <- function(data){
  # NOTE: notation: lower case => local variable
  c <- dim(data)[[2]]
  n.c <- sum(data[,1])
  
  # number of pairs of elements that should be in the 
  # same cluster
  tp.expected <- (n.c*(n.c-1)*c) / 2
  tpr <- tp(data) / tp.expected
  return(tpr)
}

#### tnr  ###############################################
# Returns the True Negative Rate for the given clustering 
# assignment
#
# INPUT:
# data - matrix of size (K.u+K.n x C) that gives the number of
#     elements of each class that are assigned to each cluster
# OUTPUT:
# tnr - True Negative Rate
tnr <- function(data){
  # NOTE: notation: lower case => local variable
  c <- dim(data)[[2]]
  n.c <- sum(data[,1])
  
  # number of pairs of elements that should be assigned to
  # different clusters
  tn.expected <- (n.c * (n.c * (c-1)) * c) / 2
  tnr <- 1 - (fn(data) / tn.expected)
  return(tnr)
}

test.helper <- function(data){
  res <- c(tp(data), tn(data), fp(data), fn(data))
  names <- c("True Positive", "True Negative", "False Positive", "False Negative")
  
  print(paste(names, "=", res, sep=" "))
  
  n.elem <- sum(data)
  
  print(paste("Total is:", sum(res)))
  print(paste("Total should be:", n.elem^2))
}

#### Data Output Helper Functions ###############################

#### print.sim ##################################################
# Print the details of the current simulation
print.sim <- function(){
  writeLines(paste("Current Simulation Information",
                   "\n# of classes:", C,
                   "\n# of useful clusters:", K.u, 
                   "\n# of noise clusters:", K.n,
                   "\nUseful class error:", epsilon1,
                   "\nNoise class error:", epsilon2,
                   "\n# of elements in each class:", n.C,
                   "\n# of simulation iterations:", n.iter,
                   sep="\t"))
}

#### get.param.str #############################################
# Take in a vector of parameter values, some of which can be NAs
# and return a string with the values for the non-NA parameters
#
# INPUT:
#     param.vals - vector of length 4 of parameter values, 
#         some of which can be NA. If none is specified take the
#         current simulation parameter values
# OTUPUT:
#     str - formatted string with non-NA param values
get.param.str <- function(param.vals=NULL){
  if(is.null(param.vals)){
    param.vals = c(K.u, K.n, epsilon1, epsilon2)
  }
  if(length(param.vals) != length(sim.params)){
    stop("Vector of parameter values is of the wrong length")
  }
  str <- paste(sim.params, param.vals, sep="=", collapse="_")
  return(str)
}

#### evaluate #################################################
# Evaluate the clustering assignment given using each metric.
# Returns a vector with the results from each clustering metric.
#
# INPUT:
# data - clustering assignment data - from generate.data()
# title - optional values for the title of the current simulation
# OUTPUT:
#  sim.res - vector with one entry for each defined clustering
#     metric that contains the result of evaluating the clustering
#     assignment using that metric
evaluate <- function(data, title=""){
  pv <- pseudo.v.metric(data,0.01)
  sim.res <- c(combinatorial.metric(1,1,data), 
               combinatorial.metric(2,1,data), 
               1 - (pv$normal / 4000), # scale these so that they can be compared to other metrics
               1 - (pv$sym / 4000),
               mcc.metric(data),
               pearson.metric(data))
  sim.names <- c("Combinatorial.Equal.Weight", 
                 "Combinatorial.Unequal.Weight", 
                 "Pseudo.V", 
                 "Pseudo.V.Sym",
                 "MCC",
                 "Pearson")
  if(title != ""){
    sim.res <- c(title, sim.res)
    names(sim.res) <- c("Name", sim.names)
  } else{
    names(sim.res) <- sim.names
  }
  return(sim.res)
}

#### SIMLUATION DEFINITION ####################################

#### create. data #############################################
# Create the matrix of clustering data based on the global
# simulation variables
# OUTPUT:
#    clusters.elem - (K.u + K.n) x C matrix whose i,j'th entry 
#         is the number of elements in cluster i from class j
#         where the first K.u clusters are the useful clusters
#         and the last K.n clusters are the noise clusters
create.data <- function(n.iter.value = NULL,
                        C.value = NULL,
                        n.C.value = NULL,
                        K.u.value = NULL,
                        K.n.value = NULL,
                        e1.value = NULL,
                        e2.value = NULL){
  # assign parameter values if given
  if(!is.null(n.iter.value)){
    n.iter <<- n.iter.value
  }
  if(!is.null(C.value)){
    C <<- C.value
  }
  if(!is.null(n.C.value)){
    n.C <<- n.C.value
  }
  if(!is.null(K.u.value)){
    K.u <<- K.u.value
  }
  if(!is.null(K.n.value)){
    K.n <<- K.n.value
  }
  if(!is.null(e1.value)){
    epsilon1 <<- e1.value
    epsilon <<- epsilon1 + epsilon2
  }
  if(!is.null(e2.value)){
    epsilon2 <<- e2.value
    epsilon <<- epsilon1 + epsilon2
  }
  
  # Safety checks for when parameters are 0
  if(K.n == 0){
    epsilon2 <<- 0
    epsilon <<- epsilon1
  }
  if(K.u == 0){
    epsilon1 <<- 0
    epsilon2 <<- 1
    epsilon <<- 1
  }
  
  # useful cluster/class pair assignments
  #     kc.pairs[i,j] = 1 => i and j are paired
  #     kc.pairs[i,j] = 0 => i and j are not paired
  # make sure that each cluster and class is paired with at 
  # least one class or cluster respectively
  kc.pairs <- matrix(rep(0, K.u*C), nrow = K.u, ncol = C)
  if(K.u < C){
    for(j in 1:C){
      i <- (j %% K.u) + 1
      kc.pairs[i,j] <- 1
    }
  } else {
    for(i in 1:K.u){
      j <- (i %% C) + 1
      kc.pairs[i,j] <- 1
    }
  }
  
  # Number of elements from each class that are assigned to
  # each cluster
  clusters.elem <- data.frame(matrix(nrow = K.u + K.n,
                                     ncol = C,
                                     # name the rows and columns according to names in V-measure paper
                                     dimnames = list(c(sapply(1:K.u,function(x){paste("Ku",x,sep="")}),
                                                       if(K.n > 0){sapply(1:K.n,function(x){paste("Kn",x,sep="")})}
                                     ),
                                     c(sapply(1:C, function(x){paste("C",x,sep="")})
                                     )
                                     )
  )
  )
  # Elements assigned to useful clusters
  clusters.elem[1:K.u,] <- apply(kc.pairs, 2, assign.useful)
  # Elements assigned to noise clusters
  if(K.n > 0){
    clusters.elem[(K.u+1):(K.u+K.n),] <- assign.noise()
  }
  return(clusters.elem)
}


#### run.sim ####################################################
# Run the simulation over the given parameter values
#
# INPUT:
#    K.u.values - values of K.u to iterate over for the simulation
#    K.n.values - values of K.n to iterate over for the simulation
#    e1.values - values of epsilon1 to iterate over for the simulation
#    e2.values - values of epsilon2 to iterate over for the simulation
# OUTPUT:
#    res - nested list of lists with each dimension corresponding to
#       one of the iterated parameters from the input

run.sim <- function(K.u.values=2:11, 
                    K.n.values=0:6, 
                    e1.values=c(0,0.033,0.066,0.1), 
                    e2.values=c(0,0.066,0.133,0.2)){
  res <- list()
  res <- lapply(K.u.values,
                function(k.u){
                  res2 <- lapply(K.n.values,
                         function(k.n){
                           res3 <- lapply(e1.values,
                                  function(e1){
                                    res4 <- lapply(e2.values,
                                           function(e2){
                                             #print.sim()
                                             sim.res.tot <- sapply(1:n.iter, 
                                                                   function(x){
                                                                     print(get.param.str())
                                                                     data <- create.data(,,,k.u,k.n,e1,e2)
                                                                     eval <- evaluate(data)
                                                                     print("...Complete")
                                                                     return(eval)
                                                                   }
                                             )
                                             gc()
                                             apply(data.frame(sim.res.tot), 1, mean)
                                             
                                           })
                                    names(res4) <- e2.values
                                    res4
                                  })
                           names(res3) <- e1.values
                           res3
                         })
                  names(res2) <- K.n.values
                  res2
                })
  names(res) <- K.u.values
  return(res)
}

#### DATA OUTPUT/ANALYSIS ###################################

#### parse.data #############################################
# Extract simulation data with the given parameter values 
# fixed. If a parameter value is not specified then all values
# of that parameter are returned.
# INPUT:
#     res - results of the simulation
#     K.u.value - value for K.u
#     K.n.value - value for K.n
#     e1.value - value for epsilon1
#     e2.value - value for epsilon2
# OUTPUT:
filter.data <- function(res,
                       K.u.value=NA, 
                       K.n.value=NA, 
                       e1.value=NA, 
                       e2.value=NA){
  if(!is.na(K.u.value)){
    res <- parse.help(res, K.u.value, "K.u")
  }
  if(!is.na(K.n.value)){
    for(val1 in names(res)){
      res[[val1]] <- parse.help(res[[val1]], K.n.value, "K.n")
    }
  }
  if(!is.na(e1.value)){
    for(val1 in names(res)){
      res1 <- res[[val1]]
      for(val2 in names(res1)){
        res1[[val2]] <- parse.help(res1[[val2]], e1.value, "epsilon1")
      }
      res[[val1]] <- res1
    }
  }
  if(!is.na(e2.value)){
    for(val1 in names(res)){
      res1 <- res[[val1]]
      for(val2 in names(res1)){
        res2 <- res1[[val2]]
        for(val3 in names(res2)){
          res2[[val3]]<- parse.help(res2[[val3]], e2.value, "epsilon2")
        }
        res1[[val2]]<- res2
      }
      res[[val1]] <- res1
    }
  }
  return(res)
}

# Helper function for parsing data
# INPUT:
#    res - nested list of lists of simulation results, 
#       may be lower dimension than original simulation results
#    value - value to look set for the given parameter
#    param - name of the parameter
# OUTPUT:
#    res - nested list of lists of simulation results including
#       only results where the given parameter has the given value
filter.help <- function(res, value, param){
  value.str <- toString(value)
  if(is.null(res[[value.str]])){
    warning(paste("Value ", value.str, " is not a valid index for ", param, 
                    "\nReturning all values of ", param, " instead"))
  } else{
    temp <- res[[value.str]]
    res <- list()
    res[[value.str]] <- temp
  }
  return(res)
}


#### parse.data ############################################
# Parse simulation results to extract the dat in a usable
# form for a scatterplot. 
#
# INPUT:
#    res - simulation results
#    params.fixed.values - vector of length 4 of parameter values that 
#       are kept fixed throughout the simulation. Any parameters 
#       that are set to be NA will be iterated over (at most one 
#       parameter can be NA at a time).
#       Ordering of parameters: K.u, K.n, epsilon1, epsilon2
#    metrics - Optional vector of names of metrics to include in
#       the plot. If not included then all metrics are considered
# OUTPUT:
#    param.data - data frame containing the simulation results
#        where all but on of the parameters are fixed

parse.data <- function(res, params.fixed.values, metrics=NA){
  # set names for easier access
  names(params.fixed.values) <- c("K.u", "K.n", "epsilon1", "epsilon2")
  # index of the unfixed parameter
  param.ind <- match(NA, params.fixed.values)
    
  if(length(param.ind) > 1){
    stop("Too many parameter values are unspecified. Enter more parameters and try again")
  }
  # name of the unfixed parameter
  param.name <- names(params.fixed.values)[param.ind]

  # get the results we care about
  res <- filter.data(res, 
                    params.fixed.values[1], 
                    params.fixed.values[2], 
                    params.fixed.values[3], 
                    params.fixed.values[4])
  
  # values for the unfixed parameter
  param.values <- switch(toString(param.ind),
                         '1' = names(res),
                         '2' = names(res[[1]]),
                         '3' = names(res[[1]][[1]]),
                         '4' = names(res[[1]][[1]][[1]])
  )
  
  # names of the different metrics being considered
  if(is.na(metrics)){
    metrics <- names(res[[1]][[1]][[1]][[1]])
  }
  
  # data frame for the ploting data
  # will be populated later
  param.data <- data.frame(matrix(nrow=length(metrics)*length(param.values),
                                  ncol=3))
  for(i in 1:length(param.values)){
    val <- param.values[[i]]
    params.values <- params.fixed.values
    params.values[[param.name]] <- val
      
    param.data[((i-1)*length(metrics)+1):(i*length(metrics)),] <- cbind(rep(as.numeric(val), length(metrics)),
                                                                      res[[params.values]][metrics],
                                                                      metrics)
    
  }
  names(param.data) <- c(param.name, "value", "metric")
  param.data[,'value'] <- as.numeric(param.data[,'value'])
  
  return(param.data)
}

#### parse.data.diff ##########################################
# Parse simulation results to extract the data in a usable
# form. Returns the difference between the metric scores for
# consecutive values of the unfixed parameter
#
# INPUT:
#    res - simulation results
#    params.fixed.values - vector of length 4 of parameter values that 
#       are kept fixed throughout the simulation. Any parameters 
#       that are set to be NA will be iterated over (at most one 
#       parameter can be NA at a time).
#       Ordering of parameters: K.u, K.n, epsilon1, epsilon2
#    metrics - Optional vector of names of metrics to include in
#       the plot. If not included then all metrics are considered
# OUTPUT:
#    param.data - data frame containing the difference in metric 
#        scores for consecutive values of the given unfixed parameter
#        where all other parameters are fixed
parse.data.diff <- function(res, params.fixed.values, metrics=NA){
  parsed.data <- parse.data(res, params.fixed.values, metrics)
  if(dim(parsed.data)[[1]] < 2){
    stop("Too few data values to perform differencing. Try using parse.data instead")
  }
  
  # number of metrics
  if(is.na(metrics)){
    metrics <- unique(parsed.data[,3]) # list of metrics being looked at
  }
  n.metrics <- length(metrics) # number of metrics
  
  param.data <- data.frame(
    matrix(nrow=(dim(parsed.data)[[1]]-n.metrics),
           ncol=dim(parsed.data)[[2]])
  )
  names(param.data) <- c(paste(names(parsed.data)[[1]], ".diff", sep=""), "value", "metric")
  
  param.data[,3] <- parsed.data[1:n.metrics,3]
    
  # a bit convoluted but relying on the structure of the parsed data
  # essentially we are subtracting consecutive "blocks" of n.metrics rows
  # from each other to get the differences, and then figuring out the names
  for(i in 1:((dim(param.data)[[1]] / n.metrics))){
    # factor level, the two parameter values that are differenced
    lvl <- paste(parsed.data[(i*n.metrics),1], 
                 "->", 
                 parsed.data[((i+1)*n.metrics),1])
    
    # differenced values
    values <- parsed.data[((i-1)*n.metrics+1):(i*n.metrics),2] - parsed.data[(i*n.metrics+1):((i+1)*n.metrics),2]
    
    param.data[((i-1)*n.metrics+1):(i*n.metrics),-3] <- cbind(lvl, values)
  }
  param.data[,2] <- as.numeric(param.data[,2])
  
  return(param.data)
}


#### plot.sim ###############################################
# Plot the changes in one or more scoring metrics as a given
# parameter changes
#
# INPUT:
#    res - simulation results
#    params.fixed.values - vector of length 4 of parameter values that 
#       are kept fixed throughout the simulation. Any parameters 
#       that are set to be NA will be iterated over (at most one 
#       parameter can be NA at a time).
#       Ordering of parameters: K.u, K.n, epsilon1, epsilon2
#    metrics - Optional vector of names of metrics to include in
#       the plot. If not included then all metrics are considered
#    diff - if true then plot the difference of the metric scores 
#       evaluated for consecutive values of the given parameter, 
#       otherwise plot the metric scores themselves  
# OUTPUT:
#    plot - plot of the given data
plot.sim <- function(res, params.fixed.values, metrics=NA, diff=FALSE, filename=NULL){
  # get the data for plotting
  if(diff){
    plot.data <- parse.data.diff(res, params.fixed.values, metrics)
      
    # need to use factors as independent variable instead of the actual 
    # metric values if using the difference
    param.formula <- paste("factor(", names(plot.data)[[1]],")", sep="")
    ylab <- "Metric Score Difference"
    
    param.name <- paste(substr(names(plot.data)[[1]], 0, nchar(names(plot.data)[[1]]) - 5),
                        "Change") # name of parameter being changed
    title <- paste("Effect of", 
                   param.name, 
                   "on Changes in Scoring Metric Value")
    xrot <- 45 # rotation of xaxis labels
  } else{
    plot.data <- parse.data(res, params.fixed.values, metrics)
      
    param.formula <- paste("factor(", names(plot.data)[[1]],")", sep="")
    ylab <- "Metric Score"
    param.name <- names(plot.data)[[1]] # name of parameter being changed
    title <- paste("Effect of", 
                   param.name, 
                   "on Scoring Metric")
    xrot <- 0 # rotation of xaxis labels
  }
  
  print("Assigned data...")
  # name of the parameter that is not being held fixed and
  # the metrics that are being considered
  data.names <- names(plot.data)
  
  # formula for scatter plot
  plot.formula <- as.formula(paste("value ~", param.formula))
  print("Calculated formula...")
  
  # details of simulation for plot
  metrics <- unique(plot.data[,3])
  levels <- unique(plot.data[,1])
  descrip <- paste("Metrics:", 
                   metrics, 
                   " Parameters", 
                   get.param.str(params.fixed.values), 
                   collapse=",")
  print("Created description...")
  
  # calculate limits for y-axis
  ymin <- min(plot.data[,2])
  ymax <- max(plot.data[,2])
  yrange <- ymax - ymin
  if(yrange != 0){
    ylim <- c(min(ymin, max(ymin - yrange/4, -0.1)), max(ymax,min(ymax + yrange/4, 1.1)))
  } else {
    ylim = c(-0.1, 1.1)
  }
  print("Calculated y limit...")
  
  # assign graphical values
  p.pch <- 19
  p.colours <- default.colours(length(metrics))
  p.cex = 1
  
  # create legend
  metric.legend <- list(
    legend = list(
      pch = p.pch,
      cex = p.cex,
      colours = p.colours,
      labels = metrics,
      title = 'Metrics',
      border = 'transparent'
    ),
    legend = list(
      pch = p.pch,
      cex = p.cex,
      #colours = default.colours(length(params.fixed.values)),
      labels = paste(sim.params,as.character(params.fixed.values), sep="="),
      title = 'Parameters',
      border = 'transparent'
      )
  );
  
  metric.legend.grob <- legend.grob(
    legends = metric.legend
  );
  print("Created legends...")
  
  # If we are graphing the differences then a few modifications need to be made
  if(diff){
    
  }
  
  create.scatterplot(plot.formula, 
                     plot.data,
                     filename = filename,
                     group = plot.data$metric,
                     pch = p.pch,
                     col = p.colours,
                     cex = p.cex,
                     legend = list(right = list(fun = metric.legend.grob)),
                     abline.h = c(0,1),
                     abline.col = 'grey',
                     xlab.label = param.name,
                     ylab.label = ylab,
                     xaxis.cex = 1,
                     yaxis.cex = 1,
                     xaxis.rot = xrot,
                     xlab.cex = 1.5,
                     ylab.cex = 1.5,
                     main = title,
                     main.cex = 1.5,
                     description = descrip,
                     ylimits = ylim,
                     )
}

# Takes in a data frame with multiple metrics in different columns
# and returns a data frame with all the data in 2 columns
make.2d.help <- function(plot.data){
  # form the data into a data frame with only two columns
  plot.data.2d <- data.frame(matrix(ncol=2))
  names(plot.data.2d) <- c(names(plot.data)[[1]], "score")
  for(i in 2:dim(plot.data)[[2]]){
    temp <- plot.data[,c(1,i)]
    names(temp) <- names(plot.data.2d)
    plot.data.2d <- rbind(plot.data.2d, temp, make.row.names=F)
  }
  plot.data.2d <- plot.data.2d[-1,]
  return(plot.data.2d)
}

#### plot.sim.save #############################################
# Plot the changes in one or more scoring metrics as a given
# parameter changes and save to given filename as a png
#
# INPUT:
#    res - simulation results
#    params.fixed.values - vector of length 4 of parameter values that 
#       are kept fixed throughout the simulation. Any parameters 
#       that are set to be NA will be iterated over (at most one 
#       parameter can be NA at a time).
#       Ordering of parameters: K.u, K.n, epsilon1, epsilon2
#    metrics - Optional vector of names of metrics to include in
#       the plot. If not included then all metrics are considered
#    diff - if true then plot the difference of the metric scores 
#       evaluated for consecutive values of the given parameter, 
#       otherwise plot the metric scores themselves  
#     filename - filename for plot to be saved under
# OUTPUT:
#    plot - plot of the given data
plot.sim.save <- function(res, params.fixed.values, metrics, diff, filename){
  png(filename)
  plot.sim(res, params.fixed.values, metrics, diff)
  dev.off()
}

#### plot.sim.batch #############################################
# Create and save plots showing the effect of modifying each parameter
# individually. Other parameters are fixed for a given plot and are
# assigned values by sampling from the results values
#
#TODO: update this
# INPUT:
#    res - simulation results
#    diff - if true then plot the difference of the metric scores 
#       evaluated for consecutive values of the given parameter, 
#       otherwise plot the metric scores themselves  
#    filename - filename for plot to be saved under
#    iters - number of plots to save for each parameter
plot.sim.batch <- function(res, diff){
  # TODO: take these from 'res'
  K.u.values=2:8 
  K.n.values=0:6
  e1.values=c(0,0.033,0.066,0.1) 
  e2.values=c(0,0.066,0.133,0.2)
  
  for(param in 1:4){
    for(i in 1:10){
      param.values <- c(sample(K.u.values,1), 
                        sample(K.n.values,1),
                        sample(e1.values, 1),
                        sample(e2.values, 1))
      print(param.values)
      param.values[param] <- NA
      
      filename <- paste("./scoring_metric_data/sim_plots/scoringmetricplot_param=",param, "_",get.params.str(param.values), ".png", sep="")
      
      plot.sim(res, param.values,,diff, filename)
    }
  }  
}

#### EVALUATION METRICS #####################################

#### pearson.metric #########################################
# Evaluates the clustering assignment based on the pearson measure
#
# INPUT:
#   data - matrix of size (K.u+K.n x C) that gives the number of
#       elements of each class that are assigned to each cluster
# OTUPUT:
#   score - score for the clustering assignment, higher is worse
pearson.metric <- function(data){
  ccms <- get.ccm(data) # co-clustering matrices, both true and predicted
  ccm.t <- ccms$ccm.t
  ccm.p <- ccms$ccm.p
  
  uppertri <- upper.tri(ccm.t)
  
  return(cor(ccm.t[uppertri], ccm.p[uppertri], method="pearson"))
}

#### mcc.metric #############################################
# Evaluates the clustering assignment based on the MCC measure
#
# INPUT:
#   data - matrix of size (K.u+K.n x C) that gives the number of
#       elements of each class that are assigned to each cluster
# OTUPUT:
#   score - score for the clustering assignment, higher is worse
mcc.metric <- function(data){
  # Find True/False Negatives/Positives
  tp <- tp(data)
  tn <- tn(data)
  fp <- fp(data)
  fn <- fn(data)
  
  return((tp*tn - fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
}

#### pseudo.v.metric ########################################
# Evaluates the clustering assingment based on the pseudo 
# v-measure, using the K-L divergence between the true matrix
# (the classes) and the predicted matrix (the clusters).
# Can give the pseudo v-measure and/or the symmetrical pseudo
# v-measure.
#
# INPUT:
#   data - matrix of size (K.u+K.n x C) that gives the number of
#       elements of each class that are assigned to each cluster
#   epsilon - small, nonzero number to be used instead of zero
#       in the K-L divergence
#   version - vector of strings corresponding to which measures to
#       calculate
#       "normal" - regular pseudo v-measure
#       "sym" - symmetrical
# OTUPUT:
#   score - score for the clustering assignment, higher is worse
pseudo.v.metric <- function(data, epsilon, version=c("normal","sym")){
  ccms <- get.ccm(data) # co-clustering matrices, both true and predicted
  
  # add epsilon to everything, to ensure no zero values
  if(output){
    print("Replacing zeros...")
  }
  ccm.t <- ccms$ccm.t + epsilon
  ccm.p <- ccms$ccm.p + epsilon
  
  
  
  # row normalize
  if(output){
    print("Normalizing rows...")
  }
  ccm.t <- row.normalize.help(ccm.t)
  ccm.p <- row.normalize.help(ccm.p)
  
  if(output){
    print("Calculating divergence...")
  }
  
  score <- list()
  score$normal <- kl.divergence(ccm.t,ccm.p)
  score$sym <- (score$normal + kl.divergence(ccm.p,ccm.t)) / 2
  
  return(score)
}


# Row normalize a matrix so all the row sums are 1
#
# INPUT:
#     mat - matrix to be normalized
# OUTPUT:
#     mat.new - normalized matrix
row.normalize.help <- function(mat){
  t(apply(mat, 1, function(row){row / sum(row)}))
}

#### kl.divergence ########################################
# Calculate the K-L divergence between two matrices
#
# INPUT:
#     p - first matrix (reference)
#     q - second matrix (new)
# OUTPUT:
#     kl - K-L divergence
kl.divergence <- function(p,q){
  sum(p * (log(p) - log(q)))
}

#### combinatorial.metric ###################################
# Evaluates the clustering assignment based on a weighted sum of
# the True Positive Rate and the True Negative Rate for
# co-clustering of elements. Returns a normalized score between
# 0 and 1.
#
# INPUT:
#    t - weight for True Positive Rate
#    f - weight for True Negative Rate
#    data - matrix of size (K.u+K.n x C) that gives the number of
#       elements of each class that are assigned to each cluster
# OUTPUT:
#    score - score for the clustering assignment, 0 = worst score,
#       1 = best score

combinatorial.metric <- function(t, f, data){
  return((t * tpr(data) + f * tnr(data)) / (t + f))
}




