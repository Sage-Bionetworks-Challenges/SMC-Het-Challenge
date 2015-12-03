#### sc3_behavior.R #############################################################
# Analyze the behavior of different scoring metrics for Sub Challenge 3

#### PREAMBLE ###################################################################
library(BoutrosLab.plotting.general)

# Directories for tsv files and plots respectively
setwd("~/Documents/SMC-Het/SMC-Het-Challenge/smc_het_eval")
tsv_dir = "scoring_metric_data/text_files/"
plot_dir = "scoring_metric_data/metric_behaviour_plots/"

# Lists and dictionaries with info on the metrics and the rankings
penalties <- c("abs", "sq", "spearman")
method.names <- list("orig", "orig_nc", "pseudoV", "pseudoV_nc", 
                  "sqrt_nc", "sym_pseudoV_nc", "pearson_nc", "aupr_nc", "mcc_nc", 
                  c('pseudoV_nc', 'mcc_nc', 'pearson_nc'),
                  c('pseudoV_nc', 'pearson_nc', 'sym_pseudoV_nc'),
                  c('aupr_nc', 'sqrt_nc', 'sym_pseudoV_nc'),
                  c('aupr_nc', 'sqrt_nc', 'sym_pseudoV_nc', 'pearson_nc'))

ordering.names <- c("Copeland", "Paul", "Nathan", "Dave", "Adriana", "Peter", "Quaid")
method.is.good.greater <- list("orig"=T, "orig_nc"=T, "pseudoV"=F, "pseudoV_nc"=F, "simpleKL_nc"=F, 
                               "sqrt_nc"=T, "sym_pseudoV_nc"=F, "pearson_nc"=T, "aupr_nc"=T, "mcc_nc"=T)
case.groups <- list("NCluster"=c("NClusterOneLineage", "NClusterTwoLineages", "NClusterCorrectLineage"),
                    "OneCluster"=c("OneCluster"),
                    "Phelogeny"=c("ParentIsNieceWithChildren", "ParentIsSiblingWithChildren", "ParentIsCousin", 
                                  "ParentIsAunt", "ParentIsGrandparent", "ParentIsSibling"),
                    "BigExtra" = c("BigExtraTop", "BigExtraMid", "BigExtraCurBot", "BigExtraNewBot"),
                    "SmallExtra" = c("SmallExtraTop", "SmallExtraMid", "SmallExtraNewBot", 'SmallExtraCurBot'),
                    "Split" = c("SplitClusterMidMultiChild", "SplitClusterMidOneChild", "SplitClusterBotSame", 
                                "SplitClusterBotDiff"),
                    "Merge" = c("MergeClusterTop&Mid", "MergeClusterMid&BotMultiChild", "MergeClusterMid&BotOneChild",
                                "MergeClusterBot"))
num.case <- 0
for(i in case.groups){num.case <- num.case + length(i)}
df.case.groups <- data.frame(matrix(nrow=num.case,ncol=2))
colnames(df.case.groups) <- c("Case", "Group")
i <- 1
for(group.name in names(case.groups)){
  group <- case.groups[group.name][[1]]
  df.case.groups[i:(i+length(group)-1), "Group"] <- group.name
  df.case.groups[i:(i+length(group)-1), "Case"] <- group
  i <- i + length(group)
}

##### plot.SC3.amit ############################################################
# SC 3 plots that are modelled after Amit's plots for SC 2
plot.SC3.amit <- function(){
  d = read.csv(file=paste(tsv_dir, "scoring3A_split_cases.tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  png(file=paste(plot_dir, "3A_Split_Cases_pseudoV.png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("3 Metric") +
      xlab("Case") + ggtitle("3A Split Cases") + coord_flip()#ylim=c(0.8,1))
  )
  dev.off()
  
  d = read.csv(file=paste(tsv_dir, "scoring3A_merge_cases.tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  png(file=paste(plot_dir, "3A_Merge_Cases_pseudoV.png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("3 Metric") +
      xlab("Case") + ggtitle("3A Merge Cases") + coord_flip()#ylim=c(0.8,1))
  )
  dev.off()
  
  d = read.csv(file=paste(tsv_dir, "scoring3A_parent_cases.tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  png(file=paste(plot_dir, "3A_Parent_Cases_pseudoV.png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("3 Metric") +
      xlab("Case") + ggtitle("3A Incorrect Parent Cases") + coord_flip()#ylim=c(0.8,1))
  )
  dev.off()
  
  d = read.csv(file=paste(tsv_dir, "scoring3A_other_cases.tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  png(file=paste(plot_dir, "3A_Other_Cases_pseudoV.png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("3 Metric") +
      xlab("Case") + ggtitle("3A Other Cases - Not Merge, Split or Incorrect Parent") + coord_flip()#ylim=c(0.75,1))
  )
  dev.off()
}

##### plot.SC3.all #############################################################
# TODO: make interactive plots?
# TODO: combine the plots
# Create a plot with the scoring metrics for all the cases for Sub Challenge 3 
# from a tsv file of results
#
# INPUT:
#     method - The scoring metric that was used to create the tsv data
#     ordering - What parameter to order the cases by. The cases will be sorted in
#         descending order according to the given parameter
#         Options:
#             Copeland - the aggregate ranking of the cases (from good to bad) using Coplend score
#             Nathan, Paul, Dave, Adriana, Peter - an individuals ranking of the cases
#             Std Dev - the standard deviation of the ranking order across all participants
#             Case - alphabetical ordering of the cases
#             Metric - the metric score of the cases
#     display - logical for whether to print the plot or not
#     pc - number of pseudo counts to use with the method
#         Options:
#           more - sqrt(n) pseudo counts where n = # of SSMs
#           less - log(n) pseudo counts
#           none - no pseudo counts
#     full - logical for whether to consider the full CCM, and AD matrix for the method (as opposed to just the upper triangular part)
# OUPUT:
#     bp - barplot created
plot.SC3.all <- function(method="pseudoV", ordering="Copeland", display=T, pc='more', m_size='full'){
  # All Cases SC3
  
  rank <- read.csv(file=paste(tsv_dir, "aggregate_scenario_rankings.csv", sep=""), sep=",", header=TRUE)
  d <- get.scoring.data(method=method, pc=pc, m_size=m_size)
  
  d <- merge(rank, d, by="Case")
  d <- merge(df.case.groups, d, by="Case")
  d <- d[order(d[,ordering]),] # order the columns based on the given value
  
  d$Case <- factor(d$Case, levels=unique(as.character(d$Case))) # change the case column to be factors
  if(method.is.good.greater[method][[1]]){ # flip data so that bigger is always worse
    upper <- ceiling(max(d$Metric))
    d$Metric <- upper - d$Metric
  }
  
  # data for barplot
  xdata <- d$Metric
  ydata <- as.factor(d$Case)
  
  
  xlims <- c(0.9*min(xdata), 1.1*max(xdata))
  
  bp <- create.barplot(
    ydata ~ xdata,
    d,
    main=paste(method,"Score for Different Clustering Mistakes"),
    main.cex=1.5,
    col=case.colours,
    yaxis.rot=30,
    yaxis.cex=0.8,
    xlab.label="Clustering Case/Mistake",
    xlab.cex=1,
    ylab.label="Metric Score",
    ylab.cex=1,
    xlimits=xlims,
    plot.horizontal=T,
    legend = list(
      inside = case.legend
      )
  )
  if(display){
    print(bp)
  }
  return(list(bp=bp, xlims=xlims))
}

# Helper function for reading in scoring data for all mistake scenarios 
# for a given method and parameter settings
# INPUT:
#     method - The scoring metric that was used to create the tsv data
#     ordering - What parameter to order the cases by. The cases will be sorted in
#         descending order according to the given parameter
#         Options:
#             Copeland - the aggregate ranking of the cases (from good to bad) using Coplend score
#             Nathan, Paul, Dave, Adriana, Peter - an individuals ranking of the cases
#             Std Dev - the standard deviation of the ranking order across all participants
#             Case - alphabetical ordering of the cases
#             Metric - the metric score of the cases
#     display - logical for whether to print the plot or not
#     pc - number of pseudo counts to use with the method
#         Options:
#           more - sqrt(n) pseudo counts where n = # of SSMs
#           less - log(n) pseudo counts
#           none - no pseudo counts
#     full - logical for whether to consider the full CCM, and AD matrix for the method (as opposed to just the upper triangular part)
# OUPUT:
get.scoring.data <- function(method="sym_pseudoV_nc", ordering="Copeland", display=T, pc='more', m_size='full'){
  # find the correct file extension
  if(pc == 'more'){ # pseudo count extension
    pc_ext = '_more_pc'
  } else if(pc == 'less'){
    pc_ext = ''
  } else{
    pc_ext = '_no_pc'
  }
  m_size_ext = paste('_', m_size, sep='')
  
  d = read.csv(file=paste(tsv_dir, "scoring3A_all_cases_", paste(method, collapse='_'), pc_ext, m_size_ext, ".tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  return(d)
}

#### plot.multi.SC3 ####################################################
# TODO: add labels to graphs
# Plot multiple metric ranking plots on the same figure to compare them
plot.multi.SC3 <- function(){
  info <- lapply(method.names, plot.SC3.all, display=F)
  plots <- lapply(info, function(x){x$bp})
  xlims <- lapply(info, function(x){as.vector(x$xlims)})
  
  odd.names <- method.names[seq(1,length(method.names), by=2)]
  even.names <- method.names[seq(2,length(method.names), by=2)]
  
  mp <- create.multiplot(
    plots,
    plot.layout=c(2, 5),
    xlab.label=c("Metric Score"),
    ylab.label= c("Cases"),
    xaxis.lab = NULL,
    yaxis.lab = NULL,
    x.relation='free',
    xlimits = xlims
  )
  print(mp)
  return(mp)
}


# Calculate the difference between the actual ordering and the ideal ordering
#
# INPUT:
#     method - string giving the scoring metric to use
#     ordering - string giving the data column to order by. This is the ordering that
#         the metric ranking will be compared to
#     penalty - penalty to use when comparing the two orderings
#         Options:
#           abs - absolute value of the difference in ranks bewteen the scoring metric
#             ordering and the desired ordering summed over all cases
#           sq - square of the difference in ranks bewteen the scoring metric ordering 
#             and the desired ordering summed over all cases and then square rooted
#           spearman - spearman rank correlation between the scoring metric ordering 
#             and the desired ordering
#     pc - number of pseudo counts to use with the method
#         Options:
#           more - sqrt(n) pseudo counts where n = # of SSMs
#           less - log(n) pseudo counts
#           none - no pseudo counts
#     full - size of matrices (for co-clustering and ancestor-descendant) to use with the method
#         Options:
#           full - full matrix
#           triu - upper triangular part of the matrix 
#     is.good.greater - logical denoting if the given metric gives good submissions 
#         larger scores (TRUE) or bad submissions larger scores (FALSE)
ordering.diff <- function(method="sym_pseudoV", ordering="Copeland", penalty="spearman", pc='more', m_size='full'){
  # All Cases SC3                   
  rank <- read.csv(file=paste(tsv_dir, "aggregate_scenario_rankings.csv", sep=""), sep=",", header=TRUE)
  d <- get.scoring.data(method=method, pc=pc, m_size=m_size)
  
  d <- merge(rank, d, by="Case")
  actual.order <- order(order(d[,"Metric"]))
  ideal.order <- order(order(d[,ordering]))# order the columns based on the given value
  
  if(penalty=="abs"){
    diff <- sum(abs(actual.order - ideal.order))
  } else if(penalty == 'sq'){
    diff <- sqrt(sum((actual.order - ideal.order)^2))
  } else if(penalty == 'spearman'){
    n <- length(actual.order)
    diff <- (6 * (sum((actual.order - ideal.order)^2) / (n * (n^2 - 1))))
  }
  return(diff)
}

# Calculate the difference between the actual ordering and the ideal ordering for
# a weighted average between 3 different metrics, using varying weights for each metric
#
# INPUT:
#     method - string giving the scoring metric to use
#     ordering - string giving the data column to order by. This is the ordering that
#         the metric ranking will be compared to
#     penalty - penalty to use when comparing the two orderings
#         Options:
#           abs - absolute value of the difference in ranks bewteen the scoring metric
#             ordering and the desired ordering summed over all cases
#           sq - square of the difference in ranks bewteen the scoring metric ordering 
#             and the desired ordering summed over all cases and then square rooted
#           spearman - spearman rank correlation between the scoring metric ordering 
#             and the desired ordering
#     is.good.greater - logical denoting if the given metric gives good submissions 
#         larger scores (TRUE) or bad submissions larger scores (FALSE)
ordering.diff.weights <- function(method=c("pseudoV_nc", "pearson_nc", "sym_pseudoV_nc"), ordering="Copeland", penalty="spearman"){
  # All Cases SC3
  rank <- read.csv(file=paste(tsv_dir, "aggregate_scenario_rankings.csv", sep=""), sep=",", header=TRUE)
  d = read.csv(file=paste(tsv_dir, "weights3A_all_cases_", paste(method, collapse="_"), ".tsv", sep=""), sep="\t",header=F)
  # make the header more readable
  header <- sapply(d[1,], as.character)
  colnames(d) <- header
  d <- d[-1,]
  
  weight.diff <- sapply(header[-1], function(col){
    d.wght <- d[,c('Case',col)]
    colnames(d.wght) <- c('Case', 'Metric')
    
    d.wght$Metric <- as.numeric(levels(d.wght$Metric))[d.wght$Metric] # change the metric column from a factor to a numeric value
    
    d.wght <- merge(rank, d.wght, by="Case")
    actual.order <- order(order(d.wght[,"Metric"]))
    ideal.order <- order(order(d.wght[,ordering]))# order the columns based on the given value
    
    if(penalty=="abs"){
      diff <- sum(abs(actual.order - ideal.order))
    } else if(penalty == 'sq'){
      diff <- sqrt(sum((actual.order - ideal.order)^2))
    } else if(penalty == 'spearman'){
      n <- length(actual.order)
      diff <- (6 * (sum((actual.order - ideal.order)^2) / (n * (n^2 - 1))))
    }
    print(as.numeric(levels(d[,col]))[d[,col]])
    print(c(col, diff))
    return(diff)
  })
  names(weight.diff) <- rev(header[-1])
  return(weight.diff)
}

#### main #########################################################################################
# Analyze the different metric behaviors
main <- function(){
  for(p in penalties){
    diff.tot <- data.frame(matrix(nrow=0, ncol=length(method.names)))
    
    for(o in ordering.names){
      diff <- sapply(method.names, ordering.diff,ordering=o, penalty=p)
      names(diff) <- method.names
      diff.tot <- rbind(diff.tot, diff)
    }
    colnames(diff.tot) <- method.names
    rownames(diff.tot) <- ordering.names
    
    plot.diff(diff.tot)
    write.table(x = diff.tot, file = paste(tsv_dir, "metric_diff_", p,".csv", sep=""), sep = ",", )
  }
}

#### plot.diff #################################################################################
# Plot the differences between the metric rankings and the desired rankings
#
# INPUT:
#   diff.tot - data frame containing the differences between a desired ranking and a metric
#       ranking for each desired ranking (each persons personal ranking and the aggregate ranking)
#       and each metric ranking.
# OUTPUT:
#   bp - barplot created from the given data
plot.diff <- function(diff.tot){
  diff.colours <- colour.gradient("red", dim(diff.tot)[2])
  diff.tot <- diff.tot[,order(diff.tot["Copeland",])]
  # data for bar plot
  xdata <- colnames(diff.tot)
  ydata <-  diff.tot["Copeland",]
  bp <- create.barplot(ydata ~ xdata, 
                 diff.tot,
                 main = "Difference Between Metric Rankings and Desired Rankings",
                 main.cex = 1.5,
                 xaxis.rot = 60,
                 xaxis.cex = 1,
                 yaxis.cex = 1,
                 xlab.label = "Metric",
                 xlab.cex = 1.5,
                 ylab.label = "Rank Difference",
                 ylab.cex = 1.5,
                 col=diff.colours,
                 sample.order = "increasing",
                 ylimits = c(0,1.1*max(ydata))
  )
  
  print(bp)
  return(bp)
}

#### plot.diff.pc #################################################################################
# Plot the differences between the metric rankings and the desired rankings for different pseudo count levels
#
# INPUT:
#     method - string giving the scoring metric to use
#     ordering - string giving the data column to order by. This is the ordering that
#         the metric ranking will be compared to
#     penalty - penalty to use when comparing the two orderings
#         Options:
#           abs - absolute value of the difference in ranks bewteen the scoring metric
#             ordering and the desired ordering summed over all cases
#           sq - square of the difference in ranks bewteen the scoring metric ordering 
#             and the desired ordering summed over all cases and then square rooted
#           spearman - spearman rank correlation between the scoring metric ordering 
#             and the desired ordering
# OUTPUT:
#   bp - barplot created from the given data
plot.diff.pc <- function(methods=method.names, ordering="Copeland", penalty="spearman"){
  pseudo_counts <- c('more', 'less', 'none') # possible pseudo count values
  m_sizes <- c('full', 'triu') # possible matrix sizes to use when evaluating SC3
  params <- expand.grid(m_sizes, pseudo_counts) # combine both into a data frame
  
  diff.tot <- lapply(methods,  
                    function(m){
                      res <- apply(params, 1,
                                   function(p){
                                     ordering.diff(m, pc=p[2], m_size=p[1])
                                   })
                      df <- data.frame(res, method=paste(m, collapse='+'), params=apply(params, 1, function(p){paste(p, collapse='.PC')}))
                      return(df)
                    })
  diff.tot <- do.call('rbind', diff.tot)

  diff.colours <- default.colours(dim(params)[1]) # barplot colours for each parameter setting
  diff.groups <- levels(unique(diff.tot$params))
  
  print(unique(diff.tot$params))
  legend <- list(
    right = list(
      fun = draw.key,
      args = list(
        key = list(
          points = list(
            col = 'black',
            pch = 22,
            cex = 2,
            fill = diff.colours
          ),
          text = list(
            lab = diff.groups
          ),
          padding.text = 3,
          cex = 1
        )
      ),
      x = 0.65,
      y = 0.95
    )
  )
  
  # data for bar plot
  xdata <- diff.tot$method
  ydata <-  diff.tot$res
  bp <- create.barplot(ydata ~ xdata, 
                       diff.tot,
                       groups=params,
                       main = "Difference Between Metric Rankings and Desired Rankings",
                       main.cex = 1.5,
                       xaxis.rot = 60,
                       xaxis.cex = 1,
                       yaxis.cex = 1,
                       xlab.label = "Metric",
                       xlab.cex = 1.5,
                       ylab.label = "Rank Difference",
                       ylab.cex = 1.5,
                       col=diff.colours,
                       legend = legend,
                       #sample.order = "increasing",
                       ylimits = c(0,1.1*max(ydata))
  )
  
  print(bp)
  return(bp)
}

#### plot.diff.weights #################################################################################
# Plot the differences between the metric rankings and the desired rankings for
# a weighted average between 3 different metrics, using varying weights for each metric
#
# INPUT:
#     methods - list of methods to consider (methods are either strings representing a single method or 
#         vectors of strings representing a combination of metrics)
#     ordering - string giving the data column to order by. This is the ordering that
#         the metric ranking will be compared to
#     penalty - penalty to use when comparing the two orderings
#         Options:
#           abs - absolute value of the difference in ranks bewteen the scoring metric
#             ordering and the desired ordering summed over all cases
#           sq - square of the difference in ranks bewteen the scoring metric ordering 
#             and the desired ordering summed over all cases and then square rooted
#           spearman - spearman rank correlation between the scoring metric ordering 
#             and the desired ordering
# OUTPUT:
#   bp - barplot created from the given data
plot.diff.weights <- function(methods=c("pseudoV_nc", "pearson_nc", "sym_pseudoV_nc"), ordering="Copeland", penalty="spearman"){
  diff.tot <- ordering.diff.weights(method=method, ordering=ordering, penalty=penalty)
  diff.colours <- colour.gradient("red", length(diff.tot))
  # data for bar plot
  xdata <- names(diff.tot)
  ydata <-  diff.tot
  d <- data.frame(xdata=xdata, ydata=ydata)
  bp <- create.barplot(ydata ~ xdata, 
                       data.frame(diff.tot),
                       main = "Difference Between Metric Rankings and Desired Rankings",
                       main.cex = 1.5,
                       xaxis.rot = 60,
                       xaxis.cex = 1,
                       yaxis.cex = 1,
                       xlab.label = "Metric",
                       xlab.cex = 1.5,
                       ylab.label = "Rank Difference",
                       ylab.cex = 1.5,
                       col=diff.colours,
                       sample.order = "increasing",
                       ylimits = c(0,1.1*max(ydata))
  )
  
  print(bp)
  return(bp)
}


