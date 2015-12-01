#### sc3_behavior.R #############################################################
# Analyze the behavior of different scoring metrics for Sub Challenge 3

#### PREAMBLE ###################################################################
library(BoutrosLab.plotting.general)

# Directories for tsv files and plots respectively
setwd("~/Documents/SMC-Het-Challenge/smc_het_eval")
tsv_dir = "scoring_metric_data/text_files/"
plot_dir = "scoring_metric_data/metric_behaviour_plots/"

# Lists and dictionaries with info on the metrics and the rankings
penalties <- c("abs", "sq")
method.names <- c("orig", "orig_no_cc", "pseudoV", "pseudoV_no_cc", "simpleKL_no_cc", 
                  "sqrt_no_cc", "sym_pseudoV_no_cc", "pearson_no_cc", "aupr_no_cc", "mcc_no_cc")

ordering.names <- c("Aggregate", "Paul", "Nathan", "Dave", "Adriana", "Peter", "Quaid")
method.is.good.greater <- list("orig"=T, "orig_no_cc"=T, "pseudoV"=F, "pseudoV_no_cc"=F, "simpleKL_no_cc"=F, 
                               "sqrt_no_cc"=T, "sym_pseudoV_no_cc"=F, "pearson_no_cc"=T, "aupr_no_cc"=T, "mcc_no_cc"=T)
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
#             Aggregate - the aggregate ranking of the cases (from good to bad)
#             Nathan, Paul, Dave, Adriana, Peter - an individuals ranking of the cases
#             Std Dev - the standard deviation of the ranking order across all participants
#             Case - alphabetical ordering of the cases
#             Metric - the metric score of the cases
#     display - logical for whether to print the plot or not
# OUPUT:
#     bp - barplot created
plot.SC3.all <- function(method="pseudoV", ordering="Aggregate", display=T){
  # All Cases SC3
  rank <- read.csv(file=paste(tsv_dir, "Rankings of phylogeny mistakes.csv", sep=""), sep=",", header=TRUE)
  d = read.csv(file=paste(tsv_dir, "scoring3A_all_cases_", method, ".tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
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
#     is.good.greater - logical denoting if the given metric gives good submissions 
#         larger scores (TRUE) or bad submissions larger scores (FALSE)
ordering.diff <- function(method="pseudoV", ordering="Aggregate", penalty="abs", is.good.greater=F){
  # All Cases SC3
  rank <- read.csv(file=paste(tsv_dir, "Rankings of phylogeny mistakes.csv", sep=""), sep=",", header=TRUE)
  d = read.csv(file=paste(tsv_dir, "scoring3A_all_cases_", method, ".tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  d <- merge(rank, d, by="Case")
  if(is.good.greater){
    actual.order <- order(d[,"Metric"])
  } else{
    actual.order <- order(-d[,"Metric"])
  }
  ideal.order <- order(d[,ordering])# order the columns based on the given value
  
  if(penalty=="abs"){
    diff <- sum(abs(actual.order - ideal.order))
  } else if(penalty == 'sq'){
    diff <- sqrt(sum((actual.order - ideal.order)^2))
  }
  return(diff)
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
  diff.tot <- diff.tot[,order(diff.tot["Aggregate",])]
  # data for bar plot
  xdata <- colnames(diff.tot)
  ydata <-  diff.tot["Aggregate",]
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
