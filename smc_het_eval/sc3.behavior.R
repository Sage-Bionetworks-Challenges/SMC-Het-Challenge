#### sc3_behavior.R #############################################################
# Analyze the behavior of different scoring metrics for Sub Challenge 3

#### PREAMBLE ###################################################################
library(BoutrosLab.plotting.general)

# Directories for tsv files and plots respectively
setwd("~/Documents/SMC-Het/SMC-Het-Challenge/smc_het_eval")
tsv.dir = "scoring_metric_data/text_files/"
plot.dir = "scoring_metric_data/metric_behaviour_plots/"

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

seeds <- read.csv(file=paste(tsv.dir,'rnd_seeds.csv'))[,1] # 1000 seeds for performing CV runs with reproducibility 

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

write.all.scores <- function(){
  res <- sapply(method.names, function(m){
      d <- get.scoring.data(m)
      out <- d$Metric
      names(out) <- d$Case
      return(out)
    })
  colnames(res) <- sapply(method.names, format.name)
  write.table(res, file=paste(tsv.dir, "scoring3A_cases_ALL.tsv", sep=""), sep="\t")
}

##### plot.SC3.amit ############################################################
# SC 3 plots that are modelled after Amit's plots for SC 2
plot.SC3.amit <- function(){
  d = read.csv(file=paste(tsv.dir, "scoring3A_split_cases.tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  png(file=paste(plot.dir, "3A_Split_Cases_pseudoV.png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("3 Metric") +
      xlab("Case") + ggtitle("3A Split Cases") + coord_flip()#ylim=c(0.8,1))
  )
  dev.off()
  
  d = read.csv(file=paste(tsv.dir, "scoring3A_merge_cases.tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  png(file=paste(plot.dir, "3A_Merge_Cases_pseudoV.png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("3 Metric") +
      xlab("Case") + ggtitle("3A Merge Cases") + coord_flip()#ylim=c(0.8,1))
  )
  dev.off()
  
  d = read.csv(file=paste(tsv.dir, "scoring3A_parent_cases.tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  png(file=paste(plot.dir, "3A_Parent_Cases_pseudoV.png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("3 Metric") +
      xlab("Case") + ggtitle("3A Incorrect Parent Cases") + coord_flip()#ylim=c(0.8,1))
  )
  dev.off()
  
  d = read.csv(file=paste(tsv.dir, "scoring3A_other_cases.tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  png(file=paste(plot.dir, "3A_Other_Cases_pseudoV.png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("3 Metric") +
      xlab("Case") + ggtitle("3A Other Cases - Not Merge, Split or Incorrect Parent") + coord_flip()#ylim=c(0.75,1))
  )
  dev.off()
}

#### format.name ###########################################
# Format the name of a metric or a vector of metric names to 
# make it more readable for the plots.
# 
# INPUT:
#     name - string representing a metric name or vector of
#       multiple metric names
# OUTPUT:
#     name.out - formatted metric name
format.name <- function(name){
  if(name == 'orig'){return('orig (with CCM)')}
  if(name == 'pseudoV'){return('pseudoV (with CCM)')}
  if(name=='default'){name <- c('psuedoV_nc', 'mcc_nc', 'pearson_nc')}
  name.out <- gsub('_nc', '', name)
  name.out <- gsub('_', ' ', name.out)
  name.out <- paste(name.out, collapse=' + ')
  return(name.out)
}

##### get.ranking #############################################################
# Calculate the desired aggregate ranking of all the mistake scenarios given the 
# individual rankings from each of the Challenge organizers.
#
# INPUT:
#   scenarios.i - indices of the mistake scenarios/rows of the csv file to be included in the ranking. 
#     Indices should all be between 1 and 27
#     If NA then include all rows
#   aggr.meth - method to be used to aggregate the ranking scenarios
#     Options:
#       Copeland - Copeland's method using pairwise victories and pairwise defeats
#       Borda - Borda count which uses pairwise victories
# OUTPUT:
#   ranking - data frame containing two columns, the mistake scenario name and it's Copeland/Borda score
get.ranking <- function(scenarios.i=NA, aggr.meth='Copeland', filename="aggregate_scenario_rankings.csv"){
  rank <- read.csv(file=paste(tsv.dir, filename, sep=""), sep=",", header=TRUE)
  if(is.na(scenarios.i)){ # if scenarios.i is NA then include all rows/mistake scenarios
    scenarios.i <- 1:dim(rank)[1]
  }
  rank <- rank[scenarios.i,!(names(rank) %in% c('Copeland', 'Borda', 'Std.Dev', 'SD'))] # consider only the mistake scenarios specified
  
  colnames(rank)[1] <- 'Case' # TODO: remove this and name the cross validation tsv files correctly
  
  num.sc <- dim(rank)[1] # number of scenarios
  
  aggr.rank <- rep(0, num.sc) # aggregate score
  for(col in 2:dim(rank)[2]){
    for(i in 1:(num.sc-1)){
      for(j in (i+1):num.sc){
        if(rank[i, col] > rank[j, col]){ # if i is better than j
          aggr.rank[i] <- aggr.rank[i] + 1
          if(aggr.meth == 'Copeland'){
            aggr.rank[j] <- aggr.rank[j] - 1
          }
        } else if(rank[i, col] < rank[j, col]){ # if j is better than i
          aggr.rank[j] <- aggr.rank[j] + 1
          if(aggr.meth == 'Copeland'){
            aggr.rank[i] <- aggr.rank[i] - 1
          }
        } else if(aggr.meth == 'Borda'){ # if i and j are equally good
          aggr.rank[i] <- aggr.rank[i] + 0.5
          aggr.rank[j] <- aggr.rank[j] + 0.5
        }
      }
    }
  }
  ranking <- data.frame(rank$Case, aggr.rank)
  colnames(ranking) <- c('Case', aggr.meth)
  return(ranking)
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
#     m_size - size of matrices (for co-clustering and ancestor-descendant) to use with the method
#         Options:
#           full - full matrix
#           triu - upper triangular part of the matrix #     scenarios.i - indices of the mistake scenarios to be included in the ranking. 
#         Indices should all be betwee 1 and 27
# OUPUT:
#     bp - barplot created
plot.SC3.all <- function(method="pseudoV", ordering="Copeland", display=T, pc='more', m_size='full', scenarios.i=1:27){
  # All Cases SC3
  
  rank <- get.ranking(scenarios.i = scenarios.i, aggr.meth = ordering)
  d <- get.scoring.data(method=method, pc=pc, m_size=m_size)
  
  d <- merge(rank, d, by="Case")
  d <- merge(df.case.groups, d, by="Case")
  d <- d[order(d[,ordering]),] # order the columns based on the given value
  
  d$Case <- factor(d$Case, levels=unique(as.character(d$Case))) # change the case column to be factors
  d$Metric <- 1 - d$Metric
  
  # data for barplot
  xdata <- d$Metric
  ydata <- as.factor(sapply(d$Case, format.name))
  
  
  xlims <- c(0.9*min(xdata), 1.1*max(xdata))
  
  bp <- create.barplot(
    ydata ~ xdata,
    d,
    main=paste(method,"Score for Different Clustering Mistakes"),
    main.cex=1.5,
    col=case.colours,
    yaxis.rot=30,
    yaxis.cex=0.8,
    xlab.label="Metric Score",
    xlab.cex=1,
    ylab.label="Clustering Case/Mistake",
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
#     full - logical for whether to consider the full CCM, and AD matrix for the method 
#         (as opposed to just the upper triangular part)
#     scenarios.i - indices of the mistake scenarios to be included in the ranking. 
#         Indices should all be betwee 1 and 27
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
  
  d = read.csv(file=paste(tsv.dir, "scoring3A_all_cases_", paste(method, collapse='_'), pc_ext, m_size_ext, ".tsv", sep=""), sep="\t",header=FALSE)
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
#     m_size - size of matrices (for co-clustering and ancestor-descendant) to use with the method
#         Options:
#           full - full matrix
#           triu - upper triangular part of the matrix 
#     scenarios.i - indices of the mistake scenarios to be included in the ranking. 
#         Indices should all be betwee 1 and 27
ordering.diff <- function(method="sym_pseudoV", ordering="Copeland", penalty="spearman", pc='more', m_size='full', scenarios.i=1:27){
  # All Cases SC3 
  rank <- get.ranking(scenarios.i=scenarios.i, aggr.meth = ordering)
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
#     scenarios.i - indices of the mistake scenarios to be included in the ranking. 
#         Indices should all be betwee 1 and 27
ordering.diff.weights <- function(method=c("pseudoV_nc", "pearson_nc", "sym_pseudoV_nc"), ordering="Copeland", penalty="spearman", scenarios.i=1:27, display=F){
  # All Cases SC3
  rank <- get.ranking(scenarios.i = scenarios.i, aggr.meth = ordering)
  d = read.csv(file=paste(tsv.dir, "weights3A_all_cases_", paste(method, collapse="_"), ".tsv", sep=""), sep="\t",header=F)
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
    return(diff)
  })
  names(weight.diff) <- header[-1]
  return(weight.diff)
}

#### plot.diff #################################################################################
# Plot the differences between the metric rankings and the desired rankings
#
# INPUT:
#     methods - list of scoring metrics to use
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
#     display - boolean for whether to display the plot
#     scenarios.i - indices of the mistake scenarios to be included in the ranking. 
#       Indices should all be betwee 1 and 27
# OUTPUT:
#   bp - barplot created from the given data
plot.diff <- function(methods=method.names, ordering="Copeland", penalty="spearman", display=T, scenarios.i=1:27, return.plot=F){
  diff.tot <- sapply(methods, ordering.diff, ordering=ordering, penalty=penalty, scenarios.i=scenarios.i)
  names(diff.tot) <- sapply(methods, format.name)
  
  if(display){
    diff.colours <- colour.gradient("deepskyblue", length(diff.tot))
    diff.tot <- diff.tot[order(diff.tot)]
    
    diff.ylab <- list('spearman'='1 - Spearman Corr', 
                      'abs'='Absolute Diff', 
                      'sq'='Squared Diff')[penalty]
    diff.main <- "Metric Rankings vs. Desired Rankings"
    
    
    # data for bar plot
    ydata <-  diff.tot
    xdata <- names(ydata)
    
    d <- data.frame(xdata=xdata, ydata=ydata)
    bp <- create.barplot(ydata ~ xdata, 
                         data.frame(),
                         #filename = paste(plot.dir, 'Rank_Diff_All.png', sep=''),
                         main = diff.main,
                         width=11,
                         main.cex = 1.5,
                         xaxis.rot = 60,
                         xaxis.cex = 1,
                         yaxis.cex = 1,
                         xaxis.lab = NULL,
                         xlab.label = NULL,
                         xlab.cex = 1.5,
                         ylab.label = diff.ylab,
                         ylab.cex = 1.5,
                         col=diff.colours,
                         sample.order = "increasing",
                         ylimits = c(0,1.1*max(ydata))
    )
    
    #print(bp)
  }
  if(display & return.plot){
    return(bp)
  }
  return(diff.tot)
}

# Matrix of labels showing combinations of different metrics used in the plot.diff barplot
#
# INPUT:
#     methods - list of scoring metrics to use
# OUTPUT:
#     dm - dotmap indicating the metrics used for each bar in the plot.diff barplot
plot.diff.labels <- function(methods=method.names, ordering="Copeland", penalty="spearman", scenarios.i=1:27){
  diff.tot <- sapply(methods, ordering.diff, ordering=ordering, penalty=penalty, scenarios.i=scenarios.i)
  names(diff.tot) <- sapply(methods, format.name)
  
  # order the metrics according to their SCC
  diff.tot.ord <- order(diff.tot)
  diff.tot <- diff.tot[diff.tot.ord]
  xdata <- names(diff.tot)
  
  # find the list of unique metrics included in this list of methods
  methods <- sapply(unique(unlist(methods)), format.name)
  
  # calculate which individual metrics are used in each potential metric
  plot.data <- sapply(xdata, function(m){
    m <- strsplit(m, '[ ]+[+][ ]+')
    m <- m[[1]]
    res <- sapply(methods, function(x){as.numeric(x %in% m)})
    print(res)
    return(res)
  })
  
  # colouring function for the dots
  spot.fn <- function(x){ifelse(x == 1, 'black', 'white')} 
  
  # create the dotmap
  dm <- create.dotmap(
    plot.data,
    xaxis.lab = ' ' ,
    yaxis.lab = methods,
    yaxis.cex = 1.5,
    spot.colour.function=spot.fn
    )
  return(dm)
}

# plot.diff.multi <- function(methods=method.names, ordering="Copeland", penalty="spearman", scenarios.i=1:27){
#   bp <- plot.diff(methods=methods, ordering=ordering, penalty=penalty, scenarios.i=scenarios.i, display=T, return.plot=T)
#   dm <- plot.diff.labels(methods=methods, ordering=ordering, penalty=penalty, scenarios.i=scenarios.i)
#   plots <- list(dm, bp)
#   
#   mp <- create.multiplot(
#     plots,
#     retrieve.plot.labels)
#   
#   return(mp)
#   
# }
plot.diff.multi <- function(bp, dm, methods=method.names){
  plots <- list(dm, bp)
  methods <- rev(sapply(unique(unlist(methods)), format.name))
  
  mp <- create.multiplot(
    plots,
    yaxis.lab = list(methods,TRUE),
    yat = list(1:length(methods), c(0.05,0.10, 0.15)),
    ylimits = list(NULL, c(0,0.15)),
    xaxis.lab=NULL,
    yaxis.cex = list(1.5,3),
    panel.heights = c(3,1.3)
    )
  
  return(mp)
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
#     display - boolean for whether to display the plot
#     scenarios.i - indices of the mistake scenarios to be included in the ranking. 
#       Indices should all be betwee 1 and 27
# OUTPUT:
#   bp - barplot created from the given data
plot.diff.pc <- function(methods=method.names, ordering="Copeland", penalty="spearman", display=T, scenarios.i=1:27){
  pseudo_counts <- c('more', 'less', 'none') # possible pseudo count values
  m_sizes <- c('full', 'triu') # possible matrix sizes to use when evaluating SC3
  params <- expand.grid(m_sizes, pseudo_counts) # combine both into a data frame
  
  diff.tot <- lapply(methods,  
                    function(m){
                      res <- apply(params, 1,
                                   function(p){
                                     ordering.diff(m, pc=p[2], m_size=p[1], scenarios.i=scenarios.i)
                                   })
                      df <- data.frame(res, method=paste(m, collapse='+'), params=apply(params, 1, function(p){paste(p, collapse='.PC')}))
                      return(df)
                    })
  diff.tot <- do.call('rbind', diff.tot)
  
  if(display){
    diff.colours <- c('magenta', 'purple', 'darkmagenta', 'greenyellow', 'limegreen', 'green4') #default.colours(dim(params)[1]) # barplot colours for each parameter setting
    diff.groups <- levels(unique(diff.tot$params))
    
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
    xdata <- sapply(diff.tot$method, format.name)
    ydata <-  diff.tot$res
    bp <- create.barplot(ydata ~ xdata, 
                         diff.tot,
                         groups=params,
                         #filename = paste(plot.dir,'sc3_pc_msize_change.png', sep=''),
                         main = "Pseudo Counts & Matrix Size Changes",
                         main.cex = 1.5,
                         width = 10,
                         xaxis.rot = 30,
                         xaxis.cex = 1,
                         yaxis.cex = 1,
                         xlab.label = "Metric",
                         xlab.cex = 1.5,
                         ylab.label = "1 - Spearman",
                         ylab.cex = 1.5,
                         col=diff.colours,
                         legend = legend,
                         #sample.order = "increasing",
                         ylimits = c(0,1.1*max(ydata))
    )
    
    print(bp)
  }
  res <- diff.tot$res
  names(res) <- paste(diff.tot$method, diff.tot$params, sep='.')
  return(res)
}

#### plot.diff.input #################################################################################
# Plot the differences between the metric rankings and the desired rankings for the mistake scenarios 
# comparing the metrics using different input matrices to calculate the metric scores
#
# INPUT:
#     methods - list of strings giving the scoring metrics to use. 
#         Only use one metric at a time, no combination metrics
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
#     display - boolean for whether to display the plot
#     scenarios.i - indices of the mistake scenarios to be included in the ranking. 
#       Indices should all be betwee 1 and 27
# OUTPUT:
#   bp - barplot created from the given data
plot.diff.input <- function(methods=method.names, ordering="Copeland", penalty="spearman", display=T, scenarios.i=1:27){
  methods <- methods[sapply(methods, function(x){length(x) == 1})] # format the methods list to remove any combination metrics
  methods <- unique(sapply(methods, function(x){gsub('_nc', '', x)})) # format the methods list to remove duplicates
  
  
  diff.groups <- c('All', 'No_CCM', 'No_ADM', 'No_ADM_transpose', 'No_Cousin_Matrix')
  diff.tot <- lapply(methods,  
                     function(m){
                       res <- c(ordering.diff(paste(m,'_all', sep=''), scenarios.i=scenarios.i), 
                                ordering.diff(paste(m,'_nc', sep=''), scenarios.i=scenarios.i), 
                                ordering.diff(paste(m,'_na', sep=''), scenarios.i=scenarios.i), 
                                ordering.diff(paste(m,'_nat', sep=''), scenarios.i=scenarios.i), 
                                ordering.diff(paste(m,'_ncous', sep=''), scenarios.i=scenarios.i)
                                )
                       df <- data.frame(res, m, diff.groups)
                       return(df)
                     })
  diff.tot <- do.call('rbind', diff.tot)
  colnames(diff.tot) <- c('res', 'method', 'inmat')
  
  
  if(display){
    diff.colours <- default.colours(12)[8:12] # barplot colours for each parameter setting
    
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
              lab = sort(diff.groups)
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
    xdata <- sapply(paste(diff.tot$method, '_nc', sep=''), format.name) # format the metric names
    # need to add '_nc' to everything to avoid weird naming with PseudoV
    ydata <-  diff.tot$res
    bp <- create.barplot(ydata ~ xdata, 
                         diff.tot,
                         groups=inmat,
                         #filename = paste(plot.dir,'sc3_input_change.png', sep=''),
                         main = "Changing the Input Matrices for SC3 Scoring",
                         main.cex = 1.5,
                         width = 10,
                         xaxis.rot = 30,
                         xaxis.cex = 1,
                         yaxis.cex = 1,
                         xlab.label = "Metric",
                         xlab.cex = 1.5,
                         ylab.label = "1 - Spearman",
                         ylab.cex = 1.5,
                         col=diff.colours,
                         legend = legend,
                         #sample.order = "increasing",
                         ylimits = c(0,1.1*max(ydata))
    )
    
    print(bp)
  }
  res <- diff.tot$res
  names(res) <- paste(diff.tot$method, diff.tot$inmat, sep='.')
  return(res)
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
plot.diff.weights <- function(methods=c("pseudoV_nc", "pearson_nc", "sym_pseudoV_nc"), ordering="Copeland", penalty="spearman", display=T, scenarios.i=1:27){
  diff.tot <- ordering.diff.weights(method=methods, ordering=ordering, penalty=penalty)
  
  if(display){
    diff.colours <- rev(colour.gradient("deepskyblue", length(diff.tot)))
    diff.colours <- diff.colours[order(order(diff.tot))]
    
    diff.ylab <- list('spearman'='1 - Spearman Corr', 
                      'abs'='Absolute Diff', 
                      'sq'='Squared Diff')[penalty]
    diff.main <- paste("Metric Rankings vs. Desired Rankings - (", 
                       format.name(methods),
                       ")", sep="")
    
    # data for bar plot
    xdata <- names(diff.tot)
    ydata <-  diff.tot
    print(ydata)
    d <- data.frame(xdata=xdata, ydata=ydata)
    bp <- create.barplot(ydata ~ xdata, 
                         data.frame(diff.tot),
                         #filename = paste(plot.dir, 'Rank_Diff_Weighted_Avg.png', sep=''),
                         main = diff.main,
                         width=11,
                         main.cex = 1.5,
                         xaxis.lab = sapply(xdata, format.name),
                         xaxis.rot = 60,
                         xaxis.cex = 1,
                         yaxis.cex = 1,
                         xlab.label = "Metric",
                         xlab.cex = 1.5,
                         ylab.label = diff.ylab,
                         ylab.cex = 1.5,
                         col=diff.colours,
                         #sample.order = "increasing",
                         ylimits = c(0,1.1*max(ydata))
    )
    
    print(bp)
  }
  return(diff.tot)
}


##### cross.validate ##########################################################
# Cross validate the results of any of the tests performed on different subsets
# of the mistake scenarios.
#
# INPUT:
#   n.folds - number of cross validation folds to use
#   seed - random seed to use in assigning scenarios to fold (for reproducibility)
#     Default is NA, meaning no seed if set.
#   inc.all - boolean indicating whether the function should also be applied to 
#     all the mistake scenarios and the results outputted for comparison. If this 
#     flag is True then the mean and SD include the results from all the mistake 
#     scenarios together
#   inc.stats - boolean indicating whether to output the mean and standard deviation 
#     across all CV iterations
#   func - function to apply the cross validation to
#   fielname - filename to save the data to. If NA then data is not saved
#   ... - additional parameters for func
# OTUPUT:
#   res - data frame with the results from applying the test to each of the CV folds
cross.validate <- function(n.folds = 5, seed=NA, inc.all=T, inc.stats=T, func=ordering.diff, filename=NA, ...){
  if(!is.na(seed)){ # set the seed for randomization, if one is specified, for reproducibility
    set.seed(seed)
  }
  folds <- sapply(1:27, function(x){x %% n.folds}) + 1 # create the folds by randomizing fold assignments
  folds <- sample(folds) 
  print(folds)
  
  n.iter <- ifelse(inc.all, n.folds+1, n.folds) # if inc.all is True then include the results of func using all mistake scenarios
  # (go to n.folds+1 in that case since the last iteration will include all scenarios)
  
  res <- sapply(1:n.iter, function(i){ 
    sc.i <- which(folds != i)
    print(sc.i)
    func(scenarios.i = sc.i, ...)
  })
  
  if(is.null(dim(res))){ # for vectors/lists
    res.names <- paste('Fold', 1:n.folds, sep='.')
    if(inc.all){
      res.names <- c(res.names, 'All')
    }
    if(inc.stats){
      res <- c(res, mean(res), sd(res)) # add the mean and standard deviation of the results to the output
      res.names <- c(res.names, 'Mean', 'SD')
    }
    names(res) <- res.names
  } else{
    res.names <- paste('Fold', 1:n.folds, sep='.')
    if(inc.all){
      res.names <- c(res.names, 'All')
    }
    if(inc.stats){
      res <- cbind(res, apply(res, 1, mean), apply(res, 1, sd)) # add the mean and standard deviation of the results to the output
      res.names <- c(res.names, 'Mean', 'SD')
    }
    colnames(res) <- res.names
  }
  
  if(!is.na(filename)){
    write.csv(file=filename, res)
  }
  return(res)
}

funcl <- list('pc'=plot.diff.pc,
              'input'=plot.diff.input,
              'weight'=ordering.diff.weights,
              'normal'=plot.diff)


##### batch.cv ##########################################################
# Perfrom a batch run of cross validation tests. Calculate and return the 
# best metrics for each run. 
#
# INPUT:
#   n.folds - number of cross validation folds to use
#   seed - random seed to use in assigning scenarios to fold (for reproducibility)
#     Default is NA, meaning no seed if set.
#   inc.all - boolean indicating whether the function should also be applied to 
#     all the mistake scenarios and the results outputted for comparison. If this 
#     flag is True then the mean and SD include the results from all the mistake 
#     scenarios together
#   inc.stats - boolean indicating whether to output the mean and standard deviation 
#     across all CV iterations
#   func.name - name of function to apply the cross validation to (from the funcl list)
#   rewrite - boolean for whether to perform the cross validation steps again if they 
#     have already been performed and saved
#   ... - additional parameters for func
# OTUPUT:
#   res - data frame with the results from applying the test to each of the CV folds
batch.cv <- function(n.iter=2, n.folds = 5, inc.all=T, inc.stats=T, func.name='ordering', rewrite=F, ...){
  func <- funcl[[func.name]]
  results <- sapply(1:n.iter, function(x){
    seed <- seeds[x]
    filename <- paste(tsv.dir, 'sc3_cv_scores_', func.name, x, '.csv', sep='')
    
    if(!file.exists(filename) || rewrite){ # perform the CV if necessary/requested
      res <- cross.validate(n.folds = n.folds, inc.all=inc.all, inc.stats=inc.stats, func=func, filename=filename, ...)
    } else{
      res <- read.csv(filename, row.names=1)
    }
    col.to.ignore <- c('SD') # columns to exclude when calculating the best scenarios 
    if(!inc.all){ 
      col.to.ignore <- c(col.to.ignore, 'All')
    }
    if(!inc.stats){
      col.to.ignore <- c(col.to.ignore, 'Mean')
    }
    res <- res[,!(colnames(res) %in% col.to.ignore)] # exclude specified columns
    
    res.min.i <- apply(res, 2, function(x){which(x == min(x))}) # best scenario for each fold
    res.best <- rep(0, dim(res)[1])
    for(i in res.min.i){
      res.best[i] <- res.best[i] + 1
    }
    names(res.best) <- rownames(res)
    
    if(func.name == 'pc'){ # for the pseudo count metric also add in a measure of total 'wins' for each metric
      metrics <- unique(sapply(rownames(res), function(x){strsplit(x, '[.]')[[1]][[1]]}))
      res.best.metric <- rep(0, length(metrics))
      for(i in res.min.i){
        i.metric <- ceiling(i / (dim(res)[1] / length(metrics))) # find metric corresponding to 'winning' parameter setting
        res.best.metric[i.metric] <- res.best.metric[i.metric] + 1
      }
      res.best <- c(res.best, res.best.metric)
      names(res.best) <- c(rownames(res), metrics)
    }
    
    res.best <- res.best / dim(res)[2]
    
    return(res.best)
  })
  
  return(results)
}

cv.all <- function(n.iter=200){
  for(f.name in c('pc')){#names(funcl)){
    res <- batch.cv(n.iter=n.iter, func.name=f.name, rewrite=T, display=F)
    write.csv(file=paste(tsv.dir, 'sc3_cv_best_scenarios_', f.name, '.csv', sep=''), apply(res, 1, function(x){sum(x) / n.iter}))
  }
}

cv.ranking <- function(f.name = 'weight'){
  res <- sapply(1:200, 
                function(x){
                  outl <- get.ranking(filename=paste('sc3_cv_scores_', f.name, x, '.csv', sep=''))
                  out <- outl$Copeland
                  if(f.name == 'weight'){names(out) <- rev(outl$Case)}
                  else{names(out) <- outl$Case}
                  
                  return(out)
                })
  print(res)
  res <- apply(res, 1, sum)
  rank <- order(order(res))
  names(rank) <- rownames(res)
  return(rank)
}
