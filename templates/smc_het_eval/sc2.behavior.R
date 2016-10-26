#### sc2.behavior.R #############################################################
# Analyze the behavior of different scoring metrics for Sub Challenge 2

# Scoring metrics analyzed include: sum of squared error (original metric), 
#   square-root, pseudoV, symmetric pseudoV, spearman, pearson, 
#   Area Under the Precision Recall curve (AUPR), Matthews Correlation Coefficient (MCC)
# **Default metric is an average of MCC, Pearson and pseudoV**

# For more details on each metric look in SMCScoring.py

#### PREAMBLE ###################################################################
library(BoutrosLab.plotting.general)

# Directories for tsv files and plots respectively
tsv.dir = "scoring_metric_data/text_files/"
plot.dir = "scoring_metric_data/metric_behaviour_plots/"
# Directory that the script is running in (should be <SMC-Het-Challenge git repo>/smc_het_eval)
script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)

# Lists and dictionaries with info on the metrics and the rankings
method.names <- c("orig",
                  "sqrt",
                  "pseudoV",
                  "sym_pseudoV",
                  "spearman",
                  "pearson",
                  "aupr",
                  "mcc")

colour.scheme <- list('pseudoV'='gold', 
                      'mcc'='black', 
                      'pearson'='purple', 
                      'sym_pseudoV'='dodgerblue', 
                      'default'='limegreen', 
                      'aupr'='brown', 
                      'orig'='red', 
                      'spearman'='pink')

# Values to use for the standard deviation of the error in the co-clustering probabilities
std.values <- c(0, 0.01,0.03,0.05,0.1,0.15,0.2)

write.all.scores <- function(){
  res <- sapply(method.names, function(m){
      d <- read.csv(file=paste(tsv.dir, "scoring2A_cases_", m, ".tsv", sep=""), sep="\t", row.names=1, header=FALSE)
      out <- d$V2
      names(out) <- rownames(d)
      return(out)
    })
  write.table(res, file=paste(tsv.dir, "scoring2A_cases_ALL.tsv", sep=""), sep="\t")
}

write.all.scores.big <- function(){
  res <- sapply(method.names, function(m){
    d <- read.csv(file=paste(tsv.dir, "scoring2A_big_cases_", m, ".tsv", sep=""), sep="\t", row.names=1, header=FALSE)
    out <- d$V2
    names(out) <- rownames(d)
    return(out)
  })
  write.table(res, file=paste(tsv.dir, "scoring2A_big_cases_ALL.tsv", sep=""), sep="\t")
}

write.all.scores.rndreassign <- function(){
  res <- sapply(method.names, function(m){
    d <- read.csv(file=paste(tsv.dir, "scoring2A_random_reassignment_", m, ".tsv", sep=""), sep="\t", header=FALSE)
    out <- d$V2
    return(out)
  })
  reassign.prob <- read.csv(file=paste(tsv.dir, "scoring2A_random_reassignment_default.tsv", sep=""), sep="\t", header=FALSE)$V1
  reassign.prob.vals <- unique(reassign.prob)
  n.iter <- length(reassign.prob) / length(reassign.prob.vals)
  
  res <- cbind(1:n.iter, reassign.prob, res)
  mean.score <- t(sapply(reassign.prob.vals, function(x){
    apply(res[which(reassign.prob == x),], 2, mean)
  }))
  mean.score[,1] <- 0
  print(mean.score)
  res <- rbind(res, mean.score)
  
  write.table(res, file=paste(tsv.dir, "scoring2A_random_reassignment_ALL.tsv", sep=""), sep="\t", row.names=F)
}

write.all.scores.prob <- function(){
  res <- sapply(method.names, function(m){
    d <- read.table(file=paste(tsv.dir, "2B_prob_scoring_", m, ".tsv", sep=""), sep=",", row.names=1, header=TRUE)
    print(d)
    out <- d
    colnames(out) <- rownames(d)
    return(out)
  })
  print(res)
  write.table(res, file=paste(tsv.dir, "2B_prob_scoring_ALL.tsv", sep=""), sep="\t")
}

##### plot.SC2.amit ##########################################################################
# For Scoring Design Paper - Supplementary
# plot all the figures that Amit made. These include: 

plot.SC2.amit <- function(method='default'){
  # For Scoring Design Paper - Main Paper
  # Scoring behavior for different 'mistake scenarios' including
  #     - splitting one cluster into two
  #     - merging two clusters into one
  #     - assigning all mutations to the same cluster
  #     - assigning each mutation to its own cluster
  #     - adding an extra cluster with mutations from each true cluster (either small or large)
  d = read.csv(file=paste(tsv.dir, "scoring2A_cases_", method, ".tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  
  # All Cases SC3
  # data for barplot
  xdata <- d$Metric
  ydata <- as.factor(d$Case)
  
  plot.col <- colour.gradient(colour.scheme[method], length(xdata))#default.colours(length(xdata))
  
  xlims <- c(0.9*min(xdata), 1.1*max(xdata))
  
  bp.mistake.sc <- create.barplot(
    ydata ~ xdata,
    d,
    #filename = paste(plot.dir, 'SC2_Plot_', method, '.png', sep=''),
    main=paste("SC2 Scores -", format.name(method)),
    main.cex=1.5,
    col=plot.col,
    yaxis.rot=30,
    yaxis.cex=0.8,
    xlab.label="Metric Score",
    xlab.cex=1,
    ylab.label="Clustering Case/Mistake",
    ylab.cex=1,
    xlimits=xlims,
    plot.horizontal=T,
  )
  
  print(bp.mistake.sc)

#   # For Scoring Design Paper - Supplementary (maybe)
#   # Scoring behavior when predicted CCM is created by adding random zero-mean beta error
#   # to the true CCM matrix entries (with different values for the concentration parameter in the beta error) 
#   d = read.csv(paste(tsv.dir, "scoring2B_beta_", method, ".tsv", sep=""),sep="\t",header=FALSE)
#   colnames(d) = c("Error","Metric")
#   print(d)
#   
#   # All Cases SC3
#   # data for barplot
#   xdata <- as.ordered(d$Error)
#   ydata <- d$Metric
#   
#   plot.col <- default.colours(length(unique(xdata)))
#   
#   #xlims <- c(0.9*min(xdata), 1.1*max(xdata))
#   
#   bp.beta.err <- create.scatterplot(
#     Metric ~ as.ordered(Error),
#     d,
#     groups=Error,
#     cex=1,
#     #filename = paste(plot.dir, 'SC2_beta_err_', method, '.png', sep=''),
#     main=paste("SC2 w/ Err - ", format.name(method)),
#     main.cex=1.5,
#     col=plot.col,
#     xaxis.rot=30,
#     xaxis.cex=0.8,
#     yaxis.cex = 0.8,
#     ylab.label="Metric Score",
#     ylab.cex=1,
#     xlab.label="Concentration of Beta Err",
#     xlab.cex=1,
#     abline.h=c(0,1),
#     abline.lty=c(1,2),
#     abline.col=c('black', 'grey'),
#     type=c('h','p')
#     #xlimits=xlims,
#     #plot.horizontal=F,
#   )
#   
#   print(bp.beta.err)
  
  # SAME AS ABOVE BUT ALL METRICS
  # For Scoring Design Paper - Supplementary (maybe)
  # Scoring behavior when predicted CCM is created by adding random zero-mean beta error
  # to the true CCM matrix entries (with different values for the concentration parameter in the beta error) 
  res <- data.frame(matrix(nrow=0, ncol=3))
  colnames(res) <- c("Error", "Score", "Metric")
  for(m in method.names){
    d <- read.csv(paste(tsv.dir, "scoring2B_beta_", m, ".tsv", sep=""),sep="\t",header=FALSE)
    colnames(d) = c("Error","Score")
    d$Score = d$Score
    d$Metric = m
    res <- rbind(res, d)
  }
  
  
  # All Cases SC3
  # data for barplot  
  plot.col <- default.colours(length(unique(res$Metric)))
  groups <- sort(unique(res$Metric))
  
  legend = list(
    text = list(
      lab = groups,
      cex = 1,
      col = 'black'
    ),
    points = list(
      pch = 19,
      col = plot.col,
      cex = 1
    ),
    x = 0.74,
    y = 0.45,
    padding.text = 2
    )  

  #xlims <- c(0.9*min(xdata), 1.1*max(xdata))
  
  bp.beta.err <- create.scatterplot(
    Score ~ as.ordered(Error),
    res,
    groups=Metric,
    cex=1,
    #filename = paste(plot.dir, 'SC2_beta_err_', method, '.png', sep=''),
    main=paste("SC2 w/ Err - ", format.name(method)),
    main.cex=1.5,
    col=plot.col,
    xaxis.rot=30,
    xaxis.cex=0.8,
    yaxis.cex = 0.8,
    ylab.label="Metric Score",
    ylab.cex=1,
    xlab.label="Concentration of Beta Err",
    xlab.cex=1,
    abline.h=c(0,1),
    abline.lty=c(1,2),
    abline.col=c('black', 'grey'),
    type=c('l','p'),
    key=legend
    #xlimits=xlims,
    #plot.horizontal=F,
  )
  
  print(bp.beta.err)
  
#   # For Scoring Design Paper - Supplementary
#   # Scoring behavior for different 'mistake scenarios' but with 10 true clusters (instead of
#   # 3 as is the case above)
#   d = read.csv(file=paste(tsv.dir, "scoring2A_big_cases_", method, ".tsv", sep=""), sep="\t",header=FALSE)
#   colnames(d) = c("Case","Metric")
#   
#   png(file=paste(plot.dir, "2A_Big_Cases_", method, ".png", sep=""))
#   print(
#     ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
#       geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
#       theme(legend.position="none") + ylab("2 Metric") +
#       xlab("Case") + ggtitle(paste("2A Cases with 10 Clusters - Metric =", method)) + coord_flip(ylim=c(0,1))
#   )
#   dev.off()
#   
#   # For Scoring Design Paper - Supplementary
#   # Scoring behavior when predicted CCM is created by reassigning mutations to another random cluster
#   # (with varying probability)
#   d = read.csv(paste(tsv.dir, "scoring2A_random_reassignment_", method, ".tsv", sep=""),sep="\t",header=FALSE)
#   colnames(d) = c("Error","Metric")
#   png(file=paste(plot.dir, "2A_random_", method, ".png", sep=""))
#   
#   print(
#     ggplot(d, aes(x=as.ordered(Error), y=as.numeric(Metric))) + 
#       geom_jitter(aes(color=as.ordered(Error)),position = position_jitter(height = 0, width=0.05)) +
#       stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.7) +
#       theme(legend.position="none") + xlab("Error Probability") + ylab("2 Metric") + 
#       ggtitle(paste("2A Random Error - Metric =", method))
#   )
#   dev.off()
#   
#   # For Scoring Design Paper - Not included
#   # Scoring behavior when predicted CCM is created by reassigning mutations to the closest cluster
#   # (with varying probability)
#   d = read.csv(paste(tsv.dir, "scoring2A_closest_reassignment_", method, ".tsv", sep=""),sep="\t",header=FALSE)
#   colnames(d) = c("Error","Metric")
#   
#   png(file=paste(plot.dir, "2A_closest_", method, ".png", sep=""))
#   print(
#     ggplot(d, aes(x=as.ordered(Error), y=as.numeric(Metric))) + 
#       geom_jitter(aes(color=as.ordered(Error)),position = position_jitter(height = 0, width=0.05)) +
#       stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.7) +
#       theme(legend.position="none") + xlab("Error Probability") + ylab("2 Metric") + 
#       ggtitle(paste("2A Closest Error - Metric =", method))
#   )
#   dev.off()
# 

  
}


##### plot.SC2.certainty.scoring #############################################################
# For Scoring Design Paper - Supplementary (if at all)
# TODO: make interactive plots?
# TODO: combine the plots
# Create a plot of the effect on using probabilty scoring for multiple scenrios for a given
# scoring method
#
# INPUT:
#     method - The scoring metric that was used to create the tsv data
#     display - logical for whether to print the plot or not
# OUPUT:
#     bp - barplot created
#     scenarios - scenarios that are being evaluated
plot.SC2.certainty.scoring <- function(method="pseudoV", display=T){
  data <- read.csv(file=paste(tsv.dir, "2B_prob_scoring_", method, ".tsv", sep=""), sep=",", header=T)
  data.reformat <- data.frame(matrix(nrow=0, ncol=3))
  datal <- apply(data, 1, function(row){
    df <- data.frame(matrix(nrow=length(row)-1, ncol=0))
    df['Score'] <- row[-1]
    df['Scenario'] <- row[1]
    df['Certainty'] <- sapply(names(row)[-1], function(x){as.numeric(substring(x,2))})
    return(df)
  })
  for(i in datal){
    data.reformat <- rbind(data.reformat, i)
  }
  
#   data.reformat <- data.reformat[order(data.reformat$Score),]
#   
#   data.nonprob <- data.reformat[data.reformat$Certainty == 1,c('Score', 'Scenario')]
#   data.nonprob <- data.nonprob[order(data.nonprob$Score),]
  
  ydata <- as.numeric(data.reformat$Score)
  xdata <- as.numeric(data.reformat$Certainty)
  
  ylimits <- c(0,1)#range(ydata)*c(0.9,1.1)
  xlimits <- range(xdata)*c(0.9,1.1)
  
  
  groups <- sort(unique(data.reformat$Scenario))
  plot.col <- colour.gradient(colour.scheme[method], length(groups))#default.colours(length(groups))
    
  legend = list(
    text = list(
      lab = groups,
      cex = 1,
      col = 'black'
      ),
    points = list(
      pch = 19,
      col = plot.col,
      cex = 1
      ),
    x = 0.04,
    y = 0.95,
    padding.text = 2)

  if(method == 'default'){
    plot.main <- 'Certainty of Clustering - SC2'
    legend$y <- 0.5
  }else if(method == 'sym_pseudoV'){
    plot.main <- 'Certainty of Clustering - SC3'
  } else{
    plot.main <- paste('Certainty of Clustering - ', format.name(method))
  }
  
  sp <- create.scatterplot(
    Score ~ Certainty,
    data=data.reformat,
    groups=data.reformat$Scenario,
    filename=paste(plot.dir, 'certainty_', method, '.png', sep=""),
    main=plot.main,
    main.cex = 2.2,
    xlab.cex = 2,
    ylab.cex = 2,
    col=plot.col,
    type='l',
    lwd=3,
    xlimits=xlimits,
    ylimits=ylimits,
    key = legend
    )
  if(display){
    print(sp)
  }
  
  return(list(sp=sp, scenarios=groups, ymin=min(unlist(ydata))))
}


##### plot.SC2.certainty.scoring.err #############################################################
# For Scoring Design Paper - Supplementary
# Create a plot of the effect on using probabilty scoring for multiple scenrios for a given
# scoring method while also adding some random error to the probabilities
#
# INPUT:
#     method - The scoring metric that was used to create the tsv data
#     display - logical for whether to print the plot or not
# OUPUT:
#     sp - scatterplot created
#     scenarios - scenarios tested
plot.SC2.certainty.scoring.err <- function(std=0.05, method="pseudoV", display=T){
  if(std == 0){
    return(plot.SC2.certainty.scoring(method=method, display=display))
  } 
  data <- read.csv(file=paste(tsv.dir, "2B_prob_scoring_with_err_", std, '_', method, ".tsv", sep=""), sep=",", header=T)
  
  data.reformat <- data.frame(matrix(nrow=0, ncol=3))
  datal <- apply(data, 1, function(row){
    df <- data.frame(matrix(nrow=length(row)-1, ncol=0))
    df['Score'] <- row[-1]
    df['Scenario'] <- row[1]
    df['Certainty'] <- sapply(names(row)[-1], function(x){as.numeric(substring(x,2))})
    return(df)
  })
  for(i in datal){
    data.reformat <- rbind(data.reformat, i)
  }
  
  ydata <- as.numeric(data.reformat$Score)
  xdata <- as.numeric(data.reformat$Certainty)
  
  ylimits <- c(0,1)#range(ydata)*c(0.9,1.1)
  xlimits <- range(xdata)*c(0.9,1.1)
  
  groups <- sort(unique(data.reformat$Scenario))
  plot.col <- default.colours(length(groups))
  
  legend = list(
    text = list(
      lab = groups,
      cex = 1,
      col = 'black'
    ),
    points = list(
      pch = 19,
      col = plot.col,
      cex = 1
    ),
    x = 0.04,
    y = 0.95,
    padding.text = 2)
  
  sp <- create.scatterplot(
    Score ~ Certainty,
    data=data.reformat,
    groups=data.reformat$Scenario,
    main=paste('Effect of Certainty of Clustering on Score using', method, 'with Normal Error, std =', std),
    main.cex = 1.5,
    col=plot.col,
    type='l',
    lwd=3,
    xlimits=xlimits,
    ylimits=ylimits,
    #key = legend
  )
  if(display){
    print(sp)
  }
  
  return(list(sp=sp, scenarios=groups, ymin=min(unlist(ydata))))
}

##### plot.SC2.certainty.multi #############################################################
# For Scoring Design Paper - Supplementary
# Create a plot of the effect of using probabilty scoring for multiple scenarios for a given
# scoring method (multiplot of plot.SC2.certainty.scoring or plot.SC2.certainty.scoring.err plots)
#
# INPUT:
#     method - The scoring metric that was used to create the tsv data
#     display - logical for whether to print the plot or not
#     err - logical for whether to plot the error or non-err
# OUPUT:
#     bp - barplot created
plot.SC2.certainty.multi <- function(method='sym_pseudoV', display=T, err=T){
  if(err){
    # retieve plots
    data <- lapply(std.values, plot.SC2.certainty.scoring.err, method=method, display=F)
    plots <- lapply(data, function(x){x$sp})
    ymin <- sapply(data, function(x){x$ymin})
    ymin <- min(ymin)
    scenarios <- data[[1]]$scenarios
    
    print(scenarios)
    
    # calculate multiplot inputs
    n.plots <- length(plots)
    plot.col <- default.colours(length(scenarios))
    
    ylimits <- c(0.9*ymin, 1)
    print(ylimits)
    
    legend <- legend.grob(
      list(
        legend = list(
          colours = plot.col,
          title = 'Scenarios',
          labels = scenarios,
          size = 3,
          title.cex=1,
          label.cex = 1
          )
        )
      )
    
    mp <- create.multiplot(
      plots,
      plot.layout = c(2, ceiling(n.plots/2.0)),
      filename = paste(plot.dir, 'SC2_', method, '_change_certainty.png'),
      main = paste('Change in', method, 'Score with Certainty of Clustering with Different Amounts of Error'),
      main.cex = 1.6,
      xaxis.cex=1,
      yaxis.cex = 1,
      xlab.label = 'Certainty',
      ylab.label = 'Score',
      ylimits = c(0,1),
      print.new.legend=T,
      retrieve.plot.labels=F,
      legend = list(right=list(fun=legend))
      )
    
    if(display){
      print(mp)
    }
    return(mp)
  }
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
  if(name=='default'){name <- c('pseudoV_nc', 'mcc_nc', 'pearson_nc')}
  name.out <- gsub('_', ' ', name)
  name.out <- paste(name.out, collapse=' + ')
  return(name.out)
}
