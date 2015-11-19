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

# Values to use for the standard deviation of the error in the co-clustering probabilities
std.values <- c(0, 0.01,0.03,0.05,0.1,0.15,0.2)

##### plot.SC2.amit ##########################################################################
# plot all the figures that Amit made. These include: 

plot.SC2.amit <- function(method='default'){
  # Scoring behavior for different 'mistake scenarios' including:
  #     - splitting one cluster into two
  #     - merging two clusters into one
  #     - assigning all mutations to the same cluster
  #     - assigning each mutation to its own cluster
  #     - adding an extra cluster with mutations from each true cluster (either small or large)
  d = read.csv(file=paste(tsv_dir, "scoring2A_cases_", method, ".tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  
  png(file=paste(plot_dir, "2A_Cases_", method, ".png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("2 Metric") +
      xlab("Case") + ggtitle(paste("2A Cases - Metric =", method)) + coord_flip(ylim=c(0,1))
  )
  dev.off()
  
  # Scoring behavior for different 'mistake scenarios' but with 10 true clusters (instead of
  # 3 as is the case above)
  d = read.csv(file=paste(tsv_dir, "scoring2A_big_cases_", method, ".tsv", sep=""), sep="\t",header=FALSE)
  colnames(d) = c("Case","Metric")
  
  png(file=paste(plot_dir, "2A_Big_Cases_", method, ".png", sep=""))
  print(
    ggplot(d,aes(y=Metric,x=as.factor(Case))) + 
      geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
      theme(legend.position="none") + ylab("2 Metric") +
      xlab("Case") + ggtitle(paste("2A Cases with 10 Clusters - Metric =", method)) + coord_flip(ylim=c(0,1))
  )
  dev.off()
  
  # Scoring behavior when predicted CCM is created by reassigning mutations to another random cluster
  # (with varying probability)
  d = read.csv(paste(tsv_dir, "scoring2A_random_reassignment_", method, ".tsv", sep=""),sep="\t",header=FALSE)
  colnames(d) = c("Error","Metric")
  png(file=paste(plot_dir, "2A_random_", method, ".png", sep=""))
  
  print(
    ggplot(d, aes(x=as.ordered(Error), y=as.numeric(Metric))) + 
      geom_jitter(aes(color=as.ordered(Error)),position = position_jitter(height = 0, width=0.05)) +
      stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.7) +
      theme(legend.position="none") + xlab("Error Probability") + ylab("2 Metric") + 
      ggtitle(paste("2A Random Error - Metric =", method))
  )
  dev.off()
  
  # Scoring behavior when predicted CCM is created by reassigning mutations to the closest cluster
  # (with varying probability)
  d = read.csv(paste(tsv_dir, "scoring2A_closest_reassignment_", method, ".tsv", sep=""),sep="\t",header=FALSE)
  colnames(d) = c("Error","Metric")
  
  png(file=paste(plot_dir, "2A_closest_", method, ".png", sep=""))
  print(
    ggplot(d, aes(x=as.ordered(Error), y=as.numeric(Metric))) + 
      geom_jitter(aes(color=as.ordered(Error)),position = position_jitter(height = 0, width=0.05)) +
      stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.7) +
      theme(legend.position="none") + xlab("Error Probability") + ylab("2 Metric") + 
      ggtitle(paste("2A Closest Error - Metric =", method))
  )
  dev.off()
  
#   # Scoring behavior when predicted CCM is created by adding random zero-mean beta error
#   # to the true CCM matrix entries (with different values for the concentration parameter in the beta error) 
#   d = read.csv(paste(tsv_dir, "scoring2B_beta_", method, ".tsv", sep=""),sep="\t",header=FALSE)
#   colnames(d) = c("Error","Metric")
#   
#   png(file=paste(plot_dir, "2B_beta_", method, ".png", sep=""))
#   print(
#     ggplot(d, aes(x=as.ordered(Error), y=as.numeric(Metric))) + 
#       geom_jitter(aes(color=as.ordered(Error)),position = position_jitter(height = 0, width=0.05)) +
#       stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.7) +
#       theme(legend.position="none") + xlab("Concentration Parameter") + ylab("2 Metric") + 
#       ggtitle(paste("2B Beta Noise - Metric =", method))
#   )
#   dev.off()
}


##### plot.SC2.certainty.scoring #############################################################
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
    main=paste('Effect of Certainty of Clustering on Score using', method),
    main.cex = 2.2,
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


##### plot.SC2.certainty.scoring.err #############################################################
# Create a plot of the effect on using probabilty scoring for multiple scenrios for a given
# scoring method
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
# Create a plot of the effect on using probabilty scoring for multiple scenrios for a given
# scoring method
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
