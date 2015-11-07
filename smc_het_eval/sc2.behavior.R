#### sc2.behavior.R #############################################################
# Analyze the behavior of different scoring metrics for Sub Challenge 3

#### PREAMBLE ###################################################################
library(BoutrosLab.plotting.general)

# Directories for tsv files and plots respectively
setwd("~/Documents/SMC-Het/SMC-Het-Challenge/smc_het_eval")
tsv.dir = "scoring_metric_data/text_files/"
plot.dir = "scoring_metric_data/metric_behaviour_plots/"

# Lists and dictionaries with info on the metrics and the rankings
method.names <- c("orig",
                  "sqrt",
                  "pseudoV",
                  "sym_pseudoV",
                  "spearman",
                  "pearson",
                  "aupr",
                  "mcc")

std.values <- c(0, 0.01,0.03,0.05,0.1,0.15,0.2)

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
    key = legend
    )
  if(display){
    print(sp)
  }
  
  return(sp)
}


##### plot.SC2.certainty.scoring.err #############################################################
# Create a plot of the effect on using probabilty scoring for multiple scenrios for a given
# scoring method
#
# INPUT:
#     method - The scoring metric that was used to create the tsv data
#     display - logical for whether to print the plot or not
# OUPUT:
#     bp - barplot created
plot.SC2.certainty.scoring.err <- function(method="pseudoV", display=T, std=0.05){
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
    key = legend
  )
  if(display){
    print(sp)
  }
  
  return(sp)
}
