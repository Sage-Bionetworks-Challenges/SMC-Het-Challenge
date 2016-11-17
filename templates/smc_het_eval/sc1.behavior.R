#library(ggplot2)
library(BoutrosLab.plotting.general)

# Directories for tsv files and plots
tsv.dir = "scoring_metric_data/text_files/"
plot.dir = "scoring_metric_data/metric_behaviour_plots/"
# Directory that the script is running in (should be <SMC-Het-Challenge git repo>/smc_het_eval)
script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)

plot.amit <- function(){
  # Show how the 1A score (using a certain scoring metric) changes based on errors in the predicted cellularity, for different true cellularity values
  metrics.1A = c('abs', 'sqr') # possible metrics for SC 1A, either using the abs error or the squared error as a penalty
  for(metric in metrics.1A){
    # read in and format data
    da = read.csv(paste(tsv.dir, "scoring1A_behavior_", metric, ".tsv", sep=""),sep="\t",header=FALSE)
    colnames(da) = c("Real","Pred","Error")
    
    # create the plot (all other sections follow the same format)
    png(file=paste(plot.dir, "1A_", metric, ".png", sep="")) # plot filepath
    print(
      ggplot(da,aes(x=Pred,y=Error,color=as.factor(Real))) + geom_point() + xlab("Predicted Cellularity") + 
        ylab("Cellularity Error") + scale_color_discrete(name="Actual Cellularity") + ggtitle(paste("1A Scoring - Error metric =", metric))
    )
    dev.off() # need this after each plot
  }
  
  # Show how the 1B score (using a certain scoring metric) changes based on errors in the predicted cellularity number of clusters, 
  # for different true cluster number values
  # possible metrics for SC 1B, either the original metric (absolute difference of true and predicted # of clusters, divided by true # of clusters)
  # or the normalized version of this metric (normalized to be between 0 and 1)
  metrics.1B = c('orig', 'normalized') 
  for(metric in metrics.1B){
    db = read.csv(paste(tsv.dir, "scoring1B_behavior_", metric, ".tsv", sep=""),sep="\t",header=FALSE)
    colnames(db) = c("Real","Pred","Error")
    
    png(file=paste(plot.dir, "1B_", metric, ".png", sep="")) # plot filepath
    print(
      ggplot(db,aes(x=Pred,y=Error,color=as.factor(Real))) + geom_point() + geom_line() + xlab("Predicted Number of Clusters") + 
        ylab("Cluster Error") + scale_color_discrete(name="Actual Number of Clusters") + ggtitle(paste("1B Scoring - Error metric =", metric))
    )
    dev.off()
  }
  
  # Various plots on the behaviour of the 1C scoring metric using either abs error penalty or squared error penalty (on the error between
  # true and predicted cellular proportion for each mutation)
  metrics.1C = c('abs', 'sqr')
  for(metric in metrics.1C){
    # Scoring behavior when the predicted submission is calculated by systematically adding error to the cellular proportion 
    # of each mutation, but keeping the size of each cluster correct.
    d = read.csv(paste(tsv.dir, "scoring1C_phi_sys_behavior_", metric, ".tsv", sep=""),sep="\t",header=FALSE)
    colnames(d) = c("Error","Metric")
    
    png(file=paste(plot.dir, "1C_phi_systematic_", metric, ".png", sep=""))
    print(
      ggplot(d,aes(x=Error,y=Metric)) + geom_point() + xlab("Systematic Phi Error") + ylab("1C Metric") + 
        ggtitle(paste("Systematic Phi Error - Error metric =", metric))
    )
    dev.off()
    
    # Scoring behavior when the predicted submission is calculated by adding beta error (centered around zero) to the cellular proportion
    # of each mutation, but keeping the size of the clusters correct. Plot shows metric behaviour using various values for the beta error
    # concentration parameter (which determines the size of the error)
    d = read.csv(paste(tsv.dir, "scoring1C_phi_ZM_behavior_", metric, ".tsv", sep=""),sep="\t",head=FALSE)
    colnames(d) = c("Concentration","PhisUsed","Metric")
    
    png(file=paste(plot.dir, "1C_phi_ZM_", metric, ".png", sep=""))
    print(
      ggplot(d, aes(x=as.ordered(Concentration), y=as.numeric(Metric))) + 
        geom_jitter(aes(color=as.ordered(Concentration)),position = position_jitter(height = 0, width=0.05)) +
        stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.7) +
        theme(legend.position="none") + xlab("Concentration Parameter") + ylab("1C Metric") + 
        ggtitle(paste("Zero-mean Phi Noise - Error metric =", metric))
    )
    dev.off()
    
    # Scoring behavior when the predicted submission is calculated by adding beta error (centered around zero) to the size of each cluster
    # of each mutation, but keeping the cellular proportions correct. Plot shows metric behaviour using various values for the beta error
    # concentration parameter (which determines the size of the error)
    d = read.csv(file=paste(tsv.dir, "scoring1C_nssm_behavior_", metric, ".tsv", sep=""), sep="\t",header=FALSE)
    colnames(d) = c("Concentration","NssmsUsed","Metric")
    
    png(file=paste(plot.dir, "1C_nssms_ZM_", metric, ".png", sep=""))
    print(
      ggplot(d, aes(x=as.ordered(Concentration), y=as.numeric(Metric))) + 
        geom_jitter(aes(color=as.ordered(Concentration)),position = position_jitter(height = 0, width=0.05)) +
        stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.7) +
        theme(legend.position="none") + xlab("Concentration Parameter") + ylab("1C Metric") +
        ggtitle(paste("Zero-mean N_ssms Noise - Error metric =", metric))
    )
    dev.off()
    
    # SUPPLEMENTARY
    # Scoring behavior for different 'mistake scnearios', i.e. all mutations in one cluster, two clusters are merged together, 
    # one cluster split into two. 
    d = read.csv(file=paste(tsv.dir, "scoring1C_cases_", metric, ".tsv", sep=""), sep="\t",header=FALSE)
    colnames(d) = c("Case","Metric")
    
    png(file=paste(plot.dir, "1C_Cases_", metric, ".png", sep=""))
    print(
      ggplot(d,aes(y=1-Metric,x=as.factor(Case))) + 
        geom_bar(aes(fill=as.factor(Case)),stat="identity",width=.6) + 
        theme(legend.position="none") + ylab("1 - 1C Metric") +
        xlab("Case") + ggtitle(paste("1C Cases - Error metric =", metric)) + coord_flip()
    )
    dev.off()
  }
}


#### sc1c.heatmap ##############################################
# Heatmap comparing the phi error and the SSM assignment error
# METHODS/SUPP (with 1C Analysis)

sc1c.heatmap <- function(method='abs'){
  data <- read.csv(file=paste(tsv.dir, 'scoring1C_interaction_behavior_', method, '.tsv', sep=''), 
                   sep='\t',
                   row.names = 1
                     )
  xaxis.lab <- sapply(rownames(data), function(x){strsplit(x,'=')[[1]][[2]]})
  yaxis.lab <- sapply(colnames(data), function(x){strsplit(x,'[.]')[[1]][[2]]})
  print(data)
  
  hm <- create.heatmap(
    x = data,
    main = '1C Error (Beta Conc)',
    xlab.label = 'Phi Error (Beta Conc)',
    ylab.label = 'SSM Assignment Error',
    xaxis.lab = xaxis.lab,
    xaxis.cex = 1.4,
    xaxis.rot=0,
    yaxis.lab = yaxis.lab,
    yaxis.cex = 1.4,
    yaxis.rot = 0,
    clustering.method='none',
    filename=paste(plot.dir, 'sc1c_heatmap.png', sep='')
    )
  return(hm)
}






