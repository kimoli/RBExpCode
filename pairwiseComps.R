# this script will load csv files and run paired comparisons on them
#   (either parametric or non-parametric, depending on whether the
#   data pass or fail a Shapiro-Wilk test)
# the script expects that csv's will be arranged in the following manner:
#   - CSV's will have 2 columns of paired data
#   - The first row of the CSV will be column headers labelling the data
#       in that row. Each header should contain no spaces.
#   - The first column will be the column that the experimenter expects
#       to be higher than the second column

rm(list=ls())
library(coin)
library(effsize)

setwd("E:/pcp2ChR2 data/rebound")

files <- list.files(path = ".", pattern = "*forPairedComp*",
           all.files = FALSE, full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

firstiter <- 1 # keep track of whether the loop is in the first iteration
for (filename in files){
  data<-read.csv(filename, header = TRUE)
  
  # check if this is a test for the BStartle parameters
  #if (grepl("bstart", colnames(data[1]), ignore.case=TRUE)) {
    # if true, do a two-sided comparison
    tailsToUse <- "two.sided"
  #} else if(grepl("mvtlatency", colnames(data[1]), ignore.case=TRUE)) {
  #  tailsToUse <- "two.sided"
  #} else { tailsToUse <- "greater"}
  
  # first, run a Shapiro-Wilk test on the differences between the pairs
  diffs <- data[[1]]-data[[2]]
  shapRes <- shapiro.test(diffs)
  
  # if the Shapiro-Wilk result comes out significant, use a Wilcoxon Signed Rank test. Otherwise, use a paired t-test
  if (shapRes[2]<0.05){
    
    # signed rank test
    if (sum(is.nan(data[[1]])) == sum(is.nan(data[[2]]))){
      ispaired<-TRUE
      testname <- "signed rank"
      
      # need to kick out nan values
      tempvals <- data[[1]]
      data[[1]]<-tempvals[is.nan(tempvals)==FALSE]
      tempvals <- data[[2]]
      data[[2]]<-tempvals[is.nan(tempvals)==FALSE]
    } else {
      ispaired<-FALSE
      testname <- "ranksum"
    }
    signrankOutput<-wilcox.test(data[[1]], data[[2]], paired=ispaired, alternative = tailsToUse) # tail specified is that first group is greater than second group
    vval <- as.numeric(signrankOutput[1])
    exactSignrankOutput<-wilcoxsign_test(data[[1]] ~ data[[2]], distribution="exact", alternative = "greater")
    zval <- as.numeric(exactSignrankOutput@statistic@teststatistic)
    pval <- as.numeric(exactSignrankOutput@distribution@pvalue(exactSignrankOutput@statistic@teststatistic))
    group1 <- colnames(data[1])
    group2 <- colnames(data[2])
    rval<-zval/sqrt(length(data[[1]]))
  
    if (firstiter==1) {
      testResults <- data.frame(groupA = unlist(group1),
                                groupB = unlist(group2),
                                statistic = vval,
                                Zordf = zval, # Z or df
                                p = pval,
                                effectsize = rval,
                                test = testname,
                                tails = tailsToUse)
      firstiter <- 0
    } else {
      tempFrame <- data.frame(groupA = unlist(group1),
                                groupB = unlist(group2),
                                statistic = vval,
                              Zordf = zval, #z or df
                                p = pval,
                                effectsize = rval,
                                test = testname,
                              tails = tailsToUse)
      testResults <- rbind(testResults, tempFrame)
    }
  } else {
    # t-test
    if (sum(is.nan(data[[1]])) == sum(is.nan(data[[2]]))){
      ispaired<-TRUE
      testname <- "paired t-test"
      
      # need to kick out nan values
      tempvals <- data[[1]]
      data[[1]]<-tempvals[is.nan(tempvals)==FALSE]
      tempvals <- data[[2]]
      data[[2]]<-tempvals[is.nan(tempvals)==FALSE]
    } else {
        ispaired<-FALSE
        testname <- "unpaired t-test"
    }
    ttestOutput<-t.test(data[[1]], data[[2]], paired = ispaired, alternative = tailsToUse) # tail specified is that first group is greater than second group
    tval <- as.numeric(ttestOutput[1])
    df <- as.numeric(ttestOutput[2])
    pval <- as.numeric(ttestOutput[3])
    group1 <- colnames(data[1])
    group2 <- colnames(data[2])
    g<-cohen.d(data[[2]], data[[1]],pooled=TRUE,paired=ispaired,
               na.rm=FALSE, mu=0, hedges.correction=TRUE,
               conf.level=0.95,noncentral=FALSE,within=TRUE, subject=NA)
    
    if (firstiter==1) {
      testResults <- data.frame(groupA = unlist(group1),
                                groupB = unlist(group2),
                                statistic = tval,
                                Zordf = df,
                                p = pval,
                                effectsize = as.numeric(g[4]),
                                test = testname,
                                tails = tailsToUse)
      firstiter <- 0
    } else {
      tempFrame <- data.frame(groupA = unlist(group1),
                              groupB = unlist(group2),
                              statistic = tval,
                              Zordf = df,
                              p = pval,
                              effectsize = as.numeric(g[4]),
                              test = testname,
                              tails = tailsToUse)
      testResults <- rbind(testResults, tempFrame)
    }
  }
}

write.csv(testResults,"PairwiseComps_20200514.csv", row.names = FALSE)
