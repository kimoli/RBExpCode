# check assumptions for repeated measures anova on the n = 7 PCP2-ChR2 rebound
# probability dataset & then substitute with a binomial test
# data collected by OAKim in JFMedina lab 2019

# install the required packages if they are not already on board
if(!require(ggpubr)){install.packages("ggpubr")}

library(ggpubr)

Data<-read.table('E:/pcp2ChR2 data/rebound/RBPRobs_20200104.csv', 
	header = TRUE, sep = ",")

###  Order factors by the order in data frame
###  Otherwise, R will alphabetize them

Data$phase = factor(Data$phase,
                         levels=unique(Data$phase))

###  Check the data frame

library(psych)

headTail(Data)

str(Data)

summary(Data)

# look at the data
ggboxplot(Data, x = "phase", y = "RB_prob", add = "jitter")
ggboxplot(Data, x = "power_mW", y = "RB_prob", add = "jitter")

### See if you can use an ANOVA

# check for normality of residulals
meanval <- mean(Data[["RB_prob"]])
residuals <- Data[["RB_prob"]] - meanval
shapiro.test(residuals) # residuals are very skewed
hist(Data[["RB_prob"]])

meanval <- mean(log(Data[["RB_prob"]]+0.00001))
residuals <- log(Data[["RB_prob"]]+0.00001) - meanval
shapiro.test(residuals) # log transforming + a constant doesn't help (tried 1, 0.01, 0.05, 0.5, 0.00001)
hist(log(Data[["RB_prob"]]+0.00001))

### set up main effect test to be binomial test
attach(Data)
RB_Observed <- RB_prob>0

power <- 15
RBObserved_pre <- sum(RB_Observed[phase=="pretest" & power_mW==power], na.rm=TRUE)
RBObserved_postTrain <- sum(RB_Observed[phase=="post_training" & power_mW==power], na.rm=TRUE)
expectedSuccess <- RBObserved_pre # let the number of observations of a rebound before training be the number of successess expected in the binomial test
numSamples <- sum(phase=="post_training" & power_mW==power, na.rm=TRUE)
# run test
binom.test(RBObserved_postTrain, numSamples, expectedSuccess/numSamples, alternative="greater") # assumes same number of samples pre and post
	# you are more likely to see rebounds after training than before training

power <- 30
RBObserved_pre <- sum(RB_Observed[phase=="pretest" & power_mW==power], na.rm=TRUE)
RBObserved_postTrain <- sum(RB_Observed[phase=="post_training" & power_mW==power], na.rm=TRUE)
expectedSuccess <- RBObserved_pre # let the number of observations of a rebound before training be the number of successess expected in the binomial test
numSamples <- sum(phase=="post_training" & power_mW==power, na.rm=TRUE)
# run test
binom.test(RBObserved_postTrain, numSamples, expectedSuccess/numSamples, alternative="greater") # assumes same number of samples pre and post
	# you are more likely to see rebounds after training than before training

power <- 60
RBObserved_pre <- sum(RB_Observed[phase=="pretest" & power_mW==power], na.rm=TRUE)
RBObserved_postTrain <- sum(RB_Observed[phase=="post_training" & power_mW==power], na.rm=TRUE)
expectedSuccess <- RBObserved_pre # let the number of observations of a rebound before training be the number of successess expected in the binomial test
numSamples <- sum(phase=="post_training" & power_mW==power, na.rm=TRUE)
# run test
binom.test(RBObserved_postTrain, numSamples, expectedSuccess/numSamples, alternative="greater") # assumes same number of samples pre and post
	# you are more likely to see rebounds after training than before training

## set up main effect test as wilcoxon signed rank
power <- 15
RBProb_pre <- RB_prob[phase=="pretest" & power_mW==power]
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
# run test
wilcox.test(RBProb_postTrain, RBProb_pre, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=6, n.s., p = 0.07446

power <- 30
RBProb_pre <- RB_prob[phase=="pretest" & power_mW==power]
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
# run test
wilcox.test(RBProb_postTrain, RBProb_pre, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=28, p = 0.01077

power <- 60
RBProb_pre <- RB_prob[phase=="pretest" & power_mW==power]
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
# run test
wilcox.test(RBProb_postTrain, RBProb_pre, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=21, p = 0.01802


detach(Data)


#### Now compare experimental group post-training with an unpaired control
Data2<-read.table('E:/pcp2ChR2 data/rebound/RBPRobs_expCont_20200104.csv', 
	header = TRUE, sep = ",")

###  Order factors by the order in data frame
###  Otherwise, R will alphabetize them

Data2$phase = factor(Data2$phase,
                         levels=unique(Data2$phase))

###  Check the data frame

headTail(Data2)

str(Data2)

summary(Data2)

# look at the data
ggboxplot(Data2, x = "phase", y = "RB_prob", add = "jitter")

### set up main effect test to be binomial test
attach(Data2)
RB_Observed <- RB_prob>0


power <- 15
RBObserved_postUnp <- sum(RB_Observed[phase=="post_unpaired" & power_mW==power], na.rm=TRUE)
RBObserved_postTrain <- sum(RB_Observed[phase=="post_training" & power_mW==power], na.rm=TRUE)
expectedSuccess <- RBObserved_postUnp # let the number of observations of a rebound before training be the number of successess expected in the binomial test
numSamples <- sum(phase=="post_training" & power_mW==power, na.rm=TRUE)
# run test
binom.test(RBObserved_postTrain, numSamples, expectedSuccess/numSamples, alternative="greater") # assumes same number of samples pre and post
	# you are more likely to see rebounds after training than before training

power <- 30
RBObserved_postUnp <- sum(RB_Observed[phase=="post_unpaired" & power_mW==power], na.rm=TRUE)
RBObserved_postTrain <- sum(RB_Observed[phase=="post_training" & power_mW==power], na.rm=TRUE)
expectedSuccess <- RBObserved_postUnp # let the number of observations of a rebound before training be the number of successess expected in the binomial test
numSamples <- sum(phase=="post_training" & power_mW==power, na.rm=TRUE)
# run test
binom.test(RBObserved_postTrain, numSamples, expectedSuccess/numSamples, alternative="greater") # assumes same number of samples pre and post
	# you are more likely to see rebounds after training than before training

power <- 60
RBObserved_postUnp <- sum(RB_Observed[phase=="post_unpaired" & power_mW==power], na.rm=TRUE)
RBObserved_postTrain <- sum(RB_Observed[phase=="post_training" & power_mW==power], na.rm=TRUE)
expectedSuccess <- RBObserved_postUnp # let the number of observations of a rebound before training be the number of successess expected in the binomial test
numSamples <- sum(phase=="post_training" & power_mW==power, na.rm=TRUE)
# run test
binom.test(RBObserved_postTrain, numSamples, expectedSuccess/numSamples, alternative="greater") # assumes same number of samples pre and post
	# you are more likely to see rebounds after training than before training

## set up main effect test as wilcoxon signed rank
power <- 15
RBProb_postUnp <- RB_prob[phase=="post_unpaired" & power_mW==power]
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
# run test
wilcox.test(RBProb_postTrain, RBProb_postUnp, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=9, n.s., p = 0.09294

power <- 30
RBProb_postUnp <- RB_prob[phase=="post_unpaired" & power_mW==power]
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
# run test
wilcox.test(RBProb_postTrain, RBProb_postUnp, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=28, p = 0.01113

power <- 60
RBProb_postUnp <- RB_prob[phase=="post_unpaired" & power_mW==power]
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
# run test
wilcox.test(RBProb_postTrain, RBProb_postUnp, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=28, p = 0.007813

detach(Data2)
