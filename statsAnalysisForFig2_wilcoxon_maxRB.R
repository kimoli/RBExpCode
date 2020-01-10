# check assumptions for repeated measures anova on the n = 7 PCP2-ChR2 rebound
# probability dataset & then substitute with a binomial test
# data collected by OAKim in JFMedina lab 2019

# install the required packages if they are not already on board
if(!require(ggpubr)){install.packages("ggpubr")}

library(ggpubr)

Data<-read.table('E:/pcp2ChR2 data/rebound/RBPRobs_max_20200108.csv', 
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
hist(Data$RB_prob)

### You can't use an ANOVA because the data are very skewed
#use the wilcoxon signed rank for the paired data

attach(Data)
## set up main effect test as wilcoxon signed rank
RBProb_pre <- RB_prob[phase=="pretest"]
RBProb_postTrain <- RB_prob[phase=="post_training"]
RBProb_postExt <- RB_prob[phase=="post_extinction"]
median(RBProb_pre)
mad(RBProb_pre)
median(RBProb_postTrain)
mad(RBProb_postTrain)
median(RBProb_postExt)
mad(RBProb_postExt)
# run tests
wilcox.test(RBProb_postTrain, RBProb_pre, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=28, p = 0.007813 (signif @ bonf-holm)
wilcox.test(RBProb_postTrain, RBProb_postExt, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=28, p = 0.01113 (signif @ bonf-holm)
wilcox.test(RBProb_postExt, RBProb_pre, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=15, p = 0.02953 (signif @ bonf-holm)

detach(Data)

### now compare the rebound amplitudes between post-training and post-extinction
	# doesn't make sense to look at rebound amplitudes at pretest because
	# median rebound probability during the pretest is 0
DataAmps<-read.table('E:/pcp2ChR2 data/rebound/RBHitAmps_max_20200108.csv', 
	header = TRUE, sep = ",")

DataAmps$phase = factor(DataAmps$phase,
                         levels=unique(DataAmps$phase))

###  Check the data frame

headTail(DataAmps)

str(DataAmps)

summary(DataAmps)

# look at the data
ggboxplot(DataAmps, x = "phase", y = "RB_hitAmp", add = "jitter")
hist(DataAmps$RB_hitAmp)

### the histogram doesn't look that bad, actually
# use shapiro.test to find out if you can do the t-test
attach(DataAmps)
RBHitAmp_T <- RB_hitAmp[phase=="post_training"]
RBHitAmp_E <- RB_hitAmp[phase=="post_extinction"]
shapiro.test(RBHitAmp_T-RBHitAmp_E)
	# W = 0.96629, p-value = 0.851, you can use the ttest
t.test(RBHitAmp_T, RBHitAmp_E, paired = TRUE, alternative = "greater")
	# t = 4.7482, df = 4, p-value = 0.004492 (comes out significant with the
	#	two-tailed, too, but for consistency with the signed rank tests I
	#	am using the more specific test here
detach(DataAmps)


#### Now compare experimental group post-training with an unpaired control
Data2<-read.table('E:/pcp2ChR2 data/rebound/RBProbs_max_expVsCont_20200108.csv', 
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
hist(Data2$RB_prob)

# data are very skewed, use nonparametric test
attach(Data2)
## set up main effect test as wilcoxon rank sum
RBProb_postUnp <- RB_prob[phase=="post_unpaired"]
RBProb_postTrain <- RB_prob[phase=="post_training"]
median(RBProb_postUnp)
mad(RBProb_postUnp)
median(RBProb_postTrain)
mad(RBProb_postTrain)
# run test
wilcox.test(RBProb_postTrain, RBProb_postUnp, paired=FALSE, alternative="greater") # assumes same number of samples pre and post
	# W=74, p = 0.002329
	#	also significant at the two-tailed test but to maintain consistency
	#	with the other tests I have used the one-tailed here

detach(Data2)
