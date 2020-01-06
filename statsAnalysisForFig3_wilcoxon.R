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

attach(Data)
## set up main effect test as wilcoxon signed rank
power <- 15
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
RBProb_postExt <- RB_prob[phase=="post_extinction" & power_mW==power]
median(RBProb_postTrain)
mad(RBProb_postTrain)
median(RBProb_postExt)
mad(RBProb_postExt)
# run test
wilcox.test(RBProb_postTrain, RBProb_postExt, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=6, n.s., p = 0.07446

power <- 30
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
RBProb_postExt <- RB_prob[phase=="post_extinction" & power_mW==power]
median(RBProb_postTrain)
mad(RBProb_postTrain)
median(RBProb_postExt)
mad(RBProb_postExt)
# run test
wilcox.test(RBProb_postTrain, RBProb_postExt, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=28, p = 0.01077

power <- 60
RBProb_postTrain <- RB_prob[phase=="post_training" & power_mW==power]
RBProb_postExt <- RB_prob[phase=="post_extinction" & power_mW==power]
median(RBProb_postTrain)
mad(RBProb_postTrain)
median(RBProb_postExt)
mad(RBProb_postExt)
# run test
wilcox.test(RBProb_postTrain, RBProb_postExt, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=19, p = 0.04675


detach(Data)
