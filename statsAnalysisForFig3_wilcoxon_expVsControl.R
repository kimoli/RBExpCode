# check assumptions for repeated measures anova on the n = 7 PCP2-ChR2 rebound
# probability dataset & then substitute with a binomial test
# data collected by OAKim in JFMedina lab 2019

# these analyses are for figure 3

# install the required packages if they are not already on board
if(!require(ggpubr)){install.packages("ggpubr")}

library(ggpubr)

Data<-read.table('E:/pcp2ChR2 data/rebound/RBProbs_expVsCont_forFig2_20200108.csv', 
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
# --> use the wilcoxon signed rank for the paired data

attach(Data)
## set up main effect test as wilcoxon signed rank
RBProb_preT <- RB_prob[phase=="pre_training"]
	median(RBProb_preT)
	mad(RBProb_preT)

RBProb_postT <- RB_prob[phase=="post_training"]
	median(RBProb_postT)
	mad(RBProb_postT)

RBProb_preU <- RB_prob[phase=="pre_unpaired"]
	median(RBProb_preU)
	mad(RBProb_preU)

RBProb_postU <- RB_prob[phase=="post_unpaired"]
	median(RBProb_postU)
	mad(RBProb_postU)

# run tests
wilcox.test(RBProb_postT, RBProb_preT, paired=TRUE)
wilcox.test(RBProb_postU, RBProb_preU, paired=TRUE)
wilcox.test(RBProb_preT, RBProb_preU, paired=FALSE)
wilcox.test(RBProb_postT, RBProb_postU, paired=FALSE)

detach(Data)

## now look at the dose response properties

Data<-read.table('E:/pcp2ChR2 data/rebound/RBData_20200108.csv', 
	header = TRUE, sep = ",")

###  Order factors by the order in data frame
###  Otherwise, R will alphabetize them

Data$phase = factor(Data$phase,
                         levels=unique(Data$phase))

###  Check the data frame

headTail(Data)

str(Data)

summary(Data)

attach(Data)
## get data
RBProb_15 <- RB_prob[phase=="post_training" & power_mW==15]
	median(RBProb_15)
	mad(RBProb_15)

RBProb_30 <- RB_prob[phase=="post_training" & power_mW==30]
	median(RBProb_30)
	mad(RBProb_30)

RBProb_60 <- RB_prob[phase=="post_training" & power_mW==60]
	median(RBProb_60)
	mad(RBProb_60)

RBHitAmp_30 <- RB_hitAmp[phase=="post_training" & power_mW==30]
	median(RBHitAmp_30)
	mad(RBHitAmp_30)

RBHitAmp_60 <- RB_hitAmp[phase=="post_training" & power_mW==60]
	median(RBHitAmp_60)
	mad(RBHitAmp_60)

# check if you can use the t-test on the amplitude data because they aren't terribly skewed
shapiro.test(RBHitAmp_60-RBHitAmp_30) # can

# run tests
wilcox.test(RBProb_15, RBProb_30, paired=TRUE, alternative="less")
wilcox.test(RBProb_15, RBProb_60, paired=TRUE, alternative="less")
wilcox.test(RBProb_60, RBProb_30, paired=TRUE, alternative="greater")
t.test(RBHitAmp_30, RBHitAmp_60, paired=TRUE, alternative="less")