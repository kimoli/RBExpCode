# check assumptions for repeated measures anova on the n = 7 PCP2-ChR2 rebound
# probability dataset & then substitute with a binomial test
# data collected by OAKim in JFMedina lab 2019

# these analyses are for figure 3

# install the required packages if they are not already on board
if(!require(ggpubr)){install.packages("ggpubr")}

library(ggpubr)

Data<-read.table('E:/pcp2ChR2 data/rebound/RBProbs_max_20200108.csv', 
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
# does extinction learning affect rebounds in terms of probability?
RBProb_pre <- RB_prob[phase=="pretest"]
	median(RBProb_pre)
	mad(RBProb_pre)

RBProb_postT <- RB_prob[phase=="post_training"]
	median(RBProb_postT)
	mad(RBProb_postT)

RBProb_postE <- RB_prob[phase=="post_extinction"]
	median(RBProb_postE)
	mad(RBProb_postE)

# run tests
wilcox.test(RBProb_postT, RBProb_postE, paired=TRUE, alternative="greater")
wilcox.test(RBProb_postE, RBProb_pre, paired=TRUE, alternative="greater")

# how does extinction affect the rebound itself?
RBHitAmp_postT <- RB_hitAmp[phase=="post_training"]
	median(RBHitAmp_postT)
	mad(RBHitAmp_postT)

RBHitAmp_postE <- RB_hitAmp[phase=="post_extinction"]
	median(RBHitAmp_postE, na.rm=TRUE)
	mad(RBHitAmp_postE, na.rm=TRUE)

# check data
hist(RBHitAmp_postT)
hist(RBHitAmp_postE)
	# data are not skewed, try checking residuals for normality
shapiro.test(RBHitAmp_postT - RBHitAmp_postE)

# run test
t.test(RBHitAmp_postT, RBHitAmp_postE, paired=TRUE)

detach(Data)
