# check assumptions for repeated measures anova on the n = 7 PCP2-ChR2 rebound
# amplitude dataset & then substitute with multiple friedman tests
# data collected by OAKim in JFMedina lab 2019

# install the required packages if they are not already on board
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(rstatix)){install.packages("rstatix")}
if(!require(e1071)){install.packages("e1071")}

library(e1071)  
library(ggpubr)

Data<-read.table('E:/pcp2ChR2 data/rebound/RBAmps_20200104.csv', 
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
ggboxplot(Data, x = "phase", y = "RB_amp", add = "jitter")
ggboxplot(Data, x = "power_mW", y = "RB_amp", add = "jitter")


# check for normality of residulals
attach(Data)

meanval <- mean(RB_amp)
residuals <- RB_amp - meanval
shapiro.test(residuals) # residuals are very skewed
hist(RB_amp)

# just do tests looking at 60 mW rebound amplitudes because there were no pretest rebounds at less than 60 mW
power <- 60
RBAmp_pre <- RB_amp[phase=="pretest" & power_mW==power]
RBAmp_postTrain <- RB_amp[phase=="post_training" & power_mW==power]
median(RBAmp_pre)
mad(RBAmp_pre)
median(RBAmp_postTrain)
mad(RBAmp_postTrain)
# run test
wilcox.test(RBAmp_postTrain, RBAmp_pre, paired=TRUE, alternative="greater") # assumes same number of samples pre and post
	# V=26, p = 0.0.0234

detach(Data)