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


## check if there is an effect of training on RBs across the different intensities
	# result summary: significant at all intensities
Data$power_mW<-factor(Data$power_mW)
Data$phase<-factor(Data$phase)
Data$mouse<-factor(Data$mouse)

attach(Data)

idx <- power_mW == 15
amp15<-RB_amp[idx]
phase15<-phase[idx]
mouse15<-mouse[idx]
friedman.test(amp15 ~ phase15|mouse15) #p = 0.01193
	# post hoc
	idx1 <- phase15 == "pretest"
	idx2 <- phase15 == "post_training"
	idx3 <- phase15 == "post_extinction"
	wilcox.test(amp15[idx1], amp15[idx2], paired=TRUE) # p = 0.01563, s @ bonf
	wilcox.test(amp15[idx1], amp15[idx3], paired=TRUE) # p = 0.07813, ns
	wilcox.test(amp15[idx2], amp15[idx3], paired=TRUE) # p = 0.07813, ns

idx <- power_mW == 30
amp30<-RB_amp[idx]
phase30<-phase[idx]
mouse30<-mouse[idx]
friedman.test(amp30 ~ phase30|mouse30) #p = 0.002149, s @ bonf-holm
	# post hoc
	idx1 <- phase30 == "pretest"
	idx2 <- phase30 == "post_training"
	idx3 <- phase30 == "post_extinction"
	wilcox.test(amp30[idx1], amp30[idx2], paired=TRUE) # p = 0.01563, s @ bonf
	wilcox.test(amp30[idx1], amp30[idx3], paired=TRUE) # p = 0.04688, s @ bonf
	wilcox.test(amp30[idx2], amp30[idx3], paired=TRUE) # p = 0.01563, s @ bonf-holm


idx <- power_mW == 60
amp60<-RB_amp[idx]
phase60<-phase[idx]
mouse60<-mouse[idx]
friedman.test(amp60 ~ phase60|mouse60) # p = 0.01193
	# post hoc
	idx1 <- phase60 == "pretest"
	idx2 <- phase60 == "post_training"
	idx3 <- phase60 == "post_extinction"
	wilcox.test(amp60[idx1], amp60[idx2], paired=TRUE) # p = 0.04688, ns @ bonf-holm
	wilcox.test(amp60[idx1], amp60[idx3], paired=TRUE) # p = 0.1094, ns
	wilcox.test(amp60[idx2], amp60[idx3], paired=TRUE) # p = 0.01563, s @ bonf
	
	# try to save it
	# the data are skewed and violate the assumptions of the wilcoxon signed rank test (symmetry), so they were log transformed, checked for normality, then run through a paired t-test
	t.test(log(amp60[idx1]), log(amp60[idx2]), paired=TRUE) # p = 0.0292, ns @ bonf-holm
	t.test(log(amp60[idx1]), log(amp60[idx3]), paired=TRUE) # p = 0.2313, ns
	t.test(log(amp60[idx2]), log(amp60[idx3]), paired=TRUE) # p = 0.008146 % s
