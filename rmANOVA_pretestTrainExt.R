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


## check if there is an effect of training on RBs across the different intensities
	# result summary: significant at all intensities
Data$power_mW<-factor(Data$power_mW)
Data$phase<-factor(Data$phase)
Data$mouse<-factor(Data$mouse)

attach(Data)

idx <- power_mW == 15
prob15<-RB_prob[idx]
phase15<-phase[idx]
mouse15<-mouse[idx]
friedman.test(prob15 ~ phase15|mouse15) #p = 0.04979, ns @ bonf holm

idx <- power_mW == 30
prob30<-RB_prob[idx]
phase30<-phase[idx]
mouse30<-mouse[idx]
friedman.test(prob30 ~ phase30|mouse30) #p = 0.001409, s @ bonf
	# post hoc
	idx1 <- phase30 == "pretest"
	idx2 <- phase30 == "post_training"
	idx3 <- phase30 == "post_extinction"
	# two-sided post-hocs
	wilcox.test(prob30[idx1], prob30[idx2], paired=TRUE) # p = 0.02154, s @ bonf
	wilcox.test(prob30[idx1], prob30[idx3], paired=TRUE) # p = 0.3711, ns
	wilcox.test(prob30[idx2], prob30[idx3], paired=TRUE) # p = 0.02154, ns
	# one-sided post-hocs
	wilcox.test(prob30[idx1], prob30[idx2], paired=TRUE, alternative="less") # p = 0.01077, s @ bonf
	wilcox.test(prob30[idx1], prob30[idx3], paired=TRUE, alternative="less") # p = 0.1855, ns
	wilcox.test(prob30[idx2], prob30[idx3], paired=TRUE, alternative="greater") # p = 0.01077, s @ bonf-holm,



idx <- power_mW == 60
prob60<-RB_prob[idx]
phase60<-phase[idx]
mouse60<-mouse[idx]
friedman.test(prob60 ~ phase60|mouse60) # p = 0.01832, s @ bonf-holm
	# post hoc
	idx1 <- phase60 == "pretest"
	idx2 <- phase60 == "post_training"
	idx3 <- phase60 == "post_extinction"
	# 2-sided post-hocs
	wilcox.test(prob60[idx1], prob60[idx2], paired=TRUE) # p = 0.03603, ns @ bonf
	wilcox.test(prob60[idx1], prob60[idx3], paired=TRUE) # p = 0.1148, ns
	wilcox.test(prob60[idx2], prob60[idx3], paired=TRUE) # p = 0.09349, ns
	# 1-sided post-hocs
	wilcox.test(prob60[idx1], prob60[idx2], paired=TRUE, alternative="less") # p = 0.01802, s @ bonf
	wilcox.test(prob60[idx1], prob60[idx3], paired=TRUE, alternative="less") # p = 0.05742, ns
	wilcox.test(prob60[idx2], prob60[idx3], paired=TRUE, alternative="greater") # p = 0.0465, ns @ bonf-holm


### set up main effect test to be binomial test
attach(Data)
RB_Observed <- RB_prob>0
pretestIdx <- phase=="pretest"
postTrainIdx <- phase=="post_training"
postExtIdx <- phase=="post_extinction"
mouseNames <- unique(mouse)

# collapse data across laser powers
RBObserved_pre<-logical()
RBObserved_postTrain<-logical()
RBObserved_postExt<-logical()
for (m in mouseNames){
	pretestIdcs <- phase=="pretest" & mouse==m;
	RBObserved_pre <- append(RBObserved_pre, sum(RB_Observed[pretestIdcs], na.rm=TRUE)>0);

	postTrainIdcs <- phase=="post_training" & mouse==m;
	RBObserved_postTrain <- append(RBObserved_postTrain, sum(RB_Observed[postTrainIdcs], na.rm=TRUE)>0);

	postExtIdcs <- phase=="post_extinction" & mouse==m;
	RBObserved_postExt <- append(RBObserved_postExt, sum(RB_Observed[postExtIdcs], na.rm=TRUE)>0);
}

## check whether there are more rebounds after training compared to before training
# set up threshold to pass: number of rebounds observed before training
expectedSuccess <- sum(RBObserved_pre) # let the number of observations of a rebound before training be the number of successess expected in the binomial test

# run test
binom.test(sum(RBObserved_postTrain, na.rm=TRUE), length(RBObserved_postTrain), expectedSuccess/length(RBObserved_pre))
	# you are more likely to see rebounds after training than before training

## check whether there more rebounds after extinction compared to before training
# use same threshold as last test, but change comparison group
binom.test(sum(RBObserved_postExt, na.rm=TRUE), length(RBObserved_postTrain), expectedSuccess/length(RBObserved_pre))
	# you are not more likely to see rebounds after extinction than before extinction

## check whether there are more rebouns after training compared with after extinction
expectedSuccess <- sum(RBObserved_postExt)
binom.test(sum(RBObserved_postTrain, na.rm=TRUE), length(RBObserved_postTrain), expectedSuccess/length(RBObserved_postExt))
	# you are not more likely to observe rebounds after training than after extinction


detach(Data)


#### Now compare experimental group post-training with an unpaired control
Data2<-read.table('E:/pcp2ChR2 data/rebound/RBPRobs_expCont_20200104.csv', 
	header = TRUE, sep = ",")

###  Order factors by the order in data frame
###  Otherwise, R will alphabetize them

Data2$phase = factor(Data2$phase,
                         levels=unique(Data2$phase))

###  Check the data frame

library(psych)

headTail(Data2)

str(Data2)

summary(Data2)

# look at the data
ggboxplot(Data2, x = "phase", y = "RB_prob", add = "jitter")

### See if you can use an ANOVA

# check for normality of residulals
meanval <- mean(Data2[["RB_prob"]])
residuals <- Data2[["RB_prob"]] - meanval
shapiro.test(residuals) # residuals are very skewed
hist(Data2[["RB_prob"]])

meanval <- mean(log(Data2[["RB_prob"]]+0.00001))
residuals <- log(Data2[["RB_prob"]]+0.00001) - meanval
shapiro.test(residuals) # log transforming + a constant doesn't help (tried 1, 0.01, 0.05, 0.5, 0.00001)
hist(log(Data2[["RB_prob"]]+0.00001))

### set up main effect test to be binomial test
attach(Data2)
RB_Observed <- RB_prob>0
postTrainIdx <- phase=="post_training"
postUnpIdx <- phase=="post_unpaired"
mouseNames <- unique(mouse)

# collapse data across laser powers
RBObserved_postTrain<-logical()
RBObserved_postUnp<-logical()
for (m in mouseNames){
	postTrainIdcs <- phase=="post_training" & mouse==m;
	if (sum(postTrainIdcs, na.rm=TRUE)>0) { 
		RBObserved_postTrain <- append(RBObserved_postTrain, sum(RB_Observed[postTrainIdcs], na.rm=TRUE)>0);
	}

	postUnpIdcs <- phase=="post_unpaired" & mouse==m;
	if (sum(postUnpIdcs, na.rm=TRUE)>0) {
		RBObserved_postUnp <- append(RBObserved_postUnp, sum(RB_Observed[postUnpIdcs], na.rm=TRUE)>0);
	}
}

## check whether there are more rebounds after training compared to before training
# set up threshold to pass: number of rebounds observed before training
expectedSuccess <- sum(RBObserved_postUnp) # let the number of observations of a rebound before training be the number of successess expected in the binomial test

# run test
binom.test(sum(RBObserved_postTrain, na.rm=TRUE), length(RBObserved_postTrain), expectedSuccess/length(RBObserved_postUnp))
	# you are more likely to see rebounds after training than before training

detach(Data2)
