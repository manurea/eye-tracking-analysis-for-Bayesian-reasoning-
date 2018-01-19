#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("qualvar")
#install.packages("lattice")
#install.packages("rasterVis")
#install.packages("gplots")
#install.packages("lme4")
#install.packages("languageR")
#install.packages("Hmisc")
#install.packages('Rmisc')
#install.packages('corrplot')
#install.packages('car')
#install.packages('mlogit')
#install.packages("gtools")
#install.packages("MASS")
#install.packages("foreign")
#install.packages("reshape2")
#install.packages("plyr")
library(tidyr)
library(dplyr)
library(ggplot2)
library(qualvar)
library(lattice)
library(rasterVis)
library(gplots)
library(lme4)
library(languageR)
library("Hmisc")
library("Rmisc")
library(corrplot)
library("car")
library("mlogit")
library("gtools")
library("MASS")
library("foreign")
library("reshape2")
library("plyr")

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#                     EXPERIMENT 2A - Bayesian reasoning and eye-tracking analysis                      #
#                     Analysis of the data from the file bayeScanpaths.csv                             #
#                                    Manuele Reani                                                     #
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#


# create the dataframe (here I use the DF without repetitions but 
# I could use the one with repetitions and it woul dbe the same)
dfANALYSIS <- read.csv(file="bayeScanpaths.csv",header=TRUE,sep=",")

# Some random plots
plot(dfANALYSIS$Age, dfANALYSIS$SNS)

boxplot(dfANALYSIS$SNS, dfANALYSIS$Correct)

plot(density(dfANALYSIS$SNS))

plot(density(dfANALYSIS$Age))

boxplot(dfANALYSIS$SNS~dfANALYSIS$Correct)
boxplot(dfANALYSIS$Age~dfANALYSIS$Correct)
boxplot(dfANALYSIS$Time~dfANALYSIS$Correct)

boxplot(dfANALYSIS$SNS~dfANALYSIS$Score)
boxplot(dfANALYSIS$Age~dfANALYSIS$Score)
boxplot(dfANALYSIS$Time~dfANALYSIS$Score)

# Some decriptive stats for Score 
mean(dfANALYSIS$Score[dfANALYSIS$Format=="tree"])
mean(dfANALYSIS$Score[dfANALYSIS$Format=="venn"])
sum(dfANALYSIS$Score[dfANALYSIS$Format == "tree"])
sum(dfANALYSIS$Score[dfANALYSIS$Format == "venn"])

# Some decriptive stats for Correct 
sum(dfANALYSIS$Correct[dfANALYSIS$Format=="tree"]) 
sum(dfANALYSIS$Correct[dfANALYSIS$Format=="venn"]) 

# plot SNS by correct 
meansSNS <- c(mean(dfANALYSIS$SNS[dfANALYSIS$Correct==1]), mean(dfANALYSIS$SNS[dfANALYSIS$Correct==0]))
SDsSNS <- c(sd(dfANALYSIS$SNS[dfANALYSIS$Correct==1]), sd(dfANALYSIS$SNS[dfANALYSIS$Correct==0]))
dfSNScorrect <- data.frame(Means = meansSNS, SD = SDsSNS, Correct = c(1,0))

ggplot(data=dfSNScorrect, aes(x=as.factor(Correct), y=Means)) + 
  geom_line(group=1, size = 1) + geom_point() + 
  geom_errorbar(aes(ymin=Means- SD, ymax=Means + SD), width=.03, position=position_dodge(0.05)) +
  scale_y_continuous(breaks = seq(0, 6, by = 0.5)) +
  expand_limits(y=seq(1, 6, by = 1)) + 
  labs(x = "Correct", y = "Numeracy") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))



#----------------------------------------------------------------------------#
#----- Analyze dwell time for the AOIs by Correctness -----------------------#
#----------------------------------------------------------------------------#

# ---------- TREE -----------------------------
# subset the dataframe and replace NA with zeros
tree1<-subset(dfANALYSIS,Format=="tree")
tree2<-tree1[,1:19]
tree2[is.na(tree2)] <- 0

# verticalize the data using the function gather (in tidyr) 

tree3 <- gather(tree2, AOI,TimeAOI, A:I)

# plot some data

ggplot(tree3, aes(x=AOI, y=TimeAOI, fill= as.factor(Correct))) + 
  geom_boxplot() +
  labs(x = "AOIs (Tree)", y = "Total Dwell Time (ms)") +
  scale_fill_grey(start = 0.1, end = 0.5, name="Correct")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        legend.key.size = unit(0.5, "cm")) 

median(tree3[tree3$Correct==1 & tree3$AOI == "C", 12])
IQR(tree3[tree3$Correct==1 & tree3$AOI == "C", 12])       
median(tree3[tree3$Correct==0 & tree3$AOI == "C", 12])
IQR(tree3[tree3$Correct==0 & tree3$AOI == "C", 12])

median(tree3[tree3$Correct==1 & tree3$AOI == "G", 12])
IQR(tree3[tree3$Correct==1 & tree3$AOI == "G", 12])       
median(tree3[tree3$Correct==0 & tree3$AOI == "G", 12])
IQR(tree3[tree3$Correct==0 & tree3$AOI == "G", 12])

median(tree3[tree3$Correct==1 & tree3$AOI == "H", 12])
IQR(tree3[tree3$Correct==1 & tree3$AOI == "H", 12])       
median(tree3[tree3$Correct==0 & tree3$AOI == "H", 12])
IQR(tree3[tree3$Correct==0 & tree3$AOI == "H", 12])
  
# tital for correct and incorrect for TREE 

median(tree3[tree3$Correct==1, 12])
IQR(tree3[tree3$Correct==1, 12])    
median(tree3[tree3$Correct==0, 12])
IQR(tree3[tree3$Correct==0, 12]) 

# --------- VENN ------------------------------
# subset the dataframe and replace NA with zeros
venn1<-subset(dfANALYSIS,Format=="venn")
venn2<-venn1[,-(11:19)]
venn2[is.na(venn2)] <- 0

# verticalize the data using the function gather (in tidyr) 

venn3 <- gather(venn2, AOI,TimeAOI, L:Q)

# plot some data

ggplot(venn3, aes(x=AOI, y=TimeAOI, fill= as.factor(Correct))) + 
  geom_boxplot() +
  labs(x = "AOIs (Venn)", y = "Total Dwell Time (ms)") +
  scale_fill_grey(start = 0.1, end = 0.5, name="Correct")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        legend.key.size = unit(0.5, "cm"))

median(venn3[venn3$Correct==1 & venn3$AOI == "Q", 12])
IQR(venn3[venn3$Correct==1 & venn3$AOI == "Q", 12])       
median(venn3[venn3$Correct==0 & venn3$AOI == "Q", 12])
IQR(venn3[venn3$Correct==0 & venn3$AOI == "Q", 12])

median(venn3[venn3$Correct==1 & venn3$AOI == "O", 12])
IQR(venn3[venn3$Correct==1 & venn3$AOI == "O", 12])       
median(venn3[venn3$Correct==0 & venn3$AOI == "O", 12])
IQR(venn3[venn3$Correct==0 & venn3$AOI == "O", 12])

median(venn3[venn3$Correct==1 & venn3$AOI == "M", 12])
IQR(venn3[venn3$Correct==1 & venn3$AOI == "M", 12])       
median(venn3[venn3$Correct==0 & venn3$AOI == "M", 12])
IQR(venn3[venn3$Correct==0 & venn3$AOI == "M", 12])

median(venn3[venn3$Correct==1 & venn3$AOI == "N", 12])
IQR(venn3[venn3$Correct==1 & venn3$AOI == "N", 12])       
median(venn3[venn3$Correct==0 & venn3$AOI == "N", 12])
IQR(venn3[venn3$Correct==0 & venn3$AOI == "N", 12])

# tital for correct and incorrect for TREE 

median(venn3[venn3$Correct==1, 12])
IQR(venn3[venn3$Correct==1, 12])    
median(venn3[venn3$Correct==0, 12])
IQR(venn3[venn3$Correct==0, 12]) 

# producing dataframes for statistics of dwelltime NBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
#TREE
tree3corr <- tree3[tree3$Correct==1,]
topogigio<- aggregate(TimeAOI ~ AOI, tree3corr, median) 
topogigio2<- aggregate(TimeAOI ~ AOI, tree3corr, IQR)
stats4treeCorrect <- cbind(topogigio, topogigio2)

tree3incorr <- tree3[tree3$Correct==0,]
topogigio<- aggregate(TimeAOI ~ AOI, tree3incorr, median) 
topogigio2<- aggregate(TimeAOI ~ AOI, tree3incorr, IQR)
stats4treeIncorrect <- cbind(topogigio, topogigio2)

#VENN
venn3corr <-venn3[venn3$Correct==1,]
topogigio<- aggregate(TimeAOI ~ AOI, venn3corr, median) 
topogigio2<- aggregate(TimeAOI ~ AOI, venn3corr, IQR)
stats4vennCorrect <- cbind(topogigio, topogigio2)

venn3incorr <- venn3[venn3$Correct==0,]
topogigio<- aggregate(TimeAOI ~ AOI, venn3incorr, median) 
topogigio2<- aggregate(TimeAOI ~ AOI, venn3incorr, IQR)
stats4vennIncorrect <- cbind(topogigio, topogigio2)


# ----------- Multilevel Log-Regression for Correct -------------------


dfANALYSIS$Format <- relevel(as.factor(dfANALYSIS$Format),"tree")
correctModel <- glmer(Correct ~ Format + SNS + Time + (1 | Ps),family = binomial("logit"), data = dfANALYSIS)
summary(correctModel)

# ----------- Multilevel Regression for Score -------------------

scoreModel <- lmer(Score ~ Format + SNS + Time + (1 | Ps), data = dfANALYSIS)
summary(scoreModel)

#..............
Anova(scoreModel)

# ----------- Multilevel Regression for Time -------------------
# Nomral distribution ?
plot(density(dfANALYSIS$Time))
ggplot(dfANALYSIS, aes(Time)) +
  geom_density()

timeModel <- lmer(Time ~ Format + SNS + (1 | Ps), data = dfANALYSIS)
summary(timeModel)

#..............
Anova(timeModel)

# Correlation between metric (Also to see if people who got it right are fatser or not)
timeModel2 <- lmer(Time ~ Correct + Score + (1 | Ps), data = dfANALYSIS)
summary(timeModel2)

# plotting the average time to solve 1 problem, by correctness 
ggplot(dfANALYSIS, aes(Time, fill = as.factor(Correct))) + 
  geom_density(alpha = 0.6) + 
  xlim(-2, 80) + 
  labs(fill = "Correct") + scale_fill_grey()


#SOME COMMENTS: 
#Format doesn't correlate with Correctness nor with Score nor with Time (speed). 
#Numeracy doesn't correlate with Score nor with Time (speed). 
#It correlates slightly however with Correctness (e.g. people with higher self-reported level of numeracy get slightly more correct answers). 
#Time (speed) doesn't correlate with Correctness nor with Score. 
#For both Tree and Venn, there is no difference in total dwell time for single AOIs for incorrect and incorrect groups.

#CONCLUSION: 
#Graphical format does not predict performance. Numeracy slightly predicts performance.
#Moreover, people’s eye-behaviour (measured only as dwell time) appears to be similar across conditions. 

#But what about Scanpath analysis?  
#I have not done this analysis yet. 
#I want to use Markov chain permutation test to see if there is a significant difference in scapaths.
#Compare this with another inferential test using permutations, but looking at histogram distance (frequency) instead of Jensen-shannon distance. 
#And then, if there is a difference (hopefully),  find the optimal scanpath with two methods:
#(1)Distance based (Sukru style)
#(2)2 Frequency based methods (using hystogram distance and SVM).




# some fursther analysis 
wilcox.test(dfANALYSIS$SNS~dfANALYSIS$Correct)
median(dfANALYSIS$SNS[dfANALYSIS$Correct== "1"])
median(dfANALYSIS$SNS[dfANALYSIS$Correct== "0"])

logit <- glm(Correct~Format+SNS, data = dfANALYSIS, family = binomial)
summary(logit)

# END REGRESSION ANALYSIS
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################












# ----------------------------------------------------------------------------------------------------------------- #
#                                                                                                                   #
#                                      Scanpath analysys                                                            #
#                                                                                                                   #
# ----------------------------------------------------------------------------------------------------------------- #




###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
#                                                                                                                                             #
#                              WE DO HERE THE ANALYSIS FOR THE DATAFRAME WITH THE REPETITIONS INCLUDED                                        #
#                                                                                                                                             #
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################


dfANALYSISrep <- read.csv(file="bayeScanpathsrep.csv",header=TRUE,sep=",")







#---------------------------------------------------------------------------------#
#                                                                                 #
#        Create the DFs to work with                                              #    
#                                                                                 #
#---------------------------------------------------------------------------------#




# DFs TREE 
# Subset by Tree only to create a "longdata" df 
GroupTree <- subset(dfANALYSISrep,Format == "tree")
dfTransitionTree <- data.frame(subject = NULL, state = NULL) 
rigaT<-1 
for(p in GroupTree$Ps){
  #print(p)
  scanna <- c(unlist(strsplit(toString(GroupTree[rigaT,10]), split="")))
  for (n in scanna){
    newRow <- data.frame(subject = p, state = n)
    dfTransitionTree <- rbind(dfTransitionTree,newRow)
  }
  rigaT <- rigaT+1
}
dfTransitionTree$state <- factor(dfTransitionTree$state, levels=c("A","B","C","D","E","F","G","H","I"))

# Subset by correctness and Tree and create 2 "longdata" dfs
correctGroupTree <- subset(dfANALYSISrep,Correct==1 & Format == "tree")
dfTransitionTreeCorrect <- data.frame(subject = NULL, state = NULL) 
rigaTC<-1 
for(p in correctGroupTree$Ps){
  scanna <- NULL
  scanna <- c(unlist(strsplit(toString(correctGroupTree[rigaTC,10]), split="")))
  for (n in scanna){
    newRow <- NULL
    newRow <- data.frame(subject = p, state = n)
    dfTransitionTreeCorrect <- rbind(dfTransitionTreeCorrect,newRow)
  }
  rigaTC <- rigaTC+1
}
dfTransitionTreeCorrect$state <- factor(dfTransitionTreeCorrect$state, levels=c("A","B","C","D","E","F","G","H","I"))

incorrectGroupTree <- subset(dfANALYSISrep,Correct==0 & Format == "tree")
dfTransitionTreeIncorrect <- data.frame(subject = NULL, state = NULL) 
rigaTI<-1 
for(p in incorrectGroupTree$Ps){
  #print(p)
  scanna <- NULL
  scanna <- c(unlist(strsplit(toString(incorrectGroupTree[rigaTI,10]), split="")))
  for (n in scanna){
    newRow <- NULL
    newRow <- data.frame(subject = p, state = n)
    dfTransitionTreeIncorrect <- rbind(dfTransitionTreeIncorrect,newRow)
  }
  rigaTI <- rigaTI+1
}
dfTransitionTreeIncorrect$state <- factor(dfTransitionTreeIncorrect$state, levels=c("A","B","C","D","E","F","G","H","I"))

# DFs VENN 
# Subset by Venn only to create a "longdata" df
GroupVenn <- subset(dfANALYSISrep,Format == "venn")
dfTransitionVenn <- data.frame(subject = NULL, state = NULL) 
rigaV<-1 
for(p in GroupVenn$Ps){
  #print(p)
  scanna <- c(unlist(strsplit(toString(GroupVenn[rigaV,10]), split="")))
  for (n in scanna){
    newRow <- data.frame(subject = p, state = n)
    dfTransitionVenn <- rbind(dfTransitionVenn,newRow)
  }
  rigaV <- rigaV+1
}
dfTransitionVenn$state <- factor(dfTransitionVenn$state, levels=c("L","M","N","O","P","Q"))

# Subset by correctness and Venn and create 2 "longdata" dfs
correctGroupVenn <- subset(dfANALYSISrep,Correct==1 & Format == "venn")
dfTransitionVennCorrect <- data.frame(subject = NULL, state = NULL) 
rigaVC<-1 # rigaTC
for(p in correctGroupVenn$Ps){ 
  scanna<-NULL
  scanna <- c(unlist(strsplit(toString(correctGroupVenn[rigaVC,10]), split="")))
  for (n in scanna){
    newRow<-NULL
    newRow <- data.frame(subject = p, state = n)
    dfTransitionVennCorrect <- rbind(dfTransitionVennCorrect,newRow)
  }
  rigaVC <- rigaVC+1
}
dfTransitionVennCorrect$state <- factor(dfTransitionVennCorrect$state, levels=c("L","M","N","O","P","Q"))

incorrectGroupVenn <- subset(dfANALYSISrep,Correct==0 & Format == "venn")
dfTransitionVennIncorrect <- data.frame(subject = NULL, state = NULL) 
rigaVI<-1  
for(p in incorrectGroupVenn$Ps){ # 
  #print(p)
  scanna <- NULL
  scanna <- c(unlist(strsplit(toString(incorrectGroupVenn[rigaVI,10]), split="")))
  for (n in scanna){
    newRow<- NULL
    newRow <- data.frame(subject = p, state = n)
    dfTransitionVennIncorrect <- rbind(dfTransitionVennIncorrect,newRow)
  }
  rigaVI <- rigaVI+1
}

dfTransitionVennIncorrect$state <- factor(dfTransitionVennIncorrect$state, levels=c("L","M","N","O","P","Q"))













#---------------------------------------------------------------------------------
# FUNCTION:     convertToMatrix()
# INPUT:        longdata DF
# OUTPUT:       a probability matrix 
# DESCRIPTION:  This function returns a probability matric for the transitions
#               of the AOIs (first order Markov Chain, for 2-grams).
# e.g. 
# convertToMatrix(dfTransitionTreeCorrect)
# convertToMatrix(dfTransitionTreeIncorrect)
#---------------------------------------------------------------------------------
convertToMatrix <- function(longdata){
  
  tmp0 <- longdata
  tmp <- tmp0 %>% group_by(subject) %>% mutate(to=lead(state))
  tmp2 <- tmp[complete.cases(tmp),]
  with(tmp2, table(state, to))
  outmat <- as.matrix(with(tmp2, table(state, to)))
  probMat <- outmat / rowSums(outmat)
  return(probMat)
  
}




#---------------------------------------------------------------------------------
# FUNCTION:     convertToPrior()
# INPUT:        probability matrix (transition matrix)
# OUTPUT:       a probability matrix where the zeros are replaced 
# DESCRIPTION:  This function returns a probability matrix for the transitions
#               of the AOIs (first order Markov Chain, for 2-grams) where the entrances 
#               with zero probabilities borrow some probability from the other entrances  
# Details:
# KL divergence on multinomials is defined when they have 
# only nonzero entries. When there are zero entries, 
# you have two choices:
# (1) Smooth the distributions in some way, for instance with a Bayesian prior, 
# or (similarly) (2) taking the convex combination of the observation with some 
# valid (nonzero) distribution.

# Bayesian prior (Dirichlet prior) method:
# For each entry:
# ni/m (ni should be an integer so we multiply it x10)
# Now replace the entrance with (ni+1)/(m+|x|)
# where m=∑i(ni) , |x|= number of outcomes (thye length of the row vector)
#---------------------------------------------------------------------------------
convertToPrior <- function(matrice){
  m <- 0
  x <- nrow(matrice) # lenght of the matrix
  for(i in 1:x){
    m <- sum(matrice[i,]*10)
    for(j in 1:x){
      matrice[i,j] <- ((matrice[i,j]*10)+1)/(m+x)
    }
  }
  return(matrice)
}
  
# Further comments:
#
# (1) Dirichlet prior proved to be numerically unstable because they assume a uniform prior.  
# see https://arxiv.org/pdf/1703.10364.pdf 
#
# (2) Bayesian methods are widely studied in 
# the context of medical decision making (7). 
# One advantage of using Bayesian methods in 
# the context of PSA for Markov models is in 
# overcoming the problem of estimating transition 
# probabilities when zero counts are observed 
# for the event of interest.
# If in the specification 
# of distributions for transition probabilities the 
# parameters are estimated uniquely from observed 
# event counts, zero probability will be associated to 
# those transitions that have never been observed. 
# With a Bayesian approach, it is natural to combine 
# data with prior information and, as a result, a non 
# null probability, even if small, will be associated 
# with the unobserved transitions.
# See http://ijphjournal.it/article/viewFile/7537/6796 
# - Armero C, Garcia-Donato G, Lopez-Quilez A. Bayesian 
# methods in cost-effectiveness studies: objectivity, 
# computation and other relevant aspects. Health Econ 
# 2010; 19: 629-43.
# - Briggs A, Ades A, Price M. Probabilistic sensitivity 
# analysis for decision trees with multiple branches: use of 
# the Dirichlet distribution in a Bayesian framework. Med 
# Decis Making 2003; 23: 341–50.




#---------------------------------------------------------------------------------
# FUNCTION:     calculateJSDistance()
# INPUT:        matrix, matrix
# OUTPUT:       JSD and JSDistance
# DESCRIPTION:  Computes correct Jensen Shannon distance between two matrices. 
#               
#---------------------------------------------------------------------------------
calculateJSDistance <- function(m1, m2){
  results <- list()
  JSDistance <- 0
  summedJSDistance <- 0
  matrix_length <- nrow(m1)
  coeff1 <- 0
  coeff2 <- 0
  for(i in 1:matrix_length){
    for(j in 1:matrix_length){
      # calculate KLD per row
      coeff1 <- coeff1 + (m1[i,j] * log(m1[i,j] / (0.5*(m1[i,j] +  m2[i,j]))))
      coeff2 <- coeff2 + (m2[i,j] * log(m2[i,j] / (0.5*(m1[i,j] +  m2[i,j])))) 
    }
    JSDistance <- sqrt(0.5*(coeff1+coeff2))
    summedJSDistance <- summedJSDistance + JSDistance
    JSDistance <- 0
    coeff1 <- 0
    coeff2 <- 0
  }
  averageDistance <- summedJSDistance/matrix_length
  return(averageDistance)
}



#---------------------------------------------------------------------------------
# FUNCTION:     GENERALrandomPermutations()
# parameters: 
# data <- e.g. dfTransitionTree
# group1_size <- number of correct (e.g. nrow(GroupTree[GroupTree$Correct==1,]), 30, length(correctGroupTree$Ps))
# distanceMeasure <- the type of distance to be used  (e.g. BhatDistance)
# perm <- 10000 (overraid the previous value of 10)
# convertPrior <- TRUE (used for JSDistance) 
#---------------------------------------------------------------------------------

GENERALrandomPermutations <- function(data, group1_size, distanceMeasure, perms = 10, convertPrior = FALSE){
  distance <- NULL
  for (i in 1:perms){
    group1_participants<- NULL
    group2_participants<- NULL
    group1<- NULL
    group2<- NULL
    m_1<- NULL
    m_2<- NULL
    distance_results<-NULL
    
    # create 2 subDF with random groups (not correct and incorrect) but the same size of correct and incorrect 
    group1_participants <- sample(unique(data$subject), group1_size)
    group2_participants <- unique(data$subject)[!(unique(data$subject) %in% group1_participants)]
    group1 <- subset(data, subject %in% group1_participants)
    group2 <- subset(data, subject %in% group2_participants)
  
    if (convertPrior){
      m_1 <- convertToPrior(convertToMatrix(group1))
      m_2 <- convertToPrior(convertToMatrix(group2))
    }
    else{
      # create 2 matrices (Markov chains)
      m_1 <- convertToMatrix(group1)
      m_2 <- convertToMatrix(group2)
    }
    
    # get distance 
    distance_results <- distanceMeasure(m_1, m_2)
    
    # store computed distance in vector and return
    distance <- c(distance, distance_results)
  }
  
  return(distance)  
  
}


#---------------------------------------------------------------------------------
# FUNCTION:     calculatePvalue()
# INPUT:        vector, vector
# OUTPUT:       void
# DESCRIPTION:  Calculate p-value (% of values > correct/incorrect value)
#               
#---------------------------------------------------------------------------------
calculatePvalue <- function(correct_and_incorrect, shuffled_distances)
{
  gtr <- length(shuffled_distances[shuffled_distances > correct_and_incorrect])
  pvalue <- gtr / length(shuffled_distances) 
  return(pvalue)
}









#---------------------------------------------------------------------------------#
#                                                                                 #
#        Create Transition Matrces 1st order                                      #    
#                                                                                 #
#---------------------------------------------------------------------------------#


# TREE
# Correct
matTC <- convertToMatrix(dfTransitionTreeCorrect) 
  
# Incorrect
matTI <- convertToMatrix(dfTransitionTreeIncorrect) 

# VENN
# Correct
matVC <- convertToMatrix(dfTransitionVennCorrect)
  
# Incorrect
matVI <- convertToMatrix(dfTransitionVennIncorrect) 




# Plot the two transition matrices for TREE for correct and incorrect--------- #



# tree_min <- min(matTC, matTI)
# tree_max <- max(matTC, matTI)

levelplot(matTC, col.regions = colorpanel(1000, "white", "grey15"),
          at = unique(c(seq(from = 0, to = 1, by = 0.001))),
          main = "Tree Correct", ylab = "AOI (destination)", xlab = "AOI (Origin)")
levelplot(matTI, col.regions = colorpanel(1000, "white", "grey15"),
          at = unique(c(seq(from = 0, to = 1, by = 0.001))),
          main = "Tree Incorrect", ylab = "AOI (destination)", xlab = "AOI (Origin)")

# Plot the two transition matrices for VENN for correct and incorrect--------- #

# venn_min <- min(matVC, matVI)
# venn_max <- max(matVC, matVI)

levelplot(matVC, col.regions = colorpanel(1000, "white", "grey15"),
          at = unique(c(seq(from = 0, to = 1, by = 0.001))),
          main = "Venn Correct", ylab = "AOI (destination)", xlab = "AOI (Origin)")
levelplot(matVI, col.regions = colorpanel(1000, "white", "grey15"),
          at = unique(c(seq(from = 0, to = 1, by = 0.001))),
          main = "Venn Incorrect", ylab = "AOI (destination)", xlab = "AOI (Origin)")
















total_permutations <- 10 # WE DEFINE THE TOTAL PERMUTATIONS FOR THE WHOLE EXPERIMENT
#---------------------------------------------------------------------------------#
#                                                                                 #
#   HERE I AM RUNNING THE MARKOV CHAIN PERMUTATION TESTS                          #    
#                                                                                 #
#---------------------------------------------------------------------------------#


#--------------- TREE

TREEnumber_correct <- length(correctGroupTree$Ps)
TREEnumber_incorrect <- length(incorrectGroupTree$Ps)

# store chains
matrix_tree_1 <- convertToPrior(matTC)  # I am using here convertToPrior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
matrix_tree_2 <- convertToPrior(matTI)  # I am using here convertToPrior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# work out distance between correct and incorrect groups
TREEcorrect_and_incorrect <- calculateJSDistance(matrix_tree_1, matrix_tree_2)
TREEylabel <- "Jensen-Shannon Distance"

TREEshuffled_distances <- GENERALrandomPermutations(dfTransitionTree, 
                                                         TREEnumber_correct,calculateJSDistance,
                                                         perms = total_permutations, convertPrior = TRUE)


# generate dot plot
TREEshuffled_distances <- append(TREEshuffled_distances, TREEcorrect_and_incorrect)
plot(TREEshuffled_distances, col = ifelse(TREEshuffled_distances == TREEcorrect_and_incorrect, 'blue', 'red'), 
     main = paste0("Tree Diagram ", " (",TREEylabel,")"), xlab = "Number of permutations", ylab = TREEylabel)
abline(h = TREEcorrect_and_incorrect, col = "purple")

# generate density plot
TREEdensity_plot <- density(TREEshuffled_distances)
plot(TREEdensity_plot, type = "n", main = paste0("Tree Diagram ", " (",TREEylabel,")"), xlab = "Jensen-Shannon Distance")
polygon(TREEdensity_plot, col = "lightgray", border = "grey")
rug(TREEshuffled_distances, col = ifelse(TREEshuffled_distances == TREEcorrect_and_incorrect, 'blue', 'red'))
abline(v = TREEcorrect_and_incorrect, col = "purple")

# results
JSDpvalueTree <- calculatePvalue(TREEcorrect_and_incorrect, TREEshuffled_distances)
JSDResultTREE <- list("The JSDistance distance is:", TREEcorrect_and_incorrect, "The p-value is:", JSDpvalueTree)
print(JSDResultTREE)





#------------ VENN                  


number_correct_venn <- length(correctGroupVenn$Ps) 
number_incorrect_venn <- length(incorrectGroupVenn$Ps)

# store chains
matrix_venn_1 <- convertToPrior(matVC)  # I am using here convertToPrior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
matrix_venn_2 <- convertToPrior(matVI)  # I am using here convertToPrior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# work out distance between correct and incorrect groups
VENNcorrect_and_incorrect <- calculateJSDistance(matrix_venn_1, matrix_venn_2) 
ylabel_venn <- "Jensen-Shannon Distance" 

VENNshuffled_distances <- GENERALrandomPermutations(dfTransitionVenn, 
                                                         number_correct_venn, calculateJSDistance,
                                                         perms = total_permutations, convertPrior = TRUE)


# generate dot plot
VENNshuffled_distances <- append(VENNshuffled_distances, VENNcorrect_and_incorrect)
plot(VENNshuffled_distances, col = ifelse(VENNshuffled_distances == VENNcorrect_and_incorrect, 'blue', 'red'), 
     main = paste0("Venn Diagram ", " (",ylabel_venn,")"), xlab = "Number of permutations", ylab = ylabel_venn)
abline(h = VENNcorrect_and_incorrect, col = "purple")

# generate density plot
density_plot_venn <- density(VENNshuffled_distances)
plot(density_plot_venn, type = "n", main = paste0("Venn Diagram ", " (",ylabel_venn,")"), xlab = "Jensen-Shannon Distance") 
polygon(density_plot_venn, col = "lightgray", border = "grey")
rug(VENNshuffled_distances, col = ifelse(VENNshuffled_distances == VENNcorrect_and_incorrect, 'blue', 'red'))
abline(v = VENNcorrect_and_incorrect, col = "purple")

# results
JSDpvalueVenn <- calculatePvalue(VENNcorrect_and_incorrect, VENNshuffled_distances)
JSDResultsVENN <- list("The JSDistance distance is:", VENNcorrect_and_incorrect, "The p-value is:", JSDpvalueVenn)
print(JSDResultsVENN)










###################################################################################
#---------------------------------------------------------------------------------#
#                                                                                 #
#                 Implementation of 2 other distance measures                     #    
#                                                                                 #
#---------------------------------------------------------------------------------#
###################################################################################

#---------------------------------------------------------------------------------
# FUNCTION:     BhatDistance()
# INPUT:        matrix, matrix
# OUTPUT:       Summed Bhattacharyya Distance
# DESCRIPTION:  Bhattacharyya distance is defined as: -ln(∑√(P(x)Q(x)))
#---------------------------------------------------------------------------------
BhatDistance <- function(m1,m2){
  length_of_matrix <- nrow(m1)
  DB <- 0
  BC <- 0
  for (i in 1:length_of_matrix){
    for (j in 1:length_of_matrix){
      BC <- BC + sqrt(m1[i, j]*m2[i, j])
    }
    DB<-DB + (-log(BC))
    BC <- 0
  }
  return(DB/length_of_matrix)
}
# Comments:
# Bhattacharyya distance measures the similarity 
# of two discrete or continuous probability distributions. 
# It is closely related to the Bhattacharyya coefficient 
# which is a measure of the amount of overlap between 
# two statistical samples or populations.
# This distance does not obey the triangle inequality, 
# but the Hellinger distance does obey the triangle inequality.

#---------------------------------------------------------------------------------
# FUNCTION:     HellDistance()
# INPUT:        matrix, matrix
# OUTPUT:       Summed Hellinger Distance
# DESCRIPTION:  Hellinger Distance is defined as: 1/sqrt(2) * sqrt(sum(suqare(sqrt(pi - squrt(qi))))
#---------------------------------------------------------------------------------
HellDistance <- function(m1,m2){
  length_of_matrix <- nrow(m1)
  HPQ <- 0
  HJ <- 0
  for (i in 1:length_of_matrix){
    for (j in 1:length_of_matrix){
      HJ <- HJ + ((sqrt(m1[i, j]) - sqrt(m2[i, j]))^2)
    }
    HPQ<-HPQ + ((1/sqrt(2))*sqrt(HJ))
    HJ <- 0
  }
  return(HPQ/length_of_matrix) 
}

# Comments:
# Hellinger distance (closely related to, although different from, the Bhattacharyya distance) 
# is used to quantify the similarity between two probability distributions.







#---------------------------------------------------------------------------------#
#                                                                                 #
#   HERE I AM RUNNING THE PERMUTATION TESTS FOR Bhattacharyya distance            #    
#                                                                                 #
#---------------------------------------------------------------------------------#

#----------------- TREE
# store chains
Bhatmatrix_tree_1 <- matTC 
Bhatmatrix_tree_2 <- matTI 

# work out distance between correct and incorrect groups
BhatDistance_resultTREE <- BhatDistance(Bhatmatrix_tree_1, Bhatmatrix_tree_2)
BhatTREEylabel <- "Bhattacharyya distance"

BhatTREEshuffled_distances <- GENERALrandomPermutations(dfTransitionTree, 
                                                     TREEnumber_correct,
                                                     BhatDistance,
                                             perms = total_permutations, convertPrior = FALSE)


# generate dot plot
BhatTREEshuffled_distances <- append(BhatTREEshuffled_distances, BhatDistance_resultTREE)
plot(BhatTREEshuffled_distances, col = ifelse(BhatTREEshuffled_distances == BhatDistance_resultTREE, 'blue', 'red'), 
     main = paste0("Tree Diagram ", " (",BhatTREEylabel,")"), xlab = "Number of permutations", ylab = BhatTREEylabel)
abline(h = BhatDistance_resultTREE, col = "purple")

# generate density plot
BhatTREEdensity_plot <- density(BhatTREEshuffled_distances)
plot(BhatTREEdensity_plot, type = "n", main = paste0("Tree Diagram ", " (",BhatTREEylabel,")"), xlab = "Bhattacharyya distance")
polygon(BhatTREEdensity_plot, col = "lightgray", border = "grey")
rug(BhatTREEshuffled_distances, col = ifelse(BhatTREEshuffled_distances == BhatDistance_resultTREE, 'blue', 'red'))
abline(v = BhatDistance_resultTREE, col = "purple")

BhatpvalueTREE <- calculatePvalue(BhatDistance_resultTREE, BhatTREEshuffled_distances)
BhatResultsTREE <- list("The Bhattacharyya distance is:", BhatDistance_resultTREE, "The p-value is:", BhatpvalueTREE)
print(BhatResultsTREE)



#---------------------- VENN
# store chains
Bhatmatrix_venn_1 <- matVC 
Bhatmatrix_venn_2 <- matVI 

# work out distance between correct and incorrect groups
BhatDistance_resultVENN <- BhatDistance(Bhatmatrix_venn_1, Bhatmatrix_venn_2)
BhatVENNylabel <- "Bhattacharyya distance"

BhatVENNshuffled_distances <- GENERALrandomPermutations(dfTransitionVenn, 
                                                     number_correct_venn,
                                                     BhatDistance,
                                                     perms = total_permutations, convertPrior = FALSE)


# generate dot plot
BhatVENNshuffled_distances <- append(BhatVENNshuffled_distances, BhatDistance_resultVENN)
plot(BhatVENNshuffled_distances, col = ifelse(BhatVENNshuffled_distances == BhatDistance_resultVENN, 'blue', 'red'), 
     main = paste0("Venn Diagram ", " (",BhatVENNylabel,")"), xlab = "Number of permutations", ylab = BhatVENNylabel)
abline(h = BhatDistance_resultVENN, col = "purple")

# generate density plot
BhatVENNdensity_plot <- density(BhatVENNshuffled_distances) 
plot(BhatVENNdensity_plot, type = "n", main = paste0("Venn Diagram ", " (",BhatVENNylabel,")"), xlab = "Bhattacharyya distance")
polygon(BhatVENNdensity_plot, col = "lightgray", border = "grey")
rug(BhatVENNshuffled_distances, col = ifelse(BhatVENNshuffled_distances == BhatDistance_resultVENN, 'blue', 'red'))
abline(v = BhatDistance_resultVENN, col = "purple")

BhatpvalueVENN <- calculatePvalue(BhatDistance_resultVENN, BhatVENNshuffled_distances)
BhatResultsVENN <- list("The Bhattacharyya distance is:", BhatDistance_resultVENN, "The p-value is:", BhatpvalueVENN)
print(BhatResultsVENN)





#---------------------------------------------------------------------------------#
#                                                                                 #
#   HERE I AM RUNNING THE PERMUTATION TESTS FOR Hellinger Distance                #    
#                                                                                 #
#---------------------------------------------------------------------------------#

#-------------------- TREE
# store chains
Hellmatrix_tree_1 <- matTC 
Hellmatrix_tree_2 <- matTI 

# work out distance between correct and incorrect groups
HellDistance_resultTREE <- HellDistance(Hellmatrix_tree_1, Hellmatrix_tree_2)
HellTREEylabel <- "Hellinger distance"

HellTREEshuffled_distances <- GENERALrandomPermutations(dfTransitionTree, 
                                                        TREEnumber_correct,
                                                        HellDistance,
                                                        perms = total_permutations, convertPrior = FALSE)


# generate dot plot
HellTREEshuffled_distances <- append(HellTREEshuffled_distances, HellDistance_resultTREE)
plot(HellTREEshuffled_distances, col = ifelse(HellTREEshuffled_distances == HellDistance_resultTREE, 'blue', 'red'), 
     main = paste0("Tree Diagram ", " (",HellTREEylabel,")"), xlab = "Number of permutations", ylab = HellTREEylabel)
abline(h = HellDistance_resultTREE, col = "purple")

# generate density plot
HellTREEdensity_plot <- density(HellTREEshuffled_distances)
plot(HellTREEdensity_plot, type = "n", main = paste0("Tree Diagram ", " (",HellTREEylabel,")"), xlab = "Hellinger distance")
polygon(HellTREEdensity_plot, col = "lightgray", border = "grey")
rug(HellTREEshuffled_distances, col = ifelse(HellTREEshuffled_distances == HellDistance_resultTREE, 'blue', 'red'))
abline(v = HellDistance_resultTREE, col = "purple")

HellpvalueTREE <- calculatePvalue(HellDistance_resultTREE, HellTREEshuffled_distances)
HellResultsTREE <- list("The Hellinger distance is:", HellDistance_resultTREE, "The p-value is:", HellpvalueTREE)
print(HellResultsTREE)


#------------------------- VENN
# store chains
Hellmatrix_VENN_1 <- matVC 
Hellmatrix_VENN_2 <- matVI 

# work out distance between correct and incorrect groups
HellDistance_resultVENN <- HellDistance(Hellmatrix_VENN_1, Hellmatrix_VENN_2)
HellVENNylabel <- "Hellinger distance"

HellVENNshuffled_distances <- GENERALrandomPermutations(dfTransitionVenn, 
                                                        number_correct_venn, 
                                                        HellDistance,
                                                        perms = total_permutations, convertPrior = FALSE)


# generate dot plot
HellVENNshuffled_distances <- append(HellVENNshuffled_distances, HellDistance_resultVENN)
plot(HellVENNshuffled_distances, col = ifelse(HellVENNshuffled_distances == HellDistance_resultVENN, 'blue', 'red'), 
     main = paste0("Venn Diagram ", " (",HellVENNylabel,")"), xlab = "Number of permutations", ylab = HellVENNylabel)
abline(h = HellDistance_resultVENN, col = "purple")

# generate density plot
HellVENNdensity_plot <- density(HellVENNshuffled_distances)
plot(HellVENNdensity_plot, type = "n", main = paste0("Venn Diagram ", " (",HellVENNylabel,")"), xlab = "Hellinger distance")
polygon(HellVENNdensity_plot, col = "lightgray", border = "grey")
rug(HellVENNshuffled_distances, col = ifelse(HellVENNshuffled_distances == HellDistance_resultVENN, 'blue', 'red'))
abline(v = HellDistance_resultVENN, col = "purple")

HellpvalueVENN <- calculatePvalue(HellDistance_resultVENN, HellVENNshuffled_distances)
HellResultsVENN <- list("The Hellinger distance is:", HellDistance_resultVENN, "The p-value is:", HellpvalueVENN)
print(HellResultsVENN)






####################################################################################################################
####################################################################################################################
#                                                                                                                  #
#  I NEED TO LOOK AT THE ANALYSIS OF TRANSITION MATRICES WITHOUT THE BACKGROUND and no repetitions)                #
#                                                                                                                  # 
####################################################################################################################
####################################################################################################################



#---------------------------------------------------------------------------------#
#                                                                                 #
#   Create the DFs with no background and  no repetitions                         #    
#                                                                                 #
#---------------------------------------------------------------------------------#

NEWdfANALYSIS <- dfANALYSISrep[,1:10] 
vettoreStringhe <- NULL
vettoreStringhe2 <- NULL
for (i in 1:length(dfANALYSISrep[,10])){
  vettoreStringhe <- c(vettoreStringhe, gsub('[AL]', '', NEWdfANALYSIS[i,10]))    
}

# This code is to eliminate repetitions 
for (i in 1:length(vettoreStringhe)){
  vettoreStringhe2 <- c(vettoreStringhe2, gsub('([[:alpha:]])\\1+', '\\1', vettoreStringhe[i]))    
}

NEWdfANALYSIS["Scanpath2"] <- vettoreStringhe2 

# rimpiazzare dfANALYSISrep con NEWdfANALYSIS & 10 (Scanpath) con 11 (Scanpath2) 

# DFs TREE 
# Subset by Tree only to create a "longdata" df 
unGroupTree <- subset(NEWdfANALYSIS,Format == "tree")                                 
undfTransitionTree <- data.frame(subject = NULL, state = NULL) 
rigaT<-1 
for(p in unGroupTree$Ps){
  #print(p)
  scanna <- c(unlist(strsplit(toString(unGroupTree[rigaT,11]), split="")))
  for (n in scanna){
    newRow <- data.frame(subject = p, state = n)
    undfTransitionTree <- rbind(undfTransitionTree,newRow)
  }
  rigaT <- rigaT+1
}
undfTransitionTree$state <- factor(undfTransitionTree$state, levels=c("B","C","D","E","F","G","H","I"))

# Subset by correctness and Tree and create 2 "longdata" dfs
uncorrectGroupTree <- subset(NEWdfANALYSIS,Correct==1 & Format == "tree")
undfTransitionTreeCorrect <- data.frame(subject = NULL, state = NULL) 
rigaTC<-1 
for(p in uncorrectGroupTree$Ps){
  scanna <- NULL
  scanna <- c(unlist(strsplit(toString(uncorrectGroupTree[rigaTC,11]), split="")))
  for (n in scanna){
    newRow <- NULL
    newRow <- data.frame(subject = p, state = n)
    undfTransitionTreeCorrect <- rbind(undfTransitionTreeCorrect,newRow)
  }
  rigaTC <- rigaTC+1
}
undfTransitionTreeCorrect$state <- factor(undfTransitionTreeCorrect$state, levels=c("B","C","D","E","F","G","H","I"))

unincorrectGroupTree <- subset(NEWdfANALYSIS,Correct==0 & Format == "tree")
undfTransitionTreeIncorrect <- data.frame(subject = NULL, state = NULL) 
rigaTI<-1 
for(p in unincorrectGroupTree$Ps){
  #print(p)
  scanna <- NULL
  scanna <- c(unlist(strsplit(toString(unincorrectGroupTree[rigaTI,11]), split="")))
  for (n in scanna){
    newRow <- NULL
    newRow <- data.frame(subject = p, state = n)
    undfTransitionTreeIncorrect <- rbind(undfTransitionTreeIncorrect,newRow)
  }
  rigaTI <- rigaTI+1
}
undfTransitionTreeIncorrect$state <- factor(undfTransitionTreeIncorrect$state, levels=c("B","C","D","E","F","G","H","I"))

# DFs VENN 
# Subset by Venn only to create a "longdata" df
unGroupVenn <- subset(NEWdfANALYSIS,Format == "venn")
undfTransitionVenn <- data.frame(subject = NULL, state = NULL) 
rigaV<-1 
for(p in unGroupVenn$Ps){
  #print(p)
  scanna <- c(unlist(strsplit(toString(unGroupVenn[rigaV,11]), split="")))
  for (n in scanna){
    newRow <- data.frame(subject = p, state = n)
    undfTransitionVenn <- rbind(undfTransitionVenn,newRow)
  }
  rigaV <- rigaV+1
}
undfTransitionVenn$state <- factor(undfTransitionVenn$state, levels=c("M","N","O","P","Q"))

# Subset by correctness and Venn and create 2 "longdata" dfs
uncorrectGroupVenn <- subset(NEWdfANALYSIS,Correct==1 & Format == "venn")
undfTransitionVennCorrect <- data.frame(subject = NULL, state = NULL) 
rigaVC<-1 # rigaTC
for(p in uncorrectGroupVenn$Ps){ 
  scanna<-NULL
  scanna <- c(unlist(strsplit(toString(uncorrectGroupVenn[rigaVC,11]), split="")))
  for (n in scanna){
    newRow<-NULL
    newRow <- data.frame(subject = p, state = n)
    undfTransitionVennCorrect <- rbind(undfTransitionVennCorrect,newRow)
  }
  rigaVC <- rigaVC+1
}
undfTransitionVennCorrect$state <- factor(undfTransitionVennCorrect$state, levels=c("M","N","O","P","Q"))

unincorrectGroupVenn <- subset(NEWdfANALYSIS,Correct==0 & Format == "venn")
undfTransitionVennIncorrect <- data.frame(subject = NULL, state = NULL) 
rigaVI<-1  
for(p in unincorrectGroupVenn$Ps){  
  #print(p)
  scanna <- NULL
  scanna <- c(unlist(strsplit(toString(unincorrectGroupVenn[rigaVI,11]), split="")))
  for (n in scanna){
    newRow<- NULL
    newRow <- data.frame(subject = p, state = n)
    undfTransitionVennIncorrect <- rbind(undfTransitionVennIncorrect,newRow)
  }
  rigaVI <- rigaVI+1
}

undfTransitionVennIncorrect$state <- factor(undfTransitionVennIncorrect$state, levels=c("M","N","O","P","Q"))













#---------------------------------------------------------------------------------#
#                                                                                 #
#   Create Transition Matrces 1st order with no Background and no repetitions     #    
#                                                                                 #
#---------------------------------------------------------------------------------#


# TREE
# Correct
NEWmatTC <- convertToMatrix(undfTransitionTreeCorrect) 

# Incorrect
NEWmatTI <- convertToMatrix(undfTransitionTreeIncorrect) 

# VENN
# Correct
NEWmatVC <- convertToMatrix(undfTransitionVennCorrect)

# Incorrect
NEWmatVI <- convertToMatrix(undfTransitionVennIncorrect)

# Plot the two transition matrices for TREE for correct and incorrect

# NEWtree_max <- max(NEWmatTC, NEWmatTI)

levelplot(NEWmatTC, col.regions = colorpanel(1000, "white", "grey15"),
          at = unique(c(seq(from = 0, to = 1, by = 0.001))),
          main = "Tree Correct", ylab = "AOI (destination)", xlab = "AOI (Origin)")
levelplot(NEWmatTI, col.regions = colorpanel(1000, "white", "grey15"),
          at = unique(c(seq(from = 0, to = 1, by = 0.001))),
          main = "Tree Incorrect", ylab = "AOI (destination)", xlab = "AOI (Origin)")


# Plot the two transition matrices for VENN for correct and incorrect--------- #


levelplot(NEWmatVC, col.regions = colorpanel(1000, "white", "grey15"),
          at = unique(c(seq(from = 0, to = 1, by = 0.001))),
          main = "Venn Correct", ylab = "AOI (destination)", xlab = "AOI (Origin)")
levelplot(NEWmatVI, col.regions = colorpanel(1000, "white", "grey15"),
          at = unique(c(seq(from = 0, to = 1, by = 0.001))),
          main = "Venn Incorrect", ylab = "AOI (destination)", xlab = "AOI (Origin)")






# ---------------------------------------------------------------------------------------#
#                                                                                        #
#                       MARKOV CHAIN TEST WITH JSDistance                                #
#                       no background and  no repetitions                                #
# ---------------------------------------------------------------------------------------#

#-------------------- TREE
# store chains
NEWmatrix_tree_1 <- convertToPrior(NEWmatTC)  # I am using here convertToPrior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NEWmatrix_tree_2 <- convertToPrior(NEWmatTI)  # I am using here convertToPrior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# work out distance between correct and incorrect groups
NEWTREEcorrect_and_incorrect <- calculateJSDistance(NEWmatrix_tree_1, NEWmatrix_tree_2)
NEWTREEylabel <- "Jensen-Shannon Distance"

NEWTREEshuffled_distances <- GENERALrandomPermutations(undfTransitionTree, 
                                             TREEnumber_correct,calculateJSDistance,
                                             perms = total_permutations, convertPrior = TRUE)


# generate dot plot
NEWTREEshuffled_distances <- append(NEWTREEshuffled_distances, NEWTREEcorrect_and_incorrect)
plot(NEWTREEshuffled_distances, col = ifelse(NEWTREEshuffled_distances == NEWTREEcorrect_and_incorrect, 'blue', 'red'), 
     main = paste0("Tree Diagram ", " (",NEWTREEylabel,")"), xlab = "Number of permutations", ylab = NEWTREEylabel)
abline(h = NEWTREEcorrect_and_incorrect, col = "purple")

# generate density plot
NEWTREEdensity_plot <- density(NEWTREEshuffled_distances)
plot(NEWTREEdensity_plot, type = "n", main = paste0("Tree Diagram ", " (",NEWTREEylabel,")"), xlab = "Jensen-Shannon Distance")
polygon(NEWTREEdensity_plot, col = "lightgray", border = "grey")
rug(NEWTREEshuffled_distances, col = ifelse(NEWTREEshuffled_distances == NEWTREEcorrect_and_incorrect, 'blue', 'red'))
abline(v = NEWTREEcorrect_and_incorrect, col = "purple")


# results
NEWJSDpvalueTree <- calculatePvalue(NEWTREEcorrect_and_incorrect, NEWTREEshuffled_distances)
NEWJSDResultTREE <- list("The JSDistance distance is:", NEWTREEcorrect_and_incorrect, "The p-value is:", NEWJSDpvalueTree)
print(NEWJSDResultTREE)



#------------------ VENN
# store chains
NEWmatrix_venn_1 <- convertToPrior(NEWmatVC)  # I am using here convertToPrior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NEWmatrix_venn_2 <- convertToPrior(NEWmatVI)  # I am using here convertToPrior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# work out distance between correct and incorrect groups
NEWVENNcorrect_and_incorrect <- calculateJSDistance(NEWmatrix_venn_1, NEWmatrix_venn_2)
NEWVENNylabel <- "Jensen-Shannon Distance"

NEWVENNshuffled_distances <- GENERALrandomPermutations(undfTransitionVenn, 
                                                number_correct_venn,calculateJSDistance,
                                                perms = total_permutations, convertPrior = TRUE)


# generate dot plot
NEWVENNshuffled_distances <- append(NEWVENNshuffled_distances, NEWVENNcorrect_and_incorrect)
plot(NEWVENNshuffled_distances, col = ifelse(NEWVENNshuffled_distances == NEWVENNcorrect_and_incorrect, 'blue', 'red'), 
     main = paste0("Venn Diagram ", " (",NEWVENNylabel,")"), xlab = "Number of permutations", ylab = NEWVENNylabel)
abline(h = NEWVENNcorrect_and_incorrect, col = "purple")

# generate density plot
NEWVENNdensity_plot <- density(NEWVENNshuffled_distances)
plot(NEWVENNdensity_plot, type = "n", main = paste0("Venn Diagram ", " (",NEWTREEylabel,")"), xlab = "Jensen-Shannon Distance")
polygon(NEWVENNdensity_plot, col = "lightgray", border = "grey")
rug(NEWVENNshuffled_distances, col = ifelse(NEWVENNshuffled_distances == NEWVENNcorrect_and_incorrect, 'blue', 'red'))
abline(v = NEWVENNcorrect_and_incorrect, col = "purple")

# results
NEWJSDpvalueVenn <- calculatePvalue(NEWVENNcorrect_and_incorrect, NEWVENNshuffled_distances)
NEWJSDResultVENN <- list("The JSDistance distance is:", NEWVENNcorrect_and_incorrect, "The p-value is:", NEWJSDpvalueVenn)
print(NEWJSDResultVENN)





#---------------------------------------------------------------------------------#
#                                                                                 #
#   HERE I AM RUNNING THE PERMUTATION TESTS FOR Bhattacharyya distance            #    
#           no background and  no repetitions                                     #
#---------------------------------------------------------------------------------#

#--------------------- TREE
# store chains
NEWBhatmatrix_tree_1 <- NEWmatTC 
NEWBhatmatrix_tree_2 <- NEWmatTI 

# work out distance between correct and incorrect groups
NEWBhatDistance_resultTREE <- BhatDistance(NEWBhatmatrix_tree_1, NEWBhatmatrix_tree_2)
BhatTREEylabel <- "Bhattacharyya distance"

NEWBhatTREEshuffled_distances <- GENERALrandomPermutations(undfTransitionTree, 
                                                        TREEnumber_correct,
                                                        BhatDistance,
                                                        perms = total_permutations, convertPrior = FALSE)


# generate dot plot
NEWBhatTREEshuffled_distances <- append(NEWBhatTREEshuffled_distances, NEWBhatDistance_resultTREE)
plot(NEWBhatTREEshuffled_distances, col = ifelse(NEWBhatTREEshuffled_distances == NEWBhatDistance_resultTREE, 'blue', 'red'), 
     main = paste0("Tree Diagram ", " (",BhatTREEylabel,")"), xlab = "Number of permutations", ylab = BhatTREEylabel)
abline(h = NEWBhatDistance_resultTREE, col = "purple")

# generate density plot
NEWBhatTREEdensity_plot <- density(NEWBhatTREEshuffled_distances)
plot(NEWBhatTREEdensity_plot, type = "n", main = paste0("Tree Diagram ", " (",BhatTREEylabel,")"), xlab = "Bhattacharyya distance")
polygon(NEWBhatTREEdensity_plot, col = "lightgray", border = "grey")
rug(NEWBhatTREEshuffled_distances, col = ifelse(NEWBhatTREEshuffled_distances == NEWBhatDistance_resultTREE, 'blue', 'red'))
abline(v = NEWBhatDistance_resultTREE, col = "purple")

NEWBhatpvalueTREE <- calculatePvalue(NEWBhatDistance_resultTREE, NEWBhatTREEshuffled_distances)
NEWBhatResultsTREE <- list("The Bhattacharyya distance is:", NEWBhatDistance_resultTREE, "The p-value is:", NEWBhatpvalueTREE)
print(NEWBhatResultsTREE)



#------------------- VENN
# store chains
NEWBhatmatrix_venn_1 <- NEWmatVC 
NEWBhatmatrix_venn_2 <- NEWmatVI 

# work out distance between correct and incorrect groups
NEWBhatDistance_resultVENN <- BhatDistance(NEWBhatmatrix_venn_1, NEWBhatmatrix_venn_2)
BhatVENNylabel <- "Bhattacharyya distance"

NEWBhatVENNshuffled_distances <- GENERALrandomPermutations(undfTransitionVenn, 
                                                        number_correct_venn,
                                                        BhatDistance,
                                                        perms = total_permutations, convertPrior = FALSE)


# generate dot plot
NEWBhatVENNshuffled_distances <- append(NEWBhatVENNshuffled_distances, NEWBhatDistance_resultVENN)
plot(NEWBhatVENNshuffled_distances, col = ifelse(NEWBhatVENNshuffled_distances == NEWBhatDistance_resultVENN, 'blue', 'red'), 
     main = paste0("Venn Diagram ", " (",BhatVENNylabel,")"), xlab = "Number of permutations", ylab = BhatVENNylabel)
abline(h = NEWBhatDistance_resultVENN, col = "purple")

# generate density plot
NEWBhatVENNdensity_plot <- density(NEWBhatVENNshuffled_distances) 
plot(NEWBhatVENNdensity_plot, type = "n", main = paste0("Venn Diagram ", " (",BhatVENNylabel,")"), xlab = "Bhattacharyya distance")
polygon(NEWBhatVENNdensity_plot, col = "lightgray", border = "grey")
rug(NEWBhatVENNshuffled_distances, col = ifelse(NEWBhatVENNshuffled_distances == NEWBhatDistance_resultVENN, 'blue', 'red'))
abline(v = NEWBhatDistance_resultVENN, col = "purple")

NEWBhatpvalueVENN <- calculatePvalue(NEWBhatDistance_resultVENN, NEWBhatVENNshuffled_distances)
NEWBhatResultsVENN <- list("The Bhattacharyya distance is:", NEWBhatDistance_resultVENN, "The p-value is:", NEWBhatpvalueVENN)
print(NEWBhatResultsVENN)





#---------------------------------------------------------------------------------#
#                                                                                 #
#   HERE I AM RUNNING THE PERMUTATION TESTS FOR Hellinger Distance                #    
#         no background and  no repetitions                                       #
#---------------------------------------------------------------------------------#

#------------------------ TREE
# store chains
NEWHellmatrix_tree_1 <- NEWmatTC 
NEWHellmatrix_tree_2 <- NEWmatTI 

# work out distance between correct and incorrect groups
NEWHellDistance_resultTREE <- HellDistance(NEWHellmatrix_tree_1, NEWHellmatrix_tree_2)
HellTREEylabel <- "Hellinger distance"

NEWHellTREEshuffled_distances <- GENERALrandomPermutations(undfTransitionTree, 
                                                        TREEnumber_correct,
                                                        HellDistance,
                                                        perms = total_permutations, convertPrior = FALSE)


# generate dot plot
NEWHellTREEshuffled_distances <- append(NEWHellTREEshuffled_distances, NEWHellDistance_resultTREE)
plot(NEWHellTREEshuffled_distances, col = ifelse(NEWHellTREEshuffled_distances == NEWHellDistance_resultTREE, 'blue', 'red'), 
     main = paste0("Tree Diagram ", " (",HellTREEylabel,")"), xlab = "Number of permutations", ylab = HellTREEylabel)
abline(h = NEWHellDistance_resultTREE, col = "purple")

# generate density plot
NEWHellTREEdensity_plot <- density(NEWHellTREEshuffled_distances)
plot(NEWHellTREEdensity_plot, type = "n", main = paste0("Tree Diagram ", " (",HellTREEylabel,")"), xlab = "Hellinger distance")
polygon(NEWHellTREEdensity_plot, col = "lightgray", border = "grey")
rug(NEWHellTREEshuffled_distances, col = ifelse(NEWHellTREEshuffled_distances == NEWHellDistance_resultTREE, 'blue', 'red'))
abline(v = NEWHellDistance_resultTREE, col = "purple")

NEWHellpvalueTREE <- calculatePvalue(NEWHellDistance_resultTREE, NEWHellTREEshuffled_distances)
NEWHellResultsTREE <- list("The Hellinger distance is:", NEWHellDistance_resultTREE, "The p-value is:", NEWHellpvalueTREE)
print(NEWHellResultsTREE)

#----------------------- VENN
# store chains
NEWHellmatrix_VENN_1 <- NEWmatVC 
NEWHellmatrix_VENN_2 <- NEWmatVI 

# work out distance between correct and incorrect groups
NEWHellDistance_resultVENN <- HellDistance(NEWHellmatrix_VENN_1, NEWHellmatrix_VENN_2)
HellVENNylabel <- "Hellinger distance"

NEWHellVENNshuffled_distances <- GENERALrandomPermutations(undfTransitionVenn, 
                                                        number_correct_venn, 
                                                        HellDistance,
                                                        perms = total_permutations, convertPrior = FALSE)


# generate dot plot
NEWHellVENNshuffled_distances <- append(NEWHellVENNshuffled_distances, NEWHellDistance_resultVENN)
plot(NEWHellVENNshuffled_distances, col = ifelse(NEWHellVENNshuffled_distances == NEWHellDistance_resultVENN, 'blue', 'red'), 
     main = paste0("Venn Diagram ", " (",HellVENNylabel,")"), xlab = "Number of permutations", ylab = HellVENNylabel)
abline(h = NEWHellDistance_resultVENN, col = "purple")

# generate density plot
NEWHellVENNdensity_plot <- density(NEWHellVENNshuffled_distances)
plot(NEWHellVENNdensity_plot, type = "n", main = paste0("Venn Diagram ", " (",HellVENNylabel,")"), xlab = "Hellinger distance")
polygon(NEWHellVENNdensity_plot, col = "lightgray", border = "grey")
rug(NEWHellVENNshuffled_distances, col = ifelse(NEWHellVENNshuffled_distances == NEWHellDistance_resultVENN, 'blue', 'red'))
abline(v = NEWHellDistance_resultVENN, col = "purple")

NEWHellpvalueVENN <- calculatePvalue(NEWHellDistance_resultVENN, NEWHellVENNshuffled_distances)
NEWHellResultsVENN <- list("The Hellinger distance is:", NEWHellDistance_resultVENN, "The p-value is:", NEWHellpvalueVENN)
print(NEWHellResultsVENN)

# ------------------------------------------------------------------------------------------------------------------











############################################################################################################################
#                                                                                                                          #
#            HERE WE LOOK AT VARIANCE IN THE DISTRIBUTION Of Transitions OF THE MARTIX (FOR CORRECT AND INCORRECT)         #
#                                                                                                                          #  
############################################################################################################################


# Calculate the average deviation between the two markov chain (correct and incorrect)
# The "qualvar" package has a function called ADA
# This takes a vector of values (in frequencies) and return an index of variation.
# from: 
# Wilcox, Allen R. ’Indices of Qualitative Variation and Political Measurement.’ The Western Political
# Quarterly 26, no. 2 (1 June 1973): 325-43. doi:10.2307/446831.


# examples
vettore1 <- c(10,1,0,0)
vettore2 <- c(4,4,4,4)
vettore3 <- c(10,0,0,0)
ADA(vettore1)
ADA(vettore2)
ADA(vettore3)


#---------------------------------------------------------------------------------
# FUNCTION:     matrixADA()
# INPUT:        matrix
# OUTPUT:       The Average Deviation Analog for the matrix
# DESCRIPTION:  Computes the the average deviation for a matrix
#---------------------------------------------------------------------------------
matrixADA <- function(m)
{
  devia <- 0
  summedDevia <- 0
  matrix_length <- nrow(m)
  for(i in 1:matrix_length){
    devia <- ADA(m[i,])
    #print(devia)
    summedDevia <- summedDevia + devia
    devia <- 0
  }
  return(summedDevia / matrix_length)           
}

#---------------------------------------------------------------------------------
# FUNCTION:     convertToMatrixFreq()
# INPUT:        longdata DF
# OUTPUT:       a Frequency matrix 
# DESCRIPTION:  This function returns a frequency matric for the transitions
#               of the AOIs (first order Markov Chain, for 2-grams).
#---------------------------------------------------------------------------------
convertToMatrixFreq <- function(longdata){
  
  tmp0 <- longdata
  tmp <- tmp0 %>% group_by(subject) %>% mutate(to=lead(state))
  tmp2 <- tmp[complete.cases(tmp),]
  with(tmp2, table(state, to))
  outmat <- as.matrix(with(tmp2, table(state, to)))
  return(outmat)
  
}

# Create the two matrices for Tree (Correct and Incorrect) from the dataframe with no repetitions

FreqmatriceTC <- convertToMatrixFreq(undfTransitionTreeCorrect) 
FreqmatriceTI <- convertToMatrixFreq(undfTransitionTreeIncorrect) 

matrixADA(FreqmatriceTC)
matrixADA(FreqmatriceTI)

FreqmatriceVC <- convertToMatrixFreq(undfTransitionVennCorrect) 
FreqmatriceVI <- convertToMatrixFreq(undfTransitionVennIncorrect) 

matrixADA(FreqmatriceVC)
matrixADA(FreqmatriceVI)







# WRITE THE DATA FILE FOR THE ANALYSIS IN THE SCRIPT CALLED bayesLCS.R 
write.csv(NEWdfANALYSIS, file = "bayeScanpathsrepNEW.csv", row.names=F)

























#TODOs:
# change per,mutation test to take into account also the kind of matrix generator that we need
# run the permutation test for Tree and Venn with and without background and repetitions for JSDistance

