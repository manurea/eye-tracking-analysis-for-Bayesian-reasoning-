#install.packages("stringi")
#install.packages("qualvar")
#install.packages("gtools")
#install.packages("ggplot2")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("plyr")
library(stringi)
library(qualvar)
library("gtools")
library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#                     EXPERIMENT 2B - Bayesian reasoning and eye-tracking analysis                      #
#                     Find the optimal subsequence (scanpath)                                          #
#                                    Manuele Reani                                                     #
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#




NEWdfANALYSIS <- read.csv(file="bayeScanpathsrepNEW.csv",header=TRUE,sep=",")

uncorrectGroupTree <- subset(NEWdfANALYSIS,Correct==1 & Format == "tree")
unincorrectGroupTree <- subset(NEWdfANALYSIS,Correct==0 & Format == "tree")

# WORK WITH 2-GRAM 

# create a vector with the number of 2gram with reduced alphabet 
#REDUCEDalphabetTree <- "BCDEFGHI" 

# Find all the 2-gram permutations with repetitions, escluding sunsequent repetitions, for the alphabet of Tree 
REDUCEDbiGramTree <- permutations(n = 8, r = 2, v = c("B", "C", "D", "E", "F", "G", "H", "I"), repeats.allowed = TRUE)

REDUCEDbiGramsT <- apply(REDUCEDbiGramTree,1,paste,collapse=" ")

REDUCEDtotalBiGrams__tree <- NULL 
for (bigram in REDUCEDbiGramsT){
  bigram <- gsub(" ", "", bigram)
  bigram <- gsub('([[:alpha:]])\\1+', '\\1', bigram)
  if(nchar(bigram) == 2){
    REDUCEDtotalBiGrams__tree <- append(REDUCEDtotalBiGrams__tree, bigram)
  }
}




# Find all the 2-gram permutations with repetitions, escluding sunsequent repetitions, for the alphabet of Venn 
REDUCEDbiGramVenn <- permutations(n = 5, r = 2, v = c("M", "N", "O", "P", "Q"), repeats.allowed = TRUE)

REDUCEDbiGramsV <- apply(REDUCEDbiGramVenn,1,paste,collapse=" ")

REDUCEDtotalBiGrams__venn <- NULL 
for (bigram in REDUCEDbiGramsV){
  bigram <- gsub(" ", "", bigram)
  bigram <- gsub('([[:alpha:]])\\1+', '\\1', bigram)
  if(nchar(bigram) == 2){
    REDUCEDtotalBiGrams__venn <- append(REDUCEDtotalBiGrams__venn, bigram)
  }
}





#---------------------------------------------------------------------------------
# FUNCTION:     subSeqfreq()
# INPUT:        ngrams = list of bigrams availables (e.g. REDUCEDtotalBiGrams__tree)
#               dataframe = the dataframe with all the scanpaths (e.g. uncorrectGroupTree)
#               columnNumber = the column number where the scanpaths are found (e.g. 11)
# OUTPUT:       A wide DF with the list of 2-gram and their frequency
# DESCRIPTION:  This function should produce a datframe where to store for each group 
# (correct and incorrect) the bi-gram and their respective frequencies
# Comments: Here I use the function "stri_count_regex" from the library "stringi"
#---------------------------------------------------------------------------------
subSeqfreq <- function(ngrams, dataframe, columnNumber){
  ngram_frquencyDF <- data.frame(n_gram = NULL, freqz = NULL)
  number_ngram <- length(ngrams)
  number_observ <- nrow(dataframe)
  freq <- 0
  for (i in 1:number_ngram){
    for (j in 1:number_observ){
      freq <- freq + stri_count_regex(dataframe[j,columnNumber], ngrams[i])
    }
    newRow <- data.frame(n_gram = ngrams[i], freqz = freq)
    ngram_frquencyDF <- rbind(ngram_frquencyDF,newRow)
    freq <- 0
  }
  return(ngram_frquencyDF) 
}



# produce the frequency data frame for subsequences for Tree (correct & incorrect)
freqSubsq_tree_correct <- subSeqfreq(REDUCEDtotalBiGrams__tree, uncorrectGroupTree, 11)
freqSubsq_tree_incorrect <- subSeqfreq(REDUCEDtotalBiGrams__tree, unincorrectGroupTree, 11)




# ---- combine the dataframe 
freqSubsq_tree_correct$correct <- rep("correct",nrow(freqSubsq_tree_correct))
freqSubsq_tree_incorrect$correct <- rep("incorrect",nrow(freqSubsq_tree_incorrect))
combineCorrectIncorrect <- rbind(freqSubsq_tree_correct,freqSubsq_tree_incorrect)

ggplot(combineCorrectIncorrect, aes(x=n_gram, y = freqz, fill=factor(correct), group=factor(correct))) + 
  geom_bar(stat="identity") + 
  labs(x = "2-gram", y = "Frequency")  
  #scale_y_continuous(breaks = seq(0, 28, by = 1), limits=c(0,28))
ggplot(combineCorrectIncorrect, aes(x = n_gram, y = freqz, fill = correct)) +
  geom_bar(stat="identity",colour="black",position = "dodge", alpha = 0.6)





# ------- freqz


# get the difference for each transition betwene corr and incorr and store it in a DF calle "subtractionFreqz"
combineCorrectIncorrectWide <- combineCorrectIncorrect %>% spread(n_gram,freqz, drop=FALSE)
subtractionFreqz <- data.frame(Trans = NULL, diff = NULL)
for (i in 2:57){
  #colnames(combineCorrectIncorrectWide)[i] <- combineCorrectIncorrectWide[!rowSums(is.na(combineCorrectIncorrectWide[i])), i]
  newRow <- data.frame(Trans = colnames(combineCorrectIncorrectWide)[i], 
                       diff =abs(diff(combineCorrectIncorrectWide[,i])))                         #(combineCorrectIncorrectWide[!rowSums(is.na(combineCorrectIncorrectWide[i])), i]) ))
  subtractionFreqz <- rbind(subtractionFreqz,newRow)
}

# print the top 5 most different transitions
print(subset(subtractionFreqz, diff %in% tail(sort(subtractionFreqz$diff),5)))
# get a vector with the names of the top 5
transizTop5 <- subset(subtractionFreqz, diff %in% tail(sort(subtractionFreqz$diff),5))$Trans
# create a df with the top 5
dfTop5 <- subset(subtractionFreqz, diff %in% tail(sort(subtractionFreqz$diff),5))

# reorder the label in crescendo
dfTop5$Trans <- factor(dfTop5$Trans, levels = dfTop5$Trans[order(dfTop5$diff)])
#plot the data
ggplot(dfTop5, aes(x = Trans, y = diff))+ geom_bar(stat="identity", width=0.8) + coord_flip() +
  xlab("2-gram") + ylab("Δ|Correct - Incorrect|")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))


#create a top5 DF fromt he original
combineCorrectIncorrectTop5 <- subset(combineCorrectIncorrect, n_gram %in% transizTop5)
# reorder the label in crescendo
combineCorrectIncorrectTop5$n_gram <- factor(combineCorrectIncorrectTop5$n_gram, levels = rev(dfTop5$Trans[order(dfTop5$diff)]))

#plot the data
ggplot(combineCorrectIncorrectTop5, aes(x = n_gram, y = freqz, fill = correct)) +
  geom_bar(stat="identity",colour="black",position = "dodge", width=0.8) +
  xlab("2-gram") + ylab("Count") +
  scale_fill_grey(start = 0.2, end = .5)+
  theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"))











# ------- probz 


# ------- create a dataaframe for probabilities transitions (probz), called combineCorrectIncorrectWide2

combineCorrectIncorrectWide2 <- combineCorrectIncorrectWide

# This code trasnform the data frequencies  into probabilities (relative frequencies), 
# over the total for each node in the chain (1 AOI) 

#combineCorrectIncorrectWide2[1,2:57] <- 0
#combineCorrectIncorrectWide2[2,2:57] <- 0

#combineCorrectIncorrectWide2[1,2:8]<-combineCorrectIncorrectWide[1,2:8]/sum(combineCorrectIncorrectWide[1,2:8])
#combineCorrectIncorrectWide2[1,9:15]<-combineCorrectIncorrectWide[1,9:15]/sum(combineCorrectIncorrectWide[1,9:15])
#combineCorrectIncorrectWide2[1,16:22]<-combineCorrectIncorrectWide[1,16:22]/sum(combineCorrectIncorrectWide[1,16:22])
#combineCorrectIncorrectWide2[1,23:29]<-combineCorrectIncorrectWide[1,23:29]/sum(combineCorrectIncorrectWide[1,23:29])
#combineCorrectIncorrectWide2[1,30:36]<-combineCorrectIncorrectWide[1,30:36]/sum(combineCorrectIncorrectWide[1,30:36])
#combineCorrectIncorrectWide2[1,37:43]<-combineCorrectIncorrectWide[1,37:43]/sum(combineCorrectIncorrectWide[1,37:43])
#combineCorrectIncorrectWide2[1,44:50]<-combineCorrectIncorrectWide[1,44:50]/sum(combineCorrectIncorrectWide[1,44:50])
#combineCorrectIncorrectWide2[1,51:57]<-combineCorrectIncorrectWide[1,51:57]/sum(combineCorrectIncorrectWide[1,51:57])

#combineCorrectIncorrectWide2[2,2:8]<-combineCorrectIncorrectWide[2,2:8]/sum(combineCorrectIncorrectWide[2,2:8])
#combineCorrectIncorrectWide2[2,9:15]<-combineCorrectIncorrectWide[2,9:15]/sum(combineCorrectIncorrectWide[2,9:15])
#combineCorrectIncorrectWide2[2,16:22]<-combineCorrectIncorrectWide[2,16:22]/sum(combineCorrectIncorrectWide[2,16:22])
#combineCorrectIncorrectWide2[2,23:29]<-combineCorrectIncorrectWide[2,23:29]/sum(combineCorrectIncorrectWide[2,23:29])
#combineCorrectIncorrectWide2[2,30:36]<-combineCorrectIncorrectWide[2,30:36]/sum(combineCorrectIncorrectWide[2,30:36])
#combineCorrectIncorrectWide2[2,37:43]<-combineCorrectIncorrectWide[2,37:43]/sum(combineCorrectIncorrectWide[2,37:43])
#combineCorrectIncorrectWide2[2,44:50]<-combineCorrectIncorrectWide[2,44:50]/sum(combineCorrectIncorrectWide[2,44:50])
#combineCorrectIncorrectWide2[2,51:57]<-combineCorrectIncorrectWide[2,51:57]/sum(combineCorrectIncorrectWide[2,51:57])


# transform the data into relative frequencies (relative to the total by group)
combineCorrectIncorrectWide2[1,2:57]<-combineCorrectIncorrectWide2[1,2:57]/sum(combineCorrectIncorrectWide2[1,2:57])
combineCorrectIncorrectWide2[2,2:57]<-combineCorrectIncorrectWide2[2,2:57]/sum(combineCorrectIncorrectWide2[2,2:57])


# get the difference for each transition betwene corr and incorr and store it in a DF calle "subtractionProbz"
subtractionProbz <- data.frame(Trans = NULL, diff = NULL)
for (i in 2:57){
  newRow <- data.frame(Trans = colnames(combineCorrectIncorrectWide2)[i], 
                       diff =abs(diff(combineCorrectIncorrectWide2[,i])))                        
  subtractionProbz <- rbind(subtractionProbz,newRow)
}

# print the top 5 most different transitions
print(subset(subtractionProbz, diff %in% tail(sort(subtractionProbz$diff),5)))
# get a vector with the names of the top 5
transizTop5Probz <- subset(subtractionProbz, diff %in% tail(sort(subtractionProbz$diff),5))$Trans 
# create a df with the top 5
dfTop5Probz <- subset(subtractionProbz, diff %in% tail(sort(subtractionProbz$diff),5)) 

# reorder the label in crescendo
dfTop5Probz$Trans <- factor(dfTop5Probz$Trans, levels = dfTop5Probz$Trans[order(dfTop5Probz$diff)]) # dfTop5
#plot the data
ggplot(dfTop5Probz, aes(x = Trans, y = diff))+ geom_bar(stat="identity", width=0.8) + coord_flip() +
  xlab("2-gram") + ylab("Δ|Correct - Incorrect|")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))


# from wide DF to narrow DF
combineCorrectIncorrect2 <- combineCorrectIncorrectWide2 %>% gather(probz, correct, convert = FALSE) # combineCorrectIncorrect
colnames(combineCorrectIncorrect2) <- c("correct","n_gram", "probz")

#create a top5 DF fromt he original
combineCorrectIncorrectTop5Probz <- subset(combineCorrectIncorrect2, n_gram %in% transizTop5Probz) 
# reorder the label in crescendo
combineCorrectIncorrectTop5Probz$n_gram <- factor(combineCorrectIncorrectTop5Probz$n_gram, levels = rev(dfTop5Probz$Trans[order(dfTop5Probz$diff)]))

#plot the data
ggplot(combineCorrectIncorrectTop5Probz, aes(x = n_gram, y = probz, fill = correct)) +
  geom_bar(stat="identity",colour="black",position = "dodge", width=0.8) +
  xlab("2-gram") + ylab("Prob") +
  scale_fill_grey(start = 0.2, end = .5)+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))


#------
# look at the 2-gram with the ghighest probabilities for both groups
correctTransitionsProbz <- combineCorrectIncorrect2[combineCorrectIncorrect2$correct=="correct",]
incorrectTransitionsProbz <- combineCorrectIncorrect2[combineCorrectIncorrect2$correct=="incorrect",]
correctTransitionsProbzTop10 <- correctTransitionsProbz[correctTransitionsProbz$probz %in% tail(sort(correctTransitionsProbz$probz), 10), ]
incorrectTransitionsProbzTop10 <- incorrectTransitionsProbz[incorrectTransitionsProbz$probz %in% tail(sort(incorrectTransitionsProbz$probz), 10), ]
correctTransitionsProbzTop10 <- correctTransitionsProbzTop10[order(-correctTransitionsProbzTop10$probz),]   
incorrectTransitionsProbzTop10 <- incorrectTransitionsProbzTop10[order(-incorrectTransitionsProbzTop10$probz),]
#------

#------
# look at the frequency for relavant 2-grams (given by Gigerenzer paper) for both groups

print(combineCorrectIncorrectWide2$BC)
print(combineCorrectIncorrectWide2$CB)

print(combineCorrectIncorrectWide2$CE)
print(combineCorrectIncorrectWide2$EC)

print(combineCorrectIncorrectWide2$DG)
print(combineCorrectIncorrectWide2$GD)

print(combineCorrectIncorrectWide2$EG)
print(combineCorrectIncorrectWide2$GE)

#------


# Look at transition distributions accross AOIs for Correct and Incorrect
# I want to see if Correct are more focused than Incorrect 
dfTransitionsLong <- gather(combineCorrectIncorrectWide2,Transitions, correct)
colnames(dfTransitionsLong)[3] <- "transition_prob"
#dfTransitionsLongABOVE <- dfTransitionsLong[dfTransitionsLong$transition_prob > 0,]
dfTransitionsLongABOVE <- dfTransitionsLong
CdfTransitionsLongABOVE <-dfTransitionsLongABOVE[dfTransitionsLongABOVE$correct== "correct",]
IdfTransitionsLongABOVE <- dfTransitionsLongABOVE[dfTransitionsLongABOVE$correct == "incorrect",] 

CdfTransitionsLongABOVE$Transitions <- factor(CdfTransitionsLongABOVE$Transitions, 
                                        levels = CdfTransitionsLongABOVE$Transitions[order(CdfTransitionsLongABOVE$transition_prob)])
IdfTransitionsLongABOVE$Transitions <- factor(IdfTransitionsLongABOVE$Transitions, 
                                              levels = IdfTransitionsLongABOVE$Transitions[order(IdfTransitionsLongABOVE$transition_prob)])

ggplot(data=CdfTransitionsLongABOVE, aes(x=Transitions, y=transition_prob)) +
  geom_bar(stat="identity")+
  ylim(0, 0.5) +
  theme(axis.text.x=element_blank()) 

ggplot(data=IdfTransitionsLongABOVE, aes(x=Transitions, y=transition_prob)) +
  geom_bar(stat="identity")+
  ylim(0, 0.5)+
  theme(axis.text.x=element_blank())

# facet_grid(. ~ correct)
# scale_fill_grey(start = 0.1, end = 0.5)

#geom_smooth(method = 'loess') 


ggplot(dfTransitionsLongABOVE, aes(transition_prob, fill = correct)) + 
  geom_density(alpha = 0.6) + labs(fill = "Correct") + scale_fill_grey()





# -------- look at the number of transitions (this is frequency of 2-grams, not transitions) ---------------------------

# number of transition for participants for the full dataframe 
# I want to know the number of transition per person for Correct 
Ntrans <- NULL
for (i in 1:nrow(uncorrectGroupTree)){
  Ntrans <- append( Ntrans, nchar(as.character(uncorrectGroupTree$Scanpath[i]))-1)
}
uncorrectGroupTree$transitions <- Ntrans

# I want to know the number of transition per person for Incorrect
Ntrans2 <- NULL
for (i in 1:nrow(unincorrectGroupTree)){
  Ntrans2 <- append( Ntrans2, nchar(as.character(unincorrectGroupTree$Scanpath[i]))-1)
}
unincorrectGroupTree$transitions <- Ntrans2

# number of transition for participants for the truncated dataframe 
# I want to know the number of TRUNCATED transitions per person for Correct 
Ntrans3 <- NULL
for (i in 1:nrow(uncorrectGroupTree)){
  Ntrans3 <- append( Ntrans3, nchar(as.character(uncorrectGroupTree$Scanpath2[i]))-1)
}
uncorrectGroupTree$transitions2 <- Ntrans3

# I want to know the number of TRUNCATED transitions per person for Incorrect
Ntrans4 <- NULL
for (i in 1:nrow(unincorrectGroupTree)){
  Ntrans4 <- append( Ntrans4, nchar(as.character(unincorrectGroupTree$Scanpath2[i]))-1)
}
unincorrectGroupTree$transitions2 <- Ntrans4

uncorrectGroupTree$correct <- rep("correct",nrow(uncorrectGroupTree))
unincorrectGroupTree$correct <- rep("incorrect",nrow(unincorrectGroupTree))
transitionCombinedDF <- rbind(uncorrectGroupTree,unincorrectGroupTree)


# to plot both COMPLETE transitions
ggplot(transitionCombinedDF, aes(correct, transitions))+geom_boxplot()
ggplot(transitionCombinedDF, aes(transitions, fill = correct)) +
  geom_density(alpha = 0.6) 
ggplot(transitionCombinedDF, aes(x = correct, y = transitions)) +
  geom_jitter(size = 10, alpha = 0.5, width = 0.25)  + 
  xlab("Time played in World Cup (minutes)")
# descriptive stats 
mean(uncorrectGroupTree$transitions)
sd(uncorrectGroupTree$transitions)
mean(unincorrectGroupTree$transitions)
sd(unincorrectGroupTree$transitions)

# to plot both TRUNCATED transitions
ggplot(transitionCombinedDF, aes(correct, transitions2))+geom_boxplot()
ggplot(transitionCombinedDF, aes(transitions2, fill = correct)) +
  geom_density(alpha = 0.6) 
ggplot(transitionCombinedDF, aes(x = correct, y = transitions2)) +
  geom_jitter(size = 10, alpha = 0.5, width = 0.25)  + 
  xlab("Time played in World Cup (minutes)")
# descriptive stats 
mean(uncorrectGroupTree$transitions2)
sd(uncorrectGroupTree$transitions2)
mean(unincorrectGroupTree$transitions2)
sd(unincorrectGroupTree$transitions2)



#----------------------- FIXATION COUNT for LONGER FIXATIONS ----------------------------------------------------------------
#
# Here we look at the fixation time for single sifxations
# we want to see the distributions of the fixation time for
# the two groups
#
#-------------------------------------------------------------------------------------------------------

# build the dataframe 

fixTimeDF <- read.csv(file="study2ET.csv",header=TRUE,sep=",")

# fixations length 

fixTREE<-subset(fixTimeDF,Format=="tree")
fixTREE<-fixTREE[,1:13]
#fixTREE[is.na(tree2)] <- 0
fixTREE$correct <- rep("correct", nrow(fixTREE))
fixTREE[fixTREE$Participant %in%  c("3","4","5","6","7","11","19","24","25","28","29","32","33","39","44","46","48","49","50"), 14] <- "incorrect"
#fixTREE <- fixTREE[-1]

ggplot(fixTREE, aes(FixationDuration, fill = correct)) + 
  geom_density(alpha = 0.6) + 
  #xlim(-2, 80) + 
  labs(fill = "Correct") + scale_fill_grey()

ggplot(fixTREE, aes(x=correct, y=FixationDuration)) + 
  geom_boxplot()+
  xlab("Group")
             

ggplot(fixTREE, aes(x = correct, y = FixationDuration)) +
  geom_jitter(size = 3, alpha = 0.4, width = 0.35)  + 
  xlab("Group")

# descriptive stats
median(fixTREE$FixationDuration[fixTREE$correct=="correct"])
quantile(fixTREE$FixationDuration[fixTREE$correct=="correct"])

median(fixTREE$FixationDuration[fixTREE$correct=="incorrect"])
quantile(fixTREE$FixationDuration[fixTREE$correct=="incorrect"])


# ---- Deeper analysis on where these long fixations are --------------------

# Create a dataframe for longer fixations over 217 millisecond
fixTREELong <- fixTREE[fixTREE$FixationDuration > 217, ]
# df with count of longer fixation
fixTREELongCount <- data.frame(AOI = names(fixTREELong)[5:13], Count = rep(0, length(names(fixTREELong)[5:13])))   


# Create a dataframe for longer fixations over 217 millisecond for correct 
fixTREELongCORRECT <- fixTREELong[fixTREELong$correct=="correct",]
# df with count of longer fixation for correct 
fixTREELongCountCORRECT <- data.frame(AOI = names(fixTREELongCORRECT)[5:13], Count = rep(0, length(names(fixTREELongCORRECT)[5:13])))  


# Create a dataframe for longer fixations over 217 millisecond for incorrect 
fixTREELongINCORRECT <- fixTREELong[fixTREELong$correct=="incorrect",]
# df with count of longer fixation for incorrect 
fixTREELongCountINCORRECT <- data.frame(AOI = names(fixTREELongCORRECT)[5:13], Count = rep(0, length(names(fixTREELongCORRECT)[5:13])))  



for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixTREELong, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:13]
    aoi <- names(aois)[aois == 1]

    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixTREELongCount[fixTREELongCount$AOI == aoi, 2] = fixTREELongCount[fixTREELongCount$AOI == aoi, 2] + 1
  }
}

for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixTREELongCORRECT, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:13]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixTREELongCountCORRECT[fixTREELongCountCORRECT$AOI == aoi, 2] = fixTREELongCountCORRECT[fixTREELongCountCORRECT$AOI == aoi, 2] + 1
  }
}

for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixTREELongINCORRECT, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:13]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixTREELongCountINCORRECT[fixTREELongCountINCORRECT$AOI == aoi, 2] = fixTREELongCountINCORRECT[fixTREELongCountINCORRECT$AOI == aoi, 2] + 1
  }
}


print(fixTREELongCount)
print(fixTREELongCountCORRECT)
print(fixTREELongCountINCORRECT)

fixTREELongCountCORRECT$Correct <- rep("correct", nrow(fixTREELongCountCORRECT))
fixTREELongCountCORRECT <- fixTREELongCountCORRECT[-1, ]
fixTREELongCountCORRECT$Count <- fixTREELongCountCORRECT$Count/sum(fixTREELongCountCORRECT$Count)

fixTREELongCountINCORRECT$Correct <- rep("incorrect", nrow(fixTREELongCountINCORRECT))
fixTREELongCountINCORRECT <- fixTREELongCountINCORRECT[-1, ]
fixTREELongCountINCORRECT$Count <- fixTREELongCountINCORRECT$Count/sum(fixTREELongCountINCORRECT$Count)


fixTREELongTOTAL <- rbind(fixTREELongCountCORRECT, fixTREELongCountINCORRECT)
fixTREELongTOTAL$AOI <- revalue(fixTREELongTOTAL$AOI, c("treetotal"="B", "treerainy"="C", "treesunny" = "D", "treerl"= "E", "treerh" = "F", "treesl"="G", "treesh" = "H", "treequestion" = "I"))

# Use position=position_dodge()
ggplot(data=fixTREELongTOTAL, aes(x=AOI, y=Count, fill=Correct)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_grey(start = 0.1, end = 0.5)





# -----------------------------------------------------------------------------------------------------------------
# Look at the fixation count for the whole dataframe (not just at long fixations)
# and produce a table for AOIs fixation frequency and a graph
# ------------------------------------------------------------------------------------------------------------------


# Find average participants fixations count by group 

FIXbyPART <- ddply(fixTimeDF,.(Participant, Format),nrow)
FIXbyPART$Correct <- rep("correct", nrow(FIXbyPART))
FIXbyPART[FIXbyPART$Participant %in%  c("3","4","5","6","7","11","19","24","25","28","29","32","33","39","44","46","48","49","50"), 4] <- "incorrect"

mean(FIXbyPART$V1[FIXbyPART$Format== "tree" & FIXbyPART$Correct == "correct"])
sd(FIXbyPART$V1[FIXbyPART$Format== "tree" & FIXbyPART$Correct == "correct"])
mean(FIXbyPART$V1[FIXbyPART$Format== "tree" & FIXbyPART$Correct == "incorrect"])
sd(FIXbyPART$V1[FIXbyPART$Format== "tree" & FIXbyPART$Correct == "incorrect"])

mean(FIXbyPART$V1[FIXbyPART$Format== "venn" & FIXbyPART$Correct == "correct"])
sd(FIXbyPART$V1[FIXbyPART$Format== "venn" & FIXbyPART$Correct == "correct"])
mean(FIXbyPART$V1[FIXbyPART$Format== "venn" & FIXbyPART$Correct == "incorrect"])
sd(FIXbyPART$V1[FIXbyPART$Format== "venn" & FIXbyPART$Correct == "incorrect"])


# look at fixations by AOIs (not by participant like before) 

#---- TREE
fixTREE_CORRECT <- fixTREE[fixTREE$correct=="correct",]
fixTREE_INCORRECT <- fixTREE[fixTREE$correct=="incorrect",]
fixTREECountCORRECT <- data.frame(AOI = names(fixTREE_CORRECT)[5:13], Count = rep(0, length(names(fixTREE_CORRECT)[5:13])))
fixTREECountINCORRECT <- data.frame(AOI = names(fixTREE_INCORRECT)[5:13], Count = rep(0, length(names(fixTREE_INCORRECT)[5:13])))

for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixTREE_CORRECT, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:13]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixTREECountCORRECT[fixTREECountCORRECT$AOI == aoi, 2] = fixTREECountCORRECT[fixTREECountCORRECT$AOI == aoi, 2] + 1
  }
}

for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixTREE_INCORRECT, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:13]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixTREECountINCORRECT[fixTREECountINCORRECT$AOI == aoi, 2] = fixTREECountINCORRECT[fixTREECountINCORRECT$AOI == aoi, 2] + 1
  }
}

# produce proportion count (probabilities)
fixTREECountCORRECT$Correct <- rep("correct", nrow(fixTREECountCORRECT))
fixTREECountCORRECT <- fixTREECountCORRECT[-1, ]
fixTREECountCORRECT$proportion <- fixTREECountCORRECT$Count/sum(fixTREECountCORRECT$Count)

fixTREECountINCORRECT$Correct <- rep("incorrect", nrow(fixTREECountINCORRECT))
fixTREECountINCORRECT <- fixTREECountINCORRECT[-1, ]
fixTREECountINCORRECT$proportion <- fixTREECountINCORRECT$Count/sum(fixTREECountINCORRECT$Count)

fixTREETOTAL <- rbind(fixTREECountCORRECT, fixTREECountINCORRECT)
fixTREETOTAL$AOI <- revalue(fixTREETOTAL$AOI, c("treetotal"="B", "treerainy"="C", "treesunny" = "D", 
                                                        "treerl"= "E", "treerh" = "F", "treesl"="G", 
                                                        "treesh" = "H", "treequestion" = "I"))

ggplot(data=fixTREETOTAL, aes(x=AOI, y=proportion, fill=Correct)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_grey(start = 0.1, end = 0.5)

#---- VENN
# I want to see the distribution for Venn as well for fixations count 
# So do the same for VENN

fixTimeDF <- read.csv(file="study2ET.csv",header=TRUE,sep=",")

  
fixVENN<-subset(fixTimeDF,Format=="venn")
fixVENN<-fixVENN[,-c(5:13)]
fixVENN$correct <- rep("correct", nrow(fixVENN))
listInorrectVenn <- NEWdfANALYSIS$Ps[NEWdfANALYSIS$Format=="venn" & NEWdfANALYSIS$Correct==0]
fixVENN[fixVENN$Participant %in%  listInorrectVenn, 11] <- "incorrect"



fixVENN_CORRECT <- fixVENN[fixVENN$correct=="correct",]
fixVENN_INCORRECT <- fixVENN[fixVENN$correct=="incorrect",]
fixVENNCountCORRECT <- data.frame(AOI = names(fixVENN_CORRECT)[5:10], Count = rep(0, length(names(fixVENN_CORRECT)[5:10])))
fixVENNCountINCORRECT <- data.frame(AOI = names(fixVENN_INCORRECT)[5:10], Count = rep(0, length(names(fixVENN_INCORRECT)[5:10])))

for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixVENN_CORRECT, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:10]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixVENNCountCORRECT[fixVENNCountCORRECT$AOI == aoi, 2] = fixVENNCountCORRECT[fixVENNCountCORRECT$AOI == aoi, 2] + 1
  }
}

for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixVENN_INCORRECT, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:10]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixVENNCountINCORRECT[fixVENNCountINCORRECT$AOI == aoi, 2] = fixVENNCountINCORRECT[fixVENNCountINCORRECT$AOI == aoi, 2] + 1
  }
}

# produce proportion count (probabilities)
fixVENNCountCORRECT$Correct <- rep("correct", nrow(fixVENNCountCORRECT))
fixVENNCountCORRECT <- fixVENNCountCORRECT[-1, ]
fixVENNCountCORRECT$proportion <- fixVENNCountCORRECT$Count/sum(fixVENNCountCORRECT$Count)

fixVENNCountINCORRECT$Correct <- rep("incorrect", nrow(fixVENNCountINCORRECT))
fixVENNCountINCORRECT <- fixVENNCountINCORRECT[-1, ]
fixVENNCountINCORRECT$proportion <- fixVENNCountINCORRECT$Count/sum(fixVENNCountINCORRECT$Count)

fixVENNTOTAL <- rbind(fixVENNCountCORRECT, fixVENNCountINCORRECT)
fixVENNTOTAL$AOI <- revalue(fixVENNTOTAL$AOI, c("vennrh"="M", "vennrl"="N", "vennsl" = "O", 
                                                "venntotal"= "P", "vennquestion" = "Q"))

ggplot(data=fixVENNTOTAL, aes(x=AOI, y=proportion, fill=Correct)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_grey(start = 0.1, end = 0.5)




# -----------------------------------------------------------------------------------------------------------------
# ADA TEST FOR FIXATION COUNT PROPORTION (not for transitions) without background 
# 
# ------------------------------------------------------------------------------------------------------------------

#distribution for fixation proportion in tree correct Vs incorrect 

distrTREE_correct <-  fixTREETOTAL$proportion[fixTREETOTAL$Correct=="correct"]

distrTREE_incorrect <-  fixTREETOTAL$proportion[fixTREETOTAL$Correct=="incorrect"]


#distribution for fixation proportion in venn correct Vs incorrect 

distrVENN_correct <-  fixVENNTOTAL$proportion[fixVENNTOTAL$Correct=="correct"]

distrVENN_incorrect <-  fixVENNTOTAL$proportion[fixVENNTOTAL$Correct=="incorrect"]



# examples
ADA(distrTREE_correct)
ADA(distrTREE_incorrect)

ADA(distrVENN_correct)
ADA(distrVENN_incorrect)






# -----------------------------------------------------------------------------------------------------------------
# ADA TEST FOR 2-gram relative FREQENCIES COUNT (not for FIXATION COUNT PROPORTION or transitions) without background 
# 
# ------------------------------------------------------------------------------------------------------------------

correctVectorTREE <- as.numeric(as.vector(combineCorrectIncorrectWide2[1,2:57]))

incorrectVectorTREE <- as.numeric(as.vector(combineCorrectIncorrectWide2[2,2:57]))


ADA(correctVectorTREE)
ADA(incorrectVectorTREE)











# --------------------------- UNNECESSARY START -------------------------------------------
# ----------------------------------------------------------------------------------------------------------------- # 
#                           LOOK AT BETWEEN FREQUENCIES (PROBABILITIES)                                             #
# ----------------------------------------------------------------------------------------------------------------- #
# further analysis for optimal scanpaths, by looking at frequencies between and 
# transforming them into probabilities 


probVectorC <- combineCorrectIncorrect$freqz[combineCorrectIncorrect$correct=="correct"]/
  sum(combineCorrectIncorrect$freqz[combineCorrectIncorrect$correct=="correct"])
probVectorI <- combineCorrectIncorrect$freqz[combineCorrectIncorrect$correct=="incorrect"]/
  sum(combineCorrectIncorrect$freqz[combineCorrectIncorrect$correct=="incorrect"])

combineCorrectIncorrect$probablity <- rep(0, nrow(combineCorrectIncorrect))


combineCorrectIncorrect$probablity[combineCorrectIncorrect$correct=="correct"] <- probVectorC
combineCorrectIncorrect$probablity[combineCorrectIncorrect$correct=="incorrect"] <- probVectorI

DFprobC <- subset(combineCorrectIncorrect, combineCorrectIncorrect$correct=="correct")
DFprobI <- subset(combineCorrectIncorrect, combineCorrectIncorrect$correct=="incorrect")

TOP5DFprobC <- DFprobC[DFprobC$probablity %in% tail(sort(DFprobC$probablity),5), ]
TOP5DFprobI <- DFprobI[DFprobI$probablity %in% tail(sort(DFprobI$probablity),5), ]

TOP5DFprobC <- TOP5DFprobC[order(-TOP5DFprobC$probablity),] 
TOP5DFprobI <- TOP5DFprobI[order(-TOP5DFprobI$probablity),]

diffFrequencyDF <- data.frame(bigram = NULL, difference = NULL, correct = NULL, incorrect = NULL)
diff <- 0
vecDiff <- NULL
for (i in 1:56){
  bgram <- combineCorrectIncorrect$n_gram[i]
  difff <- abs(combineCorrectIncorrect$probablity[i] - combineCorrectIncorrect$probablity[i+56])
  newRow <- data.frame(bigram = bgram, difference = difff, 
                       correct = combineCorrectIncorrect$probablity[i], incorrect = combineCorrectIncorrect$probablity[i+56])
  diffFrequencyDF <- rbind(diffFrequencyDF,newRow)
}

diffFrequencyDFtOP5 <- subset(diffFrequencyDF, difference %in% tail(sort(diffFrequencyDF$difference),5))
diffFrequencyDFtOP5 <- diffFrequencyDFtOP5[order(-diffFrequencyDFtOP5$difference),]
# --------------------------- UNNECESSARY END -------------------------------------------









# ----------------------------------------------------------------------------------------------------------------- # 
#                                                                                                                   #
#          PERMUTATION TEST USING HELLIGER DISTANCE OVER 2-GRAM relative FREQEUNCY                                  #
#                                                                                                                   #
# ----------------------------------------------------------------------------------------------------------------- #




#---------------------------------------------------------------------------------
# FUNCTION:     permutationTestVector()
#
# INPUT parameters: 
# DF <- subset(NEWdfANALYSIS, Format == "tree")
# LISTngrams <- REDUCEDtotalBiGrams__tree
# group1_size <- nrow(subset(NEWdfANALYSIS, Correct==1 & Format == "tree"))
# distanceMeasure <- vecHellDistance
# ncolumDATA <- 11 (where the scanpath is, the numbe rof column in the DF)
# perms <- numberPermutations (10000)
#
# OUTPUT:
# A vector of distances 
#---------------------------------------------------------------------------------

permutationTestVector <- function(DF, LISTngrams, group1_size, distanceMeasure, ncolumDATA, perms = 10){
  distance <- NULL
  for (i in 1:perms){
    GROUP1<- NULL
    GROUP2<- NULL
    freqGROUP1<- NULL
    freqGROUP2<- NULL
    distance_results<-0
    
    # create two DF groups sampling at random 
    GROUP1 <- DF[sample(1:nrow(DF), group1_size, replace=FALSE),]
    GROUP2 <- DF[ !(DF$Ps %in% GROUP1$Ps), ]
    # create a frequency dataframe (frequency of ngrams)
    freqGROUP1 <- subSeqfreq(LISTngrams, GROUP1, ncolumDATA)
    freqGROUP2 <- subSeqfreq(LISTngrams, GROUP2, ncolumDATA)
    # transform into probabilities 
    freqGROUP1$probz <- freqGROUP1$freqz/ sum(freqGROUP1$freqz)
    freqGROUP2$probz <- freqGROUP2$freqz/ sum(freqGROUP2$freqz)
    vettore1 <- as.numeric(as.vector(freqGROUP1$probz))
    vettore2 <- as.numeric(as.vector(freqGROUP2$probz))
    # caluclate distance
    distance_results <- distanceMeasure(vettore1 , vettore2)
    # add distance to fdistance vector 
    distance <- c(distance, distance_results)
  }
  return(distance) 
}

#---------------------------------------------------------------------------------
# FUNCTION:     vecHellDistance()
# INPUT:        matrix, matrix
# OUTPUT:       Summed Hellinger Distance
# DESCRIPTION:  Hellinger Distance is defined as: 1/sqrt(2) * sqrt(sum(suqare(sqrt(pi - squrt(qi))))
#---------------------------------------------------------------------------------
vecHellDistance <- function(v1,v2){
  length_of_vector <- length(v1)
  HJ <- 0
  for (i in 1:length_of_vector){
    HJ <- HJ + ((sqrt(v1[i]) - sqrt(v2[i]))^2)
  }
  return((1/sqrt(2))*sqrt(HJ)) 
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


numberPermutations <- 10
# ---------------- RUN THE PERMUTATION TEST --------------------------------------------- #

# Find the distance between correct and incorrect
distCorrectIncorrectT <- vecHellDistance(correctVectorTREE, incorrectVectorTREE)
HellTREEylabel <- "Hellinger distance"
# Run permutation test 
freqNGRAMtreeDISTANCE <- permutationTestVector(subset(NEWdfANALYSIS, Format == "tree"),
                      REDUCEDtotalBiGrams__tree,
                      nrow(subset(NEWdfANALYSIS, Correct==1 & Format == "tree")),
                      vecHellDistance,
                      11,
                      numberPermutations)




# generate dot plot
freqNGRAMtreeDISTANCE <- append(freqNGRAMtreeDISTANCE, distCorrectIncorrectT)
plot(freqNGRAMtreeDISTANCE, col = ifelse(freqNGRAMtreeDISTANCE == distCorrectIncorrectT, 'blue', 'red'), 
     main = paste0("Tree Diagram ", " (",HellTREEylabel,")"), xlab = "Number of permutations", ylab = HellTREEylabel)
abline(h = distCorrectIncorrectT, col = "purple")

# generate density plot
HellTREEdensity_plot <- density(freqNGRAMtreeDISTANCE)
plot(HellTREEdensity_plot, type = "n", main = paste0("Tree Diagram ", " (",HellTREEylabel,")"), xlab = "Hellinger distance")
polygon(HellTREEdensity_plot, col = "lightgray", border = "grey")
rug(freqNGRAMtreeDISTANCE, col = ifelse(freqNGRAMtreeDISTANCE == distCorrectIncorrectT, 'blue', 'red'))
abline(v = distCorrectIncorrectT, col = "purple")

HellpvalueTREE <- calculatePvalue(distCorrectIncorrectT, freqNGRAMtreeDISTANCE)
HellResultsTREE <- list("The Hellinger distance is:", distCorrectIncorrectT, "The p-value is:", HellpvalueTREE)
print(HellResultsTREE)






# ------ VENN

# create DFs
freqDFvennCorrect <- subSeqfreq(REDUCEDtotalBiGrams__venn, subset(NEWdfANALYSIS,Correct==1 & Format == "venn"), 11)
freqDFvennIncorrect <- subSeqfreq(REDUCEDtotalBiGrams__venn, subset(NEWdfANALYSIS,Correct==0 & Format == "venn"), 11)
# transform into probabilities 
freqDFvennCorrect$probz <- freqDFvennCorrect$freqz/ sum(freqDFvennCorrect$freqz)
freqDFvennIncorrect$probz <- freqDFvennIncorrect$freqz/ sum(freqDFvennIncorrect$freqz)
correctVectorVENN <- as.numeric(as.vector(freqDFvennCorrect$probz))
incorrectVectorVENN <- as.numeric(as.vector(freqDFvennIncorrect$probz))

# Find the distance between correct and incorrect
distCorrectIncorrectV <- vecHellDistance(correctVectorVENN, incorrectVectorVENN)
HellVENNylabel <- "Hellinger distance"
# Run permutation test 
freqNGRAMvennDISTANCE <- permutationTestVector(subset(NEWdfANALYSIS, Format == "venn"),
                                               REDUCEDtotalBiGrams__venn,
                                               nrow(subset(NEWdfANALYSIS, Correct==1 & Format == "venn")),
                                               vecHellDistance,
                                               11,
                                               numberPermutations)



# generate dot plot
freqNGRAMvennDISTANCE <- append(freqNGRAMvennDISTANCE, distCorrectIncorrectV)
plot(freqNGRAMvennDISTANCE, col = ifelse(freqNGRAMvennDISTANCE == distCorrectIncorrectV, 'blue', 'red'), 
     main = paste0("Venn Diagram ", " (",HellVENNylabel,")"), xlab = "Number of permutations", ylab = HellVENNylabel)
abline(h = distCorrectIncorrectV, col = "purple")

# generate density plot
HellVENNdensity_plot <- density(freqNGRAMvennDISTANCE)
plot(HellVENNdensity_plot, type = "n", main = paste0("Venn Diagram ", " (",HellTREEylabel,")"), xlab = "Hellinger distance")
polygon(HellVENNdensity_plot, col = "lightgray", border = "grey")
rug(freqNGRAMvennDISTANCE, col = ifelse(freqNGRAMvennDISTANCE == distCorrectIncorrectV, 'blue', 'red'))
abline(v = distCorrectIncorrectV, col = "purple")

HellpvalueVENN <- calculatePvalue(distCorrectIncorrectV, freqNGRAMvennDISTANCE)
HellResultsVENN <- list("The Hellinger distance is:", distCorrectIncorrectV, "The p-value is:", HellpvalueVENN)
print(HellResultsVENN)






# NEW ANALYSIS 28/09/2017 (AFTER NIELS' FEEDBACK)
#------------------------------------------------------------------------------------------------
#
#                   Classical statistical methods 
#
#------------------------------------------------------------------------------------------------


# create a matrix for tree
TREEfisher <- as.matrix(rbind(combineCorrectIncorrect$freqz[combineCorrectIncorrect$correct=="correct"], 
                              combineCorrectIncorrect$freqz[combineCorrectIncorrect$correct=="incorrect"]))

dimnames(TREEfisher) <- list(Correctness = c("correct", "incorrect"),
                             subsequence = c(as.vector(combineCorrectIncorrect$n_gram[combineCorrectIncorrect$correct=="correct"])))

prova<- as.data.frame(TREEfisher)

# create a matrix for venn
VENNfisher <- as.matrix(rbind(freqDFvennCorrect$freqz, 
                              freqDFvennIncorrect$freqz))

dimnames(VENNfisher) <- list(Correctness = c("correct", "incorrect"),
                             subsequence = c(as.vector(freqDFvennCorrect$n_gram)))

prova2<- as.data.frame(VENNfisher)



# ------------------------------------------------------------------------
#                               CHI-squared Test
# -------------------------------------------------------------------------
# Assumptions for Chi-square:
# the expected value in each cell is at least 5

#TREE
#chisq.test(TREEfisher)
chi1 <- chisq.test(TREEfisher, simulate.p.value = TRUE, B = 10000)
chi1$p.value

#VENN
#chisq.test(VENNfisher)
chi2 <- chisq.test(VENNfisher, simulate.p.value = TRUE, B = 10000)
chi2$p.value

# Expected values 
min(chi1$expected)
min(chi2$expected)



#-----------------------------------------------------------------
#                     FISHER's Exact Test
# -----------------------------------------------------------------
# Assumptions for Fisher's exact test:
# The marginals need to be fixed (i.e. we know in advance the totals)
# In our case the marginal are not fixed, i.e. if we re-run the exact same experiment we would 
# not obtain the same number of fixations (by row or by column)

# TREE
#fisher.test(TREEfisher)
fisher.test(TREEfisher,simulate.p.value=TRUE,B=1e5)

# VENN
#fisher.test(VENNfisher)
fisher.test(VENNfisher,simulate.p.value=TRUE,B=1e5)



#---------------------------------------------------------------------------------------------------------
#
#                       2-gram analysis for Venn (find the top 2-grams)
#
# --------------------------------------------------------------------------------------------------------

# ---- Create and combine the dataframe for VENN
freqDFvennCorrect$correct <- rep("correct",nrow(freqDFvennCorrect))
freqDFvennIncorrect$correct <- rep("incorrect",nrow(freqDFvennIncorrect))
combineCorrectIncorrectVENN <- merge(freqDFvennCorrect,freqDFvennIncorrect, by="n_gram") #rbind(freqDFvennCorrect,freqDFvennIncorrect)

# calculate abs differences
combineCorrectIncorrectVENN$diff <- abs(combineCorrectIncorrectVENN$probz.x - combineCorrectIncorrectVENN$probz.y)

# calculate odds scale, i.e. (p/(1-p)) / (q/(1-q)) 
combineCorrectIncorrectVENN$odds <- (combineCorrectIncorrectVENN$probz.x/(1-combineCorrectIncorrectVENN$probz.x)) / 
  (combineCorrectIncorrectVENN$probz.y/(1-combineCorrectIncorrectVENN$probz.y))

# find top 2-grams for frequency
top5VENNcorrect<- freqDFvennCorrect[freqDFvennCorrect$probz %in% tail(sort(freqDFvennCorrect$probz),5),]
top5VENNcorrect <- top5VENNcorrect[order(-top5VENNcorrect$probz),]
top5VENNincorrect <- freqDFvennIncorrect[freqDFvennIncorrect$probz %in% tail(sort(freqDFvennIncorrect$probz),5),]
top5VENNincorrect <- top5VENNincorrect[order(-top5VENNincorrect$probz),]
print(top5VENNcorrect)
print(top5VENNincorrect)

# find top 2-grams for differences
top5VENNdiff <- combineCorrectIncorrectVENN[combineCorrectIncorrectVENN$diff %in% tail(sort(combineCorrectIncorrectVENN$diff),5),]
top5VENNdiff<- top5VENNdiff[order(-top5VENNdiff$diff),]

# find top/bottom 2-grams for odd ratios   
top2VENNodds <- combineCorrectIncorrectVENN[combineCorrectIncorrectVENN$odds %in% tail(sort(combineCorrectIncorrectVENN$odds),2),]
top2VENNodds <- top2VENNodds[order(-top2VENNodds$odds),]
bottom2VENNodds <- combineCorrectIncorrectVENN[combineCorrectIncorrectVENN$odds %in% head(sort(combineCorrectIncorrectVENN$odds),2),] 
bottom2VENNodds <- bottom2VENNodds[order(-bottom2VENNodds$odds),]

# PLOTTING ABSOLUTE DIFFERENCES (WITH LAPLACE TRANSFORMATION)
# reorder the label in crescendo
top5VENNdiff$n_gram <- factor(top5VENNdiff$n_gram, levels = top5VENNdiff$n_gram[order(top5VENNdiff$diff)]) # dfTop5
ggplot(top5VENNdiff, aes(x = n_gram, y = diff))+ geom_bar(stat="identity", width=0.8) + coord_flip() +
  xlab("2-gram") + ylab("Δ|Correct - Incorrect|")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

# create a df for plotting 
newDFforPlot1 <- top5VENNdiff[,1:4]
colnames(newDFforPlot1)<- c("n_gram","freqz","probz","correct")
newDFforPlot2 <- top5VENNdiff[,c(1,5:7)]
colnames(newDFforPlot2)<- c("n_gram","freqz","probz","correct")
newDFforPlot3 <- rbind(newDFforPlot1, newDFforPlot2)

# reorder the label in crescendo
newDFforPlot3$n_gram <- factor(newDFforPlot3$n_gram, levels = newDFforPlot3$n_gram[order(-top5VENNdiff$diff)])
#plot the data
ggplot(newDFforPlot3, aes(x = n_gram, y = probz, fill = correct)) +
  geom_bar(stat="identity",colour="black",position = "dodge", width=0.8) +
  xlab("2-gram") + ylab("Prob") +
  scale_fill_grey(start = 0.2, end = .5)+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

# PLOTTING ODD RATIOS (WITH LAPLACE TRANSFORMATION)
# reorder the label in crescendo
oddsDF <- rbind(top2VENNodds,bottom2VENNodds)
oddsDF$n_gram <- factor(oddsDF$n_gram, levels = oddsDF$n_gram[order(oddsDF$odds)]) 
ggplot(oddsDF, aes(x = n_gram, y = odds))+ geom_bar(stat="identity", width=0.8) + coord_flip() +
  xlab("2-gram") + ylab("Odd-ratios")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

# create a df for plotting 
newDFforPlot1 <- oddsDF[,1:4]
colnames(newDFforPlot1)<- c("n_gram","freqz","probz","correct")
newDFforPlot2 <- oddsDF[,c(1,5:7)]
colnames(newDFforPlot2)<- c("n_gram","freqz","probz","correct")
newDFforPlot3 <- rbind(newDFforPlot1, newDFforPlot2)

# reorder the label in crescendo
newDFforPlot3$n_gram <- factor(newDFforPlot3$n_gram, levels = newDFforPlot3$n_gram[order(-oddsDF$odds)])
#plot the data
ggplot(newDFforPlot3, aes(x = n_gram, y = probz, fill = correct)) +
  geom_bar(stat="identity",colour="black",position = "dodge", width=0.8) +
  xlab("2-gram") + ylab("Prob") +
  scale_fill_grey(start = 0.2, end = .5)+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

#---------------------------------------------------------------------------------------------------------
#
#                       2-gram analysis for TREE (find the top 2-grams)
#
# --------------------------------------------------------------------------------------------------------

# apply "add-one" smoothing (Laplace smoothing) to the frequencies before transforming into probabilities
freqSubsq_tree_correct$probz <- (freqSubsq_tree_correct$freqz+1)/ (sum(freqSubsq_tree_correct$freqz)+56)
freqSubsq_tree_incorrect$probz <- (freqSubsq_tree_incorrect$freqz+1)/ (sum(freqSubsq_tree_incorrect$freqz)+56)
# caluclate probabilitis with no Laplace smoothing 
#freqSubsq_tree_correct$probz <- (freqSubsq_tree_correct$freqz)/ (sum(freqSubsq_tree_correct$freqz))        #-----------NB
#freqSubsq_tree_incorrect$probz <- (freqSubsq_tree_incorrect$freqz)/ (sum(freqSubsq_tree_incorrect$freqz))  #-----------NB

# ---- Create and combine the dataframe for VENN
freqSubsq_tree_correct$correct <- rep("correct",nrow(freqSubsq_tree_correct))
freqSubsq_tree_incorrect$correct <- rep("incorrect",nrow(freqSubsq_tree_incorrect))
combineCorrectIncorrectTREE <- merge(freqSubsq_tree_correct,freqSubsq_tree_incorrect, by="n_gram")

# calculate abs differences
combineCorrectIncorrectTREE$diff <- abs(combineCorrectIncorrectTREE$probz.x - combineCorrectIncorrectTREE$probz.y)

# calculate odds scale, i.e. (p/(1-p)) / (q/(1-q)) 
combineCorrectIncorrectTREE$odds <- (combineCorrectIncorrectTREE$probz.x/(1-combineCorrectIncorrectTREE$probz.x)) / 
  (combineCorrectIncorrectTREE$probz.y/(1-combineCorrectIncorrectTREE$probz.y))

# find top 2-grams for frequency
top5TREEcorrect<- freqSubsq_tree_correct[freqSubsq_tree_correct$probz %in% tail(sort(freqSubsq_tree_correct$probz),5),]
top5TREEcorrect <- top5TREEcorrect[order(-top5TREEcorrect$probz),]
top5TREEincorrect <- freqSubsq_tree_incorrect[freqSubsq_tree_incorrect$probz %in% tail(sort(freqSubsq_tree_incorrect$probz),5),]
top5TREEincorrect <- top5TREEincorrect[order(-top5TREEincorrect$probz),]
print(top5TREEcorrect)
print(top5TREEincorrect)

# find top 2-grams for differences
top5TREEdiff <- combineCorrectIncorrectTREE[combineCorrectIncorrectTREE$diff %in% tail(sort(combineCorrectIncorrectTREE$diff),5),]
top5TREEdiff<- top5TREEdiff[order(-top5TREEdiff$diff),]

# find top/bottom 2-grams for odd ratios 
top2TREEodds <- combineCorrectIncorrectTREE[combineCorrectIncorrectTREE$odds %in% tail(sort(combineCorrectIncorrectTREE$odds),2),]
top2TREEodds <- top2TREEodds[order(-top2TREEodds$odds),]
bottom2TREEodds <- combineCorrectIncorrectTREE[combineCorrectIncorrectTREE$odds %in% head(sort(combineCorrectIncorrectTREE$odds),2),] 
bottom2TREEodds <- bottom2TREEodds[order(-bottom2TREEodds$odds),]

# PLOTTING ABSOLUTE DIFFERENCES (WITH LAPLACE TRANSFORMATION)
# reorder the label in crescendo
top5TREEdiff$n_gram <- factor(top5TREEdiff$n_gram, levels = top5TREEdiff$n_gram[order(top5TREEdiff$diff)]) # dfTop5
ggplot(top5TREEdiff, aes(x = n_gram, y = diff))+ geom_bar(stat="identity", width=0.8) + coord_flip() +
  xlab("2-gram") + ylab("Δ|Correct - Incorrect|")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

# create a df for plotting 
newDFforPlot1 <- top5TREEdiff[,1:4]
colnames(newDFforPlot1)<- c("n_gram","freqz","correct","probz")
newDFforPlot2 <- top5TREEdiff[,c(1,5:7)]
colnames(newDFforPlot2)<- c("n_gram","freqz","correct","probz")
newDFforPlot3 <- rbind(newDFforPlot1, newDFforPlot2)

# reorder the label in crescendo
newDFforPlot3$n_gram <- factor(newDFforPlot3$n_gram, levels = newDFforPlot3$n_gram[order(-top5TREEdiff$diff)])
#plot the data
ggplot(newDFforPlot3, aes(x = n_gram, y = probz, fill = correct)) +
  geom_bar(stat="identity",colour="black",position = "dodge", width=0.8) +
  xlab("2-gram") + ylab("Prob") +
  scale_fill_grey(start = 0.2, end = .5)+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

# PLOTTING ODD RATIOS (WITH LAPLACE TRANSFORMATION)
# reorder the label in crescendo
oddsDF <- rbind(top2TREEodds,bottom2TREEodds)
oddsDF$n_gram <- factor(oddsDF$n_gram, levels = oddsDF$n_gram[order(oddsDF$odds)]) 
ggplot(oddsDF, aes(x = n_gram, y = odds))+ geom_bar(stat="identity", width=0.8) + coord_flip() +
  xlab("2-gram") + ylab("Odd-ratios")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

# create a df for plotting 
newDFforPlot1 <- oddsDF[,1:4]
colnames(newDFforPlot1)<- c("n_gram","freqz","correct","probz")
newDFforPlot2 <- oddsDF[,c(1,5:7)]
colnames(newDFforPlot2)<- c("n_gram","freqz","correct","probz")
newDFforPlot3 <- rbind(newDFforPlot1, newDFforPlot2)

# reorder the label in crescendo
newDFforPlot3$n_gram <- factor(newDFforPlot3$n_gram, levels = newDFforPlot3$n_gram[order(-oddsDF$odds)])
#plot the data
ggplot(newDFforPlot3, aes(x = n_gram, y = probz, fill = correct)) +
  geom_bar(stat="identity",colour="black",position = "dodge", width=0.8) +
  xlab("2-gram") + ylab("Prob") +
  scale_fill_grey(start = 0.2, end = .5)+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

# produce dataframe for the two odds in order 
VENNbyODDS <- combineCorrectIncorrectVENN[order(combineCorrectIncorrectVENN$odds),]
TREEbyODDS <- combineCorrectIncorrectTREE[order(combineCorrectIncorrectTREE$odds),]



#---------------------------------------------------------------------------------------------------------
#
#                       Long Fixations analysis for VENN
#
# --------------------------------------------------------------------------------------------------------



#----------------------- FIXATION COUNT for LONGER FIXATIONS ----------------------------------------------------------------
#
# Here we look at the fixation time for single sifxations
# we want to see the distributions of the fixation time for
# the two groups
#
#-------------------------------------------------------------------------------------------------------

ggplot(fixVENN, aes(x=correct, y=FixationDuration)) + 
  geom_boxplot()+
  xlab("Group")


ggplot(fixVENN, aes(x = correct, y = FixationDuration)) +
  geom_jitter(size = 3, alpha = 0.4, width = 0.35)  + 
  xlab("Group")

# descriptive stats
median(fixVENN$FixationDuration[fixVENN$correct=="correct"])
quantile(fixVENN$FixationDuration[fixVENN$correct=="correct"])

median(fixVENN$FixationDuration[fixVENN$correct=="incorrect"])
quantile(fixVENN$FixationDuration[fixVENN$correct=="incorrect"])


# ---- Deeper analysis on where these long fixations are --------------------

# Create a dataframe for longer fixations over 217 millisecond
fixVENNLong <- fixVENN[fixVENN$FixationDuration > 217, ]
# df with count of longer fixation
fixVENNLongCount <- data.frame(AOI = names(fixVENNLong)[5:10], Count = rep(0, length(names(fixTREELong)[5:10])))   


# Create a dataframe for longer fixations over 217 millisecond for correct 
fixVENNLongCORRECT <- fixVENNLong[fixVENNLong$correct=="correct",]
# df with count of longer fixation for correct 
fixVENNLongCountCORRECT <- data.frame(AOI = names(fixVENNLongCORRECT)[5:10], Count = rep(0, length(names(fixVENNLongCORRECT)[5:10])))  


# Create a dataframe for longer fixations over 217 millisecond for incorrect 
fixVENNLongINCORRECT <- fixVENNLong[fixVENNLong$correct=="incorrect",]
# df with count of longer fixation for incorrect 
fixVENNLongCountINCORRECT <- data.frame(AOI = names(fixVENNLongCORRECT)[5:10], Count = rep(0, length(names(fixVENNLongCORRECT)[5:10])))  



for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixVENNLong, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:10]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixVENNLongCount[fixVENNLongCount$AOI == aoi, 2] = fixVENNLongCount[fixVENNLongCount$AOI == aoi, 2] + 1
  }
}

for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixVENNLongCORRECT, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:10]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixVENNLongCountCORRECT[fixVENNLongCountCORRECT$AOI == aoi, 2] = fixVENNLongCountCORRECT[fixVENNLongCountCORRECT$AOI == aoi, 2] + 1
  }
}

for (p in 2:50){
  #participant <- p 
  # subset DF for participant
  fixa <- subset(fixVENNLongINCORRECT, Participant == p) 
  backgroundAOI <- names(fixa)[5]
  #scanpath <- NULL
  for (fixation in fixa$FixationIndex){
    aois <- fixa[fixa$FixationIndex == fixation, 6:10]
    aoi <- names(aois)[aois == 1]
    
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    fixVENNLongCountINCORRECT[fixVENNLongCountINCORRECT$AOI == aoi, 2] = fixVENNLongCountINCORRECT[fixVENNLongCountINCORRECT$AOI == aoi, 2] + 1
  }
}


print(fixVENNLongCount)
print(fixVENNLongCountCORRECT)
print(fixVENNLongCountINCORRECT)

fixVENNLongCountCORRECT$Correct <- rep("correct", nrow(fixVENNLongCountCORRECT))
fixVENNLongCountCORRECT <- fixVENNLongCountCORRECT[-1, ]
fixVENNLongCountCORRECT$Count <- fixVENNLongCountCORRECT$Count/sum(fixVENNLongCountCORRECT$Count)

fixVENNLongCountINCORRECT$Correct <- rep("incorrect", nrow(fixVENNLongCountINCORRECT))
fixVENNLongCountINCORRECT <- fixVENNLongCountINCORRECT[-1, ]
fixVENNLongCountINCORRECT$Count <- fixVENNLongCountINCORRECT$Count/sum(fixVENNLongCountINCORRECT$Count)


fixVENNLongTOTAL <- rbind(fixVENNLongCountCORRECT, fixVENNLongCountINCORRECT)
fixVENNLongTOTAL$AOI <- revalue(fixVENNLongTOTAL$AOI, c("vennrh"="M", "vennrl"="N", "vennsl" = "O", "venntotal"= "P", "vennquestion" = "Q"))

# Use position=position_dodge()
ggplot(data=fixVENNLongTOTAL, aes(x=AOI, y=Count, fill=Correct)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_grey(start = 0.1, end = 0.5)














# ----------------------------------------------------------------------------------
# ANALYSIS for AVERAGE FIXATION COUNT FOR PARTICIPANT by AOIs by Group
#--------------------------------------------------------------------------------------
fixationsDFtree <- data.frame(Participant= rep(0, 49),B=rep(0, 49),C=rep(0, 49),D=rep(0, 49),E=rep(0, 49),F=rep(0, 49),G=rep(0, 49),H=rep(0, 49),I=rep(0, 49))  

for (p in 2:50){
  vector <- c(p)
  fixa <- subset(fixTREE, Participant == p) 
  #print(sum(fixa[6]))
  for (c in 6:13){
    rowFix <- (sum(fixa[c]))
    vector <- c(vector, rowFix)
  }
  z <- p-1
  fixationsDFtree[z,]<- vector
}


fixationsDFtree$correct <- rep("correct", nrow(fixationsDFtree))
fixationsDFtree[fixationsDFtree$Participant %in%  c("3","4","5","6","7","11","19","24","25","28","29","32","33","39","44","46","48","49","50"), 10] <- "incorrect"

#for long fixations
fixationsDFtreeLong <- data.frame(Participant= rep(0, 49),B=rep(0, 49),C=rep(0, 49),D=rep(0, 49),E=rep(0, 49),F=rep(0, 49),G=rep(0, 49),H=rep(0, 49),I=rep(0, 49))  

for (p in 2:50){
  vector <- c(p)
  if (p %in% fixTREELong$Participant) {
    fixa <- subset(fixTREELong, Participant == p) 
    for (c in 6:13){
      rowFix <- (sum(fixa[c]))
      vector <- c(vector, rowFix)
    }
    z <- p-1
    fixationsDFtreeLong[z,]<- vector
  }
}

fixationsDFtreeLong$correct <- rep("correct", nrow(fixationsDFtreeLong))
fixationsDFtreeLong[fixationsDFtreeLong$Participant %in%  c("3","4","5","6","7","11","19","24","25","28","29","32","33","39","44","46","48","49","50"), 10] <- "incorrect"
fixationsDFtreeLong[21,1]<-22
#------------------------------------------------------

# VENN (DF = fixVENN)
fixationsDFvenn <- data.frame(Participant= rep(0, 49),M=rep(0, 49),N=rep(0, 49),O=rep(0, 49),P=rep(0, 49),Q=rep(0, 49))  

for (p in 2:50){
  vector <- c(p)
  fixa <- subset(fixVENN, Participant == p) 
  for (c in 6:10){
    rowFix <- (sum(fixa[c]))
    vector <- c(vector, rowFix)
  }
  z <- p-1
  fixationsDFvenn[z,]<- vector
}


fixationsDFvenn$correct <- rep("correct", nrow(fixationsDFvenn))
fixationsDFvenn[fixationsDFvenn$Participant %in%  listInorrectVenn, 7] <- "incorrect"

#for long fixations (DF = fixVENNLong)
fixationsDFVennLong <- data.frame(Participant= rep(0, 49),M=rep(0, 49),N=rep(0, 49),O=rep(0, 49),P=rep(0, 49),Q=rep(0, 49)) 

for (p in 2:50){
  vector <- c(p)
  if (p %in% fixVENNLong$Participant) {
    fixa <- subset(fixVENNLong, Participant == p) 
    for (c in 6:10){
      rowFix <- (sum(fixa[c]))
      vector <- c(vector, rowFix)
    }
    z <- p-1
    fixationsDFVennLong[z,]<- vector
  }
}

fixationsDFVennLong$correct <- rep("correct", nrow(fixationsDFVennLong))
fixationsDFVennLong[fixationsDFVennLong$Participant %in%  listInorrectVenn, 7] <- "incorrect"






#fixations means
fixationsDFtreeCORR <- fixationsDFtree[fixationsDFtree$correct=="correct",]
for (p in 2:9){
  print(mean(fixationsDFtreeCORR[,p]))
}
fixationsDFtreeINCORR <- fixationsDFtree[fixationsDFtree$correct=="incorrect",]
for (p in 2:9){
  print(mean(fixationsDFtreeINCORR[,p]))
}
fixationsDFVennCORR <- fixationsDFvenn[fixationsDFvenn$correct=="correct",]
for (p in 2:6){
  print(mean(fixationsDFVennCORR[,p]))
}
fixationsDFVennINCORR <- fixationsDFvenn[fixationsDFvenn$correct=="incorrect",]
for (p in 2:6){
  print(mean(fixationsDFVennINCORR[,p]))
}


#long fixations means
fixationsDFtreeLongCORR <- fixationsDFtreeLong[fixationsDFtreeLong$correct=="correct",]
for (p in 2:9){
  print(mean(fixationsDFtreeLongCORR[,p]))
}
fixationsDFtreeLongINCORR <- fixationsDFtreeLong[fixationsDFtreeLong$correct=="incorrect",]
for (p in 2:9){
  print(mean(fixationsDFtreeLongINCORR[,p]))
}
fixationsDFVennLongCORR <- fixationsDFVennLong[fixationsDFVennLong$correct=="correct",]
for (p in 2:6){
  print(mean(fixationsDFVennLongCORR[,p]))
}
fixationsDFVennLongINCORR <- fixationsDFVennLong[fixationsDFVennLong$correct=="incorrect",]
for (p in 2:6){
  print(mean(fixationsDFVennLongINCORR[,p]))
}




