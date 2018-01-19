#--------------------------------------------------------------------------------------------------------------------------- 
# Function: StringExtractor
# Description: This is a Algorithm to produce a DF with scanpath strings for each participants by information format.
# In the second part, the rest of the data (the file with experiment perfomance) is added to this DF to create a new DF.
# In the thrid part, this latter DF is extended including the duration for the dwell time for AOIs, and a final DF is created. 
# Author: manuele Reani 
# Date 03/07/2017
#----------------------------------------------------------------------------------------------------------------------------

library("car")

# we start with the datafile containing the eye-tracking data 
df <- read.csv(file="study2ET.csv",header=TRUE,sep=",")

# we create an empty DF where to store all the scanpaths
scanpathDF <- data.frame(Participant = NULL, Fomrat = NULL, Scanpath = NULL)

##############################################################
#------------- Find scanpaths for TREE ----------------------#
##############################################################

# subset DF for format(TREE)
treeDF<- subset(df, Format == "tree") 

# For TREE, find the scanpath for each participant (loop through Ps and Fixations)  
for (p in 2:50){
  participant <- p 
  # subset DF for participant
  tree <- subset(treeDF, Participant == p) 
  
  backgroundAOI <- names(tree)[5]
  scanpath <- NULL
  
  for (fixation in tree$FixationIndex){
    # print(fixation)
    aois <- tree[tree$FixationIndex == fixation, 6:13]
    aoi <- names(aois)[aois == 1]
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    #print(aoi)
    scanpath <- append(scanpath, aoi)
  }
  
  #------------------------------------------------------------------------------------------------------
  # This code converts a vector into a string, removing subsequent repetitions of characters for scanpath
  #-------------------------------------------------------------------------------------------------------
  # we need to rename the AOIs first
  treevars <- '"treebg"="A"; "treetotal"="B"; "treerainy"="C"; "treesunny"="D"; "treerl"="E"; "treerh"="F";"treesl"="G"; "treesh"="H"; "treequestion"="I"'
  scanpath <- car::recode(scanpath, treevars)
  
  # convert the vector into string
  scanpath <- toString(scanpath)
  
  # remove commas and white spaces from the string
  scanpath <- gsub(",", "", scanpath)
  scanpath <- gsub(" ", "", scanpath)
  
  # remove subsequent repeated characters in the string
  scanpath <- gsub('([[:alpha:]])\\1+', '\\1', scanpath)
  
  newRow <- data.frame(Participant = participant, Format = "tree",Scanpath = scanpath)
  scanpathDF <- rbind(scanpathDF,newRow)
}


##############################################################
#------------- Same for VENN --------------------------------#
##############################################################


# subset DF for format(VENN)
vennDF<- subset(df, Format == "venn") 

# For VENN, find the scanpath for each participant (loop through Ps and Fixations) 
for (p in 2:50){
  participant <- p 
  # subset DF for participant
  venn <- subset(vennDF, Participant == p) 
  
  backgroundAOI <- names(venn)[14]
  scanpath <- NULL
  
  for (fixation in venn$FixationIndex){
    # print(fixation)
    aois <- venn[fixation, 15:19]
    aoi <- names(aois)[aois == 1]
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    #print(aoi)
    scanpath <- append(scanpath, aoi)
  }
  
  #------------------------------------------------------------------------------------------------------
  # This code converts a vector into a string, removing subsequent repetitions of characters for scanpath
  #-------------------------------------------------------------------------------------------------------
  vennvars <- '"vennbg"="L"; "vennrh"="M";"vennrl"="N"; "vennsl"="O"; "venntotal"="P"; "vennquestion"="Q"'
  scanpath <- car::recode(scanpath, vennvars) 
  
  # convert the vector into string
  scanpath <- toString(scanpath)
  
  # remove commas and white spaces from the string
  scanpath <- gsub(",", "", scanpath)
  scanpath <- gsub(" ", "", scanpath)
  
  # remove subsequent repeated characters in the string
  scanpath <- gsub('([[:alpha:]])\\1+', '\\1', scanpath)
  
  newRow <- data.frame(Participant = participant, Format = "venn",Scanpath = scanpath)
  scanpathDF <- rbind(scanpathDF,newRow)
}

#################################################################################
#------ Merge scanpathDF with the dataBayesian.csv file (called df2) -----------#
#################################################################################

df2 <- read.csv(file="dataBayesian.csv",header=TRUE,sep=",")

# unfactorize the column Scanpath in both dataframes 
df2$Scanpath <- as.character(df2$Scanpath)
scanpathDF$Scanpath <- as.character(scanpathDF$Scanpath)

# mathc the scanpaths by Ps and Format, e.g.:  
# df2[which(df2$Ps== 2 & df2$Format == "t"), 10] <- scanpathDF[which(scanpathDF$Participant== 2 & scanpathDF$Format == "tree"), 3]

for (n in 2:50){
  df2[which(df2$Ps== n & df2$Format == "t"), 10] <- scanpathDF[which(scanpathDF$Participant== n & scanpathDF$Format == "tree"), 3]
  df2[which(df2$Ps== n & df2$Format == "v"), 10] <- scanpathDF[which(scanpathDF$Participant== n & scanpathDF$Format == "venn"), 3]
}

#################################################################################################### 
#----- The fixation duration is missing here, so we want to include it in the dataframe -----------#
#################################################################################################### 

# We need to create an empty DF to store that 
durationAOI <- data.frame(Participant = NULL, Format = NULL, AOI = NULL, Time = NULL)

# subset DF for format(TREE)
treeDF<- subset(df, Format == "tree")

# store AOIs with fixations durations
for (p in 2:50){
  # subset DF for participant
  tree <- subset(treeDF, Participant == p)
  backgroundAOI <- names(tree)[5]
  
  for (fixation in tree$FixationIndex){
    # print(fixation)
    aois <- tree[fixation, 6:13]
    duration <- tree[fixation, 4]
    aoi <- names(aois)[aois == 1]
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    newRow <- data.frame(Participant = p, Format = "tree", AOI = aoi, Time = duration)
    durationAOI <- rbind(durationAOI,newRow)
  }
}

# subset DF for format(VENN)
vennDF<- subset(df, Format == "venn")

# store AOIs with fixations durations
for (p in 2:50){
  # subset DF for participant
  venn <- subset(vennDF, Participant == p)
  backgroundAOI <- names(venn)[14]
  
  for (fixation in venn$FixationIndex){
    # print(fixation)
    aois <- venn[fixation, 15:19]
    duration <- venn[fixation, 4]
    aoi <- names(aois)[aois == 1]
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    newRow <- data.frame(Participant = p, Format = "venn", AOI = aoi, Time = duration)
    durationAOI <- rbind(durationAOI,newRow)
  }
}


####################################################################################
#------- Join the two dataframes (df2 and durationAOI) into one called df3 --------# 
####################################################################################


# install.packages("dplyr")
# install.packages("tidyr")
library(dplyr)
library(tidyr)

# Create a new df summing time over AOIs (by participant and by format), 
# in order to have only 1 value for each AOIs for each participant. 
# You could use: sum(durationAOI$Time[durationAOI$AOI=="treebg"])
# for each participant, but this will require an expensive for loop.
# The method below uses the dplyr library which is specifically developed for data wrangling.
# See: https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf

totalDuration <- durationAOI %>% group_by(Participant, Format, AOI) %>% summarise(totalTime=sum(Time))

#--- Comment to the code above --------------------------------------------------------------------------------------------#
# If you want, with this library, you could also have obtained the average fixation dwell time 9not just the sum) 
# for each AOI (by participant and by format) using this: 
# totalDuration <- durationAOI %>% group_by(Participant, Format, AOI) %>% summarise(totalTime=sum(Time), avtTime=mean(Time))
#--------------------------------------------------------------------------------------------------------------------------#

# change the veritical dataframe (long) into a a horizontal dataframe (wide) and store it into a new daframe. 
# You could do it using "reshape" package but we use "tidyr" package instead.
# See:https://rpubs.com/bradleyboehmke/data_wrangling
totalDurationWide <- totalDuration %>% spread(AOI, totalTime, drop=FALSE)

# in this new df, change the names of the AOIs into letters, using the recoding created before for Tree and for Venn
names(totalDurationWide) <- car::recode(names(totalDurationWide), treevars)
names(totalDurationWide) <- car::recode(names(totalDurationWide), vennvars)

# reorder the column names of this df in alphabetical order 
totalDurationWide <- totalDurationWide[,c("Participant", "Format", LETTERS[1:9], LETTERS[12:17])]

# in df2, rename the levels of the variable Format 
df2$Format <- car::recode(df2$Format, '"t"="tree";"v"="venn"')

# create a new df (called df3) merging the two dataframe (df2 and titalDurationWide) by Participant and by Format 
df3 <- dplyr::inner_join(df2,totalDurationWide, by = c("Ps"="Participant", "Format"="Format"))


###########################################################################################
#------- Save this final dataframe into a csv file ---------------------------------------#
###########################################################################################

# Create a CSV file using theis new data frame (df3) 
write.csv(df3, file = "bayeScanpaths.csv", row.names=F)


#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------#

#                               HERE WE DO THE SAME BUT WE INCLUDE REPETITIONS                            #

#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------#



##############################################################
#------------- Find scanpaths for TREE ----------------------#
##############################################################
scanpathDFrep <- data.frame(Participant = NULL, Fomrat = NULL, Scanpath = NULL) # scanpathDF                   produce a new df with repetitions ??? 
# subset DF for format(TREE)


# For TREE, find the scanpath for each participant (loop through Ps and Fixations)  
for (p in 2:50){
  participant <- p 
  # subset DF for participant
  tree <- subset(treeDF, Participant == p) 
  
  backgroundAOI <- names(tree)[5]
  scanpath <- NULL
  
  for (fixation in tree$FixationIndex){
    # print(fixation)
    aois <- tree[fixation, 6:13]
    aoi <- names(aois)[aois == 1]
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    #print(aoi)
    scanpath <- append(scanpath, aoi)
  }
  
  #------------------------------------------------------------------------------------------------------
  # This code converts a vector into a string, removing subsequent repetitions of characters for scanpath
  #-------------------------------------------------------------------------------------------------------
  # we need to rename the AOIs first
  treevars <- '"treebg"="A"; "treetotal"="B"; "treerainy"="C"; "treesunny"="D"; "treerl"="E"; "treerh"="F";"treesl"="G"; "treesh"="H"; "treequestion"="I"'
  scanpath <- car::recode(scanpath, treevars)
  
  # convert the vector into string
  scanpath <- toString(scanpath)
  
  # remove commas and white spaces from the string
  scanpath <- gsub(",", "", scanpath)
  scanpath <- gsub(" ", "", scanpath)
  
  newRow <- data.frame(Participant = participant, Format = "tree",Scanpath = scanpath)
  scanpathDFrep <- rbind(scanpathDFrep,newRow)
}


##############################################################
#------------- Same for VENN --------------------------------#
##############################################################


# subset DF for format(VENN)


# For VENN, find the scanpath for each participant (loop through Ps and Fixations) 
for (p in 2:50){
  participant <- p 
  # subset DF for participant
  venn <- subset(vennDF, Participant == p) 
  
  backgroundAOI <- names(venn)[14]
  scanpath <- NULL
  
  for (fixation in venn$FixationIndex){
    # print(fixation)
    aois <- venn[fixation, 15:19]
    aoi <- names(aois)[aois == 1]
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    #print(aoi)
    scanpath <- append(scanpath, aoi)
  }
  
  #------------------------------------------------------------------------------------------------------
  # This code converts a vector into a string, removing subsequent repetitions of characters for scanpath
  #-------------------------------------------------------------------------------------------------------
  vennvars <- '"vennbg"="L"; "vennrh"="M";"vennrl"="N"; "vennsl"="O"; "venntotal"="P"; "vennquestion"="Q"'
  scanpath <- car::recode(scanpath, vennvars) 
  
  # convert the vector into string
  scanpath <- toString(scanpath)
  
  # remove commas and white spaces from the string
  scanpath <- gsub(",", "", scanpath)
  scanpath <- gsub(" ", "", scanpath)
  
  newRow <- data.frame(Participant = participant, Format = "venn",Scanpath = scanpath)
  scanpathDFrep <- rbind(scanpathDFrep,newRow)
}

#################################################################################
#------ Merge scanpathDFrep with the dataBayesian.csv file (called df2) -----------#
#################################################################################

# unfactorize the column Scanpath in the dataframes 
scanpathDFrep$Scanpath <- as.character(scanpathDFrep$Scanpath)
df2rep <- read.csv(file="dataBayesian.csv",header=TRUE,sep=",")
for (n in 2:50){
  df2rep[which(df2rep$Ps== n & df2rep$Format == "t"), 10] <- scanpathDFrep[which(scanpathDFrep$Participant== n & scanpathDFrep$Format == "tree"), 3]
  df2rep[which(df2rep$Ps== n & df2rep$Format == "v"), 10] <- scanpathDFrep[which(scanpathDFrep$Participant== n & scanpathDFrep$Format == "venn"), 3]
}

#################################################################################################### 
#----- The fixation duration is missing here, so we want to include it in the dataframe -----------#
#################################################################################################### 

# We need to create an empty DF to store that 
durationAOIrep <- data.frame(Participant = NULL, Format = NULL, AOI = NULL, Time = NULL)

# subset DF for format(TREE)
treeDF<- subset(df, Format == "tree")

# store AOIs with fixations durations
for (p in 2:50){
  # subset DF for participant
  tree <- subset(treeDF, Participant == p)
  backgroundAOI <- names(tree)[5]
  
  for (fixation in tree$FixationIndex){
    # print(fixation)
    aois <- tree[fixation, 6:13]
    duration <- tree[fixation, 4]
    aoi <- names(aois)[aois == 1]
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    newRow <- data.frame(Participant = p, Format = "tree", AOI = aoi, Time = duration)
    durationAOIrep <- rbind(durationAOIrep,newRow)
  }
}

# subset DF for format(VENN)
vennDF<- subset(df, Format == "venn")

# store AOIs with fixations durations
for (p in 2:50){
  # subset DF for participant
  venn <- subset(vennDF, Participant == p)
  backgroundAOI <- names(venn)[14]
  
  for (fixation in venn$FixationIndex){
    # print(fixation)
    aois <- venn[fixation, 15:19]
    duration <- venn[fixation, 4]
    aoi <- names(aois)[aois == 1]
    if (length(aoi) == 0){
      aoi <- backgroundAOI
    } 
    else if (length(aoi) > 1) {
      stop("Multiple AOIs found")
    }
    newRow <- data.frame(Participant = p, Format = "venn", AOI = aoi, Time = duration)
    durationAOIrep <- rbind(durationAOIrep,newRow)
  }
}


####################################################################################
#------- Join the two dataframes (df2rep and durationAOIrep) into one called df3rep --------# 
####################################################################################


# install.packages("dplyr")
# install.packages("tidyr")
library(dplyr)
library(tidyr)

# Create a new df summing time over AOIs (by participant and by format), 
# in order to have only 1 value for each AOIs for each participant. 
# You could use: sum(durationAOIrep$Time[durationAOIrep$AOI=="treebg"])
# for each participant, but this will require an expensive for loop.
# The method below uses the dplyr library which is specifically developed for data wrangling.
# See: https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf

totalDurationrep <- durationAOIrep %>% group_by(Participant, Format, AOI) %>% summarise(totalTime=sum(Time))

#--- Comment to the code above --------------------------------------------------------------------------------------------#
# If you want, with this library, you could also have obtained the average fixation dwell time 9not just the sum) 
# for each AOI (by participant and by format) using this: 
# totalDurationrep <- durationAOIrep %>% group_by(Participant, Format, AOI) %>% summarise(totalTime=sum(Time), avtTime=mean(Time))
#--------------------------------------------------------------------------------------------------------------------------#

# change the veritical dataframe (long) into a a horizontal dataframe (wide) and store it into a new daframe. 
# You could do it using "reshape" package but we use "tidyr" package instead.
# See:https://rpubs.com/bradleyboehmke/data_wrangling
totalDurationWiderep <- totalDurationrep %>% spread(AOI, totalTime, drop=FALSE)

# in this new df, change the names of the AOIs into letters, using the recoding created before for Tree and for Venn
names(totalDurationWiderep) <- car::recode(names(totalDurationWiderep), treevars)
names(totalDurationWiderep) <- car::recode(names(totalDurationWiderep), vennvars)

# reorder the column names of this df in alphabetical order 
totalDurationWiderep <- totalDurationWiderep[,c("Participant", "Format", LETTERS[1:9], LETTERS[12:17])]

# in df2rep, rename the levels of the variable Format 
df2rep$Format <- car::recode(df2rep$Format, '"t"="tree";"v"="venn"')

# create a new df (called df3rep) merging the two dataframe (df2rep and titalDurationWide) by Participant and by Format 
df3rep <- dplyr::inner_join(df2rep,totalDurationWiderep, by = c("Ps"="Participant", "Format"="Format"))


###########################################################################################
#------- Save this final dataframe into a csv file ---------------------------------------#
###########################################################################################

# Create a CSV file using theis new data frame (df3rep) 
write.csv(df3rep, file = "bayeScanpathsrep.csv", row.names=F)
