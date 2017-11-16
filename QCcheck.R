####################
### Title: QC check after each run or sequence 
### Date: Nov. 16, 2017
### Authur: Mai Yamamoto 

## This file include codes for 
# --- 
##
##
##
##
##
##

getwd()
setwd("/home/maiski/R/DataProcess-Tools")
list.files()

library('scales')
### introduce data into the R environment 
qctest <- read.csv("IBD_T_pool_pos_extracted_171110.csv")

############*************FUNCTION*******************######################
# Function to generate RMT and RPA matrices 
relativematrix <- function(inputmatrix) {
  # Clean up the data frame to remove unnecessary stuff (eg. "RT", "SNR")
  cleanmatrix <- inputmatrix[!grepl("RT", inputmatrix[,3]),]
  samplelabels <- cleanmatrix[, 1:2] # Subset the data frame 
  ISmatrix <- cleanmatrix[, 3:4]
  analytes <- cleanmatrix[, 6:length(colnames(cleanmatrix))]
  cmpds <- length(colnames(analytes))/3 #number of compounds included in the matrix 
  x = 1
  y = 2
  ## Create empty data frames to input RMT & RPA calculation results 
  RMT <- data.frame(matrix(nrow = length(rownames(cleanmatrix)), ncol = cmpds))
  RPA <- data.frame(matrix(nrow = length(rownames(cleanmatrix)), ncol = cmpds))
  for (i in 1:cmpds) { # loop to calculate the RMT and RPA 
    RMT[, i] <- as.numeric(levels(analytes[, x]))[analytes[, x]] / as.numeric(levels(ISmatrix[, 1]))[ISmatrix[, 1]]
    RPA[, i] <- as.numeric(levels(analytes[, y]))[analytes[, y]] / as.numeric(levels(ISmatrix[, 2]))[ISmatrix[, 2]]
    colnames(RMT)[i] <- colnames(analytes)[x]
    colnames(RPA)[i] <- colnames(analytes)[x]
    x = x + 3
    y = y + 3 
  }
  ## Naming each sample as they were 
  RMTtotal <<- cbind(samplelabels, RMT)
  RPAtotal <<- cbind(samplelabels, RPA)
  ## Check cummulative %CV for RMT and RPA
  n = 3
  for (i in 3:length(colnames(RMTtotal))){
    rmtcv <- percent(sd(RMTtotal[, n])/mean(RMTtotal[, n]))
    rpacv <- percent(sd(RPAtotal[, n])/mean(RPAtotal[, n]))
    print(paste0("For ", colnames(RMTtotal[n]), "RMT and RPA %CVs are ", rmtcv, " and ", rpacv))
    n <- n + 1 
  }
}

relativematrix(qctest)

## If needed, generate csv file of RMT and RPA matrix 
write.csv(RMTtotal, "QC_RMT.csv")
write.csv(RPAtotal, "QC_RPA.csv")
