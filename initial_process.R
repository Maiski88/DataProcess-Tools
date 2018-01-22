####################
### Title: Initial data process after copy + paste from MassHunter
### Date: July 22, 2017
### Authur: Mai Yamamoto 

## This file include codes for 
# --- generating RMT and RPA matrices: relativematrix() ------line 
# --- merge replicates to generate one row per sample: merge_norm()
##
##
##
##
##
##

## call the required library for percent()
library('scales')

getwd()
setwd("your/home/directory")
# list.files() shows you what you have in the specified directory 
list.files()

### introduce data into the R environment 
rawdata <- read.csv("your_file_name.csv")

#Check data structure 
summary(rawdata)
head(rawdata)
class(rawdata)

### Support function ######################
# %CV calculator 
cvcalc <- function(input_matrix){
  QCmat <- subset(input_matrix, sample == 'QC')
  n = 3
  datalist = list()
  for (i in 3:length(colnames(QCmat))){
    percv <- round((sd(QCmat[, n], na.rm = TRUE)/mean(QCmat[, n], na.rm = TRUE))*100, digits = 2)
    datalist[[n]] <- percv # add it to your list
    n <- n + 1 
  }
  big_data = do.call(cbind, datalist)
  return(big_data)
}

############*************FUNCTION*******************######################
# Function to generate RMT and RPA matrices 
relativematrix <- function(inputmatrix, outputname) {
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
  ## Calculate %CV of QC samples 
  RMTQC <<- subset(RMTtotal, sample == 'QC')
  RPAQC <- subset(RPAtotal, sample == 'QC')
  RMTcv <- cvcalc(RMTQC)
  RPAcv <- cvcalc(RPAQC)
  cvmatrix <<- rbind(RMTcv, RPAcv) 
  dimnames(cvmatrix) <<-  list(c("RMTcv", "RPAcv"), c(colnames(RMTQC)[-(1:2)]))
  ## Output into CSV format 
  write.csv(RMTtotal, file = paste(outputname, "RMT.csv"), row.names=FALSE)
  write.csv(RPAtotal, file = paste(outputname, "RPA.csv"), row.names=FALSE)
  write.csv(cvmatrix, file = paste(outputname, "%CV.csv"))
}
## END of FUNCTION 


###***************** FUNCTION************************

# Function to remove compounds from the relative matrix when %CV of QC exceeds threshhold
irrep_rm <- function(cvmatrix, CVmax){
  for (compound in colnames(cvmatrix)){
    if (is.na(cvmatrix[2, compound]) || cvmatrix[2, compound] > CVmax) {
      RPAtotal <<- RPAtotal[ , -which(names(RPAtotal) %in% compound)]
    } 
  }
}
## END of FUNCTION 

###***************** FUNCTION************************

# Function to combine the replicates and normalize with osmolality 
## replicate tolerance applies when there is a large difference between two peaks
duplmerge <- function(RPAtotal, replicate_tolerance, outputname){
  # Sort the matrix based on "sample" column
  ordered <- RPAtotal[with(RPAtotal, order(sample)), ]
  
  # Remove QC and blank data 
  samples <- ordered[!grepl("QC|blank|F-Phe|Blank", ordered[, 'sample']),]
  
  # Take the average of replicates based on the sample name; 
  ## --- take available one if one of the value is NA
  ## --- replace the value with NA if %difference from average is > replicate_tolerance
  
  # subset the matrix
  RPAlabels <- samples[, 1:2]
  RPAsamples <- samples[, -c(1:2)]
  
  # Set up the number of rows and columns for the final matrix
  cmpds <- length(colnames(RPAsamples))
  samplenum <- length(rownames(samples))/2
  
  # Empty matrix for results to be added
  replmean <- data.frame(matrix(nrow = samplenum, ncol = cmpds))
  colnames(replmean) <- colnames(RPAsamples)
  
  i = 1
  for (i in 1:cmpds) {
    x = 1
    y = 2
    z = 1
    for (s in 1:samplenum){
      replicate1 <- RPAsamples[x, i]
      replicate2 <- RPAsamples[y, i]
      if (is.na(replicate1) && is.na(replicate2)) {
        replmean[z, i] <- "NA" #NA if both replicates are NA
      } else if (is.na(replicate1) || is.na(replicate2)) { #if one of them is integrated, it stands
        replmean[z, i] <- mean(c(replicate1, replicate2), na.rm = TRUE)
      } else if (abs(((mean(c(replicate1, replicate2)) - replicate1) / mean(c(replicate1, replicate2)))) > replicate_tolerance*0.01) {
        replmean[z, i] <- "NA" #if %diff > tolerance (eg. 30%), assign NA
      } else {
        replmean[z, i] <- mean(c(replicate1, replicate2)) #If peaks look pretty much the same, take the mean
      }
      rownames(replmean)[z] <- as.character(RPAlabels[x, 2])
      x = x + 2
      y = y + 2
      z = z + 1
    }
    i = i + 1
  }
  meansample <<- replmean
  # ### Final matrix of normalized values for each sample and columns for low volume and dilution observation
  write.csv(replmean, file = paste(outputname,".csv"))
}

## END of FUNCTION ####

###***************** FUNCTION************************
# Function to remove compounds that appears in ~ 75% of samples or less 
rmrare <- function(inputmatrix, na_perlimit, metabolite_num, outputname){
  n = 1
  highna <- NULL
  naperlist <- NULL
  rowlength <- length(rownames(inputmatrix))
  metabmatrix <- inputmatrix[, 2:(metabolite_num+1)]
  for (i in colnames(metabmatrix)){
    napercentage <- (sum(is.na(inputmatrix$i))/rowlength)*100
    if (napercentage > na_perlimit){
      highna <<- append(highna, i)
      naperlist <<- append(naperlist, napercentage)
      n <- n + 1 
    }
  }
  #finalmatrix <<- inputmatrix[, -which(names(inputmatrix) %in% c(highna))]
  #removed <<- data.frame(inputmatrix[, which(names(inputmatrix) %in% c(highna))], row.names = inputmatrix[, 'sample'])
  #removedmat <<- rbind(removed, naperlist)
  #rownames(removedmat)[rowlength+1] <- 'NApercentage'
  #write.csv(removedmat, file = paste0(outputname, "_highNA.csv"))
}

