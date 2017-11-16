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

library(stringr)
getwd()
setwd("/home/maiski/R/IBD")
list.files()

### introduce data into the R environment 
pos1_10 <- read.csv("IBD_T_pool_pos_extracted_171110.csv")

#Check data structure 
summary(pos1_10)
head(pos1_10)
class(pos1_10)

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
  RMTtotal <- cbind(samplelabels, RMT)
  RPAtotal <- cbind(samplelabels, RPA)
  ## Output into CSV format 
  write.csv(RMTtotal, file = paste(outputname, "RMT.csv"))
  write.csv(RPAtotal, file = paste(outputname, "RPA.csv"))
}

## END of FUNCTION #################################


# Implementation of the function relativematrix() in my IBD data 
relativematrix(pos1_10, "IBDQCtest")

# Apply this to other files; but first introduce them into R environment 
pos11_18 <- read.csv("IBD_extracted_11-18_pos.csv")
pos19_23 <- read.csv("IBD_extracted_19-23_pos.csv")
pos24_26 <- read.csv("IBD_extracted_24-26_pos.csv")
pos27_32 <- read.csv("IBD_extracted_27-32_pos.csv")
neg1_10 <- read.csv("IBD_extracted_1-10_neg.csv")
neg11_17 <- read.csv("IBD_extracted_11-17_neg.csv")
neg18_23 <- read.csv("IBD_extracted_18-23_neg.csv")
neg24_32 <- read.csv("IBD_extracted_24-32_neg.csv")
# Use the relativematrix function
relativematrix(pos11_18, "IBDpos11-18")
relativematrix(pos19_23, "IBDpos19-23")
relativematrix(pos24_26, "IBDpos24-26")
relativematrix(pos27_32, "IBDpos27-32")
relativematrix(neg1_10, "IBDneg1-10")
relativematrix(neg11_17, "IBDneg11-17")
relativematrix(neg18_23, "IBDneg18-23")
relativematrix(neg24_32, "IBDneg24-32")

###########-- Merge replicates of samples --- ##############
# Re-introduce the RPA data 
## POS
IBDposRPA1_10 <- read.csv("IBDpos1-10 RPA.csv")
IBDposRPA11_18 <- read.csv("IBDpos11-18 RPA.csv")
IBDposRPA19_23 <- read.csv("IBDpos19-23 RPA.csv")
IBDposRPA24_26 <- read.csv("IBDpos24-26 RPA.csv")
IBDposRPA27_32 <- read.csv("IBDpos27-32 RPA.csv")

## NEG
IBDnegRPA1_10 <- read.csv("IBDneg1-10 RPA.csv")
IBDnegRPA11_17 <- read.csv("IBDneg11-17 RPA.csv")
IBDnegRPA18_23 <- read.csv("IBDneg18-23 RPA.csv")
IBDnegRPA24_32 <- read.csv("IBDneg24-32 RPA.csv")

osm <- read.csv("IBDosm.csv")
rownames(osm) <- osm$sample
osm_clean <- osm[, -1]

# Remove compounds that have low %CV (separately identified)
## POS
clean1_10 <- IBDposRPA1_10[, -which(names(IBDposRPA1_10) %in% "X104.0706.1")]
clean11_18 <- IBDposRPA11_18[, -which(names(IBDposRPA11_18) %in% "X104.0706.1")]
clean19_23 <- IBDposRPA19_23[, -which(names(IBDposRPA19_23) %in% "X104.0706.1")]
clean24_26 <- IBDposRPA24_26[, -which(names(IBDposRPA24_26) %in% "X104.0706.1")]
clean27_32 <- IBDposRPA27_32[, -which(names(IBDposRPA27_32) %in% "X104.0706.1")]

## NEG
negclean1_10 <- IBDnegRPA1_10[, -which(names(IBDnegRPA1_10) %in% c("X117.0555", "X117.0555.1", "X221.0746", "X222.9916", "X225.0629", "X324.0936", "X326.0881", "X352.0868"))]
negclean11_17 <- IBDnegRPA11_17[, -which(names(IBDnegRPA11_17) %in% c("X222.9916", "X225.0629", "X324.0936", "X326.0881", "X352.0868"))]
negclean18_23 <- IBDnegRPA18_23[, -which(names(IBDnegRPA18_23) %in% c("X222.9916", "X225.0629", "X324.0936", "X326.0881", "X352.0868"))]
negclean24_32 <- IBDnegRPA24_32[, -which(names(IBDnegRPA24_32) %in% c("X222.9916", "X225.0629", "X324.0936", "X326.0881", "X352.0868"))]

# Sort the matrix based on "sample" column
ordered1_10 <- clean1_10[with(clean1_10, order(sample)), ]

# Remove QC data 
sample1_10 <- ordered1_10[!grepl("QC", ordered1_10[,3]),]
  
# Take the average of replicates vased on the sample name; 
## --- take available one if one of the value is NA
## --- take a larger value if %difference from average is > 30%

# subset the matrix
RPAlabels <- sample1_10[, 1:3]
RPAsamples <- sample1_10[, -c(1:3)]

# Set up the number of rows and columns for the final matrix
cmpds <- length(colnames(RPAsamples))
samplenum <- length(rownames(sample1_10))/2

# Empty matrix for results to be added
replmean <- data.frame(matrix(nrow = samplenum, ncol = cmpds))
colnames(replmean) <- colnames(RPAsamples)

#Set up Osmolality matrix
## POS
osm1_10 <- osm_clean[1:36, ]
osm11_18 <- osm_clean[37:72, ]
osm19_23 <- osm_clean[73:90, ]
osm24_26 <- osm_clean[91:102, ]
osm27_32 <- osm_clean[103:117, ]

## NEG
negosm1_10 <- osm_clean[1:36, ]
negosm11_17 <- osm_clean[37:64, ]
negosm18_23 <- osm_clean[65:90, ]
negosm24_32 <- osm_clean[91:117, ]

i = 1
for (i in 1:cmpds) {
  x = 1
  y = 2
  z = 1
  for (s in 1:samplenum){
    replicate1 <- RPAsamples[x, i]
    replicate2 <- RPAsamples[y, i]
    if (is.na(replicate1) && is.na(replicate2)) {
      replmean[z, i] <- "NA"
    } else if (is.na(replicate1) || is.na(replicate2)) {
      replmean[z, i] <- mean(c(replicate1, replicate2), na.rm = TRUE)
    } else if (abs(((mean(c(replicate1, replicate2)) - replicate1) / mean(c(replicate1, replicate2)))) > 0.3) {
      replmean[z, i] <- max(c(replicate1, replicate2))
    } else {
      replmean[z, i] <- mean(c(replicate1, replicate2))
    }
    rownames(replmean)[z] <- as.character(RPAlabels[x, 3])
    x = x + 2
    y = y + 2
    z = z + 1
  }
  i = i + 1
}
# Merge osmolality matrix with the matrix created above to ensure the match 
merged1_10 <- merge(replmean, osm1_10, by = "row.names")

# normalize each value with corresponding osmolality 
normal1_10 <- data.frame(matrix(nrow = samplenum, ncol = cmpds))
colnames(normal1_10) <- colnames(RPAsamples)
rownames(normal1_10) <- rownames(merged1_10)
for (c in 1:cmpds) {
  n = 2
  normal1_10[, c] <- merged1_10[, n] / merged1_10[, 58]
  n = n + 1
  c = c + 1
}
### Final matrix of normalized values for each sample and columns for low volume and dilution observation
final1_10 <- merge(normal1_10, osm1_10[, -1], by = "row.names")
write.csv(final1_10, file = "IBDpos_norm_1-10_model.csv")

cleanmatrix <-clean11_18
osmmatrix <- osm11_18
replicate_tolerance <-0.3
outputname <- "IBDpos_norm_11-18.csv"

###***************** FUNCTION************************#################################

# Function to combine the replicates and normalize with osmolality 
merge_norm <- function(cleanmatrix, osmmatrix, replicate_tolerance, outputname){
  # Sort the matrix based on "sample" column
  ordered <- cleanmatrix[with(cleanmatrix, order(sample)), ]
  
  # Remove QC data 
  samples <- ordered[!grepl("QC", ordered[,3]),]
  
  # Take the average of replicates vased on the sample name; 
  ## --- take available one if one of the value is NA
  ## --- take a larger value if %difference from average is > 30%
  
  # subset the matrix
  RPAlabels <- samples[, 1:3]
  RPAsamples <- samples[, -c(1:3)]
  
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
        replmean[z, i] <- "NA"
      } else if (is.na(replicate1) || is.na(replicate2)) {
        replmean[z, i] <- mean(c(replicate1, replicate2), na.rm = TRUE)
      } else if (abs(((mean(c(replicate1, replicate2)) - replicate1) / mean(c(replicate1, replicate2)))) > replicate_tolerance) {
        replmean[z, i] <- max(c(replicate1, replicate2))
      } else {
        replmean[z, i] <- mean(c(replicate1, replicate2))
      }
      rownames(replmean)[z] <- as.character(RPAlabels[x, 3])
      x = x + 2
      y = y + 2
      z = z + 1
    }
    i = i + 1
  }
  
  # Merge osmolality matrix with the matrix created above to ensure the match 
  merged <- merge(replmean, osmmatrix, by = "row.names")
  merged[, c(1:2)] <- str_split_fixed(merged$Row.names, "U", 2)
  colnames(merged)[c(1:2)] <- c("subject", "timepoint")
  # normalize each value with corresponding osmolality 
  # normalized <- data.frame(matrix(nrow = samplenum, ncol = cmpds))
  # colnames(normalized) <- colnames(RPAsamples)
  # rownames(normalized) <- rownames(replmean)
  # c = 1
  # n = 2
  # for (compd in 1:cmpds) {
  #   normalized[, c] <- merged[, n] / as.numeric(merged[, 58])
  #   n = n + 1
  #   c = c + 1
  # }
  # ### Final matrix of normalized values for each sample and columns for low volume and dilution observation
  write.csv(merged, file = outputname)
}

## END of FUNCTION ####


## Implement the merge_norm() function 
## POS
merge_norm(clean1_10, osm1_10, 0.3, "IBDpos_osm_1-10.csv")
merge_norm(clean11_18, osm11_18, 0.3,"IBDpos_osm_11-18.csv")
merge_norm(clean19_23, osm19_23, 0.3,  "IBDpos_osm_19-23.csv")
merge_norm(clean24_26, osm24_26, 0.3,  "IBDpos_osm_24-26.csv")
merge_norm(clean27_32, osm27_32, 0.3, "IBDpos_osm_27-32.csv")

## NEG
merge_norm(negclean1_10, negosm1_10, 0.3, "IBDneg_osm_1-10.csv")
merge_norm(negclean11_17, negosm11_17, 0.3,"IBDneg_osm_11-17.csv")
merge_norm(negclean18_23, negosm18_23, 0.3,  "IBDneg_osm_18-23.csv")
merge_norm(negclean24_32, negosm24_32, 0.3,  "IBDneg_osm_24-32.csv")

## Testing process ##
x <- RPAsamples[1, 1]
y <- RPAsamples[2, 1]
x <- RPAsamples[23, 21]
y <- RPAsamples[24, 21]
x <- RPAsamples[31, 21]
y <- RPAsamples[32, 21]
x <- RPAsamples[35, 21]
y <- RPAsamples[36, 21]

if (is.na(x) && is.na(y)) {
  print("NA")
} else if (is.na(x) || is.na(y)) {
  print(mean(c(x, y), na.rm = TRUE))
} else if (abs(((mean(c(x, y)) - x) / mean(c(x, y)))) > 0.3) {
  print (max(c(x,y)))
} else {
  print (mean(c(x, y)))
}
##############################

testIBD <- read.csv("IBDpos_osm_1-10.csv")
testIBD[, c(1:2)] <- str_split_fixed(testIBD$Row.names, "U", 2)
colnames(testIBD)[c(1:2)] <- c("subject", "timepoint")
