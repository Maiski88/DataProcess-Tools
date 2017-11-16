library(xcms)
getwd()
setwd("/media/maiski/Maiski's Hard Drive/QTOFdata")
#Q-TOF data from stool extract (pooled; IBD); 4GHz, HiRes mode
# This data was acquired with 1 scan / sec
ibdt_pool1 <- "171110-02.mzXML"
ibdt_pool2 <- "171110-03.mzXML"
ibdt_pool3 <- "171110-04.mzXML"
ibdt_pool4 <- "171110-05.mzXML"
ibdt_pool5 <- "171110-06.mzXML"
#TOF data from control urine sample in IBS study; 2GHz, Extended Dynamic Range mode
#mzXML_TOF1 <- "150522_05.mzXML"

#Assigning above files into xcmsRaw object; matrix object
msi <- xcmsRaw(ibdt_pool1)
msi2 <- xcmsRaw(ibdt_pool2)
msi3 <- xcmsRaw(ibdt_pool3)
msi4 <- xcmsRaw(ibdt_pool4)
msi5 <- xcmsRaw(ibdt_pool5)
#msi_qtof <- xcmsRaw(mzXML3)

##::::::: Visualize data and know your scan frequency ::::::###########
#Find out the time in seconds for each scan number 
which(msi@scantime > 250 & msi@scantime < 252)[1]
which(msi@scantime > 1760 & msi@scantime < 1762)[1]

# Check EIC in profile mode 
# mzrange gives an error if you put 3rd decimal places (when step is set as 0.05)! 
##
profStep(msi) <- 0.005
## With smoothing
profMethod(msi) <- "intlin"
plotTIC(msi)
plotChrom(msi, mzrange=c(216.04, 216.045), scanrange=c(1000,1750))
plotChrom(msi, mzrange=c(325.159, 325.162), scanrange=c(1300, 1800))
## Without smoothing
profMethod(msi) <- "bin"
plotChrom(msi, mzrange=c(216.039, 216.045), scanrange=c(240,2700))

# findPeaks.centwave method 
## peakwidth(3, 30) works well; noise = 1000 and integration = 2 results in 38802 ROI with variable peak# 4800-6500
## Narrow peaks (eg. 147.1129) cannot be picked up unless peakwidth(2, X), which miss many others. 
foundpeaks <- findPeaks.centWave(msi, ppm = 50, peakwidth=c(3.5,30), scanrange=c(450, 1750), noise=1000, integrate=2)
ordered_foundpeaks <- foundpeaks[order(foundpeaks[, "rt"]),]

foundpeaks2 <- findPeaks.centWave(msi2, ppm = 50, peakwidth=c(3.5,30), scanrange=c(450, 1750), noise=1000, integrate=2)
ordered_foundpeaks2 <- foundpeaks2[order(foundpeaks2[, "rt"]),]

foundpeaks3 <- findPeaks.centWave(msi3, ppm = 50, peakwidth=c(3.5,30), scanrange=c(450, 1750), noise=1000, integrate=2)
ordered_foundpeaks3 <- foundpeaks3[order(foundpeaks3[, "rt"]),]

foundpeaks4 <- findPeaks.centWave(msi4, ppm = 50, peakwidth=c(3.5,30), scanrange=c(450, 1750), noise=1000, integrate=2)
ordered_foundpeaks4 <- foundpeaks4[order(foundpeaks4[, "rt"]),]

foundpeaks5 <- findPeaks.centWave(msi5, ppm = 50, peakwidth=c(4,30), scanrange=c(450, 1750), noise=1000, integrate=2)
ordered_foundpeaks5 <- foundpeaks5[order(foundpeaks5[, "rt"]),]

##
## This prepares function:peakpick to remove NULL from a featureList
is.Nullprep <- function(x) is.null(x) | all(sapply(x, is.null))

## Recursively step down into list, removing all such objects 
rmNulls <- function(x) {
  x <- Filter(Negate(is.Nullprep), x)
  return(x)
}


# Function that defines IS peaks in order of migration based on m/z and time 
## Extract IS peaks and create a list of integrated peak area (ISorder)
ISmz <- 216.0422
mz <- 216.0422
MassErr <- 10
injection_num <- 7
massdev <- mz * (MassErr*10**(-6))
# IS picking with 32bit data
IS <- which(ordered_foundpeaks[, "mz"] > (mz - massdev) & ordered_foundpeaks[, "mz"] < (mz + massdev))
IStime <- ordered_foundpeaks[IS, "rt"]
peaknumber <- length(IStime)

Order = 0
for(peak in IS){
  #print(ordered_foundpeaks[peak, ])
  Order <- Order + 1 
  ISlist[[Order]] <- ordered_foundpeaks[peak, "into"] 
}

######### FUNCTION ::: peakpick() ::::::::::::::::::::::::::::::############
# Function takes a list of accurate mass and peak matrix to see if there are given number of peaks for each feature
peakpick <- function(masslist, xcmsmatrix, MassErr = 5, injection_num = 7) {
  i <- 1
  featureList <<- vector("list", length=length(masslist))
  for(mz in masslist) {
    massdev <- mz * (MassErr*10**(-6))
    feature <- which(xcmsmatrix[, "mz"] > (mz - massdev) & xcmsmatrix[, "mz"] < (mz + massdev))
    peaknumber <- length (feature)
    if(peaknumber == injection_num) {
      print(paste0(mz, " :Perfect peak picking!"))
      featureList[[i]] <<- mz
      i <- i + 1
    }else if(peaknumber > injection_num) {
      print(paste0(mz," :Large noise or isomers present; number of peaks is: ", peaknumber))
    }else {
      print(paste0(mz, " :There are only ", peaknumber, " peaks"))
    }
  }
  new_featureList <<- rmNulls(featureList)
  return(new_featureList)
}

# Test peakpick 
masslist <- list(216.0422, 72.0807, 101.0596, 104.1069, 106.0499, 120.0655, 131.1179, 132.0768, 137.0457, 147.0764, 154.0499, 156.0767, 162.1125, 170.0924, 175.119, 189.1598, 241.0312, 268.104, 407.2385, 471.2185)
result <- peakpick(masslist, ordered_foundpeaks, MassErr = 10, injection_num = 6)
result2 <- peakpick(masslist, ordered_foundpeaks2, MassErr = 10, injection_num = 6)
result3 <- peakpick(masslist, ordered_foundpeaks3, MassErr = 10, injection_num = 6)
result4 <- peakpick(masslist, ordered_foundpeaks4, MassErr = 10, injection_num = 6)
result5 <- peakpick(masslist, ordered_foundpeaks5, MassErr = 10, injection_num = 6)

#IS matrix generation
ISmz <- 216.0422
mz <- 216.0422
MassErr <- 10
injection_num <- 7
massdev <- mz * (MassErr*10**(-6))
# IS time check to see if there are numbers of peaks that are supposed to present
IS <- which(ordered_foundpeaks[, "mz"] > (mz - massdev) & ordered_foundpeaks[, "mz"] < (mz + massdev))
IStime <- ordered_foundpeaks[IS, "rt"]
# IS mt and area list generation
ISlist <- vector(mode="list", length=7)
ISrt <- vector(mode="list", length=7)
Order = 0
for(peak in IS){
  Order <- Order + 1 
  ISlist[[Order]] <- ordered_foundpeaks[peak, "into"]
  ISrt[[Order]] <- ordered_foundpeaks[peak, "rt"] 
}
## Get rid of a peak at blank position
ISlist <- ISlist[-4]
ISrt <- ISrt[-4]
## Combine the lists into one data frame and change the row names 
ISmat <- do.call(rbind, Map(data.frame, rt=ISrt, into=ISlist))
rownames(ISmat) <- c("QC1", "QC2", "QC3", "QC4", "QC5", "QC6")

# Reference matrix (refmat) prep
refrmt_raw <- read.csv("IBDQCtest RMT.csv")
refrpa_raw <- read.csv("IBDQCtest RPA.csv")
## Get rid of those compounds that couldn't be picked by foundpeak()
refrmt_raw <- refrmt_raw[-c(8,9, 13, 15, 17, 24)]
refrpa_raw <- refrpa_raw[-c(8,9, 13, 15, 17, 24)]
# Empty data frame generation to input mean of each rmt and rpa
refmat <- data.frame(matrix(vector(), length(4:22), 2,
                       dimnames=list(c(masslist), c("rt", "into"))),
                stringsAsFactors=F)
## n = starting column where the analyte matrix starts 
n = 4
i = 1
for (m in n:22){
  refmat[i, 1] <- mean(refrmt_raw[, n])
  refmat[i, 2] <- mean(refrpa_raw[, n])
  i <- i + 1 
  n <- n + 1
}

### TEST: Combine lists into data.frame and division of data frame by another
list1 <- c(1:5)
list2 <- c(6:10)
list3 <- c(1:5)
list4 <- c(6:10)
mat1 <- do.call(rbind, Map(data.frame, rt=list1, into=list2))
mat2 <- do.call(rbind, Map(data.frame, rt=list3, into=list4))
matdiv <- mat2/mat1

####### FUNCTION:::::::: qccheck() ::::::::::::::::::::::::::::::::::
# Input precalculated IS mt & peak area matrix to accommodate a case where blank is injected (ie. IS peak# = 7, analyte = 6)
qccheck <- function(masslist, ordered_foundpeaks, MassErr = 10, injection_num = 7, ISmat, refmat, ISposition){
  i <- 1
  featureList <<- vector("list", length=length(masslist))
  for(mz in masslist) {
    massdev <- mz * (MassErr*10**(-6))
    feature <- which(ordered_foundpeaks[, "mz"] > (mz - massdev) & ordered_foundpeaks[, "mz"] < (mz + massdev))
    peaknumber <- length (feature)
    featurerpa <<- vector(mode="list", length=injection_num)
    featurert <<- vector(mode="list", length=injection_num)
    # Perfect peakpicking case leads to rpa calculation using the ISlist
    if(peaknumber == injection_num) {
      print(paste0(mz, " :Perfect peak picking!"))
      featureList[[i]] <<- mz
      i <- i + 1
      # Create a list of analyte peak area
      # featurert is the list of migration time in seconds
      Order = 1
      for(p in feature){
        metapeak <- feature[Order]
        featurerpa[[Order]] <- ordered_foundpeaks[metapeak, "into"] 
        featurert[[Order]] <- ordered_foundpeaks[metapeak, "rt"] 
        Order <- Order + 1 
      }
      rawmat <- do.call(rbind, Map(data.frame, rt=featurert, into=featurerpa))
      samplemat <- rawmat / ISmat
      percentdiff_rmt <- (samplemat[ISposition, rt] - refmat[as.character(mz), rt])/refmat[as.character(mz), rt]
      percentdiff_rpa <- (samplemat[ISposition, into] - refmat[as.character(mz), into])/refmat[as.character(mz), into]
      print(paste0(mz, " RMT % difference is: ", percentdiff_rmt, " RPA % difference is: ", percentdiff_rpa))
    }else if(peaknumber > injection_num) {
      print(paste0(mz," :Large noise or isomers present; number of peaks is: ", peaknumber))
    }else {
      print(paste0(mz, " :There are only ", peaknumber, " peaks"))
    }
  }
}


qccheck(masslist, ordered_foundpeaks, injection_num = 6, ISmat, refmat, 5)


# Modified peakpick() to incorporate RPA calculation 
rpacheck <- function(masslist, ordered_foundpeaks, MassErr = 10, injection_num = 7, ISmz) {
  ISmassdev <- ISmz * (MassErr*10**(-6))
  IS <- which(ordered_foundpeaks[, "mz"] > (ISmz - ISmassdev) & ordered_foundpeaks[, "mz"] < (ISmz + ISmassdev))
  # Just in case IS peak picking fails
  if (!(length(IS) == 7)) {
    print("Adjust peakpicking parameter to detect all IS peaks")
    stop
  }
  ISlist <- vector(mode="list", length=7)
  ISrt <- vector(mode="list", length=7)
  Order = 0
  # ISlist is a "list" that contains integrated peak area of IS
  for(peak in IS){
    Order <- Order + 1 
    ISlist[[Order]] <- ordered_foundpeaks[peak, "into"] 
    ISrt[[Order]] <- ordered_foundpeaks[peak, "rt"] 
  }
  rpamatrix <<- data.frame(matrix(ncol = (length(masslist)*2), nrow = 7))
  i <- 1
  # 
  featureList <<- vector("list", length=length(masslist))
  for(mz in masslist) {
    massdev <- mz * (MassErr*10**(-6))
    feature <- which(ordered_foundpeaks[, "mz"] > (mz - massdev) & ordered_foundpeaks[, "mz"] < (mz + massdev))
    peaknumber <- length (feature)
    featurerpa <<- vector(mode="list", length=injection_num)
    featurert <<- vector(mode="list", length=injection_num)
    # Perfect peakpicking case leads to rpa calculation using the ISlist
    if(peaknumber == injection_num) {
      print(paste0(mz, " :Perfect peak picking!"))
      featureList[[i]] <<- mz
      # Create a list of analyte peak area
      # featurert is the list of migration time in seconds
      Order = 0
      for(peak in feature){
        Order <- Order + 1 
        featurerpa[[Order]] <- ordered_foundpeaks[peak, "into"] 
        featurert[[Order]] <- ordered_foundpeaks[peak, "rt"] 
      }
      rpalist <<- lapply(seq_along(ISlist),function(n)
        unlist(featurerpa[n])/unlist(ISlist[n]))
      rtlist <<- lapply(seq_along(ISrt),function(n)
        unlist(featurert[n])/unlist(ISrt[n]))
      rpamatrix[[i]] <<- unlist(rtlist)
      rpamatrix[[i+1]] <<- unlist(rpalist)
      colnames(rpamatrix)[i] <<- paste(mz, ":RMT")
      colnames(rpamatrix)[i+1] <<- paste(mz, ":RPA")
      i <- i + 2
      
    }else if(peaknumber > 7) {
      print(paste0(mz," :Large noise or isomers present; number of peaks is: ", peaknumber))
    }else {
      print(paste0(mz, " :There are only ", peaknumber, " peaks"))
    }
  }
  new_featureList <<- rmNulls(featureList)
  return(new_featureList)
  rpamatrix <<- Filter(function(x)!all(is.na(x)), rpamatrix)
  rownames(rpamatrix) <- c(1:7)
}

result2 <- rpacheck(test_masslist, ordered_foundpeaks, MassErr = 10, injection_num = 6, ISmz = 216.0422)

dt <- data.frame(10, 10)
colnames(dt)[1] <- 1
colnames(dt)[1+1] <- 2
rownames(dt) <- "test"

#xcmsSet::: This can be exported in csv format 
### For Q-TOF files, this code results in warning message about this not being centroid
#xtestMSI_selected4 <- xcmsSet(mzXML3, method="centWave", ppm=30, peakwidth=c(8,25), scanrange=c(240, 1800))
msi_ibdpool1 <- xcmsSet(ibdt_pool1, method="centWave", ppm=50, peakwidth=c(10,30), scanrange=c(450, 1800))

# Ourput csv files:::this is close to raw data; just a peak picking step is done by centWave 
setwd("/home/maiski/R/DataProcess-Tools")
msi_peaks <- peakTable(msi_ibdpool1, filebase="PeakList")
ordered_msi_peaks <- msi_peaks[order(msi_peaks[, "rt"]),]
write.csv(ordered_msi_peaks, "msi.csv") 

# Set up multi-file data process
## ---CTL not aligned 
CTL_path <- "/home/maiski/R/Data/xcms/CTL"
CTL <- list.files(CTL_path, recursive = TRUE, full.names = TRUE)
CTLfiles <- xcmsSet(CTL)

## ---CTL aligned
CTLal_path <- "/home/maiski/R/Data/xcms/CTL_aligned"
CTL_al <- list.files(CTLal_path, recursive = TRUE, full.names = TRUE)
CTLal_files <- xcmsSet(CTL_al)

### CTL aligned with centWave xcmsSet 
CTLal_path <- "/home/maiski/R/Data/xcms/CTL_aligned"
CTL_al <- list.files(CTLal_path, recursive = TRUE, full.names = TRUE)
CTLal_files_cent <- xcmsSet(CTL_al, method="centWave", ppm=50, peakwidth=c(10,30), scanrange=c(800, 2900))

## All data original
original_path <- "/home/maiski/R/Data/xcms/IBS_everything/original"
original <- list.files(original_path, recursive = TRUE, full.names = TRUE)
originalfiles <- xcmsSet(original)

## All aligned data
aligned_path <- "/home/maiski/R/Data/xcms/IBS_everything/aligned"
aligned <- list.files(aligned_path, recursive = TRUE, full.names = TRUE)
alignedfiles <- xcmsSet(aligned)

