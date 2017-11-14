library(xcms)
getwd()
setwd("/media/maiski")
setwd("/media/maiski/Maiski's Hard Drive/QTOFdata")
#Q-TOF data from stool extract (pooled; IBD); 4GHz, HiRes mode
# This data was acquired with 1 scan / sec
ibdt_pool1 <- "171110-02.mzXML"
ibdt_pool1_32 <- "171110-02_32bit.mzXML"
#TOF data from control urine sample in IBS study; 2GHz, Extended Dynamic Range mode
#mzXML_TOF1 <- "150522_05.mzXML"

setwd("/home/maiski/R/DataProcess-Tools")

#Assigning above files into xcmsRaw object; matrix object
msi <- xcmsRaw(ibdt_pool1)
msi_32 <- xcmsRaw(ibdt_pool1_32)
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
foundpeaks <- findPeaks.centWave(msi, ppm = 50, peakwidth=c(10,30), scanrange=c(450, 1750), noise=500, integrate=2)
ordered_foundpeaks <- foundpeaks[order(foundpeaks[, "rt"]),]

foundpeaks_32 <- findPeaks.centWave(msi_32, ppm = 50, peakwidth=c(10,30), scanrange=c(450, 1750), noise=500, integrate=2)
ordered_foundpeaks_32 <- foundpeaks_32[order(foundpeaks_32[, "rt"]),]

#xcmsSet::: This can be exported in csv format 
### For Q-TOF files, this code results in warning message about this not being centroid
#xtestMSI_selected4 <- xcmsSet(mzXML3, method="centWave", ppm=30, peakwidth=c(8,25), scanrange=c(240, 1800))
msi_ibdpool1 <- xcmsSet(ibdt_pool1, method="centWave", ppm=50, peakwidth=c(10,30), scanrange=c(450, 1800))

# Ourput csv files:::this is close to raw data; just a peak picking step is done by centWave 
msi_peaks <- peakTable(msi_ibdpool1, filebase="PeakList")
ordered_msi_peaks <- msi_peaks[order(msi_peaks[, "rt"]),]
write.csv(ordered_msi_peaks, "msi.csv") 

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
MassErr <- 25
injection_num <- 7
xcmsmatrix <- ordered_foundpeaks
xcmsmatrix_32 <- ordered_foundpeaks_32
massdev <- mz * (MassErr*10**(-6))
# IS picking with 64 bit data
IS <- which(xcmsmatrix[, "mz"] > (mz - massdev) & xcmsmatrix[, "mz"] < (mz + massdev))
IStime <- xcmsmatrix[IS, "rt"]
testmat <- as.data.frame(xcmsmatrix)
IStime <- as.list(IStime)
# IS picking with 32 bit data
IS_32 <- which(xcmsmatrix_32[, "mz"] > (mz - massdev) & xcmsmatrix_32[, "mz"] < (mz + massdev))
IStime_32 <- xcmsmatrix_32[IS, "rt"]


ISlist <- vector(mode="list", length=7)
if(peaknumber == 7) {
  print(paste0(mz, " :Perfect peak picking!"))
  featureList[[i]] <<- mz
  i <- i + 1
}else if(peaknumber > 7) {
  print(paste0(mz," :Large noise or isomers present; number of peaks is: ", peaknumber))
}else {
  print(paste0(mz, " :There are only ", peaknumber, " peaks"))
}
Order = 0
for(peak in IS){
  #print(ordered_foundpeaks[peak, ])
  Order <- Order + 1 
  ISlist[[Order]] <- ordered_foundpeaks[peak, "into"] 
}


# Function takes a list of accurate mass and peak matrix to see if there are given number of peaks for each feature
peakpick <- function(masslist, xcmsmatrix, MassErr = 5, injection_num = 7) {
  i <- 1
  featureList <<- vector("list", length=length(masslist))
  for(mz in masslist) {
    massdev <- mz * (MassErr*10**(-6))
    feature <- which(xcmsmatrix[, "mz"] > (mz - massdev) & xcmsmatrix[, "mz"] < (mz + massdev))
    peaknumber <- length (feature)
    if(peaknumber == 7) {
      print(paste0(mz, " :Perfect peak picking!"))
      featureList[[i]] <<- mz
      i <- i + 1
    }else if(peaknumber > 7) {
      print(paste0(mz," :Large noise or isomers present; number of peaks is: ", peaknumber))
    }else {
      print(paste0(mz, " :There are only ", peaknumber, " peaks"))
    }
  }
  new_featureList <<- rmNulls(featureList)
  return(new_featureList)
}

# Test peakpick 
test_masslist <- list(90.055, 106.0499, 138.0550, 162.113, 166.0723, 182.0809, 205.0972)
result <- peakpick(test_masslist, foundpeaks, MassErr = 20)

test <- list(1, 10, 100, 200,300, 400, 500)
test_division <- as.numeric(ISorder)/test
lapply(seq_along(ISorder),function(n)
  unlist(ISorder[n])/unlist(test[n]))


# Modified peakpick() to incorporate RPA calculation 
rpacheck <- function(masslist, xcmsmatrix, MassErr = 5, injection_num = 7, ISmz) {
  ISmassdev <- ISmz * (MassErr*10**(-6))
  IS <- which(xcmsmatrix[, "mz"] > (ISmz - ISmassdev) & xcmsmatrix[, "mz"] < (ISmz + ISmassdev))
  # Just in case IS peak picking fails
  if (!(length(IS) == 7)) {
    print("Adjust peakpicking parameter to detect all IS peaks")
    stop
  }
  ISlist <- vector(mode="list", length=7)
  ISrt <- vector(mode="list", length=7)
  Order = 0
  # ISlist contains integrated peak area of IS
  for(peak in IS){
    Order <- Order + 1 
    ISlist[[Order]] <- xcmsmatrix[peak, "into"] 
    ISrt[[Order]] <- xcmsmatrix[peak, "rt"] 
  }
  rpamatrix <<- data.frame(matrix(ncol = (length(masslist)*2), nrow = 7))
  i <- 1
  # 
  featureList <<- vector("list", length=length(masslist))
  for(mz in masslist) {
    massdev <- mz * (MassErr*10**(-6))
    feature <- which(xcmsmatrix[, "mz"] > (mz - massdev) & xcmsmatrix[, "mz"] < (mz + massdev))
    peaknumber <- length (feature)
    # Perfect peakpicking case leads to rpa calculation using the ISlist
    if(peaknumber == 7) {
      print(paste0(mz, " :Perfect peak picking!"))
      featureList[[i]] <<- mz
      # Create a list of analyte peak area
      # featurert is the list of migration time in seconds
      featurerpa <- vector(mode="list", length=7)
      featurert <<- vector(mode="list", length=7)
      Order = 0
      for(peak in feature){
        Order <- Order + 1 
        featurerpa[[Order]] <- xcmsmatrix[peak, "into"] 
        featurert[[Order]] <- xcmsmatrix[peak, "rt"] 
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

test_masslist <- list(90.055, 106.0499, 138.0550, 162.113, 166.0723, 182.0809, 205.0972)
result2 <- rpacheck(test_masslist, ordered_msi_peaks, MassErr = 20, ISmz = 216.0419)

dt <- data.frame(10, 10)
colnames(dt)[1] <- 1
colnames(dt)[1+1] <- 2
rownames(dt) <- "test"






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

