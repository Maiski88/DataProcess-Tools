library(xcms)
setwd("/home/maiski/R/Data/xcms/CTL")
list.files()

# scanrange of c(900, 2700) corresponds to 10-30 min RT 
# human urine sample in positive ion mode analyzed using CE-TOF-MS 
mzXML_TOF <- "150522_05.mzXML"

# xcmsSet(): finds peaks in mzXML file and put it into xcmsSet object 
xtestMSI_TOF <- xcmsSet(mzXML_TOF, method="centWave", ppm=50, peakwidth=c(10,30), scanrange=c(900, 2700))

## peakTable takes a xcmsSet object 
msi_peaks <- peakTable(xtestMSI_TOF, filebase="PeakList")
ordered_msi_peaks <- msi_peaks[order(msi_peaks[, "rt"]),]
write.csv(ordered_msi_peaks, "msi2.csv") 

# create a xcmsRaw object from the mzXML file 
msi <- xcmsRaw(mzXML_TOF)

############################## Check the data ########################################


# Check EIC in profile mode 
# mzrange gives an error if you put 3rd decimal places because size is too large! 
##
profStep(msi) <- 0.05
## Checking TIC with smoothing
profMethod(msi) <- "intlin"
plotTIC(msi)

# check scantime that corresponds to seconds 
which(msi@scantime > 1100 & msi@scantime < 1102)[1]
which(msi@scantime > 1800 & msi@scantime < 1802)[1]

# Check EIC of m/z 325.16
plotChrom(msi, mzrange=c(325.15, 325.17), scanrange=c(1600, 2600))

########################################################################################

# centWave method to create xcmsPeaks object 
foundpeaks <- findPeaks.centWave(msi, ppm = 50, peakwidth=c(10,30), scanrange=c(900, 2700), noise=50, integrate=1, snthresh=3)
ordered_foundpeaks <- foundpeaks[order(foundpeaks[, "rt"]),]
msi_peaks <- peakTable(xtestMSI_TOF, filebase="PeakList")
ordered_msi_peaks <- msi_peaks[order(msi_peaks[, "rt"]),]

################ SUPPORT FUNCTION #################################################
## This prepares function:peakpick to remove NULL from a featureList
is.Nullprep <- function(x) is.null(x) | all(sapply(x, is.null))

## Recursively step down into list, removing all such objects 
rmNulls <- function(x) {
  x <- Filter(Negate(is.Nullprep), x)
  return(x)
}



## == MAIN FUNCTION == #######################################################################################
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
## END of FUNCTION: peakpick() #################################################

# Test peakpick() 
test_masslist <- list(90.055, 106.0499, 138.0550, 162.113, 166.0723, 182.0809, 205.0972)
result <- peakpick(test_masslist, foundpeaks, MassErr = 20)



########################################################################################
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
## ENF of FUNCTION: rpacheck() ###############################################################

# test rpacheck() 
test_masslist <- list(90.055, 106.0499, 138.0550, 162.113, 166.0723, 182.0809, 205.0972)
result2 <- rpacheck(test_masslist, ordered_msi_peaks, MassErr = 20, ISmz = 216.0419)

### NOT to RUN #####################################################################
# BASIS of peakpick() & rpacheck()
## Extract IS peaks and create a list of integrated peak area (ISorder)
ISmz <- 216.0419
MassErr <- 10
injection_num <- 7
xcmsmatrix <- ordered_foundpeaks
massdev <- mz * (MassErr*10**(-6))
IS <- which(xcmsmatrix[, "mz"] > (mz - massdev) & xcmsmatrix[, "mz"] < (mz + massdev))
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
