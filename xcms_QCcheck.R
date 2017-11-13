library(xcms)
getwd()
setwd("/media/maiski")
setwd("/media/maiski/Maiski's Hard Drive/QTOFdata")
#Q-TOF data from stool extract (pooled; IBD); 4GHz, HiRes mode
ibdt_pool1 <- "171110-02.mzXML"
#TOF data from control urine sample in IBS study; 2GHz, Extended Dynamic Range mode
#mzXML_TOF1 <- "150522_05.mzXML"

#Assigning above files into xcmsRaw object; matrix object
msi <- xcmsRaw(ibdt_pool1)
#msi_qtof <- xcmsRaw(mzXML3)

# findPeaks.centwave method 
foundpeaks <- findPeaks.centWave(msi, ppm = 50, peakwidth=c(10,30), scanrange=c(240, 2700), noise=500, integrate=2)

#xcmsSet::: This can be exported in csv format 
### For Q-TOF files, this code results in warning message about this not being centroid
xtestMSI_selected4 <- xcmsSet(mzXML3, method="centWave", ppm=30, peakwidth=c(8,25), scanrange=c(240, 1800))
xtestMSI_TOF5 <- xcmsSet(mzXML_TOF1, method="centWave", ppm=50, peakwidth=c(10,30), scanrange=c(800, 2900))

# Ourput csv files:::this is close to raw data; just a peak picking step is done by centWave 
msi_peaks <- peakTable(xtestMSI_selected4, filebase="PeakList")
write.csv(msi_peaks, "msi.csv") 
##
msi_peaks2 <- peakTable(xtestMSI_TOF5, filebase="PeakList")
write.csv(msi_peaks2, "msi_TOFCTL.csv") 


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
plotChrom(msi, mzrange=c(216.04, 216.045), scanrange=c(1767,2360))
plotChrom(msi, mzrange=c(325.159, 325.162), scanrange=c(1300, 1800))
## Without smoothing
profMethod(msi) <- "bin"
plotChrom(msi, mzrange=c(216.039, 216.045), scanrange=c(240,2700))

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

## group() based on the density of signals 
alignedset <- group(alignedfiles, method="density")
alignedset_cent <- group(CTLal_files_cent, method="density")

## retcor() to better align data
retcored <- retcor(alignedset_cent, method="loess", plottype = "mdevden")

### second group step after retcor
secondgroup <- group(retcored, method="density")

# getEIC()
## CTL no alignment
CTLeic <- getEIC(CTLfiles, mzrange=cbind(mzmin=325.15, mzmax=325.17), rtrange=cbind(rtmin=200, rtmax=1000))
plot(CTLeic)

## CTL aligned
CTLal_eic <- getEIC(CTLal_files, mzrange=cbind(mzmin=216.04, mzmax=216.045), rtrange=cbind(rtmin=950, rtmax=1700))
plot(CTLal_eic)
### CTL aligned with centWave method
CTLal_eic_cent <- getEIC(CTLal_files_cent, mzrange=cbind(mzmin=216.04, mzmax=216.045), rtrange=cbind(rtmin=950, rtmax=1700))
grouped_CTLal <- getEIC(retcored, mzrange=cbind(mzmin=216.04, mzmax=216.045), rtrange=cbind(rtmin=950, rtmax=1700))
grouped_CTLal2 <- getEIC(secondgroup, mzrange=cbind(mzmin=216.04, mzmax=216.045), rtrange=cbind(rtmin=950, rtmax=1700))
plot(CTLal_eic_cent, title="Cl-tyr after msalign")
plot(grouped_CTLal, title="Cl-tyr after msalign and retcor")
plot(grouped_CTLal2)

## test with single file
testeic <- getEIC(xtestMSI_TOF5, mzrange=cbind(mzmin=216.04, mzmax=216.045), rtrange=cbind(rtmin=1000, rtmax=1800))
plot(testeic)

## All data original
originaleic <- getEIC(originalfiles, mzrange=cbind(mzmin=216.04, mzmax=216.045), rtrange=cbind(rtmin=950, rtmax=1700))
plot(originaleic)

## CTL aligned
allal_eic <- getEIC(alignedfiles, mzrange=cbind(mzmin=325.15, mzmax=325.17), rtrange=cbind(rtmin=500, rtmax=900))
plot(allal_eic)