#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Resummarization of 2016 Chinook Troll Mixtures ####
# Kyle Shedd Mon Aug 14 15:47:51 2017
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
# The goal of this script is to re-summarize Chinook salmon mixtures from SEAK
# commercial and sport troll harvest from 2016. Results have been provided to 
# 33 and 26 fine-scale groups. The intent here is to modify the 17RG to 18RG 
# to include NSEAK as a group (appear at >5% in sport mixtures), to summarize
# to 8 "driver stock" RGs, and 4 broad scale TBR RGs. The baseline used was
# GAPS 3.0 containing 357 populations grouped in 26RGs characterized by 13
# uSats.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Specific Objectives ####
# This script will:
# 1) Import 2016 troll objects and RGs from 2015
# 2) Resummarize BAYES results from .RGN files
#     2016 Comm Troll to 8, 18, and 4
#     Sport to 8, 18, and 4
#     Heatmaps for Comm and Sport for 8
#     Retrospective Sport to 8 for 2009-2016
#     Sample sizes
# 3) Generate plots and tables of results
# 4) Explore new figure ideas

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK15")

## Get objects
SEAK15objects <- list.files(path = "Objects", recursive = FALSE)
SEAK15objects <- SEAK15objects[!SEAK15objects %in% c("sillys_sport.txt", "sillys_sportD8&11.txt", "sillys_troll.txt")]
SEAK15objects

invisible(sapply(SEAK15objects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

# Grab important objects (i.e. groupnames and groupvecs) and dput in "Objects"
objects2dput <- c("GAPSLoci", 
                  "GAPSLoci_reordered", 
                  "GroupNames17",
                  "GroupNames18", 
                  "GroupNames26", 
                  "GroupNames26Pub", 
                  "GroupNames33", 
                  "GroupNames4", 
                  "GroupNames8",
                  "GroupNames8Pub",
                  "GroupVec17", 
                  "GroupVec18",
                  "GroupVec26RG_357", 
                  "GroupVec33RG_to26RG",
                  "GroupVec33RG_to8RG",
                  "GroupVec4",
                  "GroupVec8")

setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK16")

invisible(sapply(objects2dput, function(obj) {
  dput(x = get(obj), file = paste0("Objects/", obj, ".txt"))
})); rm(objects2dput)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK16")

# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
SEAK16objects <- list.files(path = "Objects", recursive = FALSE)
SEAK16objects <- SEAK16objects[!SEAK16objects %in% c("Sport2016GenotypesRaw.txt")]
SEAK16objects

invisible(sapply(SEAK16objects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Resummarize 2016 Comm Troll Stratified Estimates to Driver Stock RGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Define RGs

writeClipboard(GroupNames26[GroupVec33RG_to26RG])
writeClipboard(GroupNames8[GroupVec8[GroupVec33RG_to26RG]])

GroupVec33RG_to8RG <- GroupVec8[GroupVec33RG_to26RG]
dput(x = GroupVec33RG_to8RG, file = "Objects/GroupVec33RG_to8RG.txt")

WinterTrollMix2016 <- c("EWintNISISO_2016", "EWintNO_2016", "LWintNISISO_2016", "LWintNO_2016")
SpringTrollMix2016 <- c("SpringNI_2016", "SpringNO_2016", "SpringSI_2016", "SpringSO_2016")
SummerTrollMix2016 <- c("Summer1NI_2016", "Summer1NO_2016", "Summer1SI_2016", "Summer1SO_2016",
                        "Summer2NI_2016", "Summer2NO_2016", "Summer2SI_2016", "Summer2SO_2016")

invisible(sapply(c("WinterTrollMix2016", "SpringTrollMix2016", "SummerTrollMix2016"), function(obj) {
  dput(x = get(obj), file = paste0("Objects/", obj, ".txt"))
}))


### Unstratified
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.create("BAYES/Output/33RG/AllYearTroll_2016")


# 8RGs
EWintTroll2016_8RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                               maindir = "BAYES/Output/33RG/EWint_2016", 
                               mixvec = WinterTrollMix2016[1:2], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

LWintTroll2016_8RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                               maindir = "BAYES/Output/33RG/LWint_2016", 
                               mixvec = WinterTrollMix2016[3:4], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SpringTroll2016_8RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                               maindir = "BAYES/Output/33RG/Spring_2016", 
                               mixvec = SpringTrollMix2016, prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SumRet1Troll2016_8RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                               maindir = "BAYES/Output/33RG/SumRet1_2016", 
                               mixvec = SummerTrollMix2016[1:4], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SumRet2Troll2016_8RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                               maindir = "BAYES/Output/33RG/SumRet2_2016", 
                               mixvec = SummerTrollMix2016[5:8], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

TrollByQuad2016_8RG_EstimatesStats <- c(EWintTroll2016_8RG_EstimatesStats,
                                        LWintTroll2016_8RG_EstimatesStats,
                                        SpringTroll2016_8RG_EstimatesStats,
                                        SumRet1Troll2016_8RG_EstimatesStats,
                                        SumRet2Troll2016_8RG_EstimatesStats)
any(sapply(Troll2016_8RG_EstimatesStats, function(mix) {any(mix[, "GR"] > 1.2)}))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dput files
# dir.create("Estimates objects")
# Grab estimates objects and dput in "Estimates objects"
objects2dput <- c("EWintTroll2016_8RG_EstimatesStats", 
                  "LWintTroll2016_8RG_EstimatesStats", 
                  "SpringTroll2016_8RG_EstimatesStats", 
                  "SumRet1Troll2016_8RG_EstimatesStats",
                  "SumRet2Troll2016_8RG_EstimatesStats",
                  "TrollByQuad2016_8RG_EstimatesStats")

invisible(sapply(objects2dput, function(obj) {
  dput(x = get(obj), file = paste0("Estimates objects/", obj, ".txt"))
})); rm(objects2dput)



### Stratified
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8RGs
EWintTroll2016_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                          maindir="BAYES/Output/33RG/EWint_2016",
                          mixvec = c("EWintNISISO_2016", "EWintNO_2016"), catchvec = c(4216, 25147), 
                          newname = "StratifiedEWint2016_90percentCI_8RG", priorname = "", nchains = 5)

LWintTroll2016_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                          maindir="BAYES/Output/33RG/LWint_2016",
                          mixvec = c("LWintNISISO_2016", "LWintNO_2016"), catchvec = c(5248, 17680), 
                          newname = "StratifiedLWint2016_90percentCI_8RG", priorname = "", nchains = 5)

SpringTroll2016_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                          maindir="BAYES/Output/33RG/Spring_2016",
                          mixvec = c("SpringNI_2016", "SpringNO_2016", "SpringSI_2016", "SpringSO_2016"), catchvec = c(9270, 17012, 15160, 1031), 
                          newname = "StratifiedSpring2016_90percentCI_8RG", priorname = "", nchains = 5)

SumRet1Troll2016_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                          maindir="BAYES/Output/33RG/SumRet1_2016",
                          mixvec = c("Summer1NI_2016", "Summer1NO_2016", "Summer1SI_2016", "Summer1SO_2016"), catchvec = c(3805, 80323, 3618, 18888), 
                          newname = "StratifiedSumRet12016_90percentCI_8RG", priorname = "", nchains = 5)

SumRet2Troll2016_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                          maindir="BAYES/Output/33RG/SumRet2_2016",
                          mixvec = c("Summer2NI_2016", "Summer2NO_2016", "Summer2SI_2016", "Summer2SO_2016"), catchvec = c(2147, 56208, 1774, 14111), 
                          newname = "StratifiedSumRet22016_90percentCI_8RG", priorname = "", nchains = 5)


AllYearTroll2016_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to8RG, groupnames = GroupNames8,
                          maindir="BAYES/Output/33RG/AllYearTroll_2016",
                          mixvec = c("EWintNISISO_2016", "EWintNO_2016", "LWintNISISO_2016", "LWintNO_2016", 
                                     "SpringNI_2016", "SpringNO_2016", "SpringSI_2016", "SpringSO_2016", 
                                     "Summer1NI_2016", "Summer1NO_2016", "Summer1SI_2016", "Summer1SO_2016", 
                                     "Summer2NI_2016", "Summer2NO_2016", "Summer2SI_2016", "Summer2SO_2016"),
                          catchvec = c(4216, 25147, 5248, 17680,
                                       9270, 17012, 15160, 1031,
                                       3805, 80323, 3618, 18888,
                                       2147, 56208, 1774, 14111), 
                          newname = "StratifiedAllYearTroll2016_90percentCI_8RG", priorname = "", nchains = 5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dput files
# Grab estimates objects and dput in "Estimates objects"
objects2dput <- c("EWintTroll2016_8RG_StratifiedEstimatesStats", 
                  "LWintTroll2016_8RG_StratifiedEstimatesStats", 
                  "SpringTroll2016_8RG_StratifiedEstimatesStats", 
                  "SumRet1Troll2016_8RG_StratifiedEstimatesStats",
                  "SumRet2Troll2016_8RG_StratifiedEstimatesStats", 
                  "AllYearTroll2016_8RG_StratifiedEstimatesStats")

invisible(sapply(objects2dput, function(obj) {
  dput(x = get(obj)$Stats, file = paste0("Estimates objects/", obj, ".txt"))
})); rm(objects2dput)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create 2016 summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get objects
SEAK16estimatesobjects <- list.files(path = "Estimates objects", recursive = FALSE, pattern = "_8RG")
SEAK16estimatesobjects <- SEAK16estimatesobjects[-c(grep(pattern = "AllYearTroll", x = SEAK16estimatesobjects), 10)]
SEAK16estimatesobjects

invisible(sapply(SEAK16estimatesobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Estimates objects", objct, sep = "/")), pos = 1) })); beep(2)


# Dget all estimates stats
SEAK16estimatesobjects <- unlist(lapply(SEAK16estimatesobjects, function(objct) {unlist(strsplit(x = objct, split = ".txt"))}))

Troll2016_8RG_EstimatesStats <- list(
  "EWintNO_2016" = EWintTroll2016_8RG_EstimatesStats[["EWintNO_2016"]],
  "EWintAllQuad_2016" = EWintTroll2016_8RG_StratifiedEstimatesStats$Stats,
  "LWintNO_2016" = LWintTroll2016_8RG_EstimatesStats[["LWintNO_2016"]],
  "LWintAllQuad_2016" = LWintTroll2016_8RG_StratifiedEstimatesStats$Stats,
  "SpringNO_2016" = SpringTroll2016_8RG_EstimatesStats[["SpringNO_2016"]],
  "SpringSI_2016" = SpringTroll2016_8RG_EstimatesStats[["SpringSI_2016"]],
  "SpringAllQuad_2016" = SpringTroll2016_8RG_StratifiedEstimatesStats$Stats,
  "SumRet1NO_2016" = SumRet1Troll2016_8RG_EstimatesStats[["Summer1NO_2016"]],
  "SumRet1AllQuad_2016" = SumRet1Troll2016_8RG_StratifiedEstimatesStats$Stats,
  "SumRet2NO_2016" = SumRet2Troll2016_8RG_EstimatesStats[["Summer2NO_2016"]],
  "SumRet2AllQuad_2016" = SumRet2Troll2016_8RG_StratifiedEstimatesStats$Stats
)
dput(x = Troll2016_8RG_EstimatesStats, file = "Estimates objects/Troll2016_8RG_EstimatesStats.txt")

# Check GR
any(sapply(Troll2016_8RG_EstimatesStats, function(mix) {any(mix[, "GR"] > 1.2)}))


# Reformat estimates stats
Troll2016_8RG_EstimatesStats_Formatted <- sapply(Troll2016_8RG_EstimatesStats, function(yr) {
  matrix(data = yr[, 1:5], nrow = 8, ncol = 5, dimnames = list(GroupNames8Pub, c("Mean", "SD", "Median", "5%", "95%")))
}, simplify = FALSE)

Troll2016PubNames <- setNames(object = c("Northern Outside Quadrant",
                                         "All Quadrants",
                                         "Northern Outside Quadrant",
                                         "All Quadrants",
                                         "Northern Outside Quadrant",
                                         "Southern Inside Quadrant",
                                         "All Quadrants",
                                         "Northern Outside Quadrant",
                                         "All Quadrants",
                                         "Northern Outside Quadrant",
                                         "All Quadrants"), 
                              nm = names(Troll2016_8RG_EstimatesStats_Formatted))
dput(x = Troll2016PubNames, file = "Objects/Troll2016PubNames.txt")

SEAK2016Mixtures <- list.files(path = "BAYES/Mixture", full.names = FALSE)
SEAK2016Mixtures_SampSizes <- sapply(SEAK2016Mixtures, function(mix) {dim(read.table(file = paste0("BAYES/Mixture/", mix)))[1]} )
names(SEAK2016Mixtures_SampSizes) <- sapply(names(SEAK2016Mixtures_SampSizes), function(mix) {unlist(strsplit(x = mix, split = ".mix"))[1]})

Troll2016MixNames <- setNames(object = list("EWintNO_2016",
                                            WinterTrollMix2016[1:2],
                                            "LWintNO_2016",
                                            WinterTrollMix2016[3:4],
                                            "SpringNO_2016",
                                            "SpringSI_2016",
                                            SpringTrollMix2016,
                                            "Summer1NO_2016",
                                            SummerTrollMix2016[1:4],
                                            "Summer2NO_2016",
                                            SummerTrollMix2016[5:8]),
                              nm = names(Troll2016_8RG_EstimatesStats_Formatted))
dput(x = Troll2016MixNames, file = "Objects/Troll2016MixNames.txt")


Troll2016_SampleSizes <- sapply(Troll2016MixNames, function(mix) {sum(SEAK2016Mixtures_SampSizes[mix])} )

# Create fully formatted spreadsheat
EstimatesStats <- Troll2016_8RG_EstimatesStats_Formatted
SampSizes <- Troll2016_SampleSizes

# dir.create("Estimates tables")

for(mix in names(EstimatesStats)) {
  
  TableX <- matrix(data = "", nrow = 11, ncol = 7)
  TableX[1, 3] <- paste(Troll2016PubNames[mix], "(n=", SampSizes[mix], ")")
  TableX[2, 6] <- "90% CI"
  TableX[3, 2:7] <- c("Reporting Group", colnames(EstimatesStats[[mix]]))
  TableX[4:11, 1] <- 1:8
  TableX[4:11, 2] <- rownames(EstimatesStats[[mix]])
  TableX[4:11, 3:7] <- formatC(x = EstimatesStats[[mix]], digits = 3, format = "f")
  
  write.xlsx(x = TableX, file = "Estimates tables/Troll2016_8RG_StratifiedEstimatesStats_FormattedPretty.xlsx",
             col.names = FALSE, row.names = FALSE, sheetName = paste(mix, " Troll 8 Driver"), append = TRUE)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create 2016 HeatMaps ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("Figures")

# Create layout
layoutmat <- matrix(c(9,1,2,11,
                      9,3,4,11,
                      9,5,6,11,
                      9,7,8,11,
                      12,10,10,13), ncol=4,nrow=5,byrow=T)
SEAKTrollLayout <- layout(layoutmat,widths=c(0.25,1,1,0.25),heights=c(1,1,1,1,0.25))
layout.show(SEAKTrollLayout)

# Set color ramp
library('lattice')
WhiteRedColPalette <- colorRampPalette(colors=c("white","red"))
WhiteRedcol <- level.colors(x=seq(from=0,to=1,by=0.01), at = seq(from=0,to=1,by=0.01), col.regions = WhiteRedColPalette(100))

# Mixture names
mixnames <- names(EstimatesStats)[-6]

# Create list object with by RG stock comps
HeatmapEstimates <- sapply(GroupNames8Pub, function(RG) {
  matrix(data = sapply(mixnames, function(mix) {EstimatesStats[[mix]][RG, "Mean"] }),
         nrow = 2, ncol = 5, dimnames = list(c("NO", "AllQuad"), c("EWint", "LWint", "Spring", "SumRet1", "SumRet2"))
  )
}, simplify = FALSE)
zmax <- max(sapply(HeatmapEstimates, max))

Testing <- matrix(c(seq(from = 0, to = zmax, length.out = 102), seq(from = 0, to = zmax, length.out = 102)), nrow = 2, ncol = 102, byrow = T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot: Can't do a nested layout, writing out as pdf then pasting in other pdf

# pdf("Figures/2016TrollByFisheryQuadrant.pdf", family = "Times", width = 6.5, height = 6.5, title = "2016 Troll By Fishery and Quadrant")
png("Figures/2016TrollByFisheryQuadrant.png", family = "Times", width = 6.5, height = 6.5, units = "in", res = 300)
# x11(width = 6.5, height = 6.5)
par(xaxt = "n", yaxt = "n", omi = rep(0.1, 4), mar = rep(0.1, 4), family = 'serif')
layout(layoutmat,widths=c(0.3,1,1,0.25),heights=c(1,1,1,1,0.4))

## Loop through Reporting Group plots
sapply(GroupNames8Pub, function(RG) {
  image(t(HeatmapEstimates[[RG]])[, c("AllQuad", "NO")], zlim = c(0, zmax), col = WhiteRedcol, xlab = "", ylab = "", breaks = seq(from = 0, to = zmax, length.out = 102), useRaster = TRUE)
  abline(h = 0.5, lwd = 2, col = 'grey')
  abline(v = c(0.135, 0.38, 0.63, 0.875), lwd= 2 , col = 'grey')
  abline(h = c(-0.5, 1.5), v = c(-0.125, 1.125),lwd = 5, col = 'black')
  text(labels = RG, cex = 2, adj = c(0, 0.5), x = -0.1, y = 1)
})

## Plot 10 - Y-axis label
plot.new()
text(labels = "Quadrant", cex = 3, srt = 90, x = 0.3, y = 0.5, adj = c(0.5, 0))
text(labels = "NO", cex = 2, x = 0.99, y = c(0.97, 0.7, 0.43, 0.16), adj = c(1, 0.5))
text(labels = "All", cex = 2, x = 0.99, y = c(0.97, 0.7, 0.43, 0.16) - 0.135, adj = c(1, 0.5))

## Plot 11 - X-axis label
plot.new()
text(labels = "Fishery", cex = 3, adj = c(0.5, 0.5), x = 0.5, y = 0.35)
text(labels = "EW", cex = 2, adj = c(0.5, 0.5), x = c(0.02, 0.56), y = 0.8)
text(labels = "LW", cex = 2, adj = c(0.5, 0.5), x = c(0.02 + 0.115, 0.56 + 0.115), y = 0.8)
text(labels = "SP", cex = 2, adj = c(0.5, 0.5), x = c(0.02 + 0.22, 0.56 + 0.22), y = 0.8)
text(labels = "SU1", cex = 2, adj = c(0.5, 0.5), x = c(0.02 + 0.33, 0.56 + 0.33), y = 0.8)
text(labels = "SU2", cex = 2, adj = c(0.5, 0.5), x = c(0.02 + 0.43, 0.56 + 0.43), y = 0.8)

## Plot 13 - Legend
image(Testing, col = WhiteRedcol, xlab = "", ylab = "", breaks = seq(from = 0, to = zmax, length.out = 102))
text(labels = "0%", cex = 2.8, adj = c(0.5, 0.5), x = 0.5, y = 0.03)
text(labels = "50%", cex = 2.8, adj = c(0.5, 0.5), x = 0.5, y = 0.98)
abline(h = c(-0.005,  1.005),  v  =  c(-0.5,  1.5), lwd = 5, col = 'black')
dev.off()
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Resummarize 2016 Comm Troll Stratified Estimates to 18RGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Define RGs

GroupVec33RG_to18RG <- GroupVec18[GroupVec33RG_to26RG]
dput(x = GroupVec33RG_to18RG, file = "Objects/GroupVec33RG_to18RG.txt")

GroupVec33RG_to4RG <- GroupVec4[GroupVec33RG_to26RG]
dput(x = GroupVec33RG_to4RG, file = "Objects/GroupVec33RG_to4RG.txt")

### Unstratified
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 18RGs
EWintTroll2016_18RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                               maindir = "BAYES/Output/33RG/EWint_2016", 
                               mixvec = WinterTrollMix2016[1:2], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

LWintTroll2016_18RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                               maindir = "BAYES/Output/33RG/LWint_2016", 
                               mixvec = WinterTrollMix2016[3:4], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SpringTroll2016_18RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                               maindir = "BAYES/Output/33RG/Spring_2016", 
                               mixvec = SpringTrollMix2016, prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SumRet1Troll2016_18RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                               maindir = "BAYES/Output/33RG/SumRet1_2016", 
                               mixvec = SummerTrollMix2016[1:4], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SumRet2Troll2016_18RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                               maindir = "BAYES/Output/33RG/SumRet2_2016", 
                               mixvec = SummerTrollMix2016[5:8], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

TrollByQuad2016_18RG_EstimatesStats <- c(EWintTroll2016_18RG_EstimatesStats,
                                        LWintTroll2016_18RG_EstimatesStats,
                                        SpringTroll2016_18RG_EstimatesStats,
                                        SumRet1Troll2016_18RG_EstimatesStats,
                                        SumRet2Troll2016_18RG_EstimatesStats)
any(sapply(TrollByQuad2016_18RG_EstimatesStats, function(mix) {any(mix[, "GR"] > 1.2)}))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dput files
# dir.create("Estimates objects")
# Grab estimates objects and dput in "Estimates objects"
objects2dput <- c("EWintTroll2016_18RG_EstimatesStats", 
                  "LWintTroll2016_18RG_EstimatesStats", 
                  "SpringTroll2016_18RG_EstimatesStats", 
                  "SumRet1Troll2016_18RG_EstimatesStats",
                  "SumRet2Troll2016_18RG_EstimatesStats",
                  "TrollByQuad2016_18RG_EstimatesStats")

invisible(sapply(objects2dput, function(obj) {
  dput(x = get(obj), file = paste0("Estimates objects/", obj, ".txt"))
})); rm(objects2dput)



### Stratified
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 18RGs
EWintTroll2016_18RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                          maindir="BAYES/Output/33RG/EWint_2016",
                          mixvec = c("EWintNISISO_2016", "EWintNO_2016"), catchvec = c(4216, 25147), 
                          newname = "StratifiedEWint2016_90percentCI_18RG", priorname = "", nchains = 5)

LWintTroll2016_18RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                          maindir="BAYES/Output/33RG/LWint_2016",
                          mixvec = c("LWintNISISO_2016", "LWintNO_2016"), catchvec = c(5248, 17680), 
                          newname = "StratifiedLWint2016_90percentCI_18RG", priorname = "", nchains = 5)

SpringTroll2016_18RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                          maindir="BAYES/Output/33RG/Spring_2016",
                          mixvec = c("SpringNI_2016", "SpringNO_2016", "SpringSI_2016", "SpringSO_2016"), catchvec = c(9270, 17012, 15160, 1031), 
                          newname = "StratifiedSpring2016_90percentCI_18RG", priorname = "", nchains = 5)

SumRet1Troll2016_18RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                          maindir="BAYES/Output/33RG/SumRet1_2016",
                          mixvec = c("Summer1NI_2016", "Summer1NO_2016", "Summer1SI_2016", "Summer1SO_2016"), catchvec = c(3805, 80323, 3618, 18888), 
                          newname = "StratifiedSumRet12016_90percentCI_18RG", priorname = "", nchains = 5)

SumRet2Troll2016_18RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                          maindir="BAYES/Output/33RG/SumRet2_2016",
                          mixvec = c("Summer2NI_2016", "Summer2NO_2016", "Summer2SI_2016", "Summer2SO_2016"), catchvec = c(2147, 56208, 1774, 14111), 
                          newname = "StratifiedSumRet22016_90percentCI_18RG", priorname = "", nchains = 5)


AllYearTroll2016_18RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to18RG, groupnames = GroupNames18,
                          maindir="BAYES/Output/33RG/AllYearTroll_2016",
                          mixvec = c("EWintNISISO_2016", "EWintNO_2016", "LWintNISISO_2016", "LWintNO_2016", 
                                     "SpringNI_2016", "SpringNO_2016", "SpringSI_2016", "SpringSO_2016", 
                                     "Summer1NI_2016", "Summer1NO_2016", "Summer1SI_2016", "Summer1SO_2016", 
                                     "Summer2NI_2016", "Summer2NO_2016", "Summer2SI_2016", "Summer2SO_2016"),
                          catchvec = c(4216, 25147, 5248, 17680,
                                       9270, 17012, 15160, 1031,
                                       3805, 80323, 3618, 18888,
                                       2147, 56208, 1774, 14111), 
                          newname = "StratifiedAllYearTroll2016_90percentCI_18RG", priorname = "", nchains = 5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dput files
# Grab estimates objects and dput in "Estimates objects"
objects2dput <- c("EWintTroll2016_18RG_StratifiedEstimatesStats", 
                  "LWintTroll2016_18RG_StratifiedEstimatesStats", 
                  "SpringTroll2016_18RG_StratifiedEstimatesStats", 
                  "SumRet1Troll2016_18RG_StratifiedEstimatesStats",
                  "SumRet2Troll2016_18RG_StratifiedEstimatesStats", 
                  "AllYearTroll2016_18RG_StratifiedEstimatesStats")

invisible(sapply(objects2dput, function(obj) {
  dput(x = get(obj)$Stats, file = paste0("Estimates objects/", obj, ".txt"))
})); rm(objects2dput)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create 2016 summary tables for Comm Troll 18 RGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get objects
SEAK16estimatesobjects <- list.files(path = "Estimates objects", recursive = FALSE, pattern = "_18RG")
SEAK16estimatesobjects

invisible(sapply(SEAK16estimatesobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Estimates objects", objct, sep = "/")), pos = 1) })); beep(2)


# Dget all estimates stats
SEAK16estimatesobjects <- unlist(lapply(SEAK16estimatesobjects, function(objct) {unlist(strsplit(x = objct, split = ".txt"))}))

Troll2016_18RG_EstimatesStats <- list(
  "EWintNO_2016" = EWintTroll2016_18RG_EstimatesStats[["EWintNO_2016"]],
  "EWintAllQuad_2016" = EWintTroll2016_18RG_StratifiedEstimatesStats,
  "LWintNO_2016" = LWintTroll2016_18RG_EstimatesStats[["LWintNO_2016"]],
  "LWintAllQuad_2016" = LWintTroll2016_18RG_StratifiedEstimatesStats,
  "SpringNO_2016" = SpringTroll2016_18RG_EstimatesStats[["SpringNO_2016"]],
  "SpringSI_2016" = SpringTroll2016_18RG_EstimatesStats[["SpringSI_2016"]],
  "SpringAllQuad_2016" = SpringTroll2016_18RG_StratifiedEstimatesStats,
  "SumRet1NO_2016" = SumRet1Troll2016_18RG_EstimatesStats[["Summer1NO_2016"]],
  "SumRet1AllQuad_2016" = SumRet1Troll2016_18RG_StratifiedEstimatesStats,
  "SumRet2NO_2016" = SumRet2Troll2016_18RG_EstimatesStats[["Summer2NO_2016"]],
  "SumRet2AllQuad_2016" = SumRet2Troll2016_18RG_StratifiedEstimatesStats
)
dput(x = Troll2016_18RG_EstimatesStats, file = "Estimates objects/Troll2016_18RG_EstimatesStats.txt")

# Check GR
any(sapply(Troll2016_18RG_EstimatesStats, function(mix) {any(mix[, "GR"] > 1.2)}))

# Publishable names
GroupNames18Pub <- c(GroupNames26Pub[c(3:7, 9:12, 16:18, 20:21, 23:25)], "Other")
dput(x = GroupNames18Pub, file = "Objects/GroupNames18Pub.txt")


# Reformat estimates stats
Troll2016_18RG_EstimatesStats_Formatted <- sapply(Troll2016_18RG_EstimatesStats, function(yr) {
  matrix(data = yr[, 1:5], nrow = 18, ncol = 5, dimnames = list(GroupNames18Pub, c("Mean", "SD", "Median", "5%", "95%")))
}, simplify = FALSE)


SEAK2016Mixtures <- list.files(path = "BAYES/Mixture", full.names = FALSE)
SEAK2016Mixtures_SampSizes <- sapply(SEAK2016Mixtures, function(mix) {dim(read.table(file = paste0("BAYES/Mixture/", mix)))[1]} )
names(SEAK2016Mixtures_SampSizes) <- sapply(names(SEAK2016Mixtures_SampSizes), function(mix) {unlist(strsplit(x = mix, split = ".mix"))[1]})


Troll2016_SampleSizes <- sapply(Troll2016MixNames, function(mix) {sum(SEAK2016Mixtures_SampSizes[mix])} )

# Create fully formatted spreadsheat
EstimatesStats <- Troll2016_18RG_EstimatesStats_Formatted
SampSizes <- Troll2016_SampleSizes

# dir.create("Estimates tables")

for(mix in names(EstimatesStats)) {
  
  TableX <- matrix(data = "", nrow = 21, ncol = 7)
  TableX[1, 3] <- paste(Troll2016PubNames[mix], "(n=", SampSizes[mix], ")")
  TableX[2, 6] <- "90% CI"
  TableX[3, 2:7] <- c("Reporting Group", colnames(EstimatesStats[[mix]]))
  TableX[4:21, 1] <- 1:18
  TableX[4:21, 2] <- rownames(EstimatesStats[[mix]])
  TableX[4:21, 3:7] <- formatC(x = EstimatesStats[[mix]], digits = 3, format = "f")
  
  write.xlsx(x = TableX, file = "Estimates tables/Troll2016_18RG_StratifiedEstimatesStats_FormattedPretty.xlsx",
             col.names = FALSE, row.names = FALSE, sheetName = paste(mix, " Troll 18RG"), append = TRUE)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Resummarize 2016 Comm Troll Stratified Estimates to 4RGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Define RGs

GroupVec33RG_to4RG <- GroupVec4[GroupVec33RG_to26RG]
dput(x = GroupVec33RG_to4RG, file = "Objects/GroupVec33RG_to4RG.txt")

### Unstratified
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 4RGs
EWintTroll2016_4RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                               maindir = "BAYES/Output/33RG/EWint_2016", 
                               mixvec = WinterTrollMix2016[1:2], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

LWintTroll2016_4RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                               maindir = "BAYES/Output/33RG/LWint_2016", 
                               mixvec = WinterTrollMix2016[3:4], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SpringTroll2016_4RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                               maindir = "BAYES/Output/33RG/Spring_2016", 
                               mixvec = SpringTrollMix2016, prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SumRet1Troll2016_4RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                               maindir = "BAYES/Output/33RG/SumRet1_2016", 
                               mixvec = SummerTrollMix2016[1:4], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

SumRet2Troll2016_4RG_EstimatesStats <- 
  CustomCombineBAYESOutput.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                               maindir = "BAYES/Output/33RG/SumRet2_2016", 
                               mixvec = SummerTrollMix2016[5:8], prior = "",
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

TrollByQuad2016_4RG_EstimatesStats <- c(EWintTroll2016_4RG_EstimatesStats,
                                         LWintTroll2016_4RG_EstimatesStats,
                                         SpringTroll2016_4RG_EstimatesStats,
                                         SumRet1Troll2016_4RG_EstimatesStats,
                                         SumRet2Troll2016_4RG_EstimatesStats)
any(sapply(TrollByQuad2016_4RG_EstimatesStats, function(mix) {any(mix[, "GR"] > 1.2)}))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dput files
# dir.create("Estimates objects")
# Grab estimates objects and dput in "Estimates objects"
objects2dput <- c("EWintTroll2016_4RG_EstimatesStats", 
                  "LWintTroll2016_4RG_EstimatesStats", 
                  "SpringTroll2016_4RG_EstimatesStats", 
                  "SumRet1Troll2016_4RG_EstimatesStats",
                  "SumRet2Troll2016_4RG_EstimatesStats",
                  "TrollByQuad2016_4RG_EstimatesStats")

invisible(sapply(objects2dput, function(obj) {
  dput(x = get(obj), file = paste0("Estimates objects/", obj, ".txt"))
})); rm(objects2dput)



### Stratified
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4RGs
EWintTroll2016_4RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                          maindir="BAYES/Output/33RG/EWint_2016",
                          mixvec = c("EWintNISISO_2016", "EWintNO_2016"), catchvec = c(4216, 25147), 
                          newname = "StratifiedEWint2016_90percentCI_4RG", priorname = "", nchains = 5)

LWintTroll2016_4RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                          maindir="BAYES/Output/33RG/LWint_2016",
                          mixvec = c("LWintNISISO_2016", "LWintNO_2016"), catchvec = c(5248, 17680), 
                          newname = "StratifiedLWint2016_90percentCI_4RG", priorname = "", nchains = 5)

SpringTroll2016_4RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                          maindir="BAYES/Output/33RG/Spring_2016",
                          mixvec = c("SpringNI_2016", "SpringNO_2016", "SpringSI_2016", "SpringSO_2016"), catchvec = c(9270, 17012, 15160, 1031), 
                          newname = "StratifiedSpring2016_90percentCI_4RG", priorname = "", nchains = 5)

SumRet1Troll2016_4RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                          maindir="BAYES/Output/33RG/SumRet1_2016",
                          mixvec = c("Summer1NI_2016", "Summer1NO_2016", "Summer1SI_2016", "Summer1SO_2016"), catchvec = c(3805, 80323, 3618, 18888), 
                          newname = "StratifiedSumRet12016_90percentCI_4RG", priorname = "", nchains = 5)

SumRet2Troll2016_4RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                          maindir="BAYES/Output/33RG/SumRet2_2016",
                          mixvec = c("Summer2NI_2016", "Summer2NO_2016", "Summer2SI_2016", "Summer2SO_2016"), catchvec = c(2147, 56208, 1774, 14111), 
                          newname = "StratifiedSumRet22016_90percentCI_4RG", priorname = "", nchains = 5)


AllYearTroll2016_4RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec33RG_to4RG, groupnames = GroupNames4,
                          maindir="BAYES/Output/33RG/AllYearTroll_2016",
                          mixvec = c("EWintNISISO_2016", "EWintNO_2016", "LWintNISISO_2016", "LWintNO_2016", 
                                     "SpringNI_2016", "SpringNO_2016", "SpringSI_2016", "SpringSO_2016", 
                                     "Summer1NI_2016", "Summer1NO_2016", "Summer1SI_2016", "Summer1SO_2016", 
                                     "Summer2NI_2016", "Summer2NO_2016", "Summer2SI_2016", "Summer2SO_2016"),
                          catchvec = c(4216, 25147, 5248, 17680,
                                       9270, 17012, 15160, 1031,
                                       3805, 80323, 3618, 18888,
                                       2147, 56208, 1774, 14111), 
                          newname = "StratifiedAllYearTroll2016_90percentCI_4RG", priorname = "", nchains = 5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dput files
# Grab estimates objects and dput in "Estimates objects"
objects2dput <- c("EWintTroll2016_4RG_StratifiedEstimatesStats", 
                  "LWintTroll2016_4RG_StratifiedEstimatesStats", 
                  "SpringTroll2016_4RG_StratifiedEstimatesStats", 
                  "SumRet1Troll2016_4RG_StratifiedEstimatesStats",
                  "SumRet2Troll2016_4RG_StratifiedEstimatesStats", 
                  "AllYearTroll2016_4RG_StratifiedEstimatesStats")

invisible(sapply(objects2dput, function(obj) {
  dput(x = get(obj)$Stats, file = paste0("Estimates objects/", obj, ".txt"))
})); rm(objects2dput)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create 2016 summary tables for Comm Troll 4 RGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get objects
SEAK16estimatesobjects <- list.files(path = "Estimates objects", recursive = FALSE, pattern = "_4RG")
SEAK16estimatesobjects

invisible(sapply(SEAK16estimatesobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Estimates objects", objct, sep = "/")), pos = 1) })); beep(2)


# Dget all estimates stats
SEAK16estimatesobjects <- unlist(lapply(SEAK16estimatesobjects, function(objct) {unlist(strsplit(x = objct, split = ".txt"))}))

Troll2016_4RG_EstimatesStats <- list(
  "EWintNO_2016" = EWintTroll2016_4RG_EstimatesStats[["EWintNO_2016"]],
  "EWintAllQuad_2016" = EWintTroll2016_4RG_StratifiedEstimatesStats,
  "LWintNO_2016" = LWintTroll2016_4RG_EstimatesStats[["LWintNO_2016"]],
  "LWintAllQuad_2016" = LWintTroll2016_4RG_StratifiedEstimatesStats,
  "SpringNI_2016" = SpringTroll2016_4RG_EstimatesStats[["SpringNI_2016"]],
  "SpringNO_2016" = SpringTroll2016_4RG_EstimatesStats[["SpringNO_2016"]],
  "SpringSI_2016" = SpringTroll2016_4RG_EstimatesStats[["SpringSI_2016"]],
  "SpringSO_2016" = SpringTroll2016_4RG_EstimatesStats[["SpringSO_2016"]],
  "SpringAllQuad_2016" = SpringTroll2016_4RG_StratifiedEstimatesStats,
  "SumRet1NO_2016" = SumRet1Troll2016_4RG_EstimatesStats[["Summer1NO_2016"]],
  "SumRet1AllQuad_2016" = SumRet1Troll2016_4RG_StratifiedEstimatesStats,
  "SumRet2NO_2016" = SumRet2Troll2016_4RG_EstimatesStats[["Summer2NO_2016"]],
  "SumRet2AllQuad_2016" = SumRet2Troll2016_4RG_StratifiedEstimatesStats
)
dput(x = Troll2016_4RG_EstimatesStats, file = "Estimates objects/Troll2016_4RG_EstimatesStats.txt")

# Check GR
any(sapply(Troll2016_4RG_EstimatesStats, function(mix) {any(mix[, "GR"] > 1.2)}))

# Publishable names
GroupNames4Pub <- c(GroupNames4[1:3], "US South")
dput(x = GroupNames4Pub, file = "Objects/GroupNames4Pub.txt")


# Reformat estimates stats
Troll2016_4RG_EstimatesStats_Formatted <- sapply(Troll2016_4RG_EstimatesStats, function(yr) {
  matrix(data = yr[, 1:5], nrow = 4, ncol = 5, dimnames = list(GroupNames4Pub, c("Mean", "SD", "Median", "5%", "95%")))
}, simplify = FALSE)


SEAK2016Mixtures <- list.files(path = "BAYES/Mixture", full.names = FALSE)
SEAK2016Mixtures_SampSizes <- sapply(SEAK2016Mixtures, function(mix) {dim(read.table(file = paste0("BAYES/Mixture/", mix)))[1]} )
names(SEAK2016Mixtures_SampSizes) <- sapply(names(SEAK2016Mixtures_SampSizes), function(mix) {unlist(strsplit(x = mix, split = ".mix"))[1]})

Troll2016MixNames_4RG <- setNames(object = list("EWintNO_2016",
                                                WinterTrollMix2016[1:2],
                                                "LWintNO_2016",
                                                WinterTrollMix2016[3:4],
                                                "SpringNI_2016",
                                                "SpringNO_2016",
                                                "SpringSI_2016",
                                                "SpringSO_2016",
                                                SpringTrollMix2016,
                                                "Summer1NO_2016",
                                                SummerTrollMix2016[1:4],
                                                "Summer2NO_2016",
                                                SummerTrollMix2016[5:8]),
                                  nm = names(Troll2016_4RG_EstimatesStats_Formatted))
dput(x = Troll2016MixNames_4RG, file = "Objects/Troll2016MixNames_4RG.txt")

Troll2016_SampleSizes <- sapply(Troll2016MixNames_4RG, function(mix) {sum(SEAK2016Mixtures_SampSizes[mix])} )

Troll2016PubNames_4RG <- setNames(object = c("Northern Outside Quadrant",
                                             "All Quadrants",
                                             "Northern Outside Quadrant",
                                             "All Quadrants",
                                             "Northern Inside Quadrant",
                                             "Northern Outside Quadrant",
                                             "Southern Inside Quadrant",
                                             "Southern Outside Quadrant",
                                             "All Quadrants",
                                             "Northern Outside Quadrant",
                                             "All Quadrants",
                                             "Northern Outside Quadrant",
                                             "All Quadrants"), 
                                  nm = names(Troll2016_4RG_EstimatesStats_Formatted))
dput(x = Troll2016PubNames_4RG, file = "Objects/Troll2016PubNames_4RG.txt")


# Create fully formatted spreadsheat
EstimatesStats <- Troll2016_4RG_EstimatesStats_Formatted
SampSizes <- Troll2016_SampleSizes

# dir.create("Estimates tables")

for(mix in names(EstimatesStats)) {
  
  TableX <- matrix(data = "", nrow = 7, ncol = 7)
  TableX[1, 3] <- paste(Troll2016PubNames_4RG[mix], "(n=", SampSizes[mix], ")")
  TableX[2, 6] <- "90% CI"
  TableX[3, 2:7] <- c("Reporting Group", colnames(EstimatesStats[[mix]]))
  TableX[4:7, 1] <- seq(dim(EstimatesStats[[1]])[1])
  TableX[4:7, 2] <- rownames(EstimatesStats[[mix]])
  TableX[4:7, 3:7] <- formatC(x = EstimatesStats[[mix]], digits = 3, format = "f")
  
  write.xlsx(x = TableX, file = "Estimates tables/Troll2016_4RG_StratifiedEstimatesStats_FormattedPretty.xlsx",
             col.names = FALSE, row.names = FALSE, sheetName = paste(mix, " Troll 4RG"), append = TRUE)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summer Retention 2 8RG Driver Stock Resummarization 2009-2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pulled a fresh harvest report from MTA lab website to get "up to date" harvest numbers
# Copied all relevant BAYES files into one location

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2016
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK16")

SumRet2Troll2016_8RG_StratifiedEstimatesStats <- dget(file = "Estimates objects/SumRet2Troll2016_8RG_StratifiedEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015
# No 2nd retention

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK14")

# Need a new groupvec of length 27, not 26, because Chilkat was broken out of NSEAK

SumRet2Troll2014_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = c(1, GroupVec8), groupnames = GroupNames8,
                          maindir="BAYES/Output/AllYearTroll_2014",
                          mixvec = c("SumRet2NISISO.2014", "SumRet2NO.2014"),
                          catchvec = c(24365, 31288), 
                          newname = "StratifiedSumRet2Troll2014_90percentCI_8RG", priorname = "", nchains = 5)

dput(x = SumRet2Troll2014_8RG_StratifiedEstimatesStats$Stats, file = "Estimates objects/SumRet2Troll2014_8RG_StratifiedEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2013
# No 2nd retention

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2012
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK12")

SumRet2Troll2012_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec8, groupnames = GroupNames8,
                          maindir="BAYES/Output/AllYearTroll_2012",
                          mixvec = c("SumRet2NISISO.2012", "SumRet2NO.2012"),
                          catchvec = c(20056, 53914), 
                          newname = "StratifiedSumRet2Troll2012_90percentCI_8RG", priorname = "", nchains = 5)

dput(x = SumRet2Troll2012_8RG_StratifiedEstimatesStats$Stats, file = "Estimates objects/SumRet2Troll2012_8RG_StratifiedEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2011
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK11")

SumRet2Troll2011_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec8, groupnames = GroupNames8,
                          maindir="BAYES/Troll/Output/AllYearTroll_2011",
                          mixvec = c("SumRet2NISISO.2011", "SumRet2NO.2011"),
                          catchvec = c(13364, 16372), 
                          newname = "StratifiedSumRet2Troll2011_90percentCI_8RG", priorname = "", nchains = 5)

dput(x = SumRet2Troll2011_8RG_StratifiedEstimatesStats$Stats, file = "Estimates objects/SumRet2Troll2011_8RG_StratifiedEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2010
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK10")

SumRet2Troll2010_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec8, groupnames = GroupNames8,
                          maindir="BAYES/Troll/Output/AllYearTroll_2010",
                          mixvec = c("SummerR2NISISO.2010", "SummerR2NO.2010"),
                          catchvec = c(22025, 26430), 
                          newname = "StratifiedSumRet2Troll2010_90percentCI_8RG", priorname = "", nchains = 5)

dput(x = SumRet2Troll2010_8RG_StratifiedEstimatesStats$Stats, file = "Estimates objects/SumRet2Troll2010_8RG_StratifiedEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2009
setwd("V:/Analysis/1_SEAK/Chinook/Mixture/SEAK09")

SumRet2Troll2009_8RG_StratifiedEstimatesStats <- 
  StratifiedEstimator.GCL(groupvec = GroupVec8, groupnames = GroupNames8,
                          maindir="BAYES/Troll/Output/AllYearTroll_2009",
                          mixvec = c("SumRet2NISI_2009", "SumRet2NOSO_2009"),
                          catchvec = c(1796, 31216), 
                          newname = "StratifiedSumRet2Troll2009_90percentCI_8RG", priorname = "", nchains = 5)

dput(x = SumRet2Troll2009_8RG_StratifiedEstimatesStats$Stats, file = "Estimates objects/SumRet2Troll2009_8RG_StratifiedEstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cbind("AY2009" = SumRet2Troll2009_8RG_StratifiedEstimatesStats$Stats[, "mean"],
      "AY2010" = SumRet2Troll2010_8RG_StratifiedEstimatesStats$Stats[, "mean"],
      "AY2011" = SumRet2Troll2011_8RG_StratifiedEstimatesStats$Stats[, "mean"],
      "AY2012" = SumRet2Troll2012_8RG_StratifiedEstimatesStats$Stats[, "mean"],
      "AY2013" = rep(NA, 8),
      "AY2014" = SumRet2Troll2014_8RG_StratifiedEstimatesStats$Stats[, "mean"],
      "AY2015" = rep(NA, 8),
      "AY2016" = SumRet2Troll2016_8RG_StratifiedEstimatesStats[, "mean"])