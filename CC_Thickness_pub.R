####################################
#Calculate Corpus callosum (CC) Thickness on the midsagittal slice 
#and provide CC division based on Witelson 1989
#Version 1, 10/16/2019
#
#Author: Jun Yi Wang  email: jyiwang@ucdavis.edu
#
#Citation: Wang JY, Hessl D, Tassone F, Kim K, Hagerman R, Rivera SM. (2019). 
#Interaction between ventricular expansion and structural changes in the corpus  
#callosum and putamen in males with FMR1 normal and premutation alleles. 
#Neurobiol Aging doi.org/10.1016/j.neurobiolaging.2019.09.009.
#
#Input:
#csv files containing cordinates of CC masks in 3 columns: x, y, z
#One csv file for each scan with file name: id_CC_xyz
#Make sure x containing sagittal slice number, y coronal slice number, and z axial slice number 
#Note: replace id with the actual IDs
#
#CC coordinates can be generated using the matlab function: ind2sub
#CC segmentation needs to be checked for smoothness. Recommend to use fslmaths -fmedian for smoothing
#
#A midsagSliceNumber.csv file containing 5 columns: id, voxel size at x-, y-, and z-axes, 
#and midsagital slice number with the following column names: id, r.x, r.y, r.z, midsag
#
#examine the output from the console
#for scans with max. voxel distance >=2, check smoothness of the segmentation 
#and make the boundary smooth manually using ITK-snap if necessary
#
#output files
#dtSurface: assignment of CC outline to inner and outer curve with the following columns
##subID: unique for each subject and each visit/scan
##sag: sagittal coordinate
##axial: axial coordinate
##cor: coronal coordinate
##type: outline type--up, low
##min.cor: smallest coronal coordinate
##max.cor: maximum coronal coordinate
#
#dtNewSegments
##subID: unique for each subject and each visit/scan
##sag: sagittal coordinate
##axial: axial coordinate
##cor: coronal coordinate
##type: outline type--up, low, and mid
##cc.type: same as type for "low" and "up" outline, but is assigned from CC1-CC7 according to 
##         Witelson 1989 for the "mid" type
##id: new CC segment ID from 1-301, equally spaced on the lower curve
##n.section: number of CC sections determined by inflection points on the lower curve
##thickness: length of the line linking the correponding points on the lower and upper curve
#
####################################
library(data.table)
library(ggplot2)

require(graphics)
library(RColorBrewer)
library(colorspace)

rm(list = ls())
source('cc_thickness_functions.r')

#read all files in directory and process one by one
#replace the following directory with your own
filePath <- '../LimbicTrax/CC_coordinates/'
ccFileNames <- dir(filePath,pattern='_CC_')

#read csv file containing midsagittal slice number for CC and voxel size
#replace the following directory with your own
dtMidSag <- data.table(read.csv('../LimbicTrax/TraxMidsagSliceNumber.csv'))
setkey(dtMidSag,id)

#number of segments of CC
n.seg <- 300
outputCC <- getCCSurfaceThickness(ccFileNames, dtMidSag, n.seg)
#for scans with max. voxel distance >=2, check smoothness of the segmentation 
#and make the boundary smooth manually using ITK-snap if necessary before continuing
#you may also reference the plots generated below to find errors

dtSurface <- outputCC[[1]]
dtNewSegments <- outputCC[[2]]
dtThickness <- outputCC[[3]]

rm(outputCC)

write.csv(dtSurface,'cc_surface.csv')
write.csv(dtNewSegments,'cc_newsegments.csv')
write.csv(dtThickness,'cc_thickness.csv')

####################################
#Plots
#Be patient. It may take several seconds to display many plots at once
#If the exported images do not have good resolution, use the Zoom function in RStudio, 
#screen-shot, paste into paint, and save
####################################
nSub <- nrow(dtThickness)
dtThickness[,n.sub:=seq(4, nSub+3)]
dtThickness[,colNum:=n.sub%/%4]
dtThickness[,rowNum:=seq(1,length(subID)),by=colNum]
tmp <- dtThickness[,list(subID,n.sub,colNum,rowNum)]
dtSurface <- merge(dtSurface,tmp,by="subID")
dtNewSegments <- merge(dtNewSegments,tmp,by="subID")

#display plots in one row with subject #, 10 at a time
#change the value of i to display the next set
i <- 1
ggplot(dtSurface[(n.sub-4)%/%10 + 1 == i], aes(x=cor, y=axial, color=factor(type))) +  
  geom_point(size=1) + facet_grid(subID ~ ., scales="free", space="free") + 
  xlab("Coronal") + ylab("Axial")

ggplot(dtNewSegments[(n.sub-4)%/%10 + 1 == i], aes(x=cor, y=axial)) + 
  geom_line(aes(group=factor(id)),color="gray") + geom_point(aes(color=factor(cc.type)),size=0.5) + 
  scale_color_manual(values=c("red","cyan","orange","blue","yellow","purple","green","gray","gray")) + 
  facet_grid(subID ~ ., scales="free", space="free") + xlab("Coronal") + ylab("Axial")

#display plots in 4 rows but without subject #
ggplot(dtSurface, aes(x=cor, y=axial, color=factor(type))) + geom_point(size=0.5) + 
  facet_grid(vars(colNum),vars(rowNum)) + xlab("Coronal") + ylab("Axial")

ggplot(dtNewSegments, aes(x=cor, y=axial)) + 
  geom_line(aes(group=factor(id)),color="gray") + geom_point(aes(color=factor(cc.type)),size=0.5) + 
  scale_color_manual(values=c("red","cyan","orange","blue","yellow","purple","green","gray","gray")) + 
  facet_grid(vars(colNum),vars(rowNum)) + xlab("Coronal") + ylab("Axial")

#for longitudinal data, display each subject with multiple visits in a row
#add columns of visits and subject (same for all visits)
demor <- data.table(read.csv('../LimbicTrax/TraxVisits.csv'))
setkey(demor,subID)
dtSurface <- merge(dtSurface,demor,by="subID")
dtNewSegments <- merge(dtNewSegments,demor,by="subID")

#display the assignment of upper and lower boundaries
ggplot(dtSurface, aes(x=cor, y=axial, color=factor(type))) + geom_point(size=0.5) + 
  facet_grid(subject ~ visits, scales="free", space="free") + xlab("Coronal") + ylab("Axial")

#display the CC segmentation and the assignment of 7 CC sections
ggplot(dtNewSegments,aes(x=cor, y=axial)) + 
  geom_line(aes(group=factor(id)),color="gray") + geom_point(aes(color=factor(cc.type)),size=0.5) + 
  scale_color_manual(values=c("red","cyan","orange","blue","yellow","purple","green","gray","gray")) + 
  facet_grid(subject ~ visits, scales="free", space="free") + xlab("Coronal") + ylab("Axial")

#for a particular scan
ggplot(dtNewSegments[subID=="S0038_1"],aes(x=cor, y=axial)) + 
  geom_line(aes(group=factor(id)),color="gray") + geom_point(aes(color=factor(cc.type)),size=1) + 
  scale_color_manual(values=c("red","cyan","orange","blue","yellow","purple","green","gray","gray")) + 
  xlab("Coronal") + ylab("Axial")

####################################
#calculate average thickness for every 5 sections
####################################
n.seg <- 300
nPoints <- 60
segments60 <- dtNewSegments[type=='mid' &id<n]
segments60[,':='(sec=(id-1) %/% (n.seg/nPoints) + 1)]

