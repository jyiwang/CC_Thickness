getCCSurfaceThickness <- function(ccFileNames, dtMidSag, n.seg) {
  dtSurface <- data.table()
  dtSegSurface <- data.table()
  dtSegments <- data.table()
  dtNewSegments <- data.table()
  
  for (i in 1:length(ccFileNames)) {
    subID <- unlist(strsplit(ccFileNames[i],'_CC'))[1]
    print(i)
    print(paste("Processing",subID))
    #cc
    mask <- data.table(read.csv(paste(filePath,ccFileNames[i],sep=''),header=F,
                                col.names=c('sag','cor','axial')))
    setkey(mask,sag,cor,axial)
    #get inner and outer surfaces
    midSag <- dtMidSag[id==subID]$midsag
    surface <- getCCSurface(mask,midSag)
    midSagOne <- dtMidSag[id==subID]
    segSurface <- orderCCVoxels(surface,midSagOne)
    segSurface <- getCumulativeDistance(segSurface,midSagOne,midSag)
    #adjust order for those with distance >=2 and involving more than 2 voxels
    if (nrow(segSurface[d.voxel>=2])>1) {
      segSurface <- adjustOrder(segSurface)
      segSurface <- getCumulativeDistance(segSurface,midSagOne,midSag)
    }
    
    #segment the curves into equal number of segments for each sagittal slice
    segments <- segmentCCSurface(segSurface,midSagOne,n.seg)
    #find closest points on the lower surface for each point on the upper surface
    #for midsagittal slice
    closePoints <- getClosestPoints(segments,midSagOne)
    #section the CC by finding deflective points in the lower boundary: 
    #the voxels having shortest distance to >1 upper voxels
    #and the corresponding upper voxels that are in the middle
    min.connect <- 5
    newSegments <- getCCSections(closePoints,segments,segSurface,n.seg,min.connect,midSag)
    #calculate CC thickness by linking points on the upper and lower boundaries
    ccThickness <- calCCThickness(newSegments,midSagOne)
    #find midline and assign CC anatomic divisions according to Witelson 1989
    newSegments <- getCCDivisions(newSegments,ccThickness)
    
    #add to the big table
    tmp <- data.table(subID=subID,surface)
    dtSurface = rbind(dtSurface,tmp)
    tmp <- data.table(subID=subID,segSurface)
    dtSegSurface = rbind(dtSegSurface,tmp)
    tmp <- data.table(subID=subID,segments)
    dtSegments = rbind(dtSegments,tmp)
    tmp <- data.table(subID=subID,newSegments)
    dtNewSegments = rbind(dtNewSegments,tmp)
  } 
  
  #Generate CC surface and calculate thickness by sectionalize inner and outer curves
  #thickness for each section and whole
  #thickness in sections 1 and 61 may not be available
  dtNewSegments[id==1 &is.na(thickness),thickness:=0]
  secThickness <- dtNewSegments[type=='mid' &id<n,
                                list(m.thick=round(mean(thickness),3)), by=list(subID,cc.type)]
  setkey(secThickness,subID,cc.type)
  #check whether all CC are assigned with 7 sections
  secThickness[,n.sec:=length(cc.type),by=subID]
  if (nrow(secThickness[n.sec<7])>0) {
    print("following scans have CC sections less than 7")
    print(secThickness[n.sec<7])
  }
  
  #thickness for the whole CC
  ccThickness <- dtNewSegments[id<n, list(m.thick=round(mean(thickness),3)),by=list(subID)]
  setkey(ccThickness,subID)
  dtThickness <- data.table(
    subID=ccThickness$subID, t.cc=ccThickness$m.thick, t.CC1=secThickness[cc.type=='CC1']$m.thick,
    t.CC2=secThickness[cc.type=='CC2']$m.thick, t.CC3=secThickness[cc.type=='CC3']$m.thick, 
    t.CC4=secThickness[cc.type=='CC4']$m.thick, t.CC5=secThickness[cc.type=='CC5']$m.thick,
    t.CC6=secThickness[cc.type=='CC6']$m.thick, t.CC7=secThickness[cc.type=='CC7']$m.thick)
  setkey(dtThickness,subID)

  return(list(dtSurface, dtNewSegments, dtThickness))
}

#get surface coordinates for the coronal slices having one sections for the splenium half
getCCsSurfaceOneSec <- function(tmpMidMask,midCor) {
  #find lowest and most anterior tip in splenium
  minAxial <- tmpMidMask[cor<midCor,list(min.axial=min(axial)),by=sag]
  tmpMidMask <- merge(tmpMidMask,minAxial,by='sag',all.x=T)
  maxTCor <- tmpMidMask[cor<midCor &axial==min.axial &count==1,list(max.t.cor=max(cor)),by=sag]
  tmpMidMask <- merge(tmpMidMask,maxTCor,by='sag',all.x=T)
  tmpSurf <- tmpMidMask[cor>=max.t.cor & count==1,list(sag=sag,cor=cor,axial=axial,type='low')]
  tmp <- tmpMidMask[cor<max.t.cor & count==1,list(sag=sag,cor=cor,axial=axial,type='up')]
  tmpSurf <- rbind(tmpSurf,tmp)
  tmp <- tmpMidMask[count==n,list(sag=sag,cor=cor,axial=axial,type='up')]
  tmpSurf <- rbind(tmpSurf,tmp)
  return (tmpSurf)
}

#get surface coordinates for the coronal slices having one sections at the genu
getCCgSurfaceOneSec <- function(tmpMidMask,midCor) {
  #find lowest and most posterior tip in genu
  minAxial <- tmpMidMask[cor>midCor,list(min.axial=min(axial)),by=sag]
  tmpMidMask <- merge(tmpMidMask,minAxial,by='sag',all.x=T)
  minTCor <- tmpMidMask[cor>midCor &axial==min.axial &count==1,list(min.t.cor=min(cor)),by=sag]
  tmpMidMask <- merge(tmpMidMask,minTCor,by='sag',all.x=T)
  tmpSurf <- tmpMidMask[cor<min.t.cor & count==1,list(sag=sag,axial=axial,cor=cor,type='low')]
  tmp <- tmpMidMask[cor>min.t.cor & count==1,list(sag=sag,axial=axial,cor=cor,type='up')]
  tmpSurf <- rbind(tmpSurf,tmp)
  tmp <- tmpMidMask[count==n,list(sag=sag,axial=axial,cor=cor,type='up')]
  tmpSurf <- rbind(tmpSurf,tmp)
  return (tmpSurf)
}

#get inner and outer surface of the corpus callosum at the midsagittal slices
getCCSurface <- function(mask,midSag) {
  #find inner and outer surface of CC on midsagittal slice
  #inferior to superior: axial slice small to large
  #mid coronal slice
  midCor <- round((min(mask[sag==midSag]$cor) + max(mask[sag==midSag]$cor))/2)
  midMask = mask[sag==midSag,
                 list(axial=axial,n=length(axial), count=seq(1,length(axial))),by=list(sag,cor)]
  trow <- nrow(midMask)
  nextAxial <- c(midMask$axial[2:trow],NA)
  midMask[,next.axial:=nextAxial]
  midMask[count==n,next.axial:=NA]
  midMask[,':='(max.cor.cc=max(cor),min.cor.cc=min(cor)),by=sag]
  midMask[,s.cor:=2]
  midMask[cor<midCor,s.cor:=1]
  midMask[,min.axial:=min(axial),by=list(sag,s.cor)]
  surface <- data.table()

  #for coronal slices containing two separate sections of the CC
  if (nrow(midMask[axial<next.axial-1 &cor!=max.cor.cc &cor!=min.cor.cc])>0) {
    twoSection <- midMask[axial<next.axial-1 &cor!=max.cor.cc &cor!=min.cor.cc,
                          list(twoSec=T,cor=cor,sag=sag,s.cor=s.cor,axial=axial,min.axial=min.axial)]
    setkey(twoSection,sag,cor)
    #remove duplcated coronal sections that happen when >2 sections occur
    twoSection <- unique(twoSection)
    #remove two sections affecting only one coronal slice if it is not at the bottom axial slice
    twoSection[,n:=length(cor),by=list(sag,s.cor)]
    #for single slice two sections, add the voxel in the two section to lower surface
    if (nrow(twoSection[n==1])>0) {
      tmp <- twoSection[n==1,list(twoSec=twoSec,sag=sag,cor=cor,axial=axial)]
      midMask <- merge(midMask,tmp,by=c('sag','cor','axial'),all.x=T)
      #check whether it is on the outer surface
      midMask[,min.cor:=min(cor),by=list(sag,axial)]
      if (nrow(midMask[twoSec==T &cor>min.cor])>0) {
        tmp <- midMask[twoSec==T &cor>min.cor,list(sag=sag,axial=next.axial,cor=cor,type='low')]
        surface <- rbind(surface,tmp)
      }
      midMask[,':='(twoSec=NULL,min.cor=NULL)]
    }
    twoSection <- twoSection[n>1 |(n==1 &axial==min.axial)]
    twoSection[,':='(n=NULL,axial=NULL,min.axial=NULL)]
    #for those having two sections at genu or splenium, remove the outside ones
    trow <- nrow(twoSection)
    nextCor <- c(twoSection$cor[2:trow],NA)
    twoSection[,next.cor:=nextCor]
    setkey(twoSection,sag,cor)
    twoSection[,':='(n=length(cor),count=seq(1,length(cor))),by=sag]
    twoSection[n==count,next.cor:=NA]
    if (nrow(twoSection[next.cor-cor>1])>0) {
      nTwoSections <- twoSection[next.cor-cor>1]
      setkey(nTwoSections,sag,cor)
      nTwoSections[,':='(n.2s=length(cor),count.2s=seq(1,length(cor))),by=list(sag,s.cor)]
      #genu and splenium
      if (nrow(nTwoSections[n.2s>1])) {
        sec2del <- nTwoSections[n.2s>1 &count.2s==1,list(sag=sag,rm.cor=next.cor,s.cor=s.cor)]
        twoSection <- merge(twoSection,sec2del,by=c('sag','s.cor'),all.x=T)
        twoSection <- twoSection[is.na(rm.cor) |(s.cor==1 &cor>=rm.cor) |(s.cor==2 &cor<rm.cor)]
        twoSection[,rm.cor:=NULL]
      }
    }
    twoSection[,':='(next.cor=NULL,s.cor=NULL,n=NULL,count=NULL)]
    midMask[,':='(s.cor=NULL,min.axial=NULL)]
    midMask <- merge(midMask,twoSection,by=c('sag','cor'),all.x=T)
    midMask[is.na(twoSec),twoSec:=F]
    #on each sag, find smallest and largest coronal slices for having two sections
    corSlices <- midMask[twoSec==T,list(min.s.cor=min(cor),max.s.cor=max(cor)),by=sag]
    midMask <- merge(midMask,corSlices,by='sag',all.x=T)
    
    #find the tip and mark as tip: smallest coronal slice # and have no smaller axial slice #
    #on the same coronal slice
    #for genu
    if (nrow(midMask[twoSec==T &cor>midCor])>0) {
      midMask[twoSec==T &cor>midCor,min.cor:=min(cor),by=sag]
      lowerInner <- midMask[cor==min.cor &axial<next.axial-1,list(l.inn.a=min(axial)),by=sag]
      midMask <- merge(midMask,lowerInner,by="sag",all.x=T)
      tmp <- midMask[cor==min.cor &axial<=l.inn.a,list(sag=sag,axial=axial,cor=cor,type='tip')]
      midMask[,tip:=F]
      if (nrow(tmp)>0) {
        midMask[cor==min.cor &axial<=l.inn.a,tip:=T]
        surface <- rbind(surface,tmp)
      } 
      midMask[,':='(min.cor=NULL,l.inn.a=NULL)]
    }
    #for splenium
    if (nrow(midMask[twoSec==T &cor<midCor])>0) {
      midMask[twoSec==T &cor<midCor,max.cor:=max(cor),by=sag]
      lowerInner <- midMask[cor==max.cor &axial<next.axial-1,list(l.inn.a=min(axial)),by=sag]
      midMask <- merge(midMask,lowerInner,by="sag",all.x=T)
      #check if there is another voxel above the tip, does not work
      tips <- midMask[cor==max.cor &axial<=l.inn.a,list(sag=sag,axial=axial,cor=cor,type='tip')]
      if (nrow(tips)>0) {
        midMask <- merge(midMask,tips,by=c('sag','axial','cor'),all.x=T)
        midMask[type=='tip',tip:=T]
        surface <- rbind(surface,tips)
        midMask[,type:=NULL]
      } 
      midMask[,':='(max.cor=NULL,l.inn.a=NULL)]
    }
    #for coronal slices with two sections on each sagittal slice
    tmp <- midMask[twoSec==T & axial<next.axial-1 &tip==F,list(sag=sag,axial=axial,cor=cor,type='low')]
    surface <- rbind(surface,tmp)
    tmp <- midMask[twoSec==T & axial<next.axial-1,list(sag=sag,axial=next.axial,cor=cor,type='low')]
    surface <- rbind(surface,tmp)
    tmp <- midMask[twoSec==T &tip==F &((count==1 &axial==next.axial-1)|count==n),
                 list(sag=sag,axial=axial,cor=cor,type='up')]
    surface <- rbind(surface,tmp)
  
    #The two sections are at the genu
    if (nrow(midMask[min.s.cor>midCor])>0) {
      #add tip curce when the two sections are at the genu
      tmp <- midMask[twoSec==F &min.s.cor>midCor &cor>max.s.cor & (count==1|count==n),
                     list(sag=sag,axial=axial,cor=cor,type='up')]
      surface <- rbind(surface,tmp)
      #add remaining surface 
      tmpMidMask <- midMask[twoSec==F &min.s.cor>midCor &cor<min.s.cor]
      tmp <- getCCsSurfaceOneSec(tmpMidMask,midCor)               
      surface <- rbind(surface,tmp)
    } 
  
    #both at the splenium
    if (nrow(midMask[max.s.cor<midCor])>0) {
      #add tip
      tmp <- midMask[twoSec==F &max.s.cor<midCor &cor<min.s.cor & (count==1|count==n),
                     list(sag=sag,axial=axial,cor=cor,type='up')]
      surface <- rbind(surface,tmp)
      #add remaining surface 
      tmpMidMask <- midMask[twoSec==F &max.s.cor<midCor &cor>max.s.cor]
      tmp <- getCCgSurfaceOneSec(tmpMidMask,midCor)               
      surface <- rbind(surface,tmp)
    }
  
    #one in the genu and the other in the splenium
    if (nrow(midMask[min.s.cor<midCor &max.s.cor>midCor])>0) {
      tmp <- midMask[twoSec==F &min.s.cor<midCor &max.s.cor>midCor &(cor>max.s.cor|cor<min.s.cor) 
                     &(count==1|count==n),list(sag=sag,axial=axial,cor=cor,type='up')]
      surface <- rbind(surface,tmp)
      tmp <- midMask[twoSec==F &min.s.cor<midCor &max.s.cor>midCor &cor<max.s.cor &cor>min.s.cor 
                     &count==1,list(sag=sag,axial=axial,cor=cor,type='low')]
      surface <- rbind(surface,tmp)
      tmp <- midMask[twoSec==F &min.s.cor<midCor &max.s.cor>midCor &(cor<max.s.cor|cor>min.s.cor) 
                     &count==n,list(sag=sag,axial=axial,cor=cor,type='up')]
      surface <- rbind(surface,tmp)
    }
  } else {
    midMask[,twoSec:=F]
    #splenium half
    tmpMidMask <- midMask[cor<=midCor]
    tmp <- getCCsSurfaceOneSec(tmpMidMask,midCor+1)               
    surface <- rbind(surface,tmp)
    #genu half
    tmpMidMask <- midMask[cor>midCor]
    tmp <- getCCgSurfaceOneSec(tmpMidMask,midCor)               
    surface <- rbind(surface,tmp)
  }
  setkey(surface,sag,axial,cor)
  surface = unique(surface)
  
  #add voxels in the coronal gaps
  midMask[,':='(next.axial=NULL,twoSec=NULL,n=NULL,count=NULL,min.s.cor=NULL,max.s.cor=NULL)]
  setkey(midMask,sag,axial,cor)
  midMask[cor<midCor,s.cor:=1]
  midMask[cor>=midCor,s.cor:=2]
  midMask[,':='(n=length(cor),count=seq(1,length(cor))),by=list(sag,axial,s.cor)]
  trow <- nrow(midMask)
  nextCor <- c(midMask$cor[2:trow],NA)
  midMask[,next.cor:=nextCor]
  midMask[count==n,next.cor:=NA]
  
  #mark voxels below the the lower tip of the inner surface as outer surface
  ccTips <- data.table()
  midMask[,min.axial:=NA]
  if (nrow(surface[cor<midCor &type=="tip"])>0 ) {
    ccTips <- surface[cor<midCor &type=="tip",list(add.min.axial=min(axial),s.cor=1),by=sag]
  }
  if (nrow(surface[cor>midCor &type=="tip"])>0) {
    tmp <- surface[cor>midCor &type=="tip",list(add.min.axial=min(axial),s.cor=2),by=sag]
    ccTips <- rbind(ccTips,tmp)
  }
  if (nrow(ccTips)>0) {
    midMask <- merge(midMask,ccTips,by=c('sag','s.cor'),all.x=T)
    midMask[,min.axial:=add.min.axial]
    midMask[,add.min.axial:=NULL]
  } 
  
  #add voxels to outer surface
  midMask[,min.s.axial:=min(axial),by=list(sag,s.cor)]
  midMask[,max.cor:=max(cor),by=list(sag,s.cor)]
  tmp <- midMask[cor==max.cor,list(max.t.axial=max(axial)),by=list(sag,s.cor)]
  midMask <- merge(midMask,tmp,by=c('sag','s.cor'),all.x=T)
  newSurf <- midMask[((cor<midCor-1|cor>midCor) 
                      &((count==1 &s.cor==1) |(count==n &s.cor==2) |(axial<min.axial &count==n &s.cor==1) 
                        |(axial<min.axial &count==1 &s.cor==2))) 
                     |(count>1 &count<n &cor<next.cor-1 &cor>midCor 
                       &(axial<=min.s.axial+5 |axial==min.axial |axial>max.t.axial)),
                     list(sag=sag,axial=axial,cor=cor)]
  if (nrow(newSurf)>0) {
    midMask[,min.c.axial:=min(axial),by=list(sag,cor)]
    tmp <- midMask[count<n &cor<next.cor-1 &cor>midCor 
                   &(axial==min.c.axial |axial==min.axial |axial>max.t.axial),
                   list(sag=sag,axial=axial,cor=next.cor)]
    newSurf <- rbind(newSurf,tmp)
    #remove duplications
    newSurf <- merge(newSurf,surface,by=c('sag','axial','cor'),all.x=T)
    newSurf <- newSurf[is.na(type)]
    newSurf[,type:="up"]
    surface <- rbind(surface,newSurf)
    midMask[,min.c.axial:=NULL]
  }
  #add voxels to inner surface
  tmp <- midMask[(cor<midCor-1|cor>midCor) & (is.na(min.axial) |axial>min.axial) &
                  ((count==n &s.cor==1) |(count==1 &s.cor==2)),list(sag=sag,axial=axial,cor=cor)]
  tmp2 <- midMask[(cor<midCor-1|cor>midCor) & cor<next.cor-1,list(sag=sag,axial=axial,cor=cor)]
  tmp <- rbind(tmp,tmp2)
  tmp2 <- midMask[(cor<midCor-1|cor>midCor) & cor<next.cor-1,list(sag=sag,axial=axial,cor=next.cor)]
  tmp <- rbind(tmp,tmp2)
  tmp <- unique(tmp)
  #remove duplications
  if (nrow(tmp)>0) {
    tmp <- merge(tmp,surface,by=c('sag','axial','cor'),all.x=T)
    tmp <- tmp[is.na(type)]
    tmp[,type:="low"]
    surface <- rbind(surface,tmp)
  }
  
  #for those with multiple tips, mark top tip as lower boundary, remaining tip as upper boundary
  setkey(surface,type,sag,cor,axial)
  surface[type=="tip",':='(n=length(type),max.axial=max(axial)),by=list(sag,cor)]
  surface[type=="tip" &n>1 &axial==max.axial, type:="low"]
  surface[type=="tip" &n>1 &axial<max.axial, type:="up"]
  surface[,':='(max.axial=NULL,n=NULL)]
  
  return(surface)
}

#label all voxels on the surface as genu, body, or splenium based on input type
getCCType <- function(surface,ccTypeFN,filePath,i) {
  #determine cc type based on the label at the outer surface
  maskGenu <- data.table(read.csv(paste(filePath,ccTypeFN$ccgFN[i],sep=''),header=F,
                                  col.names=c('sag','cor','axial')))
  setkey(maskGenu,sag,cor,axial)
  maskGenu[,cc.type:="genu"]
  surface <- merge(surface,maskGenu,by=c('sag','cor','axial'),all.x=T)
  
  #body
  maskBody <- data.table(read.csv(paste(filePath,ccTypeFN$ccbFN[i],sep=''),header=F,
                                  col.names=c('sag','cor','axial')))
  setkey(maskBody,sag,cor,axial)
  maskBody[,t.type:="body"]
  surface <- merge(surface,maskBody,by=c('sag','cor','axial'),all.x=T)
  surface[is.na(cc.type),cc.type:=t.type]
  surface[,t.type:=NULL]
  
  #splenium
  maskSplenium <- data.table(read.csv(paste(filePath,ccTypeFN$ccsFN[i],sep=''),header=F,
                                      col.names=c('sag','cor','axial')))
  setkey(maskSplenium,sag,cor,axial)
  maskSplenium[,t.type:="splenium"]
  surface <- merge(surface,maskSplenium,by=c('sag','cor','axial'),all.x=T)
  surface[is.na(cc.type),cc.type:=t.type]
  surface[,t.type:=NULL]
  return(surface)
}

#segment surface into equal length for both upper and lower boundaries
#by first order voxels on the surfaces of the CC from the tip of splenium to the tip of genu
orderCCVoxels <- function(segSurface,midSagOne) {
  r.x <- midSagOne$r.x
  r.y <- midSagOne$r.y
  r.z <- midSagOne$r.z
  midSag <- midSagOne$midsag
  #divide the surface into 4 sections 1 upper front, 2 upper back, 3 lower front, and 4 lower back
  #find upper and lower division of the surface
  midCor <- round((min(segSurface[sag==midSag]$cor) + max(segSurface[sag==midSag]$cor))/2)
  segSurface[,s.cor:=2]
  segSurface[cor<midCor,s.cor:=1]
  segSurface[,':='(min.cor=min(cor),max.cor=max(cor)),by=list(sag,s.cor,type)]
  ccTips <- segSurface[s.cor==1 &cor==min.cor &type!="tip",
                    list(cor.tip=unique(cor),axial.tip=max(axial),s.cor=unique(s.cor)),
                    by=list(sag,type)]
  tmp <- segSurface[s.cor==2 &cor==max.cor &type!="tip",
                 list(cor.tip=unique(cor),axial.tip=max(axial),s.cor=unique(s.cor)),
                 by=list(sag,type)]
  ccTips <- rbind(ccTips,tmp)
  segSurface <- merge(segSurface,ccTips,by=c('sag','type','s.cor'),all.x=T)
  #find the most posteior coronal slice in the genu demarcating the uprise of the up boundary
  segSurface[,min.axial:=min(axial),by=list(sag,s.cor)]
  tmp <- segSurface[s.cor==2 &axial==min.axial &type=="up",list(min.u.cor=min(cor)),
                    by=list(sag,s.cor,type)]
  segSurface <- merge(segSurface,tmp,by=c('sag','s.cor','type'),all.x=T)
  ccTips <- ccTips[type=='low' &s.cor==2,list(sag=sag,axial.l.tip=axial.tip)]
  segSurface <- merge(segSurface,ccTips,by='sag',all.x=T)
  segSurface[,max.axial:=max(axial),by=list(sag,cor)]
  #change sections
  segSurface[s.cor==1 &axial<axial.tip, s.cor:=3]
  segSurface[s.cor==2 &(axial<axial.tip |(axial==axial.tip &type=='up' &cor==cor.tip)
                        |(type=='up' &cor<min.u.cor &axial<axial.l.tip)
                        |(type=='up' &axial==axial.tip &cor<cor.tip-1 &axial<max.axial)), s.cor:=4]
  #for an uprising curve in the splenium
  segSurface[s.cor==1 &cor>min.cor+2 &axial==axial.tip &type=='low',s.cor:=3]
  segSurface[,':='(min.cor=NULL,max.cor=NULL,cor.tip=NULL,axial.tip=NULL,min.axial=NULL,
                   min.u.cor=NULL,axial.l.tip=NULL,max.axial=NULL)]
  segSurface[type=="tip" &cor<midCor,s.cor:=3]
  segSurface[type=="tip" &cor>midCor,s.cor:=4]
  
  #remove the edge: the vertical at the tip of the genu 
  voxels2Add <- data.table()
  if (nrow(segSurface[s.cor==3 &(type=="low"|type=="tip")])>0) {
    tmp <- segSurface[s.cor==3 &(type=="low"|type=="tip"),list(max.cor3=max(cor)),by=sag]
    segSurface <- merge(segSurface,tmp,by=c('sag'),all.x=T)
    tmp <- segSurface[s.cor==3 &cor==max.cor3,
                      list(n.voxel.c=length(axial),min.axial3=min(axial),max.axial3=max(axial)),by=sag]
    segSurface <- merge(segSurface,tmp,by=c('sag'),all.x=T)
    #number of voxels on the horizontal line
    tmp <- segSurface[s.cor==3 &axial==max.axial3,list(n.voxel.a=length(axial)),by=sag]
    segSurface <- merge(segSurface,tmp,by=c('sag'),all.x=T)
    #find whether there is a dip in the lower surface
    tmp <- segSurface[s.cor==3 &type=="low",list(min.axial.l=min(axial)),by=list(sag,s.cor)]
    segSurface <- merge(segSurface,tmp,by=c('sag','s.cor'),all.x=T)
  } else {
    segSurface[,':='(max.cor3=NA,n.voxel.c=NA,min.axial3=NA,max.axial3=NA,n.voxel.a=NA,min.axial.l=NA)]
  }
  
  #n.voxel.a also contains 2 voxels from up and low surfaces
  setkey(segSurface,sag,s.cor,cor,axial)
  segSurface[,':='(count=as.integer(0),count.low=as.integer(0))]
  if (nrow(segSurface[s.cor==3 &axial==max.axial3])>0) {
    segSurface[s.cor==3 &axial==max.axial3,count:=seq(1,length(cor)),by=list(sag,s.cor)]
    #must have >1 low voxels
    segSurface[s.cor==3 &axial==max.axial3 &type=='low',count.low:=seq(1,length(cor)),by=list(sag,s.cor)]
    if (nrow(segSurface[count>3 &max.axial3>min.axial.l &count.low>1])>0) {
      #change the most outside one as up
      segSurface[count>3 &count==n.voxel.a &max.axial3>min.axial.l &count.low>1,type:="up"]
      #remove the middle one
      if (nrow(segSurface[count>4 &max.axial3>min.axial.l &count.low>1])>0) {
        tmp <- segSurface[count>3 &max.axial3>min.axial.l &count.low>1 &count<n.voxel.a,
                          list(sag=sag,s.cor=s.cor,cor=cor,type="edge",axial=axial)]
        voxels2Add <- rbind(voxels2Add,tmp)
        tmp[,':='(type=NULL,new.type='remove')]
        segSurface <- merge(segSurface,tmp,by=c('sag','s.cor','cor','axial'),all.x=T)
        segSurface <- segSurface[is.na(new.type)]
        segSurface[,new.type:=NULL]
      }
    }
  }

  #for small splenium curve, change tip for upper and lower boundaries 
  #to the edge at the bottom axial slice
  #change the voxel type between the old tip and new tip to low and assign to s.cor 7
  segSurface[,min.cor:=min(cor),by=list(sag,type,s.cor)]
  if (nrow(segSurface[s.cor==3 &type=='low'&max.cor3-min.cor<2])>0) {
    #mark the sagittal slice for change
    change <- segSurface[s.cor==3 &type=='low'&max.cor3-min.cor<2,list(sag=unique(sag),is.ch=T)]
    segSurface <- merge(segSurface,change,by='sag',all.x=T)
    #determine the edge at the bottom axial slice
    segSurface[s.cor==3,min.axial:=min(axial),by=list(sag)]
    newTip <- segSurface[axial==min.axial,list(tip.cor=max(cor)),by=sag]
    segSurface <- merge(segSurface,newTip,by='sag',all.x=T)
    segSurface[s.cor==3 &cor>=tip.cor &axial>min.axial 
               &(axial==min.axial3 &cor==max.cor3 |axial<min.axial3 &cor<max.cor3 ),
               ':='(type='low',s.cor=7)]
    segSurface[s.cor==3 &cor>=tip.cor &type=='up' &axial>min.axial,':='(type='low',s.cor=7)]
    segSurface[axial==min.axial &cor==tip.cor,type:='up']
    #add the tip as low if the curve is affecting >1 slice
    tip <- segSurface[axial==min.axial &cor==tip.cor]
    if (nrow(segSurface[s.cor==7])>0) {
      tmp <- segSurface[s.cor==7, list(max.axial7=max(axial)),by=sag]
      tip <- merge(tip,tmp,by='sag')
      tip <- tip[axial<max.axial7-1]
      tip[,':='(type='low',s.cor=7,max.axial7=NULL)]
      segSurface <- rbind(segSurface,tip)
    }
    segSurface[,':='(min.axial=NULL,tip.cor=NULL,is.ch=NULL)]
  }
  segSurface[,':='(min.cor=NULL)]
  
  #duplicate the tip
  if (nrow(segSurface[count==3 &count==n.voxel.a &max.axial3>min.axial.l &count.low>1])>0) {
    #If there is an existing tip in front to the voxel, mark it as up
    if (nrow(segSurface[count==3 &count==n.voxel.a &max.axial3>min.axial.l&count.low>1 
                        &s.cor==3 &type=="tip"])>0) {
      tips <- segSurface[s.cor==3 &type=="tip",list(sag=sag,cor=cor,axial=axial,n.tips=length(cor))]
      tips[,':='(n.tips=length(cor),max.cor=max(cor)),by=sag]
      tips[,newType:="up"]
      tips[cor<max.cor,newType:="low"]
      tips[,':='(n.tips=NULL,max.cor=NULL)]
      segSurface <- merge(segSurface,tips,by=c('sag','axial','cor'),all.x=T)
      segSurface[!is.na(newType),type:=newType]
      segSurface[!is.na(newType),n.voxel.a:=NA]
      #change max.cor3
      maxCor <- tips[newType=='up',list(sag=sag,max.cor3=cor)]
      segSurface[,max.cor3:=NULL]
      segSurface <- merge(segSurface,maxCor,by=c('sag'),all.x=T)
      segSurface[,newType:=NULL]
    }
    #re-calculating slice numbers
    segSurface[,':='(n.voxel.c=NULL,min.axial3=NULL,max.axial3=NULL,n.voxel.a=NULL)]
    if (nrow(segSurface[s.cor==3 &cor==max.cor3])>0) {
      tmp <- segSurface[s.cor==3 &cor==max.cor3,
                        list(n.voxel.c=length(axial),min.axial3=min(axial),max.axial3=max(axial)),by=sag]
      segSurface <- merge(segSurface,tmp,by=c('sag'),all.x=T)
      #number of voxels on the horizontal line
      tmp <- segSurface[s.cor==3 &axial==max.axial3,list(n.voxel.a=length(axial)),by=sag]
      segSurface <- merge(segSurface,tmp,by=c('sag'),all.x=T)
    }   
  }   

  #remove and mark the vertical line as edge
  segSurface[,min.cor:=min(cor),by=list(sag,type,s.cor)]
  if (nrow(segSurface[s.cor==3 &cor==max.cor3 &cor>min.cor &max.cor3-min.cor>2 &n.voxel.c>2  
                      &(n.voxel.a<4 |min.axial.l==max.axial3) & axial<max.axial3 &axial>min.axial3 
                      &n.voxel.c==max.axial3-min.axial3+1])>0) {
    tmp <- segSurface[s.cor==3 &cor==max.cor3 &cor>min.cor &max.cor3-min.cor>2 &n.voxel.c>2 
                      &(n.voxel.a<4 |min.axial.l==max.axial3) & axial<max.axial3 &axial>min.axial3  
                      &n.voxel.c==max.axial3-min.axial3+1]
    tmp <- tmp[,list(sag=sag,s.cor=s.cor,type="edge",axial=axial,cor=cor)]
    voxels2Add <- rbind(voxels2Add,tmp)
    tmp[,':='(type=NULL,edge=T)]
    segSurface <- merge(segSurface,tmp,by=c('sag','s.cor','axial','cor'),all.x=T)
    segSurface <- segSurface[is.na(edge)]
    segSurface[,edge:=NULL]
  }
  segSurface[,':='(n.voxel.c=NULL,min.axial3=NULL,n.voxel.a=NULL,min.axial.l=NULL,count=NULL,
                   count.low=NULL)]
  
  #find single line of voxels in section 3 and 4 next to the inner edge and duplicate
  if (nrow(segSurface[s.cor>2])>0) {
    #find the edge of the inner curve
    tmp <- segSurface[type=='low',list(max.i.cor=max(cor),min.i.cor=min(cor)),by=sag]
    segSurface <- merge(segSurface,tmp,by='sag')
    #find those having only one voxel on each coronal slice in the two section areas 
    if (nrow(segSurface[s.cor==3])>0) {
      segSurface[(s.cor==3 |s.cor==7) &cor>min.i.cor &cor<=max.cor3,n.axial:=length(axial),
                 by=list(sag,cor)]
    }  
    if (nrow(segSurface[s.cor==4])>0) {
      tmp <- segSurface[s.cor==4,list(min.cor4=min(cor)),by=sag]
      segSurface <- merge(segSurface,tmp,by='sag')
      segSurface[s.cor==4 &cor>=min.cor4 &cor<max.i.cor,n.axial:=length(axial),by=list(sag,cor)]
    }
    #duplicate
    if (nrow(segSurface[(s.cor==4 &n.axial==1) |(s.cor==3 &n.axial==1 &cor<max.cor3 
                                                 & max.axial3-axial<2)])>0) {
      edge <- segSurface[(s.cor==4 &n.axial==1) |(s.cor==3 &n.axial==1 &cor<max.cor3
                                                  & max.axial3-axial<2)]
      edge[,type:='low']
      segSurface[(s.cor==4 &n.axial==1) |(s.cor==3 &n.axial==1 &cor<max.cor3),type:='up']
      segSurface <- rbind(segSurface,edge)
    }
    if (nrow(segSurface[s.cor==4])>0) {
      segSurface[,min.cor4:=NULL]
    }
    segSurface[,':='(max.i.cor=NULL,min.i.cor=NULL,n.axial=NULL)]
  }
  setkey(segSurface,sag,type,s.cor,cor,axial)
  
  #remove and mark the middle point in the vertical line in the tip of genu as edge
  tmp <- segSurface[s.cor==4 &type=="up",list(min.cor4=min(cor)),by=sag]
  segSurface <- merge(segSurface,tmp,by=c('sag'),all.x=T)
  tmp <- segSurface[s.cor==4 &cor==min.cor4 &type=="up",
                    list(n.voxel.c=length(axial),min.axial4=min(axial)),by=sag]
  segSurface <- merge(segSurface,tmp,by=c('sag'),all.x=T)
  if (nrow(segSurface[s.cor==4 &type=="up" &cor==min.cor4 &n.voxel.c>1 
                      &axial>min.axial4])>0) {
    tmp <- segSurface[s.cor==4 &type=="up" &cor==min.cor4 &n.voxel.c>1 &axial>min.axial4]
    tmp <- tmp[,list(sag=sag,s.cor=s.cor,type="edge",axial=axial,cor=cor)]
    voxels2Add <- rbind(voxels2Add,tmp)
    segSurface[s.cor==4 &type=="up" &cor==min.cor4 &n.voxel.c>1 &axial>min.axial4,
               type:="edge"]
  }
  segSurface[,':='(n.voxel.c=NULL,min.axial4=NULL,min.cor4=NULL,max.axial3=NULL)]

  #duplicate tips
  if (nrow(segSurface[type=='tip'])>0) {
    tips <- segSurface[type=='tip']
    segSurface[type=='tip',type:='low']
    tips[,type:="up"]
    segSurface <- rbind(segSurface,tips)
  }
  segSurface[,max.cor3:=NULL]
  setkey(segSurface,sag,type,s.cor,cor,axial)

  #find upward curve at tips of the genu and splenium and separate them to section 5 and 6
  segSurface[cor<midCor,':='(min.axial=min(axial),max.cor=max(cor),min.cor=min(cor)),by=list(sag,type)]
  segSurface[cor>midCor,':='(min.axial=min(axial),max.cor=max(cor),min.cor=min(cor)),by=list(sag,type)]
  segSurface[,half:=2]
  segSurface[cor<midCor,half:=1]
  bottomCor <- segSurface[axial==min.axial &half==1,list(bottom.cor=max(cor)),by=list(sag,type,half)]
  tmp <- segSurface[axial==min.axial &half==2,list(bottom.cor=min(cor)),by=list(sag,type,half)]
  bottomCor <- rbind(bottomCor,tmp)
  segSurface <- merge(segSurface,bottomCor,by=c('sag','type','half'),all.x=T)
  #find the tip in the splenium
  tmp <- segSurface[s.cor==3,list(max.i.cor=max(cor)),by=list(sag,type,half)]
  segSurface <- merge(segSurface,tmp,by=c('sag','type','half'),all.x=T)
  if (nrow(segSurface[cor==max.i.cor &bottom.cor<max.i.cor])>0) {
    segSurface[cor==max.i.cor &bottom.cor<max.i.cor,min.axial3:=min(axial),by=list(sag,type,half)]
    tips <- segSurface[axial==min.axial3,list(sag=sag,type=type,half=half,max.axial3=axial)]
    segSurface <- merge(segSurface,tips,by=c('sag','type','half'),all.x=T)
    if (nrow(segSurface[axial==max.axial3+1 &min.cor<max.i.cor &cor-min.cor>1])>0) {
      tmp <- segSurface[(axial==max.axial3+1|axial==min.axial3+2) &min.cor<max.i.cor &cor-min.cor>1,
                        list(sag=sag,type=type,half=half,max.axial3=axial)]
      tips <- rbind(tips,tmp)
      tips <- tips[,list(new.max.axial3=max(max.axial3)),by=list(sag,type,half)]
      segSurface <- merge(segSurface,tips,by=c('sag','type','half'),all.x=T)
      segSurface[!is.na(new.max.axial3),max.axial3:=new.max.axial3]
      segSurface[,new.max.axial3:=NULL]
    }
    #find genu and splenium tip
    segSurface[half==1 &s.cor!=7 &axial>min.axial &axial<=max.axial3 &cor>=bottom.cor &cor-min.cor>1,
               s.cor:=5]
    segSurface[,':='(max.axial3=NULL,min.axial3=NULL)]
  }
  segSurface[s.cor==4 &axial>=min.axial &cor<=bottom.cor,s.cor:=6]
  #find whether there is a downward curve in section 6
  if (nrow(segSurface[s.cor==6]>0)) {
    segSurface[s.cor==6,':='(min.cor6=min(cor),max.axial6=max(axial)),by=list(sag,type)]
    tmp <- segSurface[type=='up' &cor==min.cor6,list(tip.axial=min(axial)),by=list(sag,type)]
    segSurface <- merge(segSurface,tmp,by=c('sag','type'),all.x=T)
    tmp <- segSurface[type=='up' &s.cor==6 &axial==max.axial6,list(top.cor=min(cor)),by=list(sag,type)]
    segSurface <- merge(segSurface,tmp,by=c('sag','type'),all.x=T)
    segSurface[s.cor==6 &type=='up'&axial<=max.axial6 &max.axial6>tip.axial &cor<=top.cor,s.cor:=8]
    segSurface[,':='(min.cor6=NULL,tip.axial=NULL,max.axial6=NULL,top.cor=NULL)]
  }
  
  #for those without section 3 check whether there is a dent on the lower surface
  if (nrow(segSurface[min.cor==bottom.cor])>0) {
    trow <- nrow(segSurface)
    nextAxial <- c(segSurface$axial[2:trow],NA)
    segSurface[,next.axial:=nextAxial]
    if (nrow(segSurface[min.cor==bottom.cor &s.cor==1 &axial<next.axial-1 &type=='low'])>0) {
      tmp <- segSurface[min.cor==bottom.cor &s.cor==1 &axial<next.axial-1 &type=='low',
                        list(axial3=max(axial)),by=sag]
      segSurface <- merge(segSurface,tmp,by='sag',all.x=T)
      segSurface[type=='low' &s.cor==1 &axial<=axial3,s.cor:=7]
      segSurface[,':='(axial3=NULL)]
    }
    segSurface[,next.axial:=NULL]
  }
  segSurface[,':='(min.axial=NULL,min.cor=NULL,max.cor=NULL,bottom.cor=NULL,max.i.cor=NULL,half=NULL)]

  
  #order voxels from the tip of splenium to the tip of the genu
  segSurface <- indexCCSurface(segSurface)
  #don't forget to add voxels2Add
  if (nrow(voxels2Add)>0) {
    voxels2Add[,':='(id.voxel=NA,n.voxel=NA)]
    segSurface <- rbind(segSurface,voxels2Add)
  }
  return(segSurface)
}

#apply index to each point in the upper and lower surfaces from splenium to genu
indexCCSurface <- function(segSurface) {
  #for s.cor 3 and 4, coronal slice should be ordered from large to small
  #for s.cor 2 and 4, axial slice should be ordered from large to small
  index <- segSurface[,list(cor=unique(cor)),by=list(sag,type,s.cor)]
  setkey(index,sag,s.cor,type,cor)
  index[,id.cor:=seq(length(cor),1),by=list(sag,type,s.cor)]
  segSurface <- merge(segSurface,index,by=c('sag','type','s.cor','cor'),all.x=T)

  index <- segSurface[,list(axial=unique(axial)),by=list(sag,type,s.cor)]
  setkey(index,sag,type,s.cor,axial)
  index[,id.axial:=seq(length(axial),1),by=list(sag,type,s.cor)]
  segSurface <- merge(segSurface,index,by=c('sag','type','s.cor','axial'),all.x=T)
  
  #order all points on the surface
  surf0 <- data.table()
  surf1 <- data.table()
  surf6 <- data.table()
  surf7 <- data.table()
  if (nrow(segSurface[s.cor==7])>0) {
    surf0 <- segSurface[s.cor==7]
    setkey(surf0,sag,s.cor,type,axial,cor)
  }
  if (nrow(segSurface[s.cor==5])>0) {
    surf1 <- segSurface[s.cor==5]
    setkey(surf1,sag,s.cor,type,id.cor,id.axial)
  }
  surf2.1 <- segSurface[s.cor==3 &type=="up"]
  setkey(surf2.1,sag,s.cor,type,id.cor,axial)
  surf2.2 <- segSurface[s.cor==3 &type=="low"]
  setkey(surf2.2,sag,s.cor,type,axial,id.cor)
  surf3.1 <- segSurface[s.cor==1 &type=="up"]
  setkey(surf3.1,sag,s.cor,type,cor,axial)
  #if voxel distance is large, check surface and change the order
  surf3.2 <- segSurface[s.cor==1 &type=="low"]
  setkey(surf3.2,sag,s.cor,type,cor,axial)
#  setkey(surf3.2,sag,s.cor,type,axial,cor)
  surf4 <- segSurface[s.cor==2]
  setkey(surf4,sag,s.cor,type,cor,id.axial)
  surf5 <- segSurface[s.cor==4]
  setkey(surf5,sag,s.cor,type,id.axial,id.cor)
  if (nrow(segSurface[s.cor==6])>0) {
    surf6 <- segSurface[s.cor==6]
    setkey(surf6,sag,s.cor,type,id.cor,axial)
  }
  if (nrow(segSurface[s.cor==8])>0) {
    surf7 <- segSurface[s.cor==8]
    setkey(surf7,sag,s.cor,type,id.axial,id.cor)
  }
  segSurface <- rbind(surf0,surf1,surf2.1,surf2.2,surf3.1,surf3.2,surf4,surf5,surf6,surf7)
  segSurface[,':='(id.cor=NULL,id.axial=NULL)]
  segSurface[,':='(id.voxel=seq(1,length(cor)),n.voxel=length(cor)),by=list(sag,type)]
  setkey(segSurface,sag,type,id.voxel)
  return(segSurface)
}  

#adjust order for those with distance >2 and involving more than 2 voxels
adjustOrder <- function(segSurface) {
    #combine s.cor 4 and 6
    segSurface[d.voxel>=2,n.bigd:=length(cor),by=list(sag,s.cor)]
    segSurface[(s.cor==4|s.cor==6) &d.voxel>=2,n.bigd:=length(cor),by=sag]
    if (nrow(segSurface[n.bigd>1])>0) {
      #find the sections need to be recordered
      disorgVoxels <- segSurface[n.bigd>1,list(sag=sag,s.cor=s.cor,type=type,
                                               id.voxel=id.voxel,axial=axial,cor=cor)]
      disorgVoxels[s.cor==6,s.cor:=4]
      setkey(disorgVoxels,sag,type,s.cor,id.voxel)
      #check id distance between voxels and axial distances
      trow <- nrow(disorgVoxels)
      nextID <- c(disorgVoxels$id.voxel[2:trow],NA)
      disorgVoxels[,id.next:=nextID]
      disorgVoxels[,':='(n=length(cor),count=seq(1,length(cor))),by=list(sag,type,s.cor)]
      disorgVoxels[n==count,id.next:=NA]
      disorgVoxels[,id.d:=id.next-id.voxel]
      disorgVoxels[,':='(min.axial=min(axial),max.axial=max(axial)),by=list(sag,s.cor,type)]
      #assign section number
      sections <- disorgVoxels[count==1 |(id.d>9 &(max.axial-min.axial>3)),
                               list(sag=sag,type=type,s.cor=s.cor,count=count,n=n)]
      rowNum <- nrow(sections)
      if (nrow(sections)>1) {
        nextCount <- c(sections$count[2:rowNum],NA)
        sections[,next.count:=nextCount]
        sections[,':='(n.sec=length(count),id.sec=seq(1,length(count))),
                 by=list(sag,s.cor,type)]
        sections[id.sec==n.sec,next.count:=as.integer(n+1)]
      } else {
        sections[,':='(next.count=as.integer(n+1),n.sec=1,id.sec=1)]
      }
      sections[,sec.d:=next.count-count]
      secNumber <- vector()
      for (j in 1:nrow(sections)) {
        tmp <- rep(sections$id.sec[j],sections$sec.d[j])
        secNumber <- c(secNumber,tmp)
      }
      disorgVoxels[,sec:=secNumber]
      disorgVoxels[,':='(min.id=min(id.voxel),max.id=max(id.voxel)),
                   by=list(sag,s.cor,type,sec)]
      #retrieve the voxels to be ordered
      disorgVoxels <- disorgVoxels[count==1]
      voxels2ret <- data.table()
      for (k in 1:nrow(disorgVoxels)) {
        tmp <- data.table(sag=disorgVoxels$sag[k],type=disorgVoxels$type[k],
                                 id.voxel=seq(disorgVoxels$min.id[k],disorgVoxels$max.id[k]+1),
                                 sec=disorgVoxels$sec[k])
        voxels2ret <- rbind(voxels2ret,tmp)
      }
      segSurface <- merge(segSurface,voxels2ret,by=c('sag','type','id.voxel'),
                          all.x=T)
      orderVoxels <- segSurface[!is.na(sec)]
      #reverse order
      index <- orderVoxels[,list(cor=unique(cor)),by=list(sag,type,s.cor,sec)]
      setkey(index,sag,type,s.cor,sec,cor)
      index[,id.cor:=seq(length(cor),1),by=list(sag,type,s.cor,sec)]
      orderVoxels <- merge(orderVoxels,index,by=c('sag','type','s.cor','sec','cor'),all.x=T)
      
      index <- orderVoxels[,list(axial=unique(axial)),by=list(sag,type,s.cor,sec)]
      setkey(index,sag,type,s.cor,sec,axial)
      index[,id.axial:=seq(length(axial),1),by=list(sag,type,s.cor,sec)]
      orderVoxels <- merge(orderVoxels,index,by=c('sag','type','s.cor','sec','axial'),all.x=T)

      #reorder according to sections
      surfAdd <- data.table()
      if (nrow(orderVoxels[s.cor==5])>0) {
        surf <- orderVoxels[s.cor==5]
        setkey(surf,sag,s.cor,type,sec,id.axial,id.cor)
        surfAdd <- rbind(surfAdd,surf)
      }
      if (nrow(orderVoxels[s.cor==3 &type=="up"])>0) {
        surf <- orderVoxels[s.cor==3 &type=="up"]
        setkey(surf,sag,s.cor,type,sec,axial,id.cor)
        surfAdd <- rbind(surfAdd,surf)
      }
      if (nrow(orderVoxels[s.cor==3 &type=="low"])>0) {
        surf <- orderVoxels[s.cor==3 &type=="low"]
        setkey(surf,sag,s.cor,type,sec,id.cor,axial)
        surfAdd <- rbind(surfAdd,surf)
      }
      if (nrow(orderVoxels[s.cor==1 &type=="up"])>0) {
        surf <- orderVoxels[s.cor==1 &type=="up"]
        setkey(surf,sag,s.cor,type,sec,axial,cor)
        surfAdd <- rbind(surfAdd,surf)
      }
      #if voxel distance is large, check surface and change the order
      if (nrow(orderVoxels[s.cor==1 &type=="low"])>0) {
        surf <- orderVoxels[s.cor==1 &type=="low"]
        setkey(surf,sag,s.cor,type,sec,axial,cor)
        surfAdd <- rbind(surfAdd,surf)
      }
      if (nrow(orderVoxels[s.cor==2])>0) {
        surf <- orderVoxels[s.cor==2]
        setkey(surf,sag,s.cor,type,sec,id.axial,cor)
        surfAdd <- rbind(surfAdd,surf)
      }
      if (nrow(orderVoxels[s.cor==4])>0) {
        surf <- orderVoxels[s.cor==4]
        setkey(surf,sag,s.cor,type,sec,id.cor,id.axial)
        surfAdd <- rbind(surfAdd,surf)
      }
      if (nrow(orderVoxels[s.cor==6])>0) {
        surf <- orderVoxels[s.cor==6]
        setkey(surf,sag,s.cor,type,sec,axial,id.cor)
        surfAdd <- rbind(surfAdd,surf)
      }
      if (nrow(orderVoxels[s.cor==7])>0) {
        surf <- orderVoxels[s.cor==7]
        setkey(surf,sag,s.cor,type,sec,id.cor,axial)
        surfAdd <- rbind(surfAdd,surf)
      }
      if (nrow(orderVoxels[s.cor==8])>0) {
        surf <- orderVoxels[s.cor==8]
        setkey(surf,sag,s.cor,type,sec,id.cor,id.axial)
        surfAdd <- rbind(surfAdd,surf)
      }
      #change voxel id
      surfAdd[,min.id:=min(id.voxel),by=list(sag,s.cor,type,sec)]
      surfAdd[,id.voxel:=as.integer(seq(1,length(cor))+min.id-1),by=list(sag,s.cor,type,sec)]
      surfAdd[,':='(id.cor=NULL,id.axial=NULL,min.id=NULL,sec=NULL)]
      segSurface <- segSurface[is.na(sec)]
      segSurface[,sec:=NULL]
      segSurface <- rbind(segSurface,surfAdd)
    }
    segSurface[,n.bigd:=NULL]
    setkey(segSurface,sag,type,id.voxel)
    
  return(segSurface)
}

#calculate cumulative distance for ordered coordinates
getCumulativeDistance <- function(segSurface,midSagOne,midSag) {
  r.x <- midSagOne$r.x
  r.y <- midSagOne$r.y
  r.z <- midSagOne$r.z
#  midSag <- midSagOne$midsag
  #remove edges
  segSurface <- segSurface[!(type=="edge")]
  #calculate distance and accumulative distance between points
  tRow <- nrow(segSurface)
  previousSag <- c(NA,segSurface$sag[1:(tRow-1)])
  previousAxial <- c(NA,segSurface$axial[1:(tRow-1)])
  previousCor <- c(NA,segSurface$cor[1:(tRow-1)])
  segSurface[,':='(p.sag=previousSag,p.axial=previousAxial,p.cor=previousCor)]
  segSurface[id.voxel==1,':='(p.sag=NA,p.axial=NA,p.cor=NA)]
  segSurface[!is.na(p.sag), d.voxel:=round(sqrt(((sag-p.sag)*r.x)^2 + ((cor-p.cor)*r.y)^2
                                                   + ((axial-p.axial)*r.z)^2),3)]
  print(paste("max. voxel distance:",max(segSurface[!is.na(d.voxel)]$d.voxel)))
  if (nrow(segSurface[!is.na(d.voxel) &d.voxel>=2])>0) {
    maxNumber <- segSurface[!is.na(d.voxel) &d.voxel>=2,list(max.n=length(cor)),by=sag]
    max.n <- max(maxNumber$max.n)
  } else {
    max.n <- 0
  }
  print(paste("max.number of voxels with distance >=2 per slice:",max.n))
  print(paste("Number of voxels with distance >=2 on midsag slice:",
              nrow(segSurface[sag==midSag &!is.na(d.voxel) &d.voxel>=2])))
  
  #add points to form equal segmentation
  #accumulative sum
  segSurface[is.na(d.voxel),d.voxel:=0]
  segSurface[,':='(d.cusum=round(cumsum(d.voxel),3),d.total=round(sum(d.voxel),3)),by=list(sag,type)]
  return(segSurface)
}

#generate new data points to segment the curve into equal length for both inner and outer curves
#and assign coordinates
segmentCCSurface <- function(segSurface,midSagOne,n.seg) {
  r.x <- midSagOne$r.x
  r.y <- midSagOne$r.y
  r.z <- midSagOne$r.z
  curveLength <- segSurface[id.voxel==1,list(sag=sag,type=type,d.total=d.total)]
  newPoints <- data.table()
  for (j in 1:nrow(curveLength)) {
    tmp <- data.table(sag=curveLength$sag[j], type=curveLength$type[j],  
                      d.voxel=round(curveLength$d.total[j]/n.seg,6), id.seg=seq(1,n.seg))
    tmp[,d.cusum:=round(cumsum(d.voxel),3)]
    newPoints <- rbind(newPoints,tmp)
  }
  newPoints[,':='(d.voxel=NULL,id.voxel=NA,n.voxel=NA)]
  segments <- segSurface[,list(sag=sag,type=type,id.seg=NA,d.cusum=d.cusum,id.voxel=id.voxel,
                               n.voxel=n.voxel)]
  segments <- rbind(segments,newPoints)
  setkey(segments,sag,type,d.cusum)
  segments[,':='(id=seq(1,length(d.cusum)),n=length(d.cusum)),by=list(sag,type)]
  
  #remove new points after last voxel
  segments <- segments[!(id==n &is.na(id.voxel))]
  segments[,n:=length(d.cusum),by=list(sag,type)]
  
  #assign voxel id to new points
  newID <- segments[!is.na(id.voxel),
                    list(sag=sag,type=type,id.voxel=id.voxel,id=id,n.voxel=n.voxel,n=n)]
  tmp <- c(newID$id[2:nrow(newID)],NA)
  newID[,next.id:=tmp]
  newID[,n.seg:=next.id-id]
  newID[id==n,n.seg:=1]
  newIDList <- rep(newID$id.voxel,newID$n.seg)
  segments[,id.voxel:=newIDList]
  setkey(segments,sag,type,id)
  coordinates <- segSurface[id.voxel>1,list(sag=sag,type=type,id.voxel=id.voxel-1,pv.axial=p.axial,
                                            pv.cor=p.cor,nv.axial=axial,nv.cor=cor,v.cusum=d.cusum)]
  lastVoxels <- segSurface[id.voxel==n.voxel,list(sag=sag,type=type,id.voxel=id.voxel,pv.axial=axial,
                                                  pv.cor=cor,nv.axial=NA,nv.cor=NA,v.cusum=d.cusum)]
  coordinates <- rbind(coordinates,lastVoxels)
  setkey(coordinates,sag,type,id.voxel)
  coordinates[,p.d.cusum:=segSurface$d.cusum] 
  segments <- merge(segments,coordinates,by=c('sag','type','id.voxel'),all.x=T)
  #calculate coordinates for new points
  segments[id.seg>0,
           ':='(axial=round(((d.cusum-p.d.cusum)*(nv.axial-pv.axial)/(v.cusum-p.d.cusum)+pv.axial),2),
                cor=round(((d.cusum-p.d.cusum)*(nv.cor-pv.cor)/(v.cusum-p.d.cusum)+pv.cor),2))]
  #add first and last points
  segments[is.na(id.seg) &(id.voxel==1 |id.voxel==n.voxel),
           ':='(axial=as.double(pv.axial),cor=as.double(pv.cor))]
  segments <- segments[!is.na(cor),list(id.voxel=id.voxel,axial=axial,cor=cor,d.cusum=d.cusum),
                       by=list(sag,type)]
  segments[,':='(n=length(cor),id=seq(1,length(cor))),by=list(sag,type)]
  setkey(segments,sag,type,id)
  return(segments)
}

#get closest points on the lower boundary for each points on the upper boundary for the sagittal slice
getClosestPoints <- function(segments,voxelSize) {
  r.x <- voxelSize$r.x
  r.y <- voxelSize$r.y
  r.z <- voxelSize$r.z
  midSag <- voxelSize$midsag
  lowAxial <- segments[type=="low"  &sag==midSag]$axial
  lowCor <- segments[type=="low"  &sag==midSag]$cor
  lowID <- segments[type=="low"  &sag==midSag]$id
  #find the points on the lower boundary that are closest to each point of the upper boundary
  closePoints <- data.table()
  for (j in 1:nrow(segments[type=="up" &sag==midSag])) {
    tmp <- data.table(u.id=segments[type=="up" &sag==midSag]$id[j],
                      u.sag=segments[type=="up" &sag==midSag]$sag[j], 
                      u.axial=segments[type=="up" &sag==midSag]$axial[j],
                      u.cor=segments[type=="up" &sag==midSag]$cor[j], 
                      l.id=lowID, l.axial=lowAxial, l.cor=lowCor)
    closePoints <- rbind(closePoints,tmp)
  }
  closePoints[,d:=round(sqrt(((u.cor-l.cor)*r.y)^2+((u.axial-l.axial)*r.z)^2),3)]
  closePoints[,min.d:=min(d),by=u.id]
  closePoints <- closePoints[d==min.d]
  return(closePoints)
}

#section the CC by lower points having shortest distance to >1 upper points
#and the corresponding upper points that are in the middle
getCCSections <- function(closePoints,segments,segSurface,n.seg,min.connect,midSag) {
  #calculate mid points for up points having shortest distance to the same lower points
  closePoints[,n.up:=length(u.id),by=l.id]
  sectionCC <- closePoints[n.up>=min.connect &l.id>=10 &l.id<=n.seg-5,
                           list(u.id=round(mean(u.id)),n.up=unique(n.up)),by=list(u.sag,l.id)]
  #segment the lower boundary according to the new sections and number of segments of upper
  lastPoint <- data.table(u.sag=unique(sectionCC$u.sag),l.id=n.seg+1,u.id=n.seg+1,n.up=1)
  sectionCC <- rbind(sectionCC,lastPoint)
  trow <- nrow(sectionCC)
  previousLID <- c(1,sectionCC$l.id[1:(trow-1)])
  previousUID <- c(1,sectionCC$u.id[1:(trow-1)])
  sectionCC[,':='(p.l.id=previousLID,p.u.id=previousUID)]
  sectionCC[,':='(n.sec=u.id-p.u.id,id.s=seq(1,length(u.sag)),sag=midSag,type='low',id=l.id)]
  #get cumulative distance for the lower points
  sectionCC <- merge(sectionCC,segments,by=c('sag','type','id'))
  #add section length
  previousCusum <- c(0,sectionCC$d.cusum[1:(trow-1)])
  sectionCC[,p.d.cusum:=previousCusum]
  sectionCC[,ln.s:=d.cusum-p.d.cusum]

  #new segments for lower boundary
  newSections <- data.table()
  for (k in 1:nrow(sectionCC)) {
    tmp <- data.table(sag=sectionCC$u.sag[k], id.seg=seq(sectionCC$p.u.id[k], (sectionCC$u.id[k]-1)),
                      n.sec=sectionCC$n.sec[k], id.s=sectionCC$id.s[k],ln.s=sectionCC$ln.s[k])
    newSections <- rbind(newSections,tmp)
  }
  newSections[,d:=round(ln.s/n.sec,3)]
  newSections[,d.cusum:=cumsum(d),by=sag]
  newSegments <- newSections[id.seg<n.seg,list(sag=sag,id.seg=id.seg,d.cusum=d.cusum,id.voxel=NA,
                                               n.voxel=NA)]
  voxels <- segSurface[sag==midSag &type=='low',
                       list(sag=sag,id.seg=NA,d.cusum=d.cusum,id.voxel=id.voxel,n.voxel=n.voxel)]
  newSegments <- rbind(newSegments,voxels)
  setkey(newSegments,sag,d.cusum)
  newSegments[,':='(id=seq(1,length(d.cusum)),n=length(d.cusum)),by=sag]
  
  #remove new points after last voxel
  newSegments <- newSegments[!(id==n &is.na(id.voxel))]
  newSegments[,n:=length(d.cusum),by=sag]
  
  #assign voxel id to new points
  newID <- newSegments[!is.na(id.voxel),list(id.voxel=id.voxel,id=id,n.voxel=n.voxel,n=n)]
  tmp <- c(newID$id[2:nrow(newID)],NA)
  newID[,next.id:=tmp]
  newID[,n.seg:=next.id-id]
  newID[id==n,n.seg:=1]
  newIDList <- rep(newID$id.voxel,newID$n.seg)
  newSegments[,id.voxel:=newIDList]
  setkey(newSegments,sag,id)
  coordinates <- segSurface[id.voxel>1 &sag==midSag &type=='low',
                            list(sag=sag,id.voxel=id.voxel-1,pv.axial=p.axial,pv.cor=p.cor,
                                 nv.axial=axial,nv.cor=cor,v.cusum=d.cusum)]
  lastVoxels <- segSurface[id.voxel==n.voxel &sag==midSag &type=='low',
                           list(sag=sag,id.voxel=id.voxel,pv.axial=axial,pv.cor=cor,
                                nv.axial=NA,nv.cor=NA,v.cusum=d.cusum)]
  coordinates <- rbind(coordinates,lastVoxels)
  setkey(coordinates,sag,id.voxel)
  prevDist <- segSurface[sag==midSag &type=='low']$d.cusum
  coordinates[,p.d.cusum:=prevDist] 
  newSegments <- merge(newSegments,coordinates,by=c('sag','id.voxel'),all.x=T)
  #calculate coordinates for new points
  newSegments[id.seg>0,
              ':='(axial=round(((d.cusum-p.d.cusum)*(nv.axial-pv.axial)/(v.cusum-p.d.cusum)+pv.axial),2),
              cor=round(((d.cusum-p.d.cusum)*(nv.cor-pv.cor)/(v.cusum-p.d.cusum)+pv.cor),2))]
  #add first and last points
  newSegments[is.na(id.seg) &(id.voxel==1 |id.voxel==n.voxel),
              ':='(axial=as.double(pv.axial),cor=as.double(pv.cor))]
  newSegments <- newSegments[!is.na(cor),list(id.voxel=id.voxel,axial=axial,cor=cor,
                                              d.cusum=d.cusum),by=sag]
  newSegments[,':='(n=length(cor),id=seq(1,length(cor)),type='low'),by=sag]
  setkey(newSegments,sag,id)
  
  #add upper boundary
  tmp <- segments[sag==midSag &type=='up',list(sag=sag,id.voxel=id.voxel,axial=axial,cor=cor,
                                               d.cusum=d.cusum,n=n,id=id,type=type)]
  newSegments <- rbind(newSegments,tmp)
  #add number of cc sections
  newSegments[,n.section:=nrow(sectionCC)]
  return(newSegments)
}

#find the midline between the two surfaces
getCCMidline <- function(segments) {
  #re-arrange and calculate the mid point for each corresponding up and low points
  midline <- data.table(id=segments[type=="up"]$id,n=segments[type=="up"]$n,
                        u.sag=segments[type=="up"]$sag,u.cor=segments[type=="up"]$cor,
                        u.axial=segments[type=="up"]$axial,l.sag=segments[type=="low"]$sag,
                        l.cor=segments[type=="low"]$cor,l.axial=segments[type=="low"]$axial)
  midline <- midline[,list(sag=u.sag,type='mid',id.voxel=NA,axial=(u.axial+l.axial)/2,
                     cor=(u.cor+l.cor)/2,id=id,n=n)]
  segments <- rbind(segments,midline)
  return(segments)
}

#calcuate thickness of a curvie structure by sectionalize the inner and outer curve 
#and calculate the distance between the corresponding points on the inner and outer curves
#calculate thickness only for the midsagittal slice
calCCThicknessMidLine <- function(segments,voxelSize) {
  r.x <- voxelSize$r.x
  r.y <- voxelSize$r.y
  r.z <- voxelSize$r.z
  midSag <- voxelSize$midSag
  upAxial=segments[type=="up" &sag==midSag]$axial
  upCor=segments[type=="up"  &sag==midSag]$cor 
  lowAxial=segments[type=="low"  &sag==midSag]$axial
  lowCor=segments[type=="low"  &sag==midSag]$cor
  #find the points on the upper and lower surfaces that are closest to each point of the midline
  newSeg <- data.table()
  for (j in 1:nrow(segments[type=="mid" &sag==midSag])) {
    tmp <- data.table(n=segments[type=="mid" &sag==midSag]$n[j], 
                      id=segments[type=="mid" &sag==midSag]$id[j],
                      sag=segments[type=="mid" &sag==midSag]$sag[j], 
                      m.axial=segments[type=="mid" &sag==midSag]$axial[j],
                      m.cor=segments[type=="mid" &sag==midSag]$cor[j], 
                      s.axial=upAxial, s.cor=upCor, type="up")
    newSeg <- rbind(newSeg,tmp)
    tmp <- data.table(n=segments[type=="mid" &sag==midSag]$n[j], 
                      id=segments[type=="mid" &sag==midSag]$id[j],
                      sag=segments[type=="mid" &sag==midSag]$sag[j], 
                      m.axial=segments[type=="mid" &sag==midSag]$axial[j],
                      m.cor=segments[type=="mid" &sag==midSag]$cor[j], 
                      s.axial=lowAxial, s.cor=lowCor, type="low")
    newSeg <- rbind(newSeg,tmp)
  }
  newSeg[,d.m:=round(sqrt(((m.cor-s.cor)*r.y)^2+((m.axial-s.axial)*r.z)^2),3)]
  newSeg[,min.d.m:=min(d.m),by=list(type,id)]
  newSeg <- newSeg[d.m==min.d.m]
  #for duplicated points, get the points in the middle
  newSeg[,':='(count.p=seq(1,length(sag)),n.p=length(sag)),by=list(id,type)]
  middleP <- unique(newSeg[n.p>1,list(n=n, sag=sag, axial=mean(s.axial), cor=mean(s.cor)),
                           by=list(type,id)])
  newSegments <- newSeg[n.p==1,list(type=type, id=id, n=n, sag=sag, axial=s.axial, cor=s.cor)]
  newSegments <- rbind(newSegments,middleP)
  setkey(newSegments,sag,type,id)
  return(newSegments)
}

#calculate CC thickness after linking points on the upper and lower boundaries
calCCThickness <- function(newSegments,voxelSize) {
  r.x <- voxelSize$r.x
  r.y <- voxelSize$r.y
  r.z <- voxelSize$r.z
  ccThickness <- newSegments[type=='up',list(sag=sag,u.axial=axial,u.cor=cor,id=id)]
  tmp <- newSegments[type=='low',list(sag=sag,l.axial=axial,l.cor=cor,id=id)]
  ccThickness <- merge(ccThickness,tmp,by=c('sag','id'))
  ccThickness[,n:=length(id),by=sag]
  ccThickness[,thickness:=round(sqrt(((u.cor-l.cor)*r.y)^2+((u.axial-l.axial)*r.z)^2),3)]
  ccThickness[(id==1 |id==n) &thickness==0,thickness:=NA]
  return(ccThickness)
}

#find midline and assign CC anatomic divisions according to Witelson 1989
getCCDivisions <- function(newSegments,ccThickness) {
  #find midline
  ccThickness[,':='(m.axial=(u.axial+l.axial)/2,m.cor=(u.cor+l.cor)/2)]
  newSegments[,':='(id.voxel=NULL,d.cusum=NULL)]
  newSegments[,n:=length(axial),by=list(sag,type)]
  midLine <- ccThickness[,list(sag=sag,axial=m.axial,cor=m.cor,n=n,id=id,type="mid")]
  midLine[,n.section:=newSegments$n.section[1]]
  newSegments <- rbind(newSegments,midLine)
  #division
  minCor <- min(newSegments[type=="up"]$cor)
  maxCor <- max(newSegments[type=="up"]$cor)
  maxCorLow <- max(newSegments[type=="low"]$cor)
  cc3 <- round(2*(maxCor-minCor)/3) + minCor
  cc4 <- round((maxCor+minCor)/2)
  cc5 <- round((maxCor-minCor)/3) + minCor
  cc6 <- round((maxCor-minCor)/5) + minCor
  newSegments[,':='(cor2=maxCorLow,cor3=cc3,cor4=cc4,cor5=cc5,cor6=cc6)]
  cc1 <- round(mean(newSegments[type=="low" &cor==maxCorLow]$axial))
  newSegments[,axial1:=cc1]
  newSegments[axial<axial1 &cor<cor2 &cor>cc4,cc.type:="CC1"]
  newSegments[cor>=cor2,cc.type:="CC2"]
  newSegments[axial>=axial1 &cor<=cor2 &cor>cc3,cc.type:="CC3"]
  newSegments[axial>=axial1 &cor<=cor3 &cor>cc4,cc.type:="CC4"]
  newSegments[cor<=cor4 &cor>cc5,cc.type:="CC5"]
  newSegments[cor<=cor5 &cor>cc6,cc.type:="CC6"]
  newSegments[cor<=cor6,cc.type:="CC7"]
  newSegments[type!="mid",cc.type:=type]
  thickness <- ccThickness[,list(id=id,thickness=thickness)]
  newSegments <- merge(newSegments,thickness,by='id')
  newSegments[,':='(cor2=NULL,cor3=NULL,cor4=NULL,cor5=NULL,cor6=NULL,axial1=NULL)]
  setkey(newSegments,sag,type,id)
  return(newSegments)
}

#calculate corpus callosum thickness by calculating the shortest distance 
#from each voxel at the outer surface to the inner surface
#input: inner surface, outer surface, and voxel size
#coronal slice number increases from the back to front of the brain
#axial slice number increases from the bottom to the top of the brain
calCCThicknessUpSurface <- function(surface,voxelSize) {
  r.x <- voxelSize$r.x
  r.y <- voxelSize$r.y
  r.z <- voxelSize$r.z
  ccThickness <- data.table()
  #from outer surface to inner surface
  for (i in 1:nrow(surface[type=="up"])) {
    tmp <- data.table(o.sag=surface[type=="up"]$sag[i],o.axial=surface[type=="up"]$axial[i],
                      o.cor=surface[type=="up"]$cor[i],i.sag=surface[type=="low"]$sag,
                      i.axial=surface[type=="low"]$axial,i.cor=surface[type=="low"]$cor)
    tmp <- tmp[o.sag==i.sag]
    tmp[,thickness:=round(sqrt(((i.sag-o.sag)*r.x)^2+((i.cor-o.cor)*r.y)^2
                               +((i.axial-o.axial)*r.z)^2),3)]
    tmp[,min.thick:=min(thickness)]
    ccThickness <- rbind(ccThickness,tmp[thickness==min.thick][1])
  }
  ccThickness[,min.thick:=NULL]
  return(ccThickness)
}

#calculate corpus callosum thickness from each voxels at the outer surface to the inner surface
#and from each voxels at the inner surface to outer surface
#input: inner surface, outer surface, and voxel size
calCCThicknessBothSurface <- function(surface,voxelSize) {
  r.x <- voxelSize$r.x
  r.y <- voxelSize$r.y
  r.z <- voxelSize$r.z
  ccThickness <- data.table()
  #from outer surface to inner surface
  for (i in 1:nrow(surface[type=="out"])) {
    tmp <- data.table(o.sag=surface[type=="out"]$sag[i],o.axial=surface[type=="out"]$axial[i],
                      o.cor=surface[type=="out"]$cor[i],i.sag=surface[type=="in"]$sag,
                      i.axial=surface[type=="in"]$axial,i.cor=surface[type=="in"]$cor)
    tmp <- tmp[o.sag==i.sag]
    tmp[,thickness:=round(sqrt(((i.sag-o.sag)*r.x)^2+((i.cor-o.cor)*r.y)^2
                               +((i.axial-o.axial)*r.z)^2),3)]
    tmp[,min.thick:=min(thickness)]
    ccThickness <- rbind(ccThickness,tmp[thickness==min.thick][1])
  }
  
  #from inner surface to out surface
  for (i in 1:nrow(surface[type=="in"])) {
    tmp <- data.table(o.sag=surface[type=="out"]$sag,o.axial=surface[type=="out"]$axial,
                      o.cor=surface[type=="out"]$cor,i.sag=surface[type=="in"]$sag[i],
                      i.axial=surface[type=="in"]$axial[i],i.cor=surface[type=="in"]$cor[i])
    tmp <- tmp[o.sag==i.sag]
    tmp[,thickness:=round(sqrt(((i.sag-o.sag)*r.x)^2+((i.cor-o.cor)*r.y)^2
                               +((i.axial-o.axial)*r.z)^2),3)]
    tmp[,min.thick:=min(thickness)]
    ccThickness <- rbind(ccThickness,tmp[thickness==min.thick][1])
  }
  ccThickness[,min.thick:=NULL]
  return(unique(ccThickness))
}

#remove the voxels at the split of the genu or splenium which is not continuous, does not work
removeCCSplits <- function(surface) {
  setkey(surface,type,sag,axial,cor)
  #assign sections for first 1/5 and last 1/5
  surface[,':='(min.cor=min(cor),max.cor=max(cor)),by=sag]
  surface[type=="low" &cor<min.cor+(max.cor-min.cor)/5, s.cor:=1]
  surface[type=="low" &cor>min.cor+4*(max.cor-min.cor)/5, s.cor:=2]
  surface[!is.na(s.cor),':='(n=length(cor),count=seq(1,length(cor))),by=list(sag,axial,s.cor)]
  trow <- nrow(surface)
  previousCor <- c(NA,surface$cor[1:(trow-1)])
  nextCor <- c(surface$cor[2:trow],NA)
  surface[,':='(prev.cor=previousCor,next.cor=nextCor)]
  surface[count==1,prev.cor:=NA]
  surface[count==n,next.cor:=NA]
  #detect continous voxels for the lower boundary
  if (nrow(surface[(s.cor==1 &cor-prev.cor>1) |(s.cor==2 &next.cor-cor>1)])>0) {
    tmp <- surface[s.cor>0 &cor-prev.cor>1,list(max.s.cor=max(cor),min.s.cor=min(cor)),by=list(sag,s.cor)]
    surface <- merge(surface,tmp,by=c('sag','s.cor'),all.x=T)
  }
  surface <- surface[is.na(n) |(s.cor==1 &n>1 &count==n) |(s.cor==2 &n>1 &count==1)]
  return(surface)
}

#find  tips in each section and duplicate them as up and lower surface
#along with the single voxel line, not used
duplicateTips <- function(segSurface) {
  if (nrow(segSurface[type=="tip"])>0) {
    #find the single voxel line with the same axial #
    tipAxialNum <- segSurface[type=='tip',list(sag=sag,s.cor=s.cor,tip.axial=axial,cor.tip=cor)]
    segSurface <- merge(segSurface,tipAxialNum,by=c('sag','s.cor'),all.x=T)
    setkey(segSurface,sag,s.cor,cor,axial)
    trow <- nrow(segSurface)
    previousAxial <- c(NA,segSurface$axial[1:(trow-1)])
    nextAxial <- c(segSurface$axial[2:trow],NA)
    segSurface[,':='(axial.prev=previousAxial,axial.next=nextAxial)]
    tmp <- segSurface[(s.cor==3 &axial==tip.axial &axial!=axial.prev+1 &axial!=axial.next-1)
                      |((s.cor==3 |s.cor==4) &axial==tip.axial &axial==axial.next)]
    #find cor # discontinuity point
    tmp[,':='(count=seq(1,length(type)),n=length(type)),by=list(sag,s.cor)]
    setkey(tmp,sag,s.cor,axial,cor)
    trow <- nrow(tmp)
    previousCor <- c(NA,tmp$cor[1:(trow-1)])
    nextCor <- c(tmp$cor[2:trow],NA)
    tmp[,':='(cor.prev=previousCor,cor.next=nextCor)]
    tmp[count==1,cor.prev:=NA]
    tmp[count==n,cor.next:=NA]
    if (nrow(tmp[cor.next-cor>1])>0) {
      maxCor <- tmp[cor.next-cor>1,list(max.cor=max(cor)),by=list(sag,s.cor)]
      tmp <- merge(tmp,maxCor,by=c('sag','s.cor'),all.x = T)
      tmp <- tmp[is.na(max.cor) |cor>max.cor]
      tmp[,max.cor:=NULL]
    }
    tips <- tmp[,list(sag=sag,s.cor=s.cor,axial=axial,cor=cor,dup=T)]
    segSurface[,':='(tip.axial=NULL,cor.tip=NULL,axial.prev=NULL,axial.next=NULL)]
    segSurface <- merge(segSurface,tips,by=c('sag','s.cor','axial','cor'),all.x=T)
    #duplicate
    tips <- segSurface[dup==T]
    segSurface[dup==T,type:="low"]
    tips[,type:="up"]
    segSurface <- rbind(segSurface,tips)
    segSurface[,dup:=NULL]
  }
}

#calculate fitted value for linear mixed models. Has one interaction and 8 coefficients
getlme2Fit <- function(fittedVal,v.lme.2,years) {
  fittedVal[,':='(b1=coef(v.lme.2)[[1]], b2=coef(v.lme.2)[[2]], b3=coef(v.lme.2)[[3]], 
                  b4=coef(v.lme.2)[[4]], b5=coef(v.lme.2)[[5]], b6=coef(v.lme.2)[[6]],
                  b7=coef(v.lme.2)[[7]], b8=coef(v.lme.2)[[8]])]
  #add age increase of 1 year
  tmp <- copy(fittedVal)
  tmp[,age.c:=age.c+years]
  fittedVal <- rbind(fittedVal,tmp)
  setkey(fittedVal,subject,age.c)
  #calculate
  fittedVal[status2==1,f.v:=round((b1 +b2*age.c +b5*age.c +b6*t1.v +b7*vScale +b8*age.c*t1.v),3)]
  fittedVal[status2==2,f.v:=round((b1 +b2*age.c +b3 +b5*age.c +b6*t1.v +b7*vScale +b8*age.c*t1.v),3)]
  fittedVal[status2==3,f.v:=round((b1 +b2*age.c +b4 +b5*age.c +b6*t1.v +b7*vScale +b8*age.c*t1.v),3)]
  return(fittedVal)
}

#calculate fitted value for linear mixed models. Has 2 two way interactions and 10 coefficients
getlme2iFit <- function(fittedVal,v.lme.2,years) {
  fittedVal[,':='(b1=coef(v.lme.2)[[1]], b2=coef(v.lme.2)[[2]], b3=coef(v.lme.2)[[3]], 
                  b4=coef(v.lme.2)[[4]], b5=coef(v.lme.2)[[5]], b6=coef(v.lme.2)[[6]],
                  b7=coef(v.lme.2)[[7]], b8=coef(v.lme.2)[[8]], b9=coef(v.lme.2)[[9]],
                  b10=coef(v.lme.2)[[10]])]
  #add age increase of 1 year
  tmp <- copy(fittedVal)
  tmp[,age.c:=age.c+years]
  fittedVal <- rbind(fittedVal,tmp)
  setkey(fittedVal,subject,age.c)
  #calculate
  fittedVal[status2==1,f.v:=round((b1 +b2*age.c +b3*age.c +b6*t1.v +b7*vScale +b10*age.c*t1.v),3)]
  fittedVal[status2==2,f.v:=round((b1 +b2*age.c +b3*age.c +b4 +b6*t1.v +b7*vScale +b8*age.c +b10*age.c*t1.v),3)]
  fittedVal[status2==3,f.v:=round((b1 +b2*age.c +b3*age.c +b5 +b6*t1.v +b7*vScale +b9*age.c +b10*age.c*t1.v),3)]
  return(fittedVal)
}

#calculate fitted value for linear mixed models. Has three way interaction and 14 coefficients
getlme3Fit <- function(fittedVal,v.lme.3,years) {
  fittedVal[,':='(b1=coef(v.lme.3)[[1]], b2=coef(v.lme.3)[[2]], b3=coef(v.lme.3)[[3]], 
                  b4=coef(v.lme.3)[[4]], b5=coef(v.lme.3)[[5]], b6=coef(v.lme.3)[[6]],
                  b7=coef(v.lme.3)[[7]], b8=coef(v.lme.3)[[8]], b9=coef(v.lme.3)[[9]],
                  b10=coef(v.lme.3)[[10]], b11=coef(v.lme.3)[[11]], b12=coef(v.lme.3)[[12]],
                  b13=coef(v.lme.3)[[13]], b14=coef(v.lme.3)[[14]])]
  #add age increase of 1 year
  tmp <- copy(fittedVal)
  tmp[,age.c:=age.c+years]
  fittedVal <- rbind(fittedVal,tmp)
  setkey(fittedVal,subject,age.c)
  #calculate
  fittedVal[status2==1,f.v:=round((b1 +b2*age.c +b3*age.c +b6*t1.v +b7*vScale +b10*age.c*t1.v),3)]
  fittedVal[status2==2,f.v:=round((b1 +b2*age.c +b3*age.c +b4 +b6*t1.v +b7*vScale +b8*age.c
                                  +b10*age.c*t1.v +b11*t1.v +b13*age.c*t1.v),3)]
  fittedVal[status2==3,f.v:=round((b1 +b2*age.c +b3*age.c +b5 +b6*t1.v +b7*vScale +b9*age.c
                                  +b10*age.c*t1.v +b12*t1.v +b14*age.c*t1.v),3)]
  return(fittedVal)
}

#get change directions 
getChangeDirection <- function(fittedVal) {
  fittedVal[,':='(n=length(age.c),count=seq(1,length(age.c))),by=subject]
  changeVal <- fittedVal[count==1,list(subject=subject,fv.1=f.v)]
  tmp <- fittedVal[count==n,list(subject=subject,fv.n=f.v)]
  changeVal <- merge(changeVal,tmp,by='subject')
  changeVal[,change:=fv.n-fv.1]
  changeVal[,':='(fv.1=NULL,fv.n=NULL)]
  fittedVal <- merge(fittedVal,changeVal,by='subject')
  fittedVal[,d.change:=-1]
  fittedVal[change>0,d.change:=1]
  return(fittedVal)  
}


#find S shaped surface 
#determined by significant cubic relationship between axial slice # and coronal slice # 
#on each sagittal slice
#significant polynomial term for order of 3 in polynomial regression
getPolynomialCubicP <- function(numAxial,numCoronal) {
  shapeModel <- lm(numCoronal ~ poly(numAxial,3))
  coefs <- data.table(coef(summary(shapeModel)))
  return(coefs$`Pr(>|t|)`[4])
}
getPolynomialQuadP <- function(numAxial,numCoronal) {
  shapeModel <- lm(numCoronal ~ poly(numAxial,3))
  coefs <- data.table(coef(summary(shapeModel)))
  return(coefs$`Pr(>|t|)`[3])
}

#combine LGI, surface area, volume, and thickness from Freesurfer to a table
composeFreesurferOutput <- function(gyri, surfarea, corVolume, thickness, corticalName) {
  #get cortical data
  matLgi <- as.matrix(gyri)[, 8:41]
  matArea <- as.matrix(surfarea)[, 10:43]
  matVol <- as.matrix(corVolume)[, 8:41]
  matTh <- as.matrix(thickness)[, 7:40]
  dtGyri <- data.table()
  for (i in 1:nrow(corticalName)) {
    tmp <- gyri[,list(subject=subject, group=group, status=status, sex=sex, age=age, hemisphere=hemisphere,
                      id.area=i, region=corticalName$region[i], name=corticalName$name[i], tcv=tcv,
                      lgi=matLgi[,i], s.area=matArea[,i], v.cor=matVol[,i], th.cor=matTh[,i])]
    dtGyri <- rbind(dtGyri,tmp)
  }
  dtGyri[,':='(lgi=as.double(lgi),s.area=as.double(s.area),v.cor=as.double(v.cor),th.cor=as.double(th.cor))]
  return(dtGyri)
}

#arrange FS data from 68 cortical areas from the two hemispheres in columns
#returns a data frame
arrangeFSData4PLS <- function(dtGyri, corticalName) {
  setkey(dtGyri,id.area,hemisphere,sex,group,subject)
  plsData <- dtGyri[id.area==1 &hemisphere=='left',
                    list(subject=subject,group=group,sex=sex,age=age)]
  #create cortical names with left and right
  newCorName <- copy(corticalName)
  newCorName[, nameLF:=paste('L.',name,sep='')]
  tmp <- copy(newCorName)
  tmp[, nameLF:=paste('R.',name,sep='')]
  newCorName <- rbind(newCorName,tmp)
  setkey(newCorName,id.area,nameLF)
  #form matrix for Fs data
  matLgi <-matrix(dtGyri$z.lgi, nrow=nrow(plsData), ncol=nrow(newCorName))
  colnames(matLgi) <- newCorName$nameLF
  #data.table does not work
  plsData <- as.data.frame(plsData)
  plsData$fs.data <- matLgi
  return(plsData)
}

#multiple linear regression on 34 freesurfer cortical areas, left and right separately
#record all models without interaction and models with significant interaction
fsLMGroupSexIntTCV <- function(dtFSData,corticalName,cutoff) {
  resultLM <- data.table()
  resultIntLM <- data.table()
  side <- c('left','right')
  for (i in 1:nrow(corticalName)) {
    #print(paste("Processing area:", i))
    for (j in 1:2) {
      fs.lm.1 <- lm(fs~age.c + f.sex + status + tcv, data=dtFSData[id.area==i &hemisphere==side[j]])
      fResult <- round(summary(fs.lm.1)$fstatistic,3)
      fPValue=round(pf(fResult[1],fResult[2],fResult[3],lower=FALSE),3)
      result <- round(summary(fs.lm.1)$coefficients,3)
      tmp <- data.table(hemisphere=side[j], id.area=i, name=corticalName[id.area==i]$name, 
                        df1=fResult[2], df2=fResult[3], f=fResult[1], f.p=fPValue,
                        type=c('age','sex','group','tcv'), coeff=result[2:5,1], se=result[2:5,2], 
                        t=result[2:5,3], p=result[2:5,4])
      resultLM <- rbind(resultLM,tmp)
      
      #with interaction
      fs.lm.2 <- lm(fs~age.c + f.sex * status + tcv, data=dtFSData[id.area==i &hemisphere==side[j]])
      fResult <- round(summary(fs.lm.2)$fstatistic,3)
      fPValue=round(pf(fResult[1],fResult[2],fResult[3],lower=FALSE),3)
      if (fPValue <= cutoff) {
        result <- round(summary(fs.lm.2)$coefficients,3)
        if (result[6,4] <= cutoff) {
          intResult <- anova(fs.lm.2,fs.lm.1)
          if (intResult[2,6] <= cutoff) {
            tmp <- data.table(hemisphere=side[j], id.area=i, name=corticalName[id.area==i]$name, 
                              df1=fResult[2], df2=fResult[3], f=fResult[1], f.p=fPValue,
                              type=c('age','sex','group','tcv','sex by group'), coeff=result[2:6,1], 
                              se=result[2:6,2], t=result[2:6,3], p=result[2:6,4])
            resultIntLM <- rbind(resultIntLM,tmp)
          }
        }
      }
    }
  }
  
  #mark significant results
  setkey(resultLM,type,id.area)
  resultLM[,fdr.f:=round(p.adjust(f.p,"BH"),3),by=type]
  resultLM[,fdr:=round(p.adjust(p,"BH"),3)]
  resultLM[,sig:=0]
  if (min(resultLM$fdr.f) <= cutoff) {
    resultLM[fdr.f<=cutoff &fdr<=cutoff,sig:=1]
  }
  setkey(resultLM,hemisphere,id.area,type)
  if (nrow(resultIntLM)>0) {
    resultIntLM[,fdr:=round(p.adjust(p,"BH"),3),by=type]
  }
  return(list(resultLM,resultIntLM))
}

#multiple linear regression on 34 freesurfer cortical areas, left and right separately
#record all models without interaction and models with significant interaction
#with brain scaling factor as a covariate
fsLMGroupSexIntBSF <- function(dtFSData,corticalName,cutoff) {
  resultLM <- data.table()
  resultIntLM <- data.table()
  side <- c('left','right')
  for (i in 1:nrow(corticalName)) {
    #print(paste("Processing area:", i))
    for (j in 1:2) {
      fs.lm.1 <- lm(fs~age.c + f.sex + status + vScale, data=dtFSData[id.area==i &hemisphere==side[j]])
      fResult <- round(summary(fs.lm.1)$fstatistic,3)
      fPValue=round(pf(fResult[1],fResult[2],fResult[3],lower=FALSE),3)
      result <- round(summary(fs.lm.1)$coefficients,3)
      tmp <- data.table(hemisphere=side[j], id.area=i, name=corticalName[id.area==i]$name, 
                        df1=fResult[2], df2=fResult[3], f=fResult[1], f.p=fPValue,
                        type=c('age','sex','group','vScale'), coeff=result[2:5,1], se=result[2:5,2], 
                        t=result[2:5,3], p=result[2:5,4])
      resultLM <- rbind(resultLM,tmp)
      
      #with interaction
      fs.lm.2 <- lm(fs~age.c + f.sex * status + vScale, data=dtFSData[id.area==i &hemisphere==side[j]])
      fResult <- round(summary(fs.lm.2)$fstatistic,3)
      fPValue=round(pf(fResult[1],fResult[2],fResult[3],lower=FALSE),3)
      if (fPValue <= cutoff) {
        result <- round(summary(fs.lm.2)$coefficients,3)
        if (result[6,4] <= cutoff) {
          intResult <- anova(fs.lm.2,fs.lm.1)
          if (intResult[2,6] <= cutoff) {
            tmp <- data.table(hemisphere=side[j], id.area=i, name=corticalName[id.area==i]$name, 
                              df1=fResult[2], df2=fResult[3], f=fResult[1], f.p=fPValue,
                              type=c('age','sex','group','vScale','sex by group'), coeff=result[2:6,1], 
                              se=result[2:6,2], t=result[2:6,3], p=result[2:6,4])
            resultIntLM <- rbind(resultIntLM,tmp)
          }
        }
      }
    }
  }
  
  #mark significant results
  setkey(resultLM,type,id.area)
  resultLM[,fdr.f:=round(p.adjust(f.p,"BH"),3),by=type]
  resultLM[,fdr:=round(p.adjust(p,"BH"),3)]
  resultLM[,sig:=0]
  if (min(resultLM$fdr.f) <= cutoff) {
    resultLM[fdr.f<=cutoff &fdr<=cutoff,sig:=1]
  }
  setkey(resultLM,hemisphere,id.area,type)
  if (nrow(resultIntLM)>0) {
    resultIntLM[,fdr:=round(p.adjust(p,"BH"),3),by=type]
  }
  return(list(resultLM,resultIntLM))
}

#separate age coefficient for male and female controls
getNCAgeCoefficients <- function(dtFSData,corticalName,cutoff) {
  resultLM <- data.table()
  side <- c('left','right')
  sex <- c(0,1)
  for (i in 1:nrow(corticalName)) {
    for (j in 1:2) {
      for (q in 1:2) {
        fs.lm.1 <- lm(fs~age.c + tcv, data=dtFSData[status==0 &f.sex==sex[q] &id.area==i &hemisphere==side[j]])
        fResult <- round(summary(fs.lm.1)$fstatistic,3)
        fPValue=round(pf(fResult[1],fResult[2],fResult[3],lower=FALSE),3)
        
        if (fPValue <= cutoff) {
          result <- round(summary(fs.lm.1)$coefficients,3)
          if (result[2,4] <= cutoff) {
            tmp <- data.table(hemisphere=side[j], id.area=i, name=corticalName[id.area==i]$name, f.sex=sex[q],
                              df1=fResult[2], df2=fResult[3], f=fResult[1], f.p=fPValue, type=c('age','tcv'),  
                              coeff=result[2:3,1], se=result[2:3,2], t=result[2:3,3], p=result[2:3,4])
            resultLM <- rbind(resultLM,tmp)
          }
        }
        
      }
    }
  }
  return(resultLM)
}

#age and sex coefficients for male and female controls for all regions no matter significant or not
getNCAgeSexCoefficients <- function(dtFSData,corticalName) {
  resultLM <- data.table()
  side <- c('left','right')
  sex <- c(0,1)
  for (i in 1:nrow(corticalName)) {
    for (j in 1:2) {
      fs.lm.1 <- lm(fs~age.c * f.sex + tcv, data=dtFSData[status==0 &id.area==i &hemisphere==side[j]])
      fResult <- round(summary(fs.lm.1)$fstatistic,3)
      fPValue=round(pf(fResult[1],fResult[2],fResult[3],lower=FALSE),3)
        
      result <- round(summary(fs.lm.1)$coefficients,3)
      tmp <- data.table(id.area=i, name=corticalName[id.area==i]$name, hemisphere=side[j], 
                        df1=fResult[2], df2=fResult[3], f=fResult[1], f.p=fPValue, 
                        age.coef=result[2,1], age.se=result[2,2], age.t=result[2,3], age.p=result[2,4],
                        sex.coef=result[3,1], sex.se=result[3,2], sex.t=result[3,3], sex.p=result[3,4],
                        tcv.coef=result[4,1], tcv.se=result[4,2], tcv.t=result[4,3], tcv.p=result[4,4],
                        int.coef=result[5,1], int.se=result[5,2], int.t=result[5,3], int.p=result[5,4])
      resultLM <- rbind(resultLM,tmp)
    }
  }
  resultLM[,':='(age.fdr=round(p.adjust(age.p,"BH"),3),sex.fdr=round(p.adjust(sex.p,"BH"),3))]
  return(resultLM)
}

#calculate z score based on control mean and SD after adjusting for age 
#separate coefficients for males and females
getZScoreFS <- function(dtFSData,ageSexCoef,corticalName) {
  tmp <- ageSexCoef[,list(id.area=id.area, hemisphere=hemisphere, age.coef=age.coef, sex.coef=sex.coef,
                          tcv.coef=tcv.coef, int.coef=int.coef)]
  dtFSData <- merge(dtFSData,tmp,by=c('id.area','hemisphere'))
  dtFSData[, adj.fs:=round((fs - age.c*age.coef - f.sex*sex.coef - age.c*f.sex*int.coef - tcv*tcv.coef), 3)]
  #control mean and SD
  meanSD <- dtFSData[status==0, list(mean.fs=round(mean(adj.fs),3), sd.fs=round(sd(adj.fs),3)),
                     by=list(f.sex,hemisphere,id.area)]
  dtFSData <- merge(dtFSData,meanSD,by=c('hemisphere','id.area','f.sex'))
  dtFSData[,z.fs:=round((adj.fs-mean.fs)/sd.fs, 3)]
  #dtFSData[,':='(coeff=NULL,mean.fs=NULL,sd.fs=NULL)]
  return(dtFSData)
}

#Bartlett's test of homogeneity of variances, more sensitive than Levene's test
getBartlettResult <- function(dtFSData,corticalName) {
  resultBartlett <- data.table()
  side <- c('left','right')
  sex <- c(0,1)
  for (i in 1:nrow(corticalName)) {
    for (j in 1:2) {
      for (s in 1:2) {
        result <- bartlett.test(adj.fs ~ fstatus, dtFSData[id.area==i &hemisphere==side[j] &f.sex==sex[s]])
        tmp <- data.table(id.area=i, name=corticalName[id.area==i]$name, hemisphere=side[j],
                          f.sex=sex[s], KSqrd=round(result$statistic,3), p=round(result$p.value,3))
        resultBartlett <- rbind(resultBartlett,tmp)
      }
    }
  }
  resultBartlett[,':='(fdr=round(p.adjust(p,"BH"),3))]
  setkey(resultBartlett,p,f.sex,id.area)
  return(resultBartlett)
}

#Levene's test of equality of variances
getLeveneResult <- function(dtFSData,corticalName) {
  resultLevene <- data.table()
  side <- c('left','right')
  sex <- c(0,1)
  for (i in 1:nrow(corticalName)) {
    for (j in 1:2) {
      for (s in 1:2) {
        result <- leveneTest(dtFSData[id.area==i &hemisphere==side[j] &f.sex==sex[s]]$adj.fs,
                             dtFSData[id.area==i &hemisphere==side[j] &f.sex==sex[s]]$fstatus,
                             location="median", correction.method="zero.correction")
        tmp <- data.table(id.area=i, name=corticalName[id.area==i]$name, hemisphere=side[j],
                          f.sex=sex[s], f=round(result$F[1],4), p=round(result$Pr[1],4))
        resultLevene <- rbind(resultLevene,tmp)
      }
    }
  }
  resultLevene[,':='(fdr=round(p.adjust(p,"BH"),4))]
  return(resultLevene)
}

#perform FDR only for areas showing aberrant gyrification
getPearsonFS <- function(dt4Pearson,corticalName,aberAreas) {
  side <- c('left','right')
  grp <- c('NC','PM')
  sexx <- c('m','f')
  resultPearson <-data.table()
  for (i in 1:nrow(corticalName)) {
    for (j in 1:2) {
      for (k in 1:2) {
        for (s in 1:2) {
          oneArea <- dt4Pearson[id.area==i &hemisphere==side[j] &group==grp[k] &sex==sexx[s]]
          pearsonTest <- cor.test(x=oneArea$v2, y=abs(oneArea$v1), method='pearson')
          tmp <- data.table(id.area=i, hemisphere=side[j], name=corticalName[id.area==i]$name, group=grp[k],
                            sex=sexx[s], r=round(pearsonTest$estimate,3), p=round(pearsonTest$p.value,4))
          resultPearson <- rbind(resultPearson,tmp)
        }
      }
    }
  }
  resultPearson <- merge(resultPearson,aberAreas,by=c('id.area','name','hemisphere','sex','group'),all.x=T)
  resultPearson[n.hyper>0 |n.hypo>0, bh:=round(p.adjust(p,'BH'),3), by=list(group)]
  setkey(resultPearson,group,sex,p)
  return(resultPearson)
}

getPearsonOrigFS <- function(dt4Pearson,corticalName,aberAreas) {
  side <- c('left','right')
  grp <- c('NC','PM')
  sexx <- c('m','f')
  resultPearson <-data.table()
  for (i in 1:nrow(corticalName)) {
    for (j in 1:2) {
      for (k in 1:2) {
        for (s in 1:2) {
          oneArea <- dt4Pearson[id.area==i &hemisphere==side[j] &group==grp[k] &sex==sexx[s]]
          pearsonTest <- cor.test(x=oneArea$v2, y=oneArea$v1, method='pearson')
          tmp <- data.table(id.area=i, hemisphere=side[j], name=corticalName[id.area==i]$name, group=grp[k],
                            sex=sexx[s], r=round(pearsonTest$estimate,3), p=round(pearsonTest$p.value,4))
          resultPearson <- rbind(resultPearson,tmp)
        }
      }
    }
  }
  resultPearson <- merge(resultPearson,aberAreas,by=c('id.area','name','hemisphere','sex','group'),all.x=T)
  resultPearson[n.hyper>0 |n.hypo>0, bh:=round(p.adjust(p,'BH'),3), by=list(group,sex)]
  setkey(resultPearson,group,sex,p)
  return(resultPearson)
}

getSpearmanOrigFS <- function(dt4Pearson,corticalName,aberAreas) {
  side <- c('left','right')
  grp <- c('NC','PM')
  sexx <- c('m','f')
  resultPearson <-data.table()
  for (i in 1:nrow(corticalName)) {
    for (j in 1:2) {
      for (k in 1:2) {
        for (s in 1:2) {
          oneArea <- dt4Pearson[id.area==i &hemisphere==side[j] &group==grp[k] &sex==sexx[s]]
          pearsonTest <- cor.test(x=oneArea$v2, y=oneArea$v1, method='spearman')
          tmp <- data.table(id.area=i, hemisphere=side[j], name=corticalName[id.area==i]$name, group=grp[k],
                            sex=sexx[s], r=round(pearsonTest$estimate,3), p=round(pearsonTest$p.value,4))
          resultPearson <- rbind(resultPearson,tmp)
        }
      }
    }
  }
  resultPearson[,bh:=round(p.adjust(p,'BH'),3),by=list(group,sex)]
  setkey(resultPearson,group,sex,p)
  return(resultPearson)
}

#t test comparing two groups for Freesurfer data
#v1 numeric, v2 binary
getTTestResult <- function(dt4TTest, corticalName) {
  side <- c('left','right')
  resultTTest <- data.table()
  for (i in 1:nrow(corticalName)) {
    #print(paste("Processing area:", i))
    for (j in 1:2) {
      oneArea <- dt4TTest[id.area==i &hemisphere==side[j]]
      tresult <- t.test(oneArea$v1~oneArea$v2)
      tmp <- data.table(id.area=i, hemisphere=side[j], name=corticalName[id.area==i]$name,
                        mean.0=round(tresult$estimate[1],3),  mean.1=round(tresult$estimate[2],3),
                        t=round(tresult$statistic,3), CIl=round(tresult$conf.int[1],3), 
                        CIu=round(tresult$conf.int[2],3), p=round(tresult$p.value,4))
      resultTTest <- rbind(resultTTest,tmp)
    }
  }
  resultTTest[,bh:=round(p.adjust(p,'BH'),3)]
  setkey(resultTTest,p)
  return(resultTTest)
}

#multiple linear regression between freesurfer data, FMR1, and IQ
#CGG and mRNA as independent variable and IQ scores ad dependent variable
#CGG and mRNA should be the 4th and 5th variable
getLMResult <- function(data4LM) {
  fs.lm.1 <- lm(dvar~age.c + ivar + tcv, data=data4LM)
  fResult <- round(summary(fs.lm.1)$fstatistic,3)
  fPValue=round(pf(fResult[1],fResult[2],fResult[3],lower=FALSE),3)
  result <- round(summary(fs.lm.1)$coefficients,3)
  lmResult <- data.table(df1=fResult[2], df2=fResult[3], f=fResult[1], f.p=fPValue, 
                            age.coef=result[2,1], age.se=result[2,2], age.t=result[2,3], age.p=result[2,4],
                            tcv.coef=result[4,1], tcv.se=result[4,2], tcv.t=result[4,3], tcv.p=result[4,4],
                            ivar.coef=result[3,1], ivar.se=result[3,2], ivar.t=result[3,3], ivar.p=result[3,4])
  return(lmResult)
}

#perform multiple linear regression separately for each group and sex
getLMResultGroupSex <- function(fsData4LM, IQ.cutoff) {
  fsName <- unique(fsData4LM$var.name)
  varGroup <- unique(fsData4LM$status)
  varSex <- unique(fsData4LM$sex)
  fsMLGResult <- data.table()
  for (i in 1:length(fsName)) {
    for (j in 1:length(varGroup)) {
      for (k in 1:length(varSex)) {
        #CGG for premutation group only
        if (varGroup[j]==1) {
          data4LM <- fsData4LM[var.name==fsName[i] &status==varGroup[j] &sex==varSex[k]  &!is.na(CGG),
                               list(subject,age.c,tcv,ivar=CGG,dvar=fs)]
          tmpResult <- getLMResult(data4LM)
          tmpResult[,':='(ivar='CGG', dvar=fsName[i], status=varGroup[j], sex=varSex[k])]
          fsMLGResult <- rbind(fsMLGResult,tmpResult)
        }
        
        #mRNA
        data4LM <- fsData4LM[var.name==fsName[i] &status==varGroup[j] &sex==varSex[k] &!is.na(mRNA),
                             list(subject,age.c,tcv,ivar=mRNA,dvar=fs)]
        tmpResult <- getLMResult(data4LM)
        tmpResult[,':='(ivar='mRNA', dvar=fsName[i], status=varGroup[j], sex=varSex[k])]
        fsMLGResult <- rbind(fsMLGResult,tmpResult)

        #VCI
        data4LM <- fsData4LM[var.name==fsName[i] &status==varGroup[j] &sex==varSex[k] &VCI>=IQ.cutoff,
                             list(subject,age.c,tcv,ivar=fs,dvar=VCI)]
        tmpResult <- getLMResult(data4LM)
        tmpResult[,':='(ivar=fsName[i], dvar='VCI', status=varGroup[j], sex=varSex[k])]
        fsMLGResult <- rbind(fsMLGResult,tmpResult)
        
        #PRI
        data4LM <- fsData4LM[var.name==fsName[i] &status==varGroup[j] &sex==varSex[k] &PRI>=IQ.cutoff,
                             list(subject,age.c,tcv,ivar=fs,dvar=PRI)]
        tmpResult <- getLMResult(data4LM)
        tmpResult[,':='(ivar=fsName[i], dvar='PRI', status=varGroup[j], sex=varSex[k])]
        fsMLGResult <- rbind(fsMLGResult,tmpResult)
        
        #WMI
        data4LM <- fsData4LM[var.name==fsName[i] &status==varGroup[j] &sex==varSex[k] &WMI>=IQ.cutoff,
                             list(subject,age.c,tcv,ivar=fs,dvar=WMI)]
        tmpResult <- getLMResult(data4LM)
        tmpResult[,':='(ivar=fsName[i], dvar='WMI', status=varGroup[j], sex=varSex[k])]
        fsMLGResult <- rbind(fsMLGResult,tmpResult)
        
        #PSI
        data4LM <- fsData4LM[var.name==fsName[i] &status==varGroup[j] &sex==varSex[k] &PSI>=IQ.cutoff,
                             list(subject,age.c,tcv,ivar=fs,dvar=PSI)]
        tmpResult <- getLMResult(data4LM)
        tmpResult[,':='(ivar=fsName[i], dvar='PSI', status=varGroup[j], sex=varSex[k])]
        fsMLGResult <- rbind(fsMLGResult,tmpResult)
      }
    }
  }
  
  return(fsMLGResult)
}
