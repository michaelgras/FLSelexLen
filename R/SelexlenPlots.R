#{{{
# plotselexlen 
#
#' plots mean selectivity at length Sl across selected years 
#'
#' @param sel selexpars FLPars(s) or Sl FLquants or output from fitselexlen()  Sl
#' @param Sl observed selectivity-at-length (FLQuant) or FLStock 
#' @param obs show observations if TRUE
#' @param compounds option to show selex compounds
#' @param legend.title option for customized legend.title for discrete scenarios only
#' @param max.discrete determines max selectivities after which plot legend becomes continous 
#' 
#' @return plot  
#' @export
plotselexlen<- function(sel,Sl=NULL,obs=NULL,compounds=FALSE,colours=NULL,legend.title="S50",max.discrete=8){
  if(is.null(obs) & class(sel)=="FLQuants") obs=FALSE
  if(is.null(obs)) obs=TRUE
  if(class(sel)=="FLQuant") sel = FLQuants(sel)  
  if(class(Sl)=="FLStockLen") Sl = sellen(Sl)
  if(class(sel)=="FLQuants" & is.null(Sl)) Sl = sel[[1]]
  flq = ifelse(class(sel)=="FLQuants",TRUE,FALSE)
  object = sel
  if(is.null(colours)){colf = r4col} else {colf = colours}
  if(class(object)=="list"){
    pars = object$par
    Sl = object$fits$observed
  } else {
    pars=object  
  }
  
  if(flq){
    seldat=as.data.frame(pars)
    max.discrete=100
    compounds=FALSE
    p = ggplot(as.data.frame(seldat))+
      geom_line(aes(x=len,y=data,colour=qname))+geom_hline(yintercept = 0.5,linetype="dotted")+
      scale_colour_discrete(legend.title)+
      ylab("Selectivity")+xlab("Length")+
      #scale_x_continuous(breaks = 1:100)+scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
      scale_x_continuous(expand = c(0, 0),breaks = 1:100) + scale_y_continuous(expand = c(0, 0),limits=c(0,1.03),breaks = seq(0, 1, by = 0.25))
    
  } else {
    
    if(class(pars)=="FLPar"){
      pars = FLPars(pars)
      pars@names = paste0(round(pars[[1]][[1]],2))
    }
    if(is.null(Sl)){
      Sl = FLQuant(c(0.01,0.5,rep(1,ceiling(pars[[1]][[2]]*2)-2)),dimnames=list(len=1:ceiling(pars[[1]][[2]]*2)))  
      obs=FALSE
    }
    # predict
    pdat = FLQuant(0.5,dimnames=list(len=seq(dims(Sl)$min,dims(Sl)$max,0.05))) 
    pred = lapply(pars,function(x){
      selexlen(pdat,x)
    })
    
    if(length(pred)<max.discrete){
      
      if(compounds==TRUE & length(pred)==1){
        seldat = as.data.frame(pred[[1]][c(3:5,2)]) 
        cols=c(rainbow(3),"black")
      }
      if(compounds==FALSE& length(pred)==1){
        seldat = as.data.frame(pred[[1]][2])  
        cols=c("black")
      }
      if(length(pred)>1){
        seldat = FLQuants(lapply(pred,function(x){
          x[["fitted"]]}))
        compounds=FALSE
      }
      # Plot
      p = ggplot(as.data.frame(seldat))+
        geom_line(aes(x=len,y=data,colour=qname))+geom_hline(yintercept = 0.5,linetype="dotted")
      
      if(length(pred)==1){
        p=p + scale_color_manual("Selex",values=cols)
      } else {
        p = p +scale_colour_discrete(legend.title)
      }
      if(obs & length(pred)==1) p = p+geom_point(data=as.data.frame(Sl),aes(x=len,y=data), fill="white",shape=21,size=2)
      if(obs & length(pred)>1) p = p+geom_line(data=as.data.frame(Sl),aes(x=len,y=data),linetype="dashed")
      
      p = p +ylab("Selectivity")+xlab("Length")+
        #scale_x_continuous(breaks = 1:100)+scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
        scale_x_continuous(expand = c(0, 0),breaks = 1:100)+ scale_y_continuous(expand = c(0, 0),limits=c(0,1.03),breaks = seq(0, 1, by = 0.25))
      
    }
    
    
    if(length(pred)>=max.discrete){
      seldat = FLQuants(lapply(pred,function(x){
        x[["fitted"]]}))
      compounds=FALSE
      dat = as.data.frame(seldat)
      dat$S50 = an(as.character(dat$qname))
      if(obs){
        Sobs = as.data.frame(Sl)
        dat$ao = c(Sobs$len,rep(NA,nrow(dat)-nrow(Sobs)))
        dat$so = c(Sobs$data,rep(NA,nrow(dat)-nrow(Sobs)))
      }
      p = ggplot(data=dat,aes(x=len,y=data,group=S50))+    
        geom_line(aes(color=S50))+
        scale_color_gradientn(colours=rev(colf(20)))+
        ylab("Selectivity")+xlab("Length")+
        scale_x_continuous(expand = c(0, 0),breaks=1:100) + scale_y_continuous(expand = c(0, 0),limits=c(0,1.03))+
        
        theme(legend.key.size = unit(1, 'cm'), #change legend key size
              legend.key.height = unit(1, 'cm'),
              legend.text = element_text(size=7),
              legend.key.width = unit(0.6, 'cm'),
              legend.title=element_text(size=9)
        )
      if(obs) p = p+geom_line(aes(x=ao,y=so),linetype="dashed", na.rm=TRUE)
      
      
    }
  } # end of non flq loop
  return(p)  
}

# }}}


#{{{
# plotsellen
#
#' plots mean selectivity at length Sl across selected years 
#'
#' @param stock Input FLStockLen object.
#' @param nyears numbers of last years to compute selectivity
#' @param year specific years (will overwrite nyears)
#' @return Plot
#' @export

plotsellen<- function(stock,nyears=5,year=NULL){
  if(is.null(year)){yr= (range(stock)["maxyear"]-nyears+1):range(stock)["maxyear"]} else {
    yr = year } 
  Sl = as.data.frame(sellen(stock,nyears=nyears,year=year))
  p = ggplot(data=(as.data.frame(catch.selLen(stock[,ac(yr)]))),aes(x=len,y=data))+
    geom_line(aes(color = factor(year)))+ theme(legend.title=element_blank())+
    ylab("Selectivity")+xlab("Length")+geom_line(data=Sl,aes(x=len,y=data),size=1)
  return(p)  
}
# }}}


# to replace catch.sel that is defined only 
catch.selLen <- function(object) {
  return(harvest(object) %/% (apply(harvest(object), 2:6, max) + 1e-32))
}


#{{{
#' ploteqselex_2()
#
#' Isopleth plot of Selectivity vs F
#'
#' @param brps output from brp.selex() 
#' @param fit output from fitselex() 
#' @param Fmax upper possible limit of  F range
#' @param panels choice of plots 1:4
#' \itemize{
#'   \item 1  F over Relative Yield
#'   \item 2  F over SRP or SSB (requires ssr)
#'   \item 3  Isopleth Yield
#'   \item 4  Isopleth SPR or SSB (requires ssr)
#' }    
#' @param ncol number of columns
#' @param colours optional, e.g. terrain.col, rainbow, etc.
#' @param Ftrg  option to dynamic Ftrg=c("none","fmsy","f0.1")
#' \itemize{
#'   \item "none"  
#'   \item "msy"  
#'   \item "f0.1"
#' }
#' @return ggplot   
#' @export
ploteqselex_2 = function(brps, fit, Fmax=2., ageBased=TRUE, len50gear=NULL, 
                         vBPar, panels=NULL, ncol=NULL, colours=NULL, 
                         Ftrg=c("none","msy","f0.1"), stk_name = NULL){
  
  # Colour function
  if(is.null(colours)){colf = r4col} else {colf = colours}
  if(is.null(panels)) panels=1:4
  if(is.null(ncol)){
    if(length(panels)%in%c(1,3)){ncol=length(panels)} else {ncol=2}
  }
  # instead of observed we choose the brp that corresponds to fitted s50
  obsS <- names(brps)[an(names(brps)) >= round(an(fit$par[1]),1)][2]
  
  # add an extra point in the graph of isopleths (by adding an observed L50 of population
  # in the function)
  age50gear <- invVB(len50gear, vBPar)
  
  obsL50 <- rep(NA, length(len50gear))
  for (i in 1:length(age50gear))
  {
  obsL50[i] <- names(brps)[an(names(brps)) >= round(age50gear[i],1)][2]
  }

  # Check range 
  # WARNING: Also include na.rm = TRUE in order to get the limits (lim) in case the fmle did not converged all the times
  if(paste(brps[[1]]@model)[3]%in%c("a + ssb/ssb - 1")){
    pr = TRUE
    lim = min(Fmax,max(2*refpts(brps[[1]])["f0.1","harvest"],refpts(brps[[1]])["Fref","harvest"]*1.5,dims(brps[[1]])[["min"]]), na.rm = TRUE)
    quants = c("YPR","SPR")
    labs = c("Relative YPR",expression(SPR/SPR[0]))
    
  } else {
    pr = FALSE
    lim = min(Fmax,max(2*refpts(brps[[3]])["msy","harvest"],refpts(brps[[1]])["Fref","harvest"]*1.05,dims(brps[[1]])[["min"]]), na.rm = TRUE)
    quants = c("Yield","SSB")
    labs = c("Relative Yield",expression(SBB/SSB[0]))
  }
  
  # Prep some data for plotting
  fbar(brps[[1]]) = seq(0,lim,lim/101)[1:101]
  obs = data.frame(obs="obs",model.frame(metrics(brps[[1]],list(ssb=ssb, harvest=fbar, rec=rec, yield=landings)),drop=FALSE))
  obs[,8:11][obs[,8:11]<0] <- 0
  S50 = as.list(an(brps[-1]@names))
  isodat = do.call(rbind,Map(function(x,y){
    fbar(x) = seq(0,lim,lim/101)[1:101]
    mf =  model.frame(metrics(x,list(ssb=ssb, harvest=fbar, rec=rec, yield=landings)),drop=FALSE)
    data.frame(S50=y,as.data.frame(mf))  
  },brps[-1],S50))
  isodat$S50l <- round(vonB(isodat$S50, vBPar), 0)
  isodat$yield = isodat$yield#/max(isodat$yield)
  isodat[,8:11][isodat[,8:11]<0] <- 0
  Fobs = an(refpts(brps[[obsS]])["Fref","harvest"])
  Yobs = an(refpts(brps[[obsS]])["Fref","yield"])
  Sobs = an(refpts(brps[[obsS]])["Fref","ssb"])#/an(refpts(brps[[1]])["virgin","ssb"])
  # Sa=brps[[1]]@landings.sel/max(brps[[1]]@landings.sel)
  # S50obs = s50(Sa)
  S50obs = an(fit$par[1])
  S50lobs <- vonB(S50obs, vBPar)
  
  # Adding an extra point to the graph providing len50gear
  Flen50gear = an(refpts(brps[[obsL50[1]]])["Fref","harvest"])
  Ylen50gear = an(refpts(brps[[obsL50[1]]])["Fref","yield"])
  Slen50gear = an(refpts(brps[[obsL50[1]]])["Fref","ssb"])
  ## ================================================== ##
  isodat$Fo = c(obs$harvest,rep(NA,nrow(isodat)-nrow(obs)))
  isodat$Yo = c(obs$yield,rep(NA,nrow(isodat)-nrow(obs)))/max(isodat$yield)
  isodat$So = c(obs$ssb,rep(NA,nrow(isodat)-nrow(obs)))/max(isodat$ssb)
  
  if(Ftrg[1]!="none"){
    ftrg = do.call(rbind,Map(function(x,y){
      frp= c(refpts(x)[Ftrg,"harvest"])
      data.frame(s50=y,ftrg=frp)
    },brps[-1],S50))
    ftrg = ftrg[ftrg$ftrg<lim,]
  }
  ftrg$s50l <- vonB(ftrg$s50, vBPar)
  
  # F vs Yield
  if (ageBased) {
    P1 = ggplot(data=isodat,aes(x=harvest,y=yield/max(yield),group=S50))+
      geom_line(aes(color=S50))+geom_line(aes(x=Fo,y=Yo),size=0.7,linetype="dashed", na.rm=TRUE)+
      scale_color_gradientn(colours=rev(colf(20)))+ylab(labs[1])+
      geom_segment(aes(x = Fobs, xend = Fobs, y = 0, yend = Yobs/max(yield)), colour = "black",size=0.3,linetype="dotted")+
      geom_segment(aes(x = 0, xend = Fobs, y = Yobs/max(yield), yend = Yobs/max(yield)), colour = "black",size=0.3,linetype="dotted")+
      geom_point(aes(x=Fobs,y=Yobs/max(yield)),size=2)+
      xlab("Fishing Mortality")+
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.text = element_text(size=7),
            legend.key.width = unit(0.6, 'cm'),
            axis.title=element_text(size=10),
            legend.title=element_text(size=9)
      )+
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0),limits=c(0,1))
  } else {
    P1 = ggplot(data=isodat,aes(x=harvest,y=yield/max(yield),group=S50l))+
      geom_line(aes(color=S50l))+geom_line(aes(x=Fo,y=Yo),size=0.7,linetype="dashed", na.rm=TRUE)+
      scale_color_gradientn(colours=rev(colf(20)))+ylab(labs[1])+
      geom_segment(aes(x = Fobs, xend = Fobs, y = 0, yend = Yobs/max(yield)), colour = "black",size=0.3,linetype="dotted")+
      geom_segment(aes(x = 0, xend = Fobs, y = Yobs/max(yield), yend = Yobs/max(yield)), colour = "black",size=0.3,linetype="dotted")+
      geom_point(aes(x=Fobs,y=Yobs/max(yield)),size=2)+
      xlab("Fishing Mortality")+
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.text = element_text(size=7),
            legend.key.width = unit(0.6, 'cm'),
            axis.title=element_text(size=10),
            legend.title=element_text(size=9)
      )+
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0),limits=c(0,1))  
  }
  
  # F vs SSB  
  if (ageBased) {
    P2 = ggplot(data=isodat,aes(x=harvest,y=ssb/max(ssb),group=S50))+
      geom_line(aes(color=S50))+geom_line(aes(x=Fo,y=So),size=0.7,linetype="dashed", na.rm=TRUE)+
      scale_color_gradientn(colours=rev(colf(20)))+
      geom_segment(aes(x = Fobs, xend = Fobs, y = 0, yend = Sobs/max(ssb)), colour = "black",size=0.3,linetype="dotted")+
      geom_segment(aes(x = 0, xend = Fobs, y = Sobs/max(ssb), yend = Sobs/max(ssb)), colour = "black",size=0.3,linetype="dotted")+
      geom_point(aes(x=Fobs,y=Sobs/max(ssb)),size=2)+
      ylab(labs[2])+xlab("Fishing Mortality")+
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.key.width = unit(0.6, 'cm'),
            legend.text = element_text(size=7),
            axis.title=element_text(size=10),
            legend.title=element_text(size=9)
      )+ 
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0),limits=c(0,1))
  } else{
    P2 = ggplot(data=isodat,aes(x=harvest,y=ssb/max(ssb),group=S50l))+
      geom_line(aes(color=S50l))+geom_line(aes(x=Fo,y=So),size=0.7,linetype="dashed", na.rm=TRUE)+
      scale_color_gradientn(colours=rev(colf(20)))+
      geom_segment(aes(x = Fobs, xend = Fobs, y = 0, yend = Sobs/max(ssb)), colour = "black",size=0.3,linetype="dotted")+
      geom_segment(aes(x = 0, xend = Fobs, y = Sobs/max(ssb), yend = Sobs/max(ssb)), colour = "black",size=0.3,linetype="dotted")+
      geom_point(aes(x=Fobs,y=Sobs/max(ssb)),size=2)+
      ylab(labs[2])+xlab("Fishing Mortality")+
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.key.width = unit(0.6, 'cm'),
            legend.text = element_text(size=7),
            axis.title=element_text(size=10),
            legend.title=element_text(size=9)
      )+ 
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0),limits=c(0,1))
  }  
  # Isopleth plot Yield
  colbr = c(seq(0,0.6,0.2),seq(0.7,0.9,0.1),0.95,1)
  conbr = c(0,0.2,0.4,seq(0.5,0.9,0.1),0.95,1)
  nbr = length(colbr)
  
  if (ageBased) { 
    P3=ggplot(isodat, aes(x=harvest,y=S50))+
      geom_raster(aes(fill = yield/max(yield)), 
                  interpolate = T, hjust = 0.5, vjust = 0.5)+ 
      metR::geom_contour2(aes(z=yield/max(yield)),color = grey(0.4,1),breaks =conbr )+ 
      metR::geom_text_contour(aes(z=yield/max(yield)),stroke = 0.2,size=3,skip=0,breaks = conbr)+
      scale_fill_gradientn(colours=rev(colf(nbr+3))[-c(10:11,13)],limits=c(-0.03,1), breaks=colbr, name=paste(quants[1]))+
      geom_point(aes(x=Fobs,y=S50obs),size=2.5)+
      geom_segment(aes(x = Fobs, xend = Fobs, y = min(S50), yend = S50obs), colour = "black",size=0.3,linetype="dotted")+
      geom_segment(aes(x = 0, xend = Fobs, y = S50obs, yend = S50obs), colour = "black",size=0.3,linetype="dotted")+
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.key.width = unit(0.6, 'cm'),
            legend.text = element_text(size=7),
            axis.title=element_text(size=10),
            legend.title=element_text(size=9)
      )+
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(-0.03, 0))+
      ylab("Age-at-50%-Selectivity")+xlab("Fishing Mortality")
  } else
  {
    P3=ggplot(isodat, aes(x=harvest,y=S50l))+
      geom_raster(aes(fill = yield/max(yield)), 
                  interpolate = T, hjust = 0.5, vjust = 0.5)+ 
      metR::geom_contour2(aes(z=yield/max(yield)),color = grey(0.4,1),breaks =conbr )+ 
      metR::geom_text_contour(aes(z=yield/max(yield)),stroke = 0.2,size=3,skip=0,breaks = conbr)+
      scale_fill_gradientn(colours=rev(colf(nbr+3))[-c(10:11,13)],limits=c(-0.03,1), breaks=colbr, name=paste(quants[1]))+
      geom_point(aes(x=Fobs,y=S50lobs),size=2.5)+
      geom_segment(aes(x = Fobs, xend = Fobs, y = min(S50l), yend = S50lobs), colour = "black",size=0.3,linetype="dotted")+
      geom_segment(aes(x = 0, xend = Fobs, y = S50lobs, yend = S50lobs), colour = "black",size=0.3,linetype="dotted")+
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.key.width = unit(0.6, 'cm'),
            legend.text = element_text(size=7),
            axis.title=element_text(size=10),
            legend.title=element_text(size=9)
      )+
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(-0.03, 0))+
      ylab("Length-at-50%-Selectivity")+xlab("Fishing Mortality")
  }
  
  # Isopleth SSB
  colbr = c(0.05,seq(0,1,0.1))
  conbr = c(seq(0.05,0.4,0.05),seq(0.5,0.6,1),1)
  nbr = length(colbr)
  if (ageBased) {
    P4 = ggplot(isodat, aes(x=harvest,y=S50))+
      geom_raster(aes(fill = ssb/max(ssb)), 
                  interpolate = T, hjust = 0.5, vjust = 0.5)+ 
      metR::geom_contour2(aes(z=ssb/max(ssb)),color = grey(0.4,1),breaks =conbr )+
      metR::geom_text_contour(aes(z=ssb/max(ssb)),stroke = 0.2,size=3,skip=0,breaks = conbr)+
      scale_fill_gradientn(colours=rev(colf(nbr+4))[-c(4:7)],limits=c(-0.03,1), breaks=colbr, name=quants[2])+
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.key.width = unit(0.6, 'cm'),
            legend.text = element_text(size=7),
            axis.title=element_text(size=10),
            legend.title=element_text(size=9)
      )+
      geom_point(aes(x=Fobs,y=S50obs),size=2.5)+
      geom_segment(aes(x = Fobs, xend = Fobs, y = min(S50), yend = S50obs), colour = "black",size=0.3,linetype="dotted")+
      geom_segment(aes(x = 0, xend = Fobs, y = S50obs, yend = S50obs), colour = "black",size=0.3,linetype="dotted")+
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(-0.03, 0))+
      ylab("Age-at-50%-Selectivity")+xlab("Fishing Mortality")
  } else {
    P4 = ggplot(isodat, aes(x=harvest,y=S50l))+
      geom_raster(aes(fill = ssb/max(ssb)), 
                  interpolate = T, hjust = 0.5, vjust = 0.5)+ 
      metR::geom_contour2(aes(z=ssb/max(ssb)),color = grey(0.4,1),breaks =conbr )+
      metR::geom_text_contour(aes(z=ssb/max(ssb)),stroke = 0.2,size=3,skip=0,breaks = conbr)+
      scale_fill_gradientn(colours=rev(colf(nbr+4))[-c(4:7)],limits=c(-0.03,1), breaks=colbr, name=quants[2])+
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.key.width = unit(0.6, 'cm'),
            legend.text = element_text(size=7),
            axis.title=element_text(size=10),
            legend.title=element_text(size=9)
      )+
      geom_point(aes(x=Fobs,y=S50lobs),size=2.5)+
      geom_segment(aes(x = Fobs, xend = Fobs, y = min(S50l), yend = S50lobs), colour = "black",size=0.3,linetype="dotted")+
      geom_segment(aes(x = 0, xend = Fobs, y = S50lobs, yend = S50lobs), colour = "black",size=0.3,linetype="dotted")+
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(-0.03, 0))+
      ylab("Length-at-50%-Selectivity")+xlab("Fishing Mortality")
  }
  
  
  if(Ftrg[1]!="none"){
    if (ageBased) {
      P3 = P3+geom_point(data=ftrg,aes(x=ftrg,y=s50),shape = 21, colour = "black",size=1.1, fill = "white")
    } else {
      P3 = P3+geom_point(data=ftrg,aes(x=ftrg,y=s50l),shape = 21, colour = "black",size=1.1, fill = "white")
    }
    if (ageBased) {
      P4= P4+geom_point(data=ftrg,aes(x=ftrg,y=s50),shape = 21, colour = "black",size=1.1, fill = "white")
    } else {
      P4= P4+geom_point(data=ftrg,aes(x=ftrg,y=s50l),shape = 21, colour = "black",size=1.1, fill = "white")
    }
  }
  
  ## Adding the extra point to the graph
  if(!is.null(len50gear)){
    
    P3 = P3 + geom_point(data = data.frame(x=rep(Flen50gear, length(len50gear)), y=len50gear), aes(x = x, y = y), shape = 23, fill="white", size = 2)
    P4 = P4 + geom_point(data = data.frame(x=rep(Flen50gear, length(len50gear)), y=len50gear), aes(x = x, y = y), shape = 23, fill="white", size = 2)
  }
  plots <- list(P1=P1,P2=P2,P3=P3,P4=P4)
  
  if(length(panels)>1) return(gridExtra::grid.arrange(grobs =  plots[panels], ncol = ncol, top = stk_name))  
  if(length(panels)==1) return(plots[[panels]])  
}
