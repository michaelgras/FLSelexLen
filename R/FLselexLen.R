# {{{
# sellen 
#
#' coverts harvest() and catch.sel() selectivity at length Sl  
#'
#' @param stock Input FLStockLen object or FLQuant.
#' @param nyears numbers of last years to compute selectivity
#' @param year option to specify year range, overwrites nyears
#' @param maturity if TRUE sellen() extract mat_at_length
#' @return FLQuant of Selectivity-at-length Sl  
#' @export
sellen <- function (stock, nyears=3,year=NULL,maturity=FALSE){

nyears = min(dims(stock)$maxyear-dims(stock)$minyear+1,nyears)  
  
if(is.null(year)){yr= (dims(stock)$"maxyear"-nyears+1):dims(stock)$"maxyear"} else {
yr = year  
}
if(class(stock)=="FLStockLen"){
  if(!maturity){flq = harvest(stock)[,ac(yr)]} else {flq =mat(stock)[,ac(yr)]}
} else { flq = stock[,ac(yr)]}

Sl = apply(flq,1,mean,na.rm=T)/max(apply(flq,1,mean,na.rm=T),na.rm=T)

Sl@units = "NA"
return(Sl)
}


#{{{
# FLStock2FLStockLen() 
#
#' Transforms an FLStock object into FLStockLen object using von Bertalanffy growth curve parameters
#'
#' @param object Input FLStock object.
#' @param vBPar Parameters of the von Bertalanffy growth curve
#' @param ndec number of decimal required in the FLStockLen object
#' @return an FLStockLen object
#' @export

FLStock2FLStockLen <- function (object, vBPar, ndec=2)
{
  ages <- object@range["min"]:object@range["max"]
    lengths <- vonB(object@range["min"]:object@range["max"], vBPar)

stkL <- FLStockLen(
  FLQuant(NA, 
    dimnames=list(
      len=round(vonB(object@range["min"]:object@range["max"], vBPar),ndec),
      year=object@range["minyear"]:object@range["maxyear"]), 
    units="cm")
  ) 

stkL@catch[] <- object@catch
stkL@catch.n[] <- object@catch.n
stkL@catch.wt[] <- object@catch.wt
stkL@discards[] <- object@discards
stkL@discards.n[] <- object@discards.n
stkL@discards.wt[] <- object@discards.wt
stkL@landings[] <- object@landings
stkL@landings.n[] <- object@landings.n
stkL@landings.wt[] <- object@landings.wt
stkL@stock[] <- object@stock
stkL@stock.n[] <- object@stock.n
stkL@stock.wt[] <-  object@stock.wt
stkL@m[] <- object@m
stkL@mat[] <- object@mat
stkL@harvest[] <- object@harvest
stkL@harvest.spwn[] <- object@harvest.spwn
stkL@m.spwn[] <- object@m.spwn
stkL@name[] <- object@name
stkL@desc[] <- object@desc
stkL@range[] <- stkL@range[c("min", "max", "minyear", "maxyear", "minfbar", "maxfbar")]
return(stkL)
}


#{{{
# s50len
#
#' approximates length-at-50-selectivity from Sl 
#'
#' @param Sl selectivity at length Sl = sellen(FLStockLen)
#' @return s50len 
#' @export 
#ok
s50len = function(Sl){ 
  Sl = as.data.frame(Sl)
  slo =which(Sl$data>0.5)[1]
  if(slo>1){
    bin = Sl[(slo-1):slo,]
    reg = an(lm(data~len,bin)$coef)
    S50len = max((0.5-reg[1])/reg[2],0.5)
  } else {
    S50len = Sl[1,"data"]/2  
  }
  return(S50len)
}


#{{{
# selexlenpars
#
#' computes initial values for selex pars 
#'
#' @param Sl selectivity at length Sl = sellen(stock)
#' @param S50 length-at-50-selectivty
#' @param S95 length-at-95-selectivty
#' @param Smax length at peak of halfnormal or top of descending slop
#' @param Dcv CV of halfnormal determining the steepness of the descending slope
#' @param Dmin height of the descending slop at maximum age
#'
#' @return vector of selectivity pars 
#' @export 
# 
selexlenpars <- function(Sl,S50=NULL,S95=NULL,Smax=NULL,Dcv=NULL,Dmin=NULL){
  S50proxy =s50len(Sl)
  sl = as.data.frame(Sl)
  
  len = sl$len
  sel = sl$data
  
  #selex.dat = data.frame(len=c((0:3)[which(!0:3 %in% len)], len),
  #                       sel=c(rep(0.01,3)[which(!0:3 %in% len)], sel))
  selex.dat = data.frame(len=len,
                         sel=sel)
  peak <- which(selex.dat[,2]>0.85)[1]
  peak <- (peak+1):dim(selex.dat)[1]
  neg=sel[peak]
  dcv = ifelse(length(neg)>1, abs(0.5+2*quantile(neg[-1]-neg[-length(neg)],0.2)[[1]]),0.3)
   
  if(is.null(S50))  S50=  S50proxy
  if(is.null(S95))  S95=an(quantile(c(S50proxy,selex.dat[peak[1],1]),0.9))
  if(is.null(Smax))  Smax=selex.dat[max(which(selex.dat[,2]>0.85)[1]),1]*1.1
  if(is.null(Dcv))  Dcv=0.2 #(1-sel[nrow(sa)])#/length(age[age>peak])
  if(is.null(Dmin))  Dmin=min(sel[nrow(sl)],0.7)
  
  pars = FLPar(S50=S50,S95=S95,Smax=Smax,Dcv=Dcv,Dmin=Dmin)
  return(pars)
}
# }}}



#{{{
# selexlen
#
#' computes selex curves 
#'
#' @param Sl selectivity at length Sl = sellen(stock)
#' @param pars selexpars S50, S95, Smax, Dcv, Dmin 
#' @return FLquants selex predictions
#' @export 

selexlen <- function(Sl,pars){
  sl = as.data.frame(Sl)
  len = sl$len
  sel = sl$data
  S50 = pars[[1]]
  S95 = pars[[2]]
  Smax =pars[[3]]
  Dcv =pars[[4]]
  Dmin =pars[[5]]
  
  selex.dat = data.frame(len=c((0:3)[which(!0:3 %in% len)], len),
                         sel=c(rep(0.01,3)[which(!0:3 %in% len)], sel))
  subs = which(selex.dat$len %in% len)
  
  
  psel_a = 1/(1+exp(-log(19)*(selex.dat$len-S50)/(S95-S50)))
  psel_b = dnorm(selex.dat$len, Smax,Dcv * Smax)/max(dnorm(selex.dat$len, Smax, Dcv * Smax))
  psel_c = 1+(Dmin-1)*(psel_b-1)/-1
  psel = ifelse(selex.dat$len>=max(Smax),psel_c,psel_a)
  #psel = ifelse(psel_a<0.95,psel_a,psel_c)
  #psel = psel/max(psel)
  #psel = pmin(psel,0.999)
  
  #resids = log(selex.dat$sel)-log(psel)
  
  #fits=data.frame(age=selex.dat$age[subs],obs=selex.dat$sel[subs],fit=psel[subs],logis=psel_a[subs],halfnorm=psel_b[subs],height=psel_c[subs])
  observed = fitted = logis = hnorm = height = Sl
  fitted[]=matrix(psel[subs])
  logis[]=matrix(psel_a[subs])
  hnorm[]=matrix(psel_b[subs])
  height[]=matrix(psel_c[subs])
  
  pred = FLQuants(observed,fitted,logis,hnorm,height)
  pred@names = c("observed","fitted","logis","hnorm","height")  
  
  
 
  return(pred)
}
# }}}




#{{{
#' fitselexlen
#'
#' fits selex selectivity function to L_at_age from FLStockLen
#'
#' @param Sl selectivity at length Sl = sellen(stock)
#' @param S50 init value length-at-50-selectivty
#' @param S95 init value length-at-95-selectivty
#' @param Smax init value length at peak of halfnormal or top of descending slop
#' @param Dcv init value CV of halfnormal determining the steepness of the descending slope
#' @param Dmin init value height of the descending slop at maximum length
#' @param CVlim sets upper CV bound
#' @param ogive fit logistic 
#' @return list with fit and FLquants selex predictions
#' @export 
fitselexlen <- function(Sl,S50=NULL,S95=NULL,Smax=NULL,Dcv=NULL,Dmin=NULL,CVlim=0.5,logistic=FALSE){
  
  pars =selexlenpars(Sl=Sl,S50=S50,S95=S95,Smax=Smax,Dcv=Dcv,Dmin=Dmin)
  
 
  
  imp = c(pars)
  imp[4] = CVlim/3
  lower = imp*0.3
  upper = imp*2
  upper[3] = max(as.data.frame(Sl)$len)*2
  lower[3] = min(as.data.frame(Sl)$len[which(as.data.frame(Sl)$data==1)])
  #upper[5] = max(min(1,imp[5]),0.2)
  upper[4] = CVlim
  lower[4] = 0.05
  lower[5]= 0.0001
  upper[5]= 0.75
  if(logistic){
    lower[3]= 100
    upper[3]= 200
    imp[3] = 150
  }
  # Likelihood
  jsel.ll = function(par=imp,data=Sl){
    Sl=data
    flimp = FLPar(S50=par[1],S95=par[2],Smax=par[3],Dcv=par[4],Dmin=par[5])
    pred= selexlen(Sl=Sl,pars=flimp)
   return(sum(((pred$observed)-(pred$fitted))^2))
  }
  
  
  fit = optim(par=imp, fn = jsel.ll,method="L-BFGS-B",lower=lower,upper=upper, data=Sl, hessian = TRUE)
  fit$par = FLPar(S50=fit$par[1],S95=fit$par[2],Smax=fit$par[3],Dcv=fit$par[4],Dmin=fit$par[5])
  fit$name
  fit$fits = selexlen(Sl,fit$par)
  return(fit)
}
# }}}


#{{{
#' lopt()
#'
#' Function to compute lopt, the length where an unfished cohort attains maximum biomass  
#' @param stock class FLStock
#' @return FLQuant with annual spr0y  
#' @export
#' @author Michael Gras & Henning Winker
lopt<-function(stock, nyears=3){
  object=stock
  len <- as.numeric(dimnames(stkL)$len) # dims(object)[["min"]]:dims(object)[["max"]]
  survivors <- exp(-apply(m(object),2,cumsum))
  survivors[-1] <- survivors[-dim(survivors)[1]]
  survivors[1] <- 1
  expZ <- exp(-object@m[dim(object@m)[1]]) # exp(-m(object[dim(m(object))[1]]))
  if (!is.na(range(object)["plusgroup"]))
  survivors[dim(m(object))[1]] <- survivors[dim(m(object))[1]]*(-1.0/(expZ-1.0))
  bl <- yearMeans(tail((stock.wt(object)*survivors)[-dims(object)[["max"]],],nyears))
  lopt = len[which(bl==max(bl))[1]]
  # Condition that at aopt fish had spawned 1 or more times on average
  lopt =  max(lopt,(which(yearMeans(tail(object@mat,nyears))>0.5)+1)[1])
  # ensure that aopt <= maxage
  lopt = min(lopt,dims(object)[["max"]]-1)
  
  return(lopt)
}
# }}}


#{{{
# varselexlen
#
#' function to dynamically vary selex parameters 
#'
#' @param pars 5 selex parameters of class FLPar 
#' @param stock optional stock object for tuning of age range  
#' @param step step size of change in one or several pars
#' @param amin start of S50
#' @param amax end of S50, required if stock = NULL
#' @param amax end of S50, required if stock = NULL
#' @param nyears end years for computing aopt() as automated amax limit
#' @param type option of selectivity change "crank","shift" or "dynamic"
#' @param return type of returned object FLPars or FLQuants 
#' @return selex parameters (FLPars) 
#' @export 
varselexlen <- function(pars, stock, step=0.1, lmin=NULL, lmax=NULL, nyears=3, type=c("crank","shift","dynamic"),return=c("Pars","Sl")) {

selpar=pars  
type = type[1]
return = return[1]
if(is.null(lmin)) lmin = round(selpar[[1]]*0.7,1)
if(type=="crank"){
  if(is.null(lmax)) lmax = round(selpar[[2]]*0.9,1)
}
if(type%in%c("shift","dynamic","selective")){
  if(is.null(lmax)) lmax = min(max(lopt(stock,nyears),selpar[[1]]),selpar[[1]]+6) 
  }
seqi = seq(lmin, lmax, step)
diff = seqi-selpar[[1]]
if(type=="crank"){
pars = FLPars(lapply(as.list(diff),function(x){
  out = selpar
  out[1]=out[1]+x
  out
}))  
}  
if(type=="dynamic"){
  pars = FLPars(lapply(as.list(diff),function(x){
    out = selpar
    ds = selpar[3]-selpar[2]
    out[1]= out[1]+x
    out[2] = max(selpar[2],out[1]+0.55)# *1.2
    out[3] = max(out[2]+ds,selpar[3])
    out
  }))  
}  
if(type=="shift"){
  pars = FLPars(lapply(as.list(diff),function(x){
    out = selpar
    out[1:3]=out[1:3]+x
    out
  }))  
}  
pars@names = paste0(seqi)

if(return=="Pars"){ 
  rtn = pars} else {
  rtn = FLQuants(lapply(pars,function(x){
    selexlen(sellen(stock),x)$fitted}))
  
}
return(rtn)
}
# }}}


#' par2sl()
#' 
#' Function to convert selex pars to Sa (FLQuants)  
#' @param pars selexpars 
#' @param object FLStock or FLQuant of Sa 
#' @param nyears number of end years for reference
#' @return FLStock
#' @export 
par2sl <- function(pars,object,nyears=3){
  if(class(object)=="FLStocklen") object=sellen(object,nyears)
  if(class(pars)=="FLPar") pars = FLPars(pars)
  out = FLQuants(lapply(pars,function(x)selex(object,x)$fitted))
  return(out)
}

