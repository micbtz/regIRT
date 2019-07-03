
# functions for simulations

sim_reg_nom<-function(r,n,param0,D,lambda,itmsel,cat0,pen,adaptive)
{
  requireNamespace("mirt")
  requireNamespace("regIRT")
  requireNamespace("MASS")
  
  set.seed(r)
  abilities<-MASS::mvrnorm(n=n,mu=rep(0,D),Sigma=diag(1,D))
  dataset<-regIRT::simdatanom(param0,abilities,D=D)
  # verify that there are all categories
  allcat<-sum(abs(apply(dataset,2,FUN=function(x) length(unique(x)))-cat0))==0
  if (allcat)
  {
    mod<-mirt::mirt(dataset, D, 'nominal')
    par<-mirt::coef(mod,rotation="varimax")
    param<-lapply(par,FUN=repar,D=D)
    param<-param[1:(length(param)-1)]
    parammat<-lapply(param,FUN=function(x) matrix(x,ncol=D+1))
    # multiply slopes for -1 if correlation with true par is negative (for each dimension)
    correl<- sign(diag(cor(do.call("rbind",param0), do.call("rbind",parammat))))
    param<-lapply(param, FUN=function(x) x*rep(correl,each=length(x)/(D+1)))
    par<-unlist(param)
    
    # small lasso penalty
    mod_lasso1<-nominalmod(data=dataset,D=D,parini=par,lambda=10^-3,pen="lasso",adaptive=FALSE)
    par_lasso1<-mod_lasso1$par

    mods<-nominalmod(data=dataset,D=D,parini=par_lasso1,parW=par_lasso1,lambda=lambda,pen=pen,adaptive=adaptive)
    out<-mods$par
    ll<-mods$lik
    conv<-mods$convergence

    nitems<-ncol(dataset)
    ncat<-apply(dataset,2,max)+1
    npar<-ncat-1
    ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))
    
    # number of unique parameters
    k<-rep(0,ncol(out)) 
    for (j in 1:nitems)
    {
      outj<-out[ind==j,]
      if (nrow(outj)>(D+1)) {
        ncatj<-ncat[j]-1
        outj<-outj[1:(ncatj*D),]
        seld<-rep(1:D,each=ncatj)
        np<-0
        for (d in 1:D) {
          np<-np+apply(outj[seld==d,],2,FUN=function(x) length(unique(round(x,3)[round(x,3)!=0])))
        }
        np<-np+ncatj #number of intercepts
      }
      else np<-rep(D+1,ncol(out))
      k<-k+np
    }
    
    bic<-log(n)*k+2*ll
    aic<-2*k+2*ll
    
    # cross-validation
    cvres<-nominalCV(mods,K=5)
  }
  else out<-ll<-conv<-bic<-aic<-cvres<-NULL
  # results
  list(out=out,ll=ll,conv=conv,bic=bic,aic=aic,cvres=cvres,allcat=allcat)
  
}

