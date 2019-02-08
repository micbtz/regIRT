
# functions for simulations

sim_reg_nom<-function(r,n,param0,lambda,itmsel,cat0,nodes,weights,pen,adaptive)
{
  requireNamespace("mirt")
  requireNamespace("regIRT")
  requireNamespace("Hmisc")
  
  set.seed(r)
  abilities<-rnorm(n)
  dataset<-simdata(param0,abilities)
  # verify that there are all categories
  allcat<-sum(abs(apply(dataset,2,FUN=function(x) length(unique(x)))-cat0))==0
  if (allcat)
  {
    mod<-mirt::mirt(dataset, 1, 'nominal')
    par<-mirt::coef(mod)
    param<-lapply(par,FUN=repar)
    param<-param[1:4]
    par<-unlist(param)
    
    # collapse patterns
    rdata<-reduce_data(dataset)
    datared<-rdata$data
    numpatt<-rdata$numpatt
    
    # weights for adaptive lasso
    nitems<-ncol(dataset)
    if (adaptive)
    {
      w<-list()
      for(j in 1:nitems) {
        m<-length(param[[j]])/2
        alpha<-c(0,param[[j]][1:m])
        w[[j]]<-1/abs(outer(alpha,alpha,"-"))
      }
      
    }
    else
    {
      w<-list()
      for(j in 1:nitems) {
        m<-length(param[[j]])/2
        w[[j]]<-matrix(1,m+1,m+1)
      }
    }
    # estimation at different levels of lambda
    if (pen=="lasso")
    {
      opt<-optim(par=par,fn=nominallik_fpenRcppA,gr=gradnominallik_fpenRcppA,data=datared,method="BFGS",lambda=lambda[1],nodes=nodes,weights=weights,numpatt=numpatt,itemsselect=itmsel-1,control=list(maxit=500),eps=0.01,w=w)
      for (i in 3:9) opt<-optim(par=opt$par,fn=nominallik_fpenRcppA,gr=gradnominallik_fpenRcppA,data=datared,method="BFGS",lambda=lambda[1],nodes=nodes,weights=weights,numpatt=numpatt,itemsselect=itmsel-1,eps=10^-i,w=w)
      out_approx<-matrix(opt$par)
      for (l in 2:length(lambda)) {
        #print(l)
        opt<-optim(par=out_approx[,l-1],fn=nominallik_fpenRcppA,gr=gradnominallik_fpenRcppA,data=datared,method="BFGS",lambda=lambda[l],nodes=nodes,weights=weights,numpatt=numpatt,itemsselect=itmsel-1,eps=0.01,w=w)
        for (i in 3:10) opt<-optim(par=opt$par,fn=nominallik_fpenRcppA,gr=gradnominallik_fpenRcppA,data=datared,method="BFGS",lambda=lambda[l],nodes=nodes,weights=weights,numpatt=numpatt,itemsselect=itmsel-1,eps=10^-i,w=w)
        out_approx<-cbind(out_approx,opt$par)
      }
    out<-out_approx
    }
    if (pen=="ridge")
    {
      opt<-optim(par=par,fn=nominallikRcppA,gr=gradnominallikRcppA,data=datared,method="BFGS",nodes=nodes,weights=weights,lambda=lambda[1],numpatt=numpatt,itemsselect=itmsel-1,control=list(maxit=500))
      out<-matrix(opt$par)
      for (l in 2:length(lambda)) {
        print(l)
        opt<-optim(par=par,fn=nominallikRcppA,gr=gradnominallikRcppA,data=datared,method="BFGS",nodes=nodes,weights=weights,lambda=lambda[l],numpatt=numpatt,itemsselect=itmsel-1,control=list(maxit=500))
        out<-cbind(out,opt$par)
      }
    }
    
    ncat<-apply(dataset,2,max)+1
    npar<-ncat-1
    ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*2),ncat=ncat))
    
    # number of unique parameters
    k<-rep(0,ncol(out)) 
    for (j in 1:ncol(dataset))
    {
      outj<-out[ind==j,]
      if (nrow(outj)>2) {
        outj<-outj[1:(nrow(outj)/2),]
        np<-apply(outj,2,FUN=function(x) length(unique(round(x,3)[round(x,3)!=0])))
        np<-np+nrow(outj)
      }
      else np<-rep(2,ncol(out))
      k<-k+np
    }
    
    # likelihood
    ll<-c()
    for (l in 1:length(lambda)) {
      if (pen=="lasso") ll<-c(ll,nominallik_fpenRcppA(par=out[,l],data=datared,lambda=lambda[l],nodes=nodes,weights=weights,numpatt=numpatt,itemsselect=itmsel-1,eps=0.00,w=w))
      if (pen=="ridge") ll<-c(ll,nominallikRcppA(par=out[,l],data=datared,nodes=nodes,weights=weights,lambda=lambda[l],numpatt=numpatt,itemsselect=itmsel-1))
    }
    
    bic<-log(n)*k+2*ll
    aic<-2*k+2*ll
    
    # cross-validation
    cvres<-nominalCV(data=dataset,K=10,par=par,lambda=lambda,nodes=nodes,weights=weights,items.select=itmsel,pen=pen,w=w)
  }
  else out<-ll<-bic<-aic<-cvres<-NULL
  # results
  list(out=out,ll=ll,bic=bic,aic=aic,cvres=cvres,allcat=allcat)
  
}




sim_proxgrad<-function(r,n,param0,lambda,itmsel,cat0,nodes,weights,adaptive)
{
  requireNamespace("mirt")
  requireNamespace("regIRT")
  requireNamespace("Hmisc")
  requireNamespace("genlasso")
  requireNamespace("Matrix")
  
  set.seed(r)
  abilities<-rnorm(n)
  dataset<-simdata(param0,abilities)
  # verify that there are all categories
  allcat<-sum(abs(apply(dataset,2,FUN=function(x) length(unique(x)))-cat0))==0
  if (allcat)
  {
    mod<-mirt::mirt(dataset, 1, 'nominal')
    par<-mirt::coef(mod)
    param<-lapply(par,FUN=repar)
    param<-param[1:4]
    par<-unlist(param)
    
    # collapse patterns
    rdata<-reduce_data(dataset)
    datared<-rdata$data
    numpatt<-rdata$numpatt
    
    # weights for adaptive lasso
    nitems<-ncol(dataset)
    if (adaptive)
    {
      w<-list()
      for(j in 1:nitems) {
        m<-length(param[[j]])/2
        alpha<-c(0,param[[j]][1:m])
        w[[j]]<-1/abs(outer(alpha,alpha,"-"))
      }
      
    }
    else
    {
      w<-list()
      for(j in 1:nitems) {
        m<-length(param[[j]])/2
        w[[j]]<-matrix(1,m+1,m+1)
      }
    }
  
    parpg<-nominal_proxgr(par=par,datared=datared,weights=weights,nodes=nodes,numpatt=numpatt,lambda=lambda[1],s=0.01,w=w)
    out<-matrix(parpg)
    for (l in 2:length(lambda))
    {
      parpg<-nominal_proxgr(par=out[,l-1],datared=datared,weights=weights,nodes=nodes,numpatt=numpatt,lambda=lambda[l],s=0.01,w=w)
      out<-cbind(out,parpg)
    }
    
    # likelihood
    ll<-c()
    for (l in 1:length(lambda))
      ll<-c(ll,nominallik_fpenRcppA(par=out[,l],data=datared,lambda=lambda[l],nodes=nodes,weights=weights,numpatt=numpatt,itemsselect=itmsel-1,eps=0.00,w=w))
  
  }
  else out<-ll<-bic<-aic<-cvres<-NULL
  # results
  list(out=out,ll=ll,allcat=allcat)
}


