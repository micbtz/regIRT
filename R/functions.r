
# probabilities of a nominal model
nomprobs<-function(param,nodes,D)
{
  #m<-length(param)
  param<-matrix(param,ncol=D+1)
  param<-rbind(0,param)
  alpha<-subset(param,select=1:D) # discriminations
  beta<-param[,D+1] # thresholds
  pr<-exp(alpha%*%t(nodes)+beta)
  #pr<-exp(outer(alpha,nodes)+beta)
  #pr<-rbind(1,pr)
  t(t(pr)/colSums(pr))
}


par2list<-function(par,data,D)
{
  nitems<-ncol(data)
  ncat<-apply(data,2,max,na.rm=TRUE)+1
  ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))
  param<-split(par,ind)
  names(param)<-colnames(data)
  param
}

# data responses should start from zero and have consecutive numbers
# likelihood of a nominal response model
# responses start from zero
# if D>1 (more than 1 dimension) nodes and weights are matrices obtained from expand.grid; weights contains the product of the elements
nominallik<-function(par,data,D=1,lambda=0,nodes,weights)
{
  patterns<-apply(data,1,paste, collapse ="_")
  tab<-table(patterns)
  first_patt<-!duplicated(patterns)
  datared<-data[first_patt,]
  patterns_unique<-patterns[first_patt]
  numpatt<-as.vector(tab[patterns_unique])

  n<-nrow(datared)
  #nq<-length(nodes)
  param<-par2list(par=par,data=datared,D=D)
  nitems<-ncol(datared)
  # ncat<-apply(datared,2,max,na.rm=TRUE)+1
  # ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))
  # param<-split(par,ind)

  probs<- lapply(param, FUN=nomprobs, nodes=nodes, D=D)
  logprobs<-lapply(probs,log)
  logprodj <- matrix(0, n, nrow(nodes))
  for (j in 1:nitems) {
    logprobsj <- logprobs[[j]]
    xj <- datared[, j]
    na.ind <- is.na(xj)
    logprobsj <- logprobsj[xj+1, ]
    if (any(na.ind))
      logprobsj[na.ind, ] <- 0
    logprodj <- logprodj + logprobsj
  }
  prodj<-exp(logprodj)
  sumq <- (prodj %*% weights)
  mlik <- -sum(log(sumq)*numpatt) # minus log-likelihood

  pen<-0
  if (lambda>0)
    pen<-sum(sapply(param,FUN=ridgepennominal,D=D)) # ridge-type penalization
  out <- mlik + lambda * pen
  #print(out)
  out
}


# ridge penalty
ridgepennominal<-function(param,D)
{
  param<-matrix(param,ncol=D+1)
  param<-rbind(0,param)
  alpha<-subset(param,select=1:D) # discriminations
  out<-0
  m<-nrow(param)
  for (d in 1:D)
  {
    for (k in 1:(m-1))
    {
      for (h in (k+1):m)
      {
        out<-out+(alpha[k,d]-alpha[h,d])^2
      }
    }
  }
  out
}


# derivative of ridge penalty
derridgepennominal<-function(param,D)
{
  param<-matrix(param,ncol=D+1)
  param<-rbind(0,param)
  alpha<-subset(param,select=1:D) # discriminations
  m<-nrow(param)
  out<-matrix(0,m,D+1)
  for (d in 1:D)
  {
    for (k in 1:m)
    {
      for (h in 1:m)
      {
        out[k,d]<- out[k,d]+2*(alpha[k,d]-alpha[h,d])
      }
    }
  }
  matrix(out[-1,],ncol=1)
}


# gradient of the likelihood of a nominal response model
gradnominallik<-function(par,data,D=1,lambda=0,nodes,weights)
{
  patterns<-apply(data,1,paste, collapse ="_")
  tab<-table(patterns)
  first_patt<-!duplicated(patterns)
  datared<-data[first_patt,]
  patterns_unique<-patterns[first_patt]
  numpatt<-as.vector(tab[patterns_unique])
  
  n<-nrow(datared)
  nq<-nrow(as.matrix(nodes))
  nitems<-ncol(datared)
  ncat<-apply(datared,2,max,na.rm=TRUE)+1
  # ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))
  # param<-split(par,ind)
  param<-par2list(par=par,data=datared,D=D)
  probs<- lapply(param, FUN=nomprobs, nodes=nodes, D=D)
  logprobs<-lapply(probs,log)

  logpr <- array(0, c(n, nq, nitems))
  for (j in 1:nitems) {
    logprobsj <- logprobs[[j]]
    xj <- datared[, j]
    na.ind <- is.na(xj)
    logprobsj <- logprobsj[xj+1, ]
    if (any(na.ind)) logprobsj[na.ind, ] <- 0
    logpr[,,j] <- logprobsj
  }
  prodj<-exp(apply(logpr,c(1,2),sum))
  sumq <- (prodj %*% weights)
  der<-matrix(0,length(par),1)
  for (j in 1:nitems) {
    logprbis<-logpr
    logprbis[,,j]<-0
    prodmj<-exp(apply(logprbis,c(1,2),sum)) #prod over items without item j
    xj <- datared[, j]
    na.ind <- is.na(xj)
    for (k in 1:(ncat[j]-1)) { # loop through threshold parameters
      probs_xj<-probs[[j]][xj+1,]
      derprob<-matrix(NA,n,nq)
      derprob[xj==k & !is.na(xj),]<-  probs_xj[xj==k & !is.na(xj),]*(1-probs_xj[xj==k & !is.na(xj),])
      derprob[xj!=k & !is.na(xj),]<- -probs_xj[xj!=k & !is.na(xj),]*rep(probs[[j]][k+1,],each=sum(xj!=k & !is.na(xj)))
      #derprob[xj!=k,]<- -probs_xj[xj!=k,]*matrix(probs[[j]][k+1,],nrow=1)[rep(1,sum(xj!=k)),]
      if (any(na.ind)) derprob[na.ind, ] <- 0
      sumq_numerator <- ((prodmj*derprob) %*% weights)
      der[k+(ncat[j]-1)*D+sum(ncat[0:(j-1)]-1)*(D+1)]<-sum(sumq_numerator/sumq*numpatt) #thresholds

      for (d in 1:D)
      {
        derprob_alphad<-  derprob*rep(nodes[,d],each=n)
        sumq_numerator <- ((prodmj*derprob_alphad) %*% weights)
        der[k+(ncat[j]-1)*(d-1)+sum(ncat[0:(j-1)]-1)*(D+1)]<-sum(sumq_numerator/sumq*numpatt) #discriminations
      }
    }
  }

  # penalization
  derpen<-0
  if (lambda>0)
    derpen<-unlist(lapply(param,FUN=derridgepennominal,D=D)) # ridge-type penalization
  out <- der - lambda * derpen
  -out # derivatives of minus log-likelihood
}




# likelihood of a nominal response model with fused lasso penalty on slope parameters
# responses start from zero
# lasso penalty approximated with quadratic approximation
nominallik_fpen<-function(par,data,D=1,lambda=0,nodes,weights,items.select=1:ncol(data),eps=0.001,alphaW)
{
  patterns<-apply(data,1,paste, collapse ="_")
  tab<-table(patterns)
  first_patt<-!duplicated(patterns)
  datared<-data[first_patt,]
  patterns_unique<-patterns[first_patt]
  numpatt<-as.vector(tab[patterns_unique])
  
  n<-nrow(datared)
  #nq<-length(nodes)
  nitems<-ncol(datared)
  # ncat<-apply(datared,2,max,na.rm=TRUE)+1
  # ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))
  # param<-split(par,ind)
  param<-par2list(par=par,data=datared,D=D)
  
  probs<- lapply(param, FUN=nomprobs, nodes=nodes, D=D)
  logprobs<-lapply(probs,log)
  logprodj <- matrix(0, n, nrow(nodes))
  for (j in 1:nitems) {
    logprobsj <- logprobs[[j]]
    xj <- datared[, j]
    na.ind <- is.na(xj)
    logprobsj <- logprobsj[xj+1, ]
    if (any(na.ind))
      logprobsj[na.ind, ] <- 0
    logprodj <- logprodj + logprobsj
  }
  prodj<-exp(logprodj)
  sumq <- (prodj %*% weights)
  mlik <- -sum(log(sumq)*numpatt) # minus log-likelihood
  
  pen<-0
  if (lambda>0)
    pen<-sum(mapply(fpen,param=param[items.select],alphaW=alphaW[items.select],eps=eps,D=D))
  out <- mlik + lambda * pen
  #print(out)
  out
}


# fused lasso penalty on slope parameters of the nominal model
fpen<-function(param,eps,alphaW,D)
{
  param<-matrix(param,ncol=D+1)
  param<-rbind(0,param)
  alpha<-subset(param,select=1:D) # discriminations
  out<-0
  m<-nrow(param)
  for (d in 1:D)
  {
    for (k in 1:(m-1))
    {
      for (h in (k+1):m)
      {
        if (is.null(alphaW)) out<-out+approxabs(alpha[k,d]-alpha[h,d],eps=eps)
        else out<-out+approxabs(alpha[k,d]-alpha[h,d],eps=eps)/abs(alphaW[k,d]-alphaW[h,d])
      }
    }
  }
  out
}

# derivative of fused penalization
derfpen<-function(param,eps,alphaW,D)
{
  param<-matrix(param,ncol=D+1)
  param<-rbind(0,param)
  alpha<-subset(param,select=1:D) # discriminations
  m<-nrow(param)
  ind<-1:m
  out<-matrix(0,m,D+1)
  for (d in 1:D)
  {
    for (g in 1:m) {
      for (k in ind[ind<g]) {
        absa<-approxabs(alpha[k,d]-alpha[g,d],eps)
        if (is.null(alphaW)) w<-1
        else w<-abs(alphaW[k,d]-alphaW[g,d])
        out[g,d]<-out[g,d]-(alpha[k,d]-alpha[g,d])/absa/w
      }
      for (h in ind[ind>g]) {
        absa<-approxabs(alpha[g,d]-alpha[h,d],eps)
        if (is.null(alphaW)) w<-1
        else w<-abs(alphaW[g,d]-alphaW[h,d])
        out[g,d]<-out[g,d]+(alpha[g,d]-alpha[h,d])/absa/w
      }
    }
  }
  matrix(out[-1,],ncol=1)
}


# gradient of the likelihood of a nominal response model with fused lasso penalty on slope parameters
gradnominallik_fpen<-function(par,data,D=1,lambda=0,nodes,weights,items.select=1:ncol(data),eps=0.001,alphaW)
{
  patterns<-apply(data,1,paste, collapse ="_")
  tab<-table(patterns)
  first_patt<-!duplicated(patterns)
  datared<-data[first_patt,]
  patterns_unique<-patterns[first_patt]
  numpatt<-as.vector(tab[patterns_unique])
  
  n<-nrow(datared)
  nq<-nrow(as.matrix(nodes))
  nitems<-ncol(datared)
  ncat<-apply(datared,2,max,na.rm=TRUE)+1
  # ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))
  # param<-split(par,ind)
  param<-par2list(par=par,data=datared,D=D)
  probs<- lapply(param, FUN=nomprobs, nodes=nodes, D=D)
  logprobs<-lapply(probs,log)
  
  logpr <- array(0, c(n, nq, nitems))
  for (j in 1:nitems) {
    logprobsj <- logprobs[[j]]
    xj <- datared[, j]
    na.ind <- is.na(xj)
    logprobsj <- logprobsj[xj+1, ]
    if (any(na.ind)) logprobsj[na.ind, ] <- 0
    logpr[,,j] <- logprobsj
  }
  prodj<-exp(apply(logpr,c(1,2),sum))
  sumq <- (prodj %*% weights)
  der<-matrix(0,length(par),1)
  for (j in 1:nitems) {
    logprbis<-logpr
    logprbis[,,j]<-0
    prodmj<-exp(apply(logprbis,c(1,2),sum)) #prod over items without item j
    xj <- datared[, j]
    na.ind <- is.na(xj)
    for (k in 1:(ncat[j]-1)) { # loop through threshold parameters
      probs_xj<-probs[[j]][xj+1,]
      derprob<-matrix(NA,n,nq)
      derprob[xj==k & !is.na(xj),]<-  probs_xj[xj==k & !is.na(xj),]*(1-probs_xj[xj==k & !is.na(xj),])
      derprob[xj!=k & !is.na(xj),]<- -probs_xj[xj!=k & !is.na(xj),]*rep(probs[[j]][k+1,],each=sum(xj!=k & !is.na(xj)))
      #derprob[xj!=k,]<- -probs_xj[xj!=k,]*matrix(probs[[j]][k+1,],nrow=1)[rep(1,sum(xj!=k)),]
      if (any(na.ind)) derprob[na.ind, ] <- 0
      sumq_numerator <- ((prodmj*derprob) %*% weights)
      der[k+(ncat[j]-1)*D+sum(ncat[0:(j-1)]-1)*(D+1)]<-sum(sumq_numerator/sumq*numpatt) #thresholds
      
      for (d in 1:D)
      {
        derprob_alphad<-  derprob*rep(nodes[,d],each=n)
        sumq_numerator <- ((prodmj*derprob_alphad) %*% weights)
        der[k+(ncat[j]-1)*(d-1)+sum(ncat[0:(j-1)]-1)*(D+1)]<-sum(sumq_numerator/sumq*numpatt) #discriminations
      }
    }
  }
  
  # penalization
  if (lambda>0)
    for (j in items.select) {
      derpen<-derfpen(param[[j]],eps=eps,alphaW=alphaW[[j]],D)
      der[(sum((ncat[0:(j-1)]-1)*(D+1))+1) : (sum((ncat[0:j]-1)*(D+1))) ]<-der[(sum((ncat[0:(j-1)]-1)*(D+1))+1) : (sum((ncat[0:j]-1)*(D+1))) ] - lambda * derpen
    }
  -der # derivatives of minus log-likelihood
}




# approximation of abs function
approxabs<-function(x,eps) sqrt(x^2+eps)




# responses start from zero
# likelihood of a nominal response model + penalization
# for ADMM algorithm
nominallikaug<-function(par,data,lambda=0,nodes,weights,gamma,v,c,items.select=1:ncol(data))
{
  n<-nrow(data)
  nq<-length(nodes)
  nitems<-ncol(data)
  # ncat<-apply(data,2,max,na.rm=TRUE)+1
  # ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*2),ncat=ncat))
  # param<-split(par,ind)
  param<-par2list(par=par,data=data,D=1)
  
  probs<- lapply(param, FUN=nomprobs, nodes=nodes)
  logprobs<-lapply(probs,log)
  logprodj <- matrix(0, n, nq)
  for (j in 1:nitems) {
    logprobsj <- logprobs[[j]]
    xj <- data[, j]
    na.ind <- is.na(xj)
    logprobsj <- logprobsj[xj+1, ]
    if (any(na.ind))
      logprobsj[na.ind, ] <- 0
    logprodj <- logprodj + logprobsj
  }
  prodj<-exp(logprodj)
  sumq <- (prodj %*% weights)
  mlik <- -sum(log(sumq)) # minus log-likelihood
  pen<-0
  if (lambda>0)
    pen<-sum(mapply(aug_lagr_pen,param[items.select],gamma[items.select],v[items.select],c))
  out <- mlik + lambda * pen
  out
}



# function for ADMM algorithm
aug_lagr_pen<-function(param,gamma,v,c)
{
  m<-length(param)
  sel<-upper.tri(matrix(0,m/2+1,m/2+1))
  alpha<-c(0,param[1:(m/2)]) # thresholds
  alphadiff<-outer(alpha,alpha,"-")
  out<-0
  out<-out+sum((v*(gamma-alphadiff))[sel])
  out<-out+sum((c/2*(gamma-alphadiff)^2)[sel])
  out
}


gradnominallikaug<-function(par,data,lambda=0,nodes,weights,gamma,v,c,items.select=1:ncol(data))
{
  n<-nrow(data)
  nq<-length(nodes)
  nitems<-ncol(data)
  ncat<-apply(data,2,max,na.rm=TRUE)+1
  # ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*2),ncat=ncat))
  # param<-split(par,ind)
  param<-par2list(par=par,data=data,D=1)
  probs<- lapply(param, FUN=nomprobs, nodes=nodes)
  logprobs<-lapply(probs,log)

  logpr <- array(0, c(n, nq, nitems))
  for (j in 1:nitems) {
    logprobsj <- logprobs[[j]]
    xj <- data[, j]
    na.ind <- is.na(xj)
    logprobsj <- logprobsj[xj+1, ]
    if (any(na.ind)) logprobsj[na.ind, ] <- 0
    logpr[,,j] <- logprobsj
  }
  prodj<-exp(apply(logpr,c(1,2),sum))
  sumq <- (prodj %*% weights)
  der<-matrix(0,length(par),1)
  for (j in 1:nitems) {
    logprbis<-logpr
    logprbis[,,j]<-0
    prodmj<-exp(apply(logprbis,c(1,2),sum)) #prod over items without item j
    xj <- data[, j]
    na.ind <- is.na(xj)
    for (k in 1:(ncat[j]-1)) { # loop through threshold parameters
      probs_xj<-probs[[j]][xj+1,]
      derprob<-matrix(NA,n,nq)
      derprob[xj==k,]<-  probs_xj[xj==k,]*(1-probs_xj[xj==k,])
      derprob[xj!=k,]<- -probs_xj[xj!=k,]*rep(probs[[j]][k+1,],each=sum(xj!=k))
      if (any(na.ind)) derprob[na.ind, ] <- 0
      sumq_numerator <- ((prodmj*derprob) %*% weights)
      der[k+ncat[j]-1+sum(ncat[0:(j-1)]-1)*2]<-sum(sumq_numerator/sumq) #thresholds

      derprob<-  derprob*rep(nodes,each=n)
      sumq_numerator <- ((prodmj*derprob) %*% weights)
      der[k+sum(ncat[0:(j-1)]-1)*2]<-sum(sumq_numerator/sumq) #discriminations
    }
  }

  # penalization
  if (lambda>0)
    for (j in items.select)
    {
      derpen<-der_aug_lagr_pen(param[[j]],gamma[[j]],v[[j]],c)
      der[(sum((ncat[0:(j-1)]-1)*2)+1) : (sum((ncat[0:j]-1)*2)) ]<-der[(sum((ncat[0:(j-1)]-1)*2)+1) : (sum((ncat[0:j]-1)*2)) ] - lambda * derpen
      
    }
  -der # derivatives of minus log-likelihood
}



der_aug_lagr_pen<-function(param,gamma,v,c)
{
  m<-length(param)
  sel<-upper.tri(matrix(0,m/2,m/2))
  alpha<-c(0,param[1:(m/2)]) # thresholds
  # beta<-c(0,param[(m/2+1):m]) # discriminations
  alphadiff<-outer(alpha,alpha,"-")
  # betadiff<-outer(beta,beta,"-")
  out<-rep(0,m)
  for (k in 1:(m/2))
  {
    # out[k]<- -sum(v[k,-(1:k)])+sum(v[-(k:(m/2)),k])-c*sum((gamma-alphadiff)[k,-(1:k)])+c*sum((gamma-alphadiff)[-(k:(m/2)),k])
    # out[k+m/2]<- -sum(u[k,-(1:k)])+sum(u[-(k:(m/2)),k])-c*sum((phi-betadiff)[k,-(1:k)])+c*sum((phi-betadiff)[-(k:(m/2)),k])
    out[k]<- -sum(v[k+1,-(1:(k+1))])+sum(v[-((k+1):(m/2+1)),k+1])-c*sum((gamma-alphadiff)[k+1,-(1:(k+1))])+c*sum((gamma-alphadiff)[-((k+1):(m/2+1)),k+1])
    # out[k+m/2]<- -sum(u[k+1,-(1:(k+1))])+sum(u[-((k+1):(m/2+1)),k+1])-c*sum((phi-betadiff)[k+1,-(1:(k+1))])+c*sum((phi-betadiff)[-((k+1):(m/2+1)),k+1])
  }
  out
}






# estimation of the parameters of a nominal model
# with optional penalty on the slope parameters
nominalmod<-function(data,D=1,parini,parW=NULL,lambda=0,pen=NULL,adaptive=NULL,items.select=1:ncol(data),nq=NULL)
{
  # if (is.vector(parW)) parW<-matrix(parW)
  # if (pen=="lasso" & adaptive)
  #   if (ncol(parW)!=length(lambda)) stop("the number of columns of parW should be equal to the length of lambda")
  
  nitems<-ncol(data)

  if (any(lambda>0) & is.null(pen)) stop("specify argument pen if lambda>0")
  if (!is.null(pen)) if (pen=="lasso" & is.null(adaptive)) stop("speficy argument adaptive if pen='lasso'")   
  
  # check that the lowest category is zero
  mincat<-apply(data,2,min,na.rm=TRUE)
  if (!all(mincat==0))
  { 
    warning("categories have been rescaled to start from zero")
    for (j in 1:nitems) data[,j]<-data[,j]-mincat[j]
  }
  # check that categories have consecutive numbers
  maxcat<-apply(data,2,max,na.rm=TRUE)
  ncat<-apply(data,2,FUN=function(x) sum(!is.na(unique(x))))
  for (j in 1:nitems)
  {
    if(maxcat[j]+1 != ncat[j])
    {
      warning("catetories have been rescaled to have consecutive numbers")
      tmp<-as.factor(data[,j])
      levels(tmp)<-0:(ncat[j]-1)
      data[,j]<-as.numeric(as.character(tmp))
    }
  }
  
  if (is.null(nq)) nq<-switch(as.character(D), '1'=61, '2'=31, '3'=15, '4'=9, '5'=7, 3)
  
  gq<-statmod::gauss.quad.prob(n=nq,dist="normal")
  nodes<-matrix(gq$nodes)
  weights<-gq$weights
  nodes_list<-c()
  for (d in 1:D) nodes_list[[d]]<-nodes
  nodes<-as.matrix(expand.grid(nodes_list))
  weights_list<-c()
  for (d in 1:D) weights_list[[d]]<-weights
  weights<-as.matrix(expand.grid(weights_list))
  weights<-apply(weights,1,prod)
  
  # ncat<-apply(data,2,FUN=function(x) sum(!is.na(unique(x))))
  # ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))

  rdata<-reduce_data(data)
  datared<-rdata$data
  numpatt<-rdata$numpatt
  datared[is.na(datared)]<- -999
  
  if (!is.null(pen))
  {
    if (pen=="lasso")
    {
      if (adaptive)
      {
        paramW<-par2list(par=parW,data=datared,D=D)
        #paramW<-split(parW,ind)
        alphaW<-vector("list",nitems)
        for(j in 1:nitems) {
          paramj<-matrix(paramW[[j]],ncol=D+1)
          paramj<-rbind(0,paramj)
          alphaj<-subset(paramj,select=1:D) # discriminations
          alphaW[[j]]<-alphaj
        }
      }
      else {
        alphaW<-vector("list",nitems)
        for (i in 1:nitems) alphaW[[i]]<-matrix(1,1,1)
      }
    }
  }
  
  parout<-c()
  lik<-c()
  conv<-c()
  nlambda<-length(lambda)
  for (l in 1:nlambda)
  {
    cat("lambda =", lambda[l],"\n")
    if (l==1) ini<-parini
    else ini<-parout[,l-1]
    if (lambda[l]==0)
      opt<-optim(par=ini,fn=nominallikRcppA,gr=gradnominallikRcppA,method="BFGS",
                 data=datared,D=D,nodes=nodes,weights=weights,lambda=0,numpatt=numpatt,
                 itemsselect=items.select-1,control=list(maxit=500))
    
    if (lambda[l]>0)
    {
      if (pen=="ridge")
        opt<-optim(par=ini,fn=nominallikRcppA,gr=gradnominallikRcppA,method="BFGS",
                   data=datared,D=D,nodes=nodes,weights=weights,lambda=lambda[l],numpatt=numpatt,
                   itemsselect=items.select-1,control=list(maxit=500))
  
      if (pen=="lasso")
      {
        opt<-optim(par=ini,fn=nominallik_fpenRcppA,gr=gradnominallik_fpenRcppA,data=datared,D=D,
                   method="BFGS",lambda=lambda[l],nodes=nodes,weights=weights,numpatt=numpatt,
                   itemsselect=items.select-1,control=list(maxit=500),eps=0.01,alphaW=alphaW,adaptive=adaptive)
        for (i in 3:9) 
          opt<-optim(par=ini,fn=nominallik_fpenRcppA,gr=gradnominallik_fpenRcppA,data=datared,D=D,
                     method="BFGS",lambda=lambda[l],nodes=nodes,weights=weights,numpatt=numpatt,
                     itemsselect=items.select-1,control=list(maxit=500),eps=10^-i,alphaW=alphaW,adaptive=adaptive)
      }
    }
    parout<-cbind(parout,opt$par)
    lik<-c(lik,-opt$value)
    conv<-c(conv,opt$convergence)
  }

  return(list(data=data, D=D, parini=parini, parW=parW, lambda=lambda, pen=pen, adaptive=adaptive, items.select=items.select, nq=nq, par=parout, lik=lik, convergence=conv))

}

# cross-validation
nominalCV<-function(object,K,trace=FALSE)
{
  data<-object$data
  D<-object$D
  parini<-object$parini
  parW<-object$parW
  lambda<-object$lambda
  pen<-object$pen
  adaptive<-object$adaptive
  items.select<-object$items.select
  nq<-object$nq
  
  n<-nrow(data)
  nitems<-ncol(data)
  categ<-apply(data,2,FUN=function(x) sort(unique(x)))
  
  if (is.null(nq)) nq<-switch(as.character(D), '1'=61, '2'=31, '3'=15, '4'=9, '5'=7, 3)

  gq<-statmod::gauss.quad.prob(n=nq,dist="normal")
  nodes<-matrix(gq$nodes)
  weights<-gq$weights
  nodes_list<-c()
  for (d in 1:D) nodes_list[[d]]<-nodes
  nodes<-as.matrix(expand.grid(nodes_list))
  weights_list<-c()
  for (d in 1:D) weights_list[[d]]<-weights
  weights<-as.matrix(expand.grid(weights_list))
  weights<-apply(weights,1,prod)
  
  # ncat<-apply(data,2,FUN=function(x) sum(!is.na(unique(x))))
  # ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))

  if (pen=="lasso") {
    if (adaptive)
    {
      paramW<-par2list(par=parW,data=data,D=D)
      # paramW<-split(parW,ind)
      alphaW<-vector("list",nitems)
      for(j in 1:nitems) {
        paramj<-matrix(paramW[[j]],ncol=D+1)
        paramj<-rbind(0,paramj)
        alphaj<-subset(paramj,select=1:D) # discriminations
        alphaW[[j]]<-alphaj
      }
    }
    else {
      alphaW<-vector("list",nitems)
      for (i in 1:nitems) alphaW[[i]]<-matrix(1,1,1)
    }
  }
  
  cond<-TRUE
  count<-0
  while (cond){ # repeat until all groups have all categories
    gr <- split(sample(n, n, replace=FALSE), as.factor(1:K)) # generation of subsets for CROSS-VALIDATION
    datared<-vector("list",K)
    cond1k<-rep(TRUE,K)
    for (k in 1:K) # k = subset
    {
      data_k<-data[-gr[[k]],] # training set
      rdata_k<-reduce_data(data_k)
      rdata_k$data[is.na(rdata_k$data)]<- -999
      datared[[k]]$training<-rdata_k
      datak<-data[gr[[k]],] # validation set
      rdatak<-reduce_data(datak)
      rdatak$data[is.na(rdatak$data)]<- -999
      datared[[k]]$validation<-rdatak
      categ_trk<-apply(datared[[k]]$training$data,2,FUN=function(x) sort(unique(x[x!= -999])))
      # categ_valk<-apply(datared[[k]]$validation$data,2,FUN=function(x) sort(unique(x)))
      if (identical(categ,categ_trk)) cond1k[k]<-FALSE
    }
    cond<-any(cond1k)
    count<-count+1
    if (count==100) cond<-FALSE
  }
  if(count<100)
  {
    # ===============================
    # not penalized
    # ===============================
    
    est_nopen<-matrix(NA,K,length(parini))
    lik_nopen<-rep(NA,K)
    
    for (k in 1:K) # k = subset
    {
      if (trace) print(k)
      o<-optim(par=parini,fn=nominallikRcppA,gr=gradnominallikRcppA,method="BFGS",data=datared[[k]]$training$data,D=D,nodes=nodes,weights=weights,numpatt=datared[[k]]$training$numpatt,itemsselect=items.select-1,control=list(maxit=500))
      est_nopen[k,]<-o$par
      lik_nopen[k]<-nominallikRcppA(par=o$par,data=datared[[k]]$validation$data,D=D,nodes=nodes,weights=weights,numpatt=datared[[k]]$validation$numpatt,itemsselect=items.select-1)
    }
    
    # ===============================
    # penalized
    # ===============================
    
    est_pen<-vector("list",length(lambda)) # estimates
    lik_pen<-vector("list",length(lambda)) # likelihood
    
    for (i in 1:length(lambda)) {
      est_pen[[i]]<-matrix(NA,K,length(parini))
      lik_pen[[i]]<-rep(NA,K)
      for (k in 1:K) # k = subset
      {
        if (trace) print(c(i,k))
        if (i>1) ini<-est_pen[[i-1]][k,]
        if (i==1) ini<-est_nopen[k,]
        if (pen=="lasso")
        {
          opt<-optim(par=ini,fn=nominallik_fpenRcppA,gr=gradnominallik_fpenRcppA,data=datared[[k]]$training$data,D=D,method="BFGS",lambda=lambda[i],nodes=nodes,weights=weights,numpatt=datared[[k]]$training$numpatt,itemsselect=items.select-1,control=list(maxit=500),eps=0.01,alphaW=alphaW,adaptive=adaptive)
          for (e in 3:10) opt<-optim(par=opt$par, fn=nominallik_fpenRcppA,gr=gradnominallik_fpenRcppA,data=datared[[k]]$training$data,D=D,method="BFGS",lambda=lambda[i],nodes=nodes,weights=weights,numpatt=datared[[k]]$training$numpatt,itemsselect=items.select-1,control=list(maxit=500),eps=10^-e,alphaW=alphaW,adaptive=adaptive)
        }
        if (pen=="ridge")
        {
          opt<-optim(par=ini,fn=nominallikRcppA,gr=gradnominallikRcppA,data=datared[[k]]$training$data,D=D,method="BFGS",nodes=nodes,weights=weights,lambda=lambda[i],numpatt=datared[[k]]$training$numpatt,itemsselect=items.select-1,control=list(maxit=500))
        }
        est_pen[[i]][k,]<-opt$par
        lik_pen[[i]][k]<-nominallikRcppA(par=opt$par,data=datared[[k]]$validation$data,D=D,nodes=nodes,weights=weights,numpatt=datared[[k]]$validation$numpatt,itemsselect=items.select-1)
      }
    }
    # select the maximum likelihood
    sel<-which.min(sapply(lik_pen,mean))
  }
  else est_pen<-lik_pen<-sel<-NULL
  
  return(list(est=est_pen,lik=lik_pen,lambda=lambda,sel=sel,lambdasel=lambda[sel],par=object$par,data=data,D=D))
}




# initial values (for ADMM algorithm)
vinit<-function(gamma,lambda)
{
  -lambda*sign(gamma)
}


# returns unique patterns of responses and relative frequencies
reduce_data<-function(data)
{
  patterns<-apply(data,1,paste, collapse ="_")
  tab<-table(patterns)
  first_patt<-!duplicated(patterns)
  datared<-as.matrix(data[first_patt,])
  patterns_unique<-patterns[first_patt]
  numpatt<-as.matrix(tab[patterns_unique])
  return(list(data=datared,numpatt=numpatt))
}


# optimization through augmented lagrangian (ADMM algorithm) 
nominal_grfused<-function(par,data,maxiter=1000,maxiter_innerloop=100,v=NULL,cinit,lambda,items.select=1:ncol(data),nodes,weights)
{
  rdata<-reduce_data(data)
  datared<-rdata$data
  numpatt<-rdata$numpatt
  datared[is.na(datared)]<- -999
  
  nitems<-ncol(datared)
  ncat<-apply(datared,2,max,na.rm=TRUE)+1
  npar<-ncat-1
  ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*2),ncat=ncat))
  param<-split(par,ind)
  if (is.null (v)) {
    #initial values for v
    v<-vector("list",nitems)
    for(j in 1:nitems) {
      m<-npar[j]
      v[[j]]<-matrix(0,m+1,m+1)
      alpha<-c(0,param[[j]][1:m])
      for(k in 1:(npar[j])) {
        for(h in (k+1):(npar[j]+1)) {
          v[[j]][k,h]<-vinit(alpha[k]-alpha[h],lambda)
        }
      }
    }
  }
  gamma<-vector("list",nitems)
  for (j in 1:nitems)
  {
    m<-npar[j]
    gamma[[j]]<- matrix(NA,m+1,m+1)
  }
  parold<-0 
  ck<-cinit
  par_gamma_old<-0
  crit1old<-0
  iter<-0
  noconv<-TRUE
  while (noconv) {
    #stp<-min(100,ck)
    stp<-ck
    iter_innerloop<-0
    noconv_innerloop<-TRUE
    while (noconv_innerloop) {
      
      #estimate gamma
      for(j in 1:nitems) {
        m<-npar[j]
        alpha<-c(0,param[[j]][1:m])
        for(k in 1:(npar[j])) {
          for(h in (k+1):(npar[j]+1)) {
            tmp <- -v[[j]][k,h]-ck*(-alpha[k]+alpha[h])
            if (abs(tmp)<lambda) gamma[[j]][k,h]<-0
            if (tmp< -lambda) gamma[[j]][k,h] <- -v[[j]][k,h]/ck+alpha[k]-alpha[h]+lambda/ck
            if (tmp>lambda) gamma[[j]][k,h] <- -v[[j]][k,h]/ck+alpha[k]-alpha[h]-lambda/ck
          }
        }
      }

      # opt_par<-optim(par=par,fn=nominallikaugRcppA,gr=gradnominallikaugRcppA,data=datared,
      #                method="BFGS",nodes=nodes,weights=weights,lambda=lambda,gamma=gamma,
      #                v=v,c=ck,numpatt=numpatt,items.select-1,control=list(reltol=1e-100/ck))
      # deriv<-gradnominallikaugRcppA(opt_par$par,data=datared,nodes=nodes,weights=weights,
      #                               lambda=lambda,gamma=gamma,v=v,c=ck,
      #                               numpatt=numpatt,itemsselect=items.select-1)
      # cat("deriv",sum(abs(deriv)),"\n")
      
      solut<-nleqslv::nleqslv(x=par,fn=gradnominallikaugRcppA,data=datared,method="Newton",nodes=nodes,
                     weights=weights,lambda=lambda,gamma=gamma,v=v,c=ck,numpatt=numpatt,
                     itemsselect=items.select-1,control=list(ftol=1e-8/(iter+1)))
      cat("deriv",sum(abs(solut$fvec)),"\n")
      
      
      #plot(par,solut$x)
      #crit<-mean(abs(gradnominallikaugRcppA(par=solut$x,data=data,nodes=nodes,weights=weights,lambda=lambda,gamma=gamma,phi=phi,v=v,u=u,c=c)))
      #if (crit<0.02) noconv<-FALSE
      param<-split(solut$x,ind)
      #print(crit)
      #par<-solut$x
      #noconv_innerloop<-FALSE
      iter_innerloop<-iter_innerloop+1
      cat("iter_innerloop",iter_innerloop,"\n")
      par_gamma_new<-unlist(c(solut$x,sapply(gamma,FUN=function(x) x[upper.tri(x)])))
      #tmp<-cbind(tmp,par_gamma_phi_new)
      #tmp<<-tmp
      crit_inner_loop<-max(abs(par_gamma_new-par_gamma_old))
      cat("crit_inner_loop",crit_inner_loop,"\n")
      if(crit_inner_loop<0.01) noconv_innerloop<-FALSE
      if (iter_innerloop==maxiter_innerloop) noconv_innerloop<-FALSE
      cat("noconv_innerloop",noconv_innerloop,"\n")
      par_gamma_old<-par_gamma_new
    }
    constr<-vector("list",nitems)
    # update v and u
    for (j in 1:nitems)
    {
      m<-npar[j]
      alpha<-c(0,param[[j]][1:m])
      alphadiff<-outer(alpha,alpha,"-")
      constr[[j]]<-(gamma[[j]]-alphadiff)[upper.tri(alphadiff)]
      v[[j]]<-v[[j]]+stp*(gamma[[j]]-alphadiff)
    }
    iter<-iter+1
    if (iter==maxiter) noconv<-FALSE
    parnew<-unlist(c(solut$x,sapply(gamma,FUN=function(x) x[upper.tri(x)]),
                     sapply(v,FUN=function(x) x[upper.tri(x)])))
    crit<-max(abs(parnew-parold))
    cat("crit",crit,"\n\n")
    if(crit<0.005 & iter>3) noconv<-FALSE
    parold<-parnew
    crit1<-sum(sapply(constr,FUN=function(x) sum(abs(x))))
    if (crit1>0.25*crit1old) ck<-ck*1.05
    crit1old<-crit1
    cat("ck",ck,"\n")
  }
  return(list(par=solut$x,v=v))
}


# reparameterization of the parameters of the nominal model from mirt to regIRT package
repar<-function(param,D=1)
{
  discrm<-param[1:D]
  param<-param[-(1:D)]
  m<-length(param)/2
  c(as.vector(outer(param[2:m],discrm)),param[2:m+m])
}


# optimization through proximal gradient
nominal_proxgr<-function(par,datared,nodes,weights,numpatt,lambda,s,w)
{
  nitems<-ncol(datared)
  ncat<-apply(datared,2,max,na.rm=TRUE)+1
  npar<-ncat-1
  ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*2),ncat=ncat))
  param<-split(par,ind)
  D<-list()
  for(j in 1:nitems) {
    m<-length(param[[j]])/2
    Dj<-c()
    for (k in 1:m)
    {
      Dk<-diag(1,m+1)
      for (i in 1:(m+1-k)) Dk[i,i+k]<- -1
      Dk<-Dk[1:(m+1-k),]
      Dj<-rbind(Dj,Dk)
    }
    Dj<-subset(Dj,select= -1)
    zeros<-matrix(0,nrow(Dj),m)
    Dj<-cbind(Dj,zeros)
    D[[j]]<-Dj
  }
  D<-as.matrix(Matrix::bdiag(D))
  
  W<-as.matrix(Matrix::bdiag(lapply(w,FUN=w2W)))
  D<-W%*%D

  par_t<-par
  noconv<-TRUE
  count<-0
  while(noconv)
  {
    grcpp<-gradnominallikRcppA(par=par_t,data=datared,nodes=nodes,weights=weights,numpatt=numpatt,itemsselect=1:ncol(datared)-1)
    par_t1<-par_t-s*grcpp
    out<-genlasso::genlasso(par_t1, D=D)
    par_t2<-coef(out,lambda=s*lambda)$beta
    par_t3<-par_t2+count/(count+3)*(par_t2-par_t)
    if (max(abs(par_t-par_t3))<0.001) noconv<-FALSE
    par_t<-par_t3
    plot(par,par_t2)
    count<-count+1
    print(count)
    print(par_t2[17])
    if (count==500) noconv<-FALSE
  }
  par_t2
}


# function for proximal gradient
w2W<-function(w)
{
  out<-c()
  m<-nrow(w)-1
  for (k in 1:m)
  {
    for (i in 1:(m+1-k)) out<-c(out,w[i,i+k])
  }
  diag(out)
}




# generate data from a nominal model
simdatanom<-function(param,abilities,D)
{
  probs<- lapply(param, FUN=nomprobs, nodes=abilities,D=D)
  sapply(probs,FUN=function(x) Hmisc::rMultinom(probs=t(x), m=1)-1)
}


# regularization path
regPath<-function(cvres)
{
  data<-cvres$data
  D<-cvres$D
  lambda<-cvres$lambda
  nitems<-ncol(data)
  ncat<-apply(data,2,max)+1
  npar<-ncat-1
  ind<-unlist(lapply(as.list(1:nitems),FUN=function(x,ncat) rep(x,each=(ncat[x]-1)*(D+1)),ncat=ncat))
  for (j in 1:nitems)
  {
    outj<-cvres$par[ind==j,]
    cl<-rep(1:(D+1),each=nrow(outj)/(D+1))
    plot(lambda,outj[1,],type="l",ylim=c(min(outj),max(outj)),xlab=expression(lambda),ylab="",main=colnames(data)[j])
    abline(v=cvres$lambdasel)
    for (i in 1:sum(ind==j))
      lines(lambda,outj[i,],col=cl[i])
  }
}





