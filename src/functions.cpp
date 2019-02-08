// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat nomprobsC(arma::vec param, arma::vec nodes)
{
  int m=param.size();
  int nq=nodes.size();
  arma::vec alpha=param.subvec(0,m/2-1);
  arma::vec beta=param.subvec(m/2,m-1);
  arma::mat pr=alpha*nodes.t();
  pr.each_col() += beta;
  pr=exp(pr);
  arma::vec uno(nq);
  uno.ones();
  arma::mat uno_pr = join_cols(uno.t(),pr);//rbind
  arma::mat prsum = sum(uno_pr,0);  
  uno_pr.each_row() /= prsum;
  return uno_pr;
}



double ridgepennominal(arma::vec param)
{
  int m=param.size();
  arma::vec alpha=arma::zeros(m/2+1);
  alpha.subvec(1,m/2)=param.subvec(0,m/2-1);
  // arma::vec beta=arma::zeros(m/2+1);
  // beta.subvec(1,m/2)=param.subvec(m/2,m-1);
  double out=0;
  for (int k=0;k<(m/2);k++)
    for (int h=(k+1);h<(m/2+1);h++)
      out += pow(alpha(k)-alpha(h),2);
  return out;
}


double fpenC(arma::vec param, double eps, arma::mat w)
{
  int m=param.size();
  arma::vec alpha=arma::zeros(m/2+1);
  alpha.subvec(1,m/2)=param.subvec(0,m/2-1);
  // arma::vec beta=arma::zeros(m/2+1);
  // beta.subvec(1,m/2)=param.subvec(m/2,m-1);
  double out=0;
  for (int k=0;k<(m/2);k++) {
    for (int h=(k+1);h<(m/2+1);h++) {
      out += pow(pow(alpha(k)-alpha(h),2)+eps,0.5)*w(k,h);
    }
  }
  return out;
}

arma::vec derfpen(arma::vec param, double eps, arma::mat w)
{
  int m=param.size();
  arma::vec alpha=arma::zeros(m/2+1);
  alpha.subvec(1,m/2)=param.subvec(0,m/2-1);
  // arma::vec beta=arma::zeros(m/2+1);
  // beta.subvec(1,m/2)=param.subvec(m/2,m-1);
  arma::vec out=arma::zeros(m);
  for (int g=1;g<(m/2+1);g++)
  {
    for (int k=0;k<g;k++)
    {
      double absa=pow(pow(alpha(k)-alpha(g),2)+eps,0.5);
      out(g-1) -= (alpha(k)-alpha(g))/absa*w(k,g);
    }
    for (int h=(g+1);h<(m/2+1);h++)
    {
      double absa=pow(pow(alpha(g)-alpha(h),2)+eps,0.5);
      out(g-1) += (alpha(g)-alpha(h))/absa*w(g,h);
    }
  }
  return out;
}


arma::vec derridgepennominal(arma::vec param)
 {
   int m=param.size();
   arma::vec alpha=param.subvec(0,m/2-1);
   arma::vec out=arma::zeros(m);
   for (int k=0;k<(m/2);k++)
   {
     for (int h=0;h<(m/2);h++)
     {
       out(k) += 2*(alpha(k)-alpha(h));
     }
   }
   return out;
 }

double aug_lagr_pen(arma::vec param, arma::mat gamma, arma::mat v, double c) 
{
  int m=param.size();
  arma::vec alpha=arma::zeros(m/2+1);
  alpha.subvec(1,m/2)=param.subvec(0,m/2-1);
  //arma::vec beta=arma::zeros(m/2+1);
  //beta.subvec(1,m/2)=param.subvec(m/2,m-1);
  double out=0;
  for (int k=0;k<(m/2);k++)
  {
    for (int h=(k+1);h<(m/2+1);h++)
    {
      out += v(k,h)*(gamma(k,h)-alpha(k)+alpha(h));
      //out += u(k,h)*(phi(k,h)-beta(k)+beta(h));
      out += c/2*pow(gamma(k,h)-alpha(k)+alpha(h),2);
      //out += c/2*pow(phi(k,h)-beta(k)+beta(h),2);
    }
  }
  return out;
}

arma::vec der_aug_lagr_pen(arma::vec param, arma::mat gamma, arma::mat v, double c)
{
  int m=param.size();
  arma::vec alpha=arma::zeros(m/2+1);
  alpha.subvec(1,m/2)=param.subvec(0,m/2-1);
  // arma::vec beta=arma::zeros(m/2+1);
  // beta.subvec(1,m/2)=param.subvec(m/2,m-1);
  arma::vec out=arma::zeros(m);
  for (int g=1;g<(m/2+1);g++)
  {
    for (int k=0;k<g;k++)
    {
      out(g-1) += v(k,g) + c* (gamma(k,g)-alpha(k)+alpha(g));
    }
    for (int h=(g+1);h<(m/2+1);h++)
    {
      out(g-1) -= v(g,h) + c* (gamma(g,h)-alpha(g)+alpha(h));
    }
  }
  return out;
}




// yes
// [[Rcpp::export]]
arma::mat nominallikRcppA(arma::vec par, arma::mat data, arma::vec nodes, 
                          arma::vec weights, arma::mat numpatt,
                          arma::uvec itemsselect, double lambda=0)
{
  int j;
  int nitems = data.n_cols;
  arma::rowvec ncatm1 = max(data,0); 
  arma::umat ncatm1u = arma::conv_to< arma::uvec >::from( ncatm1 );
  arma::umat datau = arma::conv_to< arma::umat >::from( data );
  arma::uvec cumncat = cumsum(ncatm1u); 
  arma::uvec z = arma::zeros<arma::uvec>(1);
  arma::uvec cumncat0 = join_cols(z,cumncat);
  int nq = nodes.size();
  int n=data.n_rows;
  List miss(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec xxj=data.col(j);
    arma::uvec missval = find(xxj == -999.00);
    miss(j) = missval;
  }
  List probs(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
    probs(j) = nomprobsC(parj,nodes);
  }
  arma::mat  logprodj = arma::zeros<arma::mat>(n,nq);
  for (j=0;j<nitems;j++) {
    arma::mat probsj=probs(j);
    arma::mat logprobsj=log(probsj);
    arma::uvec xj=datau.col(j);
    arma::mat logprobsj_exp = logprobsj.rows(xj);
    arma::uvec missj=miss(j);
    if (missj.size()>0) logprobsj_exp.rows(missj).zeros();
    //Rcout << "logprobsj_exp : " << logprobsj_exp << "\n";
    logprodj += logprobsj_exp;
  }
  arma::mat prodj=exp(logprodj);
  arma::mat logsumq = log(prodj * weights) ;
  arma::mat mlik =  - sum(logsumq%numpatt);
  
  double pen=0;
  if (lambda>0) {
    for (auto j : itemsselect) {
      arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
      pen += ridgepennominal(parj);
    }
  }
  arma::mat out = mlik + lambda * pen;
  return out;
}


// yes
// [[Rcpp::export]]
arma::mat gradnominallikRcppA(arma::vec par, arma::mat data, arma::vec nodes, 
                              arma::vec weights, arma::mat numpatt, 
                              arma::uvec itemsselect, double lambda=0)
{
  //Rcout << "par : " << par << "\n";
  int j;
  int nitems = data.n_cols;
  arma::rowvec ncatm1 = max(data,0); 
  arma::umat ncatm1u = arma::conv_to< arma::uvec >::from( ncatm1 );
  arma::umat datau = arma::conv_to< arma::umat >::from( data );
  arma::uvec cumncat = cumsum(ncatm1u); 
  arma::uvec z = arma::zeros<arma::uvec>(1);
  arma::uvec cumncat0 = join_cols(z,cumncat);
  int nq = nodes.size();
  int npar = par.size();
  int n=data.n_rows;
  List miss(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec xxj=data.col(j);
    arma::uvec missval = find(xxj == -999.00);
    miss(j) = missval;
  }
  List probs(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
    probs(j) = nomprobsC(parj,nodes);
  }
  arma::cube logpr = arma::zeros<arma::cube>(n,nq,nitems);
  for (j=0;j<nitems;j++) {
    arma::mat probsj=probs(j);
    arma::mat logprobsj=log(probsj);
    arma::uvec xj=datau.col(j);
    arma::mat logprobs_xj = logprobsj.rows(xj);
    arma::uvec missj=miss(j);
    if (missj.size()>0) logprobs_xj.rows(missj).zeros();
    logpr.slice(j) = logprobs_xj;
  }
  arma::mat prodj = exp(sum(logpr,2));
  arma::mat sumq = prodj * weights ;
  arma::mat der(npar,1);
  for (j=0;j<nitems;j++) {
    arma::cube logprbis=logpr;
    logprbis.slice(j).zeros();
    arma::mat prodmj = exp(sum(logprbis,2));
    arma::mat probsj=probs(j);
    arma::uvec xj=datau.col(j);
    arma::uvec missj=miss(j);
    for (arma::uword k=1;k<=ncatm1u(j);k++) {
      arma::mat probs_xj = exp(logpr.slice(j)); //probsj.rows(xj);
      arma::mat derprob(n,nq);
      arma::uvec ids = find(xj == k);
      arma::uvec ids1 = find(xj != k);
      arma::mat replace = probs_xj.rows(ids)%(1-probs_xj.rows(ids));
      derprob.rows(ids) = replace;
      arma::mat mat1=repmat(probsj.row(k),ids1.size(),1);
      arma::mat replace1 = -probs_xj.rows(ids1)%mat1;
      derprob.rows(ids1) = replace1;
      
      if (missj.size()>0) derprob.rows(missj).zeros();
      
      arma::mat sumq_numerator = (prodmj%derprob)*weights;
      der(k-1+ncatm1u(j)+cumncat0(j)*2) = arma::as_scalar(sum(sumq_numerator/sumq%numpatt)); //thresholds
      derprob.each_row() %= nodes.t();
      arma::mat sumq_numerator_dscr = (prodmj%derprob)*weights;
      der(k-1+cumncat0(j)*2) = arma::as_scalar(sum(sumq_numerator_dscr/sumq%numpatt)); //discriminations
    }
  }
  if (lambda>0) {
    arma::vec derpen=arma::zeros<arma::vec>(npar);
    for (auto j : itemsselect) {
      arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
      derpen.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1) = derridgepennominal(parj);
    }
    der -= lambda*derpen;
  }
  return -der;
}


// yes
// [[Rcpp::export]]
arma::mat nominallikaugRcppA(arma::vec par, arma::mat data, arma::vec nodes, 
                             arma::vec weights, List gamma, List v, 
                             double c, arma::mat numpatt, 
                             arma::uvec itemsselect, double lambda=0)
{
  //Rcout << "par : " << par << "\n";
  int j;
  int nitems = data.n_cols;
  arma::rowvec ncatm1 = max(data,0); 
  
  arma::umat ncatm1u = arma::conv_to< arma::uvec >::from( ncatm1 );
  arma::umat datau = arma::conv_to< arma::umat >::from( data );
  arma::uvec cumncat = cumsum(ncatm1u); 
  arma::uvec z = arma::zeros<arma::uvec>(1);
  arma::uvec cumncat0 = join_cols(z,cumncat);
  int nq = nodes.size();
  int n=data.n_rows;
  List miss(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec xxj=data.col(j);
    arma::uvec missval = find(xxj == -999.00);
    miss(j) = missval;
  }
  List probs(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
    probs(j) = nomprobsC(parj,nodes);
  }
  
  arma::mat  logprodj = arma::zeros<arma::mat>(n,nq);
  for (j=0;j<nitems;j++) {
    arma::mat probsj=probs(j);
    arma::mat logprobsj=log(probsj);
    arma::uvec xj=datau.col(j);
    arma::mat logprobsj_exp = logprobsj.rows(xj);
    arma::uvec missj=miss(j);
    if (missj.size()>0) logprobsj_exp.rows(missj).zeros();
    logprodj += logprobsj_exp;
  }
  arma::mat prodj=exp(logprodj);
  arma::mat logsumq = log(prodj * weights) ;
  arma::mat mlik =  - sum(logsumq%numpatt);
  double pen=0;
  if (lambda>0) {
    for (auto j : itemsselect) {
      arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
      pen += aug_lagr_pen(parj,Rcpp::as<arma::mat>(gamma(j)),Rcpp::as<arma::mat>(v(j)),c);
    }
  }
  arma::mat out = mlik + lambda * pen;
  return out;
}


// yes
// [[Rcpp::export]]
arma::mat gradnominallikaugRcppA(arma::vec par, arma::mat data, arma::vec nodes, arma::vec weights, 
                                 List gamma, List v, double c, arma::mat numpatt, 
                                 arma::uvec itemsselect, double lambda=0)
{
  int j;
  int nitems = data.n_cols;
  arma::rowvec ncatm1 = max(data,0); 
  
  arma::umat ncatm1u = arma::conv_to< arma::uvec >::from( ncatm1 );
  arma::umat datau = arma::conv_to< arma::umat >::from( data );
  arma::uvec cumncat = cumsum(ncatm1u); 
  arma::uvec z = arma::zeros<arma::uvec>(1);
  arma::uvec cumncat0 = join_cols(z,cumncat);
  int nq = nodes.size();
  int npar = par.size();
  int n=data.n_rows;
  List miss(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec xxj=data.col(j);
    arma::uvec missval = find(xxj == -999.00);
    miss(j) = missval;
  }
  List probs(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
    probs(j) = nomprobsC(parj,nodes);
  }

  arma::cube logpr = arma::zeros<arma::cube>(n,nq,nitems);
  for (j=0;j<nitems;j++) {
    arma::mat probsj=probs(j);
    arma::mat logprobsj=log(probsj);
    arma::uvec xj=datau.col(j);
    arma::mat logprobs_xj = logprobsj.rows(xj);
    arma::uvec missj=miss(j);
    if (missj.size()>0) logprobs_xj.rows(missj).zeros();
    logpr.slice(j) = logprobs_xj;
  }
  arma::mat prodj = exp(sum(logpr,2));
  arma::mat sumq = prodj * weights ;
  arma::mat der(npar,1);
  for (j=0;j<nitems;j++) {
    arma::cube logprbis=logpr;
    logprbis.slice(j).zeros();
    arma::mat prodmj = exp(sum(logprbis,2));
    arma::mat probsj=probs(j);
    arma::uvec xj=datau.col(j);
    arma::uvec missj=miss(j);
    for (arma::uword k=1;k<=ncatm1u(j);k++) {
      arma::mat probs_xj = exp(logpr.slice(j)); //probsj.rows(xj);
      arma::mat derprob(n,nq);
      arma::uvec ids = find(xj == k);
      arma::uvec ids1 = find(xj != k);
      arma::mat replace = probs_xj.rows(ids)%(1-probs_xj.rows(ids));
      derprob.rows(ids) = replace;
      arma::mat mat1=repmat(probsj.row(k),ids1.size(),1);
      arma::mat replace1 = -probs_xj.rows(ids1)%mat1;
      derprob.rows(ids1) = replace1;
      
      if (missj.size()>0) derprob.rows(missj).zeros();
      
      arma::mat sumq_numerator = (prodmj%derprob)*weights;
      der(k-1+ncatm1u(j)+cumncat0(j)*2) = arma::as_scalar(sum(sumq_numerator/sumq%numpatt)); //thresholds
      derprob.each_row() %= nodes.t();
      arma::mat sumq_numerator_dscr = (prodmj%derprob)*weights;
      der(k-1+cumncat0(j)*2) = arma::as_scalar(sum(sumq_numerator_dscr/sumq%numpatt)); //discriminations
    }
  }
  if (lambda>0) {
    arma::vec derpen=arma::zeros<arma::vec>(npar);
    for (auto j : itemsselect) {
      arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
      derpen.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1) = der_aug_lagr_pen(parj,Rcpp::as<arma::mat>(gamma(j)),
                    Rcpp::as<arma::mat>(v(j)),c);
    }
    der -= lambda*derpen;
  }
  
  return -der;
}


//yes
// [[Rcpp::export]]
arma::mat nominallik_fpenRcppA(arma::vec par, arma::mat data, arma::vec nodes, 
                               arma::vec weights, arma::mat numpatt, 
                               arma::uvec itemsselect, List w, 
                               double lambda=0, double eps=0.001)
{
  int j;
  int nitems = data.n_cols;
  arma::rowvec ncatm1 = max(data,0); 
  
  arma::umat ncatm1u = arma::conv_to< arma::uvec >::from( ncatm1 );
  arma::umat datau = arma::conv_to< arma::umat >::from( data );
  arma::uvec cumncat = cumsum(ncatm1u); 
  arma::uvec z = arma::zeros<arma::uvec>(1);
  arma::uvec cumncat0 = join_cols(z,cumncat);
  int nq = nodes.size();
  int n=data.n_rows;
  List miss(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec xxj=data.col(j);
    arma::uvec missval = find(xxj == -999.00);
    miss(j) = missval;
  }
  List probs(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
    probs(j) = nomprobsC(parj,nodes);
  }
  arma::mat  logprodj = arma::zeros<arma::mat>(n,nq);
  for (j=0;j<nitems;j++) {
    arma::mat probsj=probs(j);
    arma::mat logprobsj=log(probsj);
    arma::uvec xj=datau.col(j);
    arma::mat logprobsj_exp = logprobsj.rows(xj);
    arma::uvec missj=miss(j);
    if (missj.size()>0) logprobsj_exp.rows(missj).zeros();
    logprodj += logprobsj_exp;
  }
  arma::mat prodj=exp(logprodj);
  arma::mat logsumq = log(prodj * weights) ;
  arma::mat mlik =  - sum(logsumq%numpatt);
  
  double pen=0;
  if (lambda>0) {
    // itemsselect.print();
    for (auto j : itemsselect) {
      //Rcout << "j : " << j << "\n";
      arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
      pen += fpenC(parj,eps,Rcpp::as<arma::mat>(w(j)));
    }
  }
  arma::mat out = mlik + lambda * pen;
  return out;
}


// yes
// [[Rcpp::export]]
arma::mat gradnominallik_fpenRcppA(arma::vec par, arma::mat data, arma::vec nodes, 
                              arma::vec weights, arma::mat numpatt, 
                              arma::uvec itemsselect, List w, double lambda=0, 
                              double eps=0.001)
{
  //Rcout << "par : " << par << "\n";
  int j;
  int nitems = data.n_cols;
  arma::rowvec ncatm1 = max(data,0); 
  arma::umat ncatm1u = arma::conv_to< arma::uvec >::from( ncatm1 );
  arma::umat datau = arma::conv_to< arma::umat >::from( data );
  arma::uvec cumncat = cumsum(ncatm1u); 
  arma::uvec z = arma::zeros<arma::uvec>(1);
  arma::uvec cumncat0 = join_cols(z,cumncat);
  int nq = nodes.size();
  int npar = par.size();
  int n=data.n_rows;
  List miss(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec xxj=data.col(j);
    arma::uvec missval = find(xxj == -999.00);
    miss(j) = missval;
  }
  List probs(nitems);
  for (j=0;j<nitems;j++) {
    arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
    probs(j) = nomprobsC(parj,nodes);
  }
  arma::cube logpr = arma::zeros<arma::cube>(n,nq,nitems);
  for (j=0;j<nitems;j++) {
    arma::mat probsj=probs(j);
    arma::mat logprobsj=log(probsj);
    arma::uvec xj=datau.col(j);
    arma::mat logprobs_xj = logprobsj.rows(xj);
    arma::uvec missj=miss(j);
    if (missj.size()>0) logprobs_xj.rows(missj).zeros();
    logpr.slice(j) = logprobs_xj;
  }
  arma::mat prodj = exp(sum(logpr,2));
  arma::mat sumq = prodj * weights ;
  arma::mat der(npar,1);
  for (j=0;j<nitems;j++) {
    arma::cube logprbis=logpr;
    logprbis.slice(j).zeros();
    arma::mat prodmj = exp(sum(logprbis,2));
    arma::mat probsj=probs(j);
    arma::uvec xj=datau.col(j);
    arma::uvec missj=miss(j);
    for (arma::uword k=1;k<=ncatm1u(j);k++) {
      arma::mat probs_xj = exp(logpr.slice(j)); //probsj.rows(xj);
      arma::mat derprob(n,nq);
      arma::uvec ids = find(xj == k);
      arma::uvec ids1 = find(xj != k);
      arma::mat replace = probs_xj.rows(ids)%(1-probs_xj.rows(ids));
      derprob.rows(ids) = replace;
      arma::mat mat1=repmat(probsj.row(k),ids1.size(),1);
      arma::mat replace1 = -probs_xj.rows(ids1)%mat1;
      derprob.rows(ids1) = replace1;
      
      if (missj.size()>0) derprob.rows(missj).zeros();
      
      arma::mat sumq_numerator = (prodmj%derprob)*weights;
      der(k-1+ncatm1u(j)+cumncat0(j)*2) = arma::as_scalar(sum(sumq_numerator/sumq%numpatt)); //thresholds
      derprob.each_row() %= nodes.t();
      arma::mat sumq_numerator_dscr = (prodmj%derprob)*weights;
      der(k-1+cumncat0(j)*2) = arma::as_scalar(sum(sumq_numerator_dscr/sumq%numpatt)); //discriminations
    }
  }
  if (lambda>0) {
    arma::vec derpen=arma::zeros<arma::vec>(npar);
    for (auto j : itemsselect) {
      arma::vec parj = par.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1);
      derpen.subvec(cumncat0(j)*2,cumncat0(j+1)*2-1) = derfpen(parj,eps,Rcpp::as<arma::mat>(w(j)));
    }
    der -= lambda*derpen;
  }
  return -der;
}



