### utils


Ltrans<-function(X,d=TRUE){ X[upper.tri(X,d)]  }

Ltrinv<-function(x,V,d=TRUE){ Y = matrix(0,ncol = V,nrow = V);
Y[upper.tri(Y,d)]=x;return(Y + t(Y) -d*diag(diag(Y)))  }


soft_thre<-function(rho,z,lambda){
  if (z==0){
    return(0)
  }
  if (rho<=-lambda){
    u=(rho+lambda)/z
  }else if(rho<lambda){
    u=0
  }else{
    u=(rho-lambda)/z
  }
  return(u)
}

normalize<-function(x){
  if (all(x==0)){
    x=x
  }else{
    x=x/norm(x,type='2')
  }
  return(x)
}



#

psdmat_inverse<-function(mat)
{
  # Originally defined for PSD matrices
  p = dim(mat)[1]
  eigendecomp = eigen(mat)


  if( min(eigendecomp$values)<=0.0001 )
  {
    message("Matrix has nearly zero eigenvalue, perturbation is added. ")

    perturb <- max(max(eigendecomp$values) - p * min(eigendecomp$values), 0)/(p - 1)
  }else
  {
    perturb = 0
  }

  return( (eigendecomp$vectors)%*%diag(1/(perturb+eigendecomp$values))%*%t(eigendecomp$vectors) )
}

### conditional mean and covariance of A

Sig_AconY_func<-function(S_cross_S,p,L,sigma2,phi2){

  Sig1= S_cross_S * (1/phi2)
  return(solve(Sig1 +diag(c(1/sigma2),nrow=L)))
}

mu_AconY_func<-function(Sig_AconY,SL,XFmat,dat,sigma2,phi2,p){

  return(t(Sig_AconY%*%((1/sigma2)*XFmat+t(SL)%*%(rep(1/phi2,p)*t(dat)))))

}


### mean effect of A

Fmat_func<-function(X,A,ni){
  t(A)%*%X%*%solve(t(X)%*%X)
}

XFmat_func<-function(X,Fmat){
  Fmat%*%t(X)
}

### variance of A and error

sigma2_func<-function(Sig_AconY,XFmat,A,H){
  ni=nrow(A)
  return(diag(Sig_AconY)+apply((A-t(XFmat))^2,2,function(x) mean(x/(rep(1,ni)-H))))

}

phi2_func<-function(dat,A,S,Sig_AconY,S_cross_S,p){
  n=nrow(dat)
  phi2=(sum(dat*(dat-2*A%*%S))/(n*p)+sum(t(A)%*%A*S_cross_S)/(n*p))+sum(S_cross_S*Sig_AconY)/(p)

  return(phi2)
}

### marginal likelihoood

inv_margin_cov_func=function(sigma2,phi2,S,S_cross_S){
  p=ncol(S)
  inv_margin_cov=diag(1/phi2,p)-1/(phi2^2)*t(S)%*%solve(diag(1/c(sigma2),nrow=length(c(sigma2)))+1/phi2*S_cross_S)%*%S

  return(inv_margin_cov)
}

LS_func<-function(dat,margin_mean,inv_margin_cov){
  sum(crossprod(dat-margin_mean)*inv_margin_cov)
}


margin_LS_func<-function(XFmat,S,inv_margin_cov,datGrp){
  margin_mean=lapply(XFmat, function(x) t(x)%*%S)

  LS= mapply(LS_func, dat=datGrp, margin_mean=margin_mean,inv_margin_cov=inv_margin_cov,SIMPLIFY = TRUE)


  return(LS)

}


det_cov_func<-function(sigma2,phi2,S,S_cross_S){
  p=ncol(S)
  return(unlist(1/2*(determinant(diag(1/c(sigma2),nrow=length(c(sigma2)))+S_cross_S/phi2)$modulus+log(phi2)*p+sum(log(sigma2)))))

}

loss_func<-function(sigma2,phi2,S,S_cross_S,XFmat,datGrp,M,ni){

  p=ncol(datGrp[[1]]);n=sum(ni)

  inv_margin_cov=mapply(inv_margin_cov_func,sigma2=sigma2,phi2=phi2,MoreArgs = list(S=S,S_cross_S=S_cross_S),SIMPLIFY = FALSE)

  det_cov=mapply(det_cov_func,sigma2=sigma2,phi2=phi2,MoreArgs = list(S=S,S_cross_S=S_cross_S),SIMPLIFY = TRUE)
  margin_LS=margin_LS_func(XFmat=XFmat,S=S,inv_margin_cov=inv_margin_cov,datGrp=datGrp)

  return(((Reduce('+',ni*det_cov))+sum(margin_LS)/2)/(n*p))
}



#### initial

eigen_initial<- function(dat,L,V){
  mean_dat<-apply(dat, 2, mean)
  mean_dat_mat<-Ltrinv(mean_dat,V)

  eigen_decomp<-eigen(mean_dat_mat)
  U_ini=eigen_decomp$vectors[,1:L,drop=FALSE]
  S_ini_mat= sapply(1:L,function(l) tcrossprod(U_ini[,l]))
  S_ini=t(apply(S_ini_mat,2, function(x) Ltrans(matrix(x,V,V))))
  A_ini = dat%*%t(S_ini)%*%solve(S_ini%*%t(S_ini))

  return(list(U=U_ini,S=S_ini,A=A_ini))

}


### update


update_func <- function(dat,A,U,batch,XFmat,X,XGrp,sigma2,phi2,HGrp,lambda_ch = 0.5,p,V,L,order.temp.batch,tau=0.3 * sqrt(log(L*V)/n),datGrp,n)
{



  rmat = matrix(rep(1:V,V),ncol=V)
  cmat = t(rmat)
  Umat<-U
  Lcoorsym = matrix(0,ncol=2,nrow=V*(V+1)/2)
  Lcoorsym[,1] = Ltrans(rmat,T)
  Lcoorsym[,2] = Ltrans(cmat,T)


  newS = array(dim=c(L,p))

  b1<-sapply(1:nrow(A),function(i) Umat%*%diag(A[i,],ncol=L))
  b2<-t(b1)
  b3<-matrix(b2,ncol=L)
  phi2_vec=rep(unlist(phi2)[batch],V)
  yvpen=sapply(1:V,function(v)c(dat[,which(Lcoorsym[,1] == v | Lcoorsym[, 2] == v)])             # dat
  )

  M <- b3 %*% t(Umat)
  rho=matrix(0,nrow=V,ncol=L)
  z=colMeans(b3^2/phi2_vec)

  if (is.null(tau)){
    lambda=lambda_ch
    for (l in 1:L) {
      rho[,l] <- colMeans((b3[,l]*(yvpen-(M - outer( b3[, l], t(Umat)[l, ]))))/phi2_vec)
      Umat[,l]=sapply(1:V,function(v) soft_thre(rho[v,l],z[l],lambda))
    }
  }else{
    lambda<-lambda_ch/tau*ifelse(abs(Umat)<=tau,1,0)
    for (l in 1:L) {
      rho[,l] <- colMeans(b3[,l]*(yvpen-(M - outer( b3[, l], t(Umat)[l, ]))))
      Umat[,l]=sapply(1:V,function(v) soft_thre(rho[v,l],z[l],lambda[v,l]))
    }
  }


  Umat=apply(Umat,2,normalize)


  newS=t(sapply(1:L,function(l) Ltrans(Umat[,l]%*%t(Umat[,l]))))

  S_cross_S=tcrossprod(newS)

  # Update A
  Sig_AconY=mapply(Sig_AconY_func,sigma2=sigma2,phi2=phi2,MoreArgs=list(S_cross_S=S_cross_S,p=p,L=L),SIMPLIFY = FALSE)

  newA_Grp =  mapply(mu_AconY_func,Sig_AconY= Sig_AconY,XFmat=XFmat,sigma2=sigma2,dat=datGrp,phi2=phi2,MoreArgs=list(SL=t(newS),p=p),SIMPLIFY = FALSE)

  newA=matrix(0, nrow = n, ncol = L)

  for(i in seq_along(order.temp.batch)){
    newA[order.temp.batch[[i]], ] =  newA[order.temp.batch[[i]], ] + newA_Grp[[i]]
  }

  Fmat_new=t(newA)%*%X%*%solve(t(X)%*%X)
  XFmat_new=mapply(XFmat_func,XGrp,MoreArgs = list(Fmat=Fmat_new), SIMPLIFY = FALSE)
  sigma2_new=mapply(sigma2_func,Sig_AconY=Sig_AconY,A=newA_Grp,XFmat=XFmat_new,H=HGrp,SIMPLIFY = FALSE)
  phi2_new=mapply(phi2_func,dat=datGrp,A=newA_Grp,Sig_AconY= Sig_AconY,MoreArgs = list(S=newS,S_cross_S=S_cross_S,p=p), SIMPLIFY = FALSE)

  return(list(A=newA,U=Umat,S=newS,S_cross_S=S_cross_S,Fmat=Fmat_new,XFmat=XFmat_new,sigma2=sigma2_new,phi2=phi2_new))

}

### reparametrize mean effect

alpha_gamma_func<-function(Fmat,M,L,ni){

  q=ncol(Fmat)-M

  if ((M>1)&(q>0)){

    Zmat=Fmat[,c(1,(q+2):(q+M)),drop=FALSE]
    beta=Fmat[,2:(q+1),drop=FALSE]
    gamma=matrix(ncol=M,nrow=L)

    gamma[,1]=rowSums(((rep(1,L)%*%t(-ni[2:M]))*Zmat[,2:M,drop=FALSE]))/(sum(ni))
    alpha=Zmat[,1]-gamma[,1]
    gamma[,2:M]=Zmat[,2:M]+gamma[,1]
  }else if((M==1)&(q>0)){
    alpha=Fmat[,1]
    beta=Fmat[,2:(q+1),drop=FALSE]
    gamma=NULL
  }else{
    Zmat=Fmat[,c(1,(q+2):(q+M)),drop=FALSE]
    beta=NULL
    gamma=matrix(ncol=M,nrow=L)

    gamma[,1]=rowSums(((rep(1,L)%*%t(-ni[2:M]))*Zmat[,2:M,drop=FALSE]))/(sum(ni))
    alpha=Zmat[,1]-gamma[,1]
    gamma[,2:M]=Zmat[,2:M]+gamma[,1]
  }


  return(list(gamma=gamma,alpha=alpha,beta=beta))

}

### harmonization function

harmonize_func<-function(estimates,covariate,ni,batch,dat,order.temp.batch,L,standarize,center,scale){

  n=sum(ni);p=ncol(dat)
  sigma2=estimates$sigma2
  phi2=estimates$phi2
  alpha=estimates$alpha
  beta=estimates$beta
  Fmat=estimates$Fmat
  XFmat=estimates$XFmat

  datGrp= split.data.frame(dat,batch)

  gamma=estimates$gamma

  S=estimates$S




  A=estimates$A
  if (!is.null(beta)){
    modbeta=covariate%*%t(beta)
  }else{
    modbeta=0
  }

  alpha_mat=rep(1,n)%*%t(alpha)


  sigma2_mat=do.call('rbind',sigma2)
  sigma2_h=apply(sigma2_mat,2,function(x) weighted.mean(x,w=ni))

  sigma_weight=(rep(1,n)%*%t(sqrt(sigma2_h)))/sqrt((sigma2_mat[batch,]))
  phi2_vec=unlist(phi2)

  phi2_h=weighted.mean(phi2_vec,w=ni)

  phi_weight=(sqrt(phi2_h)*matrix(1,ncol=p,nrow=n))/(sqrt(phi2_vec[batch])%*%t(rep(1,p)))
  E=dat-A%*%S

  XFmat_cb=matrix(0, nrow = n, ncol = L)

  for(i in seq_along(order.temp.batch)){
    XFmat_cb[order.temp.batch[[i]], ] =  XFmat_cb[order.temp.batch[[i]], ] + t(XFmat[[i]])
  }


  dat_h=(sigma_weight*(A-XFmat_cb)+alpha_mat+modbeta)%*%S+phi_weight*E

  if (standarize){
    dat_h= scale*dat_h+rep(1,n)%*%t(center)
  }

  return(dat_h)
}

