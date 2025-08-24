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
  #p=nrow(SL);L=ncol(SL)
  Sig1= S_cross_S * (1/phi2)
  return(solve(Sig1 +diag(c(1/sigma2),nrow=L)))
}

mu_AconY_func<-function(Sig_AconY,SL,xFmat,Y,sigma2,phi2,p){
  # p=nrow(SL)
  return(t(Sig_AconY%*%((1/sigma2)*xFmat+t(SL)%*%(rep(1/phi2,p)*t(Y)))))

}

### mean effect of A

Fmat_func<-function(x,A,ni){
  t(A)%*%x%*%solve(t(x)%*%x)
}

xFmat_func<-function(x,Fmat){
  Fmat%*%t(x)
}

### variance of A and error

sigma2_func<-function(Sig_AconY,xFmat,A,H){
  ni=nrow(A)
  return(diag(Sig_AconY)+apply((A-t(xFmat))^2,2,function(x) mean(x/(rep(1,ni)-H))))

}

phi2_func<-function(Y,A,S,Sig_AconY,S_cross_S,p){

  # p=ncol(Y);n=nrow(Y);V= (sqrt(1+8*p)-1)/2;L=nrow(S)
  n=nrow(Y)
  phi2=(sum(Y*(Y-2*A%*%S))/(n*p)+sum(t(A)%*%A*S_cross_S)/(n*p))+sum(S_cross_S*Sig_AconY)/(p)

  return(phi2)
}

### marginal likelihoood

inv_margin_cov_func=function(sigma2,phi2,S,S_cross_S){
  p=ncol(S)
  inv_margin_cov=diag(1/phi2,p)-1/(phi2^2)*t(S)%*%solve(diag(1/c(sigma2),nrow=length(c(sigma2)))+1/phi2*S_cross_S)%*%S
  # eigenvalues <- eigen(inv_margin_cov)$values
  # if (any(eigenvalues <=0)){
  #   inv_margin_cov= nearPD(inv_margin_cov)$mat
  # }
  return(inv_margin_cov)
}

LS_func<-function(Y,margin_mean,inv_margin_cov){
  sum(crossprod(Y-margin_mean)*inv_margin_cov)
}


margin_LS_func<-function(xFmat,S,inv_margin_cov,YGrp){
  margin_mean=lapply(xFmat, function(x) t(x)%*%S)

  LS= mapply(LS_func, Y=YGrp, margin_mean=margin_mean,inv_margin_cov=inv_margin_cov,SIMPLIFY = TRUE)


  return(LS)

}


det_cov_func<-function(sigma2,phi2,S,S_cross_S){
  p=ncol(S)
  return(unlist(1/2*(determinant(diag(1/c(sigma2),nrow=length(c(sigma2)))+S_cross_S/phi2)$modulus+log(phi2)*p+sum(log(sigma2)))))

}

loss_func<-function(sigma2,phi2,S,S_cross_S,xFmat,YGrp,M,ni){

  p=ncol(YGrp[[1]]);n=sum(ni)

  inv_margin_cov=mapply(inv_margin_cov_func,sigma2=sigma2,phi2=phi2,MoreArgs = list(S=S,S_cross_S=S_cross_S),SIMPLIFY = FALSE)
  det_cov=mapply(det_cov_func,sigma2=sigma2,phi2=phi2,MoreArgs = list(S=S,S_cross_S=S_cross_S),SIMPLIFY = TRUE)
  margin_LS=margin_LS_func(xFmat=xFmat,S=S,inv_margin_cov=inv_margin_cov,YGrp=YGrp)

  return(((Reduce('+',ni*det_cov))+sum(margin_LS)/2)/(n*p))
}



#### initial

SCH_eigen_initial<- function(Y,L,V){
  mean_Y<-apply(Y, 2, mean)
  mean_Y_mat<-Ltrinv(mean_Y,V)

  eigen_decomp<-eigen(mean_Y_mat)
  U_ini=eigen_decomp$vectors[,1:L,drop=FALSE]
  #D_ini=eigen_decomp$values[1:L]
  #S_ini_mat= sapply(1:L,function(l) quad.tform(D_ini[l],U_ini[,l]))
  S_ini_mat= sapply(1:L,function(l) tcrossprod(U_ini[,l]))
  S_ini=t(apply(S_ini_mat,2, function(x) Ltrans(matrix(x,V,V))))
  A_ini = Y%*%t(S_ini)%*%solve(S_ini%*%t(S_ini))

  return(list(U=U_ini,S=S_ini,A=A_ini))

}

SCH_warmL_initial<- function(Y,L,V,U,A,S){
  L_ini=ncol(U)

  mean_Y<-apply(Y, 2, mean)
  mean_Y_mat<-Ltrinv(mean_Y,V)

  eigen_decomp<-eigen(mean_Y_mat)
  U_ini=eigen_decomp$vectors[,1:L,drop=FALSE]
  U_ini_2=apply(U_ini-U%*%solve(t(U)%*%U)%*%t(U)%*%U_ini,2,normalize)
  U_upd=svd(U_ini_2)$u[,1:(L-L_ini)]


  S_upd_mat= sapply(1:(L-L_ini),function(l) tcrossprod(U_upd[,l]))
  S_upd=t(apply(S_upd_mat,2, function(x) Ltrans(matrix(x,V,V))))
  # S_upd=S_ini[(L_ini+1):L,]
  A_upd = (Y-A%*%S)%*%t(S_upd)%*%psdmat_inverse(S_upd%*%t(S_upd))
  # A_ini=cbind(A,A_upd)
  # S_ini=rbind(S,S_upd)
  return(list(U=U_upd,S=S_upd,A=A_upd))

}


SCH_warmL_random_initial<- function(Y,L,V,U,A,S){
  L_ini=ncol(U)
  U_upd=apply(matrix(rnorm(V*(L-L_ini)),ncol=(L-L_ini)),2,normalize)

  U_upd=apply(U_upd-U%*%psdmat_inverse(t(U)%*%U)%*%t(U)%*%U_upd,2,normalize)


  S_upd_mat= sapply(1:(L-L_ini),function(l) tcrossprod(U_upd[,l]))
  S_upd=t(apply(S_upd_mat,2, function(x) Ltrans(matrix(x,V,V))))

  A_upd = (Y-A%*%S)%*%t(S_upd)%*%psdmat_inverse(S_upd%*%t(S_upd))

  return(list(U=U_upd,S=S_upd,A=A_upd))

}


SCH_random_initial<- function(YGrp,L,V,xGrp,inval1_sigma2,inval2_sigma2,ni,M){
  p=ncol(YGrp[[1]]);q=ncol(xGrp[[1]])
  U_ini=apply(matrix(rnorm(V*L),ncol=L),2,normalize)
  #D_ini=eigen_decomp$values[1:L]
  #S_ini_mat= sapply(1:L,function(l) quad.tform(D_ini[l],U_ini[,l]))
  S_ini_mat= sapply(1:L,function(l) tcrossprod(U_ini[,l]))
  S_ini=t(apply(S_ini_mat,2, function(x) Ltrans(matrix(x,V,V))))

  Fmat=matrix(rnorm(n=L*q),nrow=L)

  xFmat=mapply(xFmat_func,xGrp,MoreArgs = list(Fmat=Fmat), SIMPLIFY = FALSE)

  sigma2_ini=lapply(1:M, function(i) runif(L,inval1_sigma2,inval2_sigma2))

  AGrp <- lapply(1:M, function(j) matrix(rnorm(ni[j]*L,mean=c(t(xFmat[[j]])),sd=rep(sqrt(sigma2_ini[[j]]),each=ni[j])),nrow=ni[j],ncol=L))

  A_ini=do.call('rbind',AGrp)


  phi2_ini=lapply(1:M,function(i) sum((YGrp[[i]]- AGrp[[i]]%*%S_ini)^2)/(ni[i]*p))



  return(list(U=U_ini,S=S_ini,A=A_ini,Fmat=Fmat,sigma2=sigma2_ini,phi2=phi2_ini))

}



### update


SCH_update_OLS <- function(Y,A,U,grp,xFmat,x=x,sigma2,phi2,Sig_AconY,HGrp, penalt = NULL,lambda_ch = 0.5,tau=0.3 * sqrt(log(L*V)/n),silent = FALSE,MaxIteration,espli,MaxIteration_U=50)
{

  if(is.null(penalt)){
    if(!silent)
      cat("SCH without penalty.")
  }else{
    if(!silent)
      cat(paste("SCH with", penalt,"penalty."))
  }



  YGrp <- split.data.frame(Y,grp)
  xGrp=split.data.frame(x,grp)


  p = dim(Y)[2]         # Number of edges
  V = (sqrt(1+8*p)-1)/2
  N = dim(Y)[1]
  L = ncol(U)

  rmat = matrix(rep(1:V,V),ncol=V)
  cmat = t(rmat)
  Umat<-U
  Lcoorsym = matrix(0,ncol=2,nrow=V*(V+1)/2)
  Lcoorsym[,1] = Ltrans(rmat,T)
  Lcoorsym[,2] = Ltrans(cmat,T)


  lambda=lambda_ch

  newS = array(dim=c(L,p))


  # ### ols
  b1<-sapply(1:nrow(A),function(i) Umat%*%diag(A[i,],ncol=L))
  b2<-t(b1)
  b3<-matrix(b2,ncol=L)

  b4=matrix(sapply(1:nrow(A),function(i) crossprod(U)*(tcrossprod(A[i,])+(Sig_AconY[[grp[i]]]))),nrow=L^2)

  b5=solve(matrix(rowSums(b4),L,L))

  yvpen=sapply(1:V,function(v)c(Y[,which(Lcoorsym[,1] == v | Lcoorsym[, 2] == v)])             # dat
  )

  Umat=t(b5%*%t(b3)%*%yvpen)

  ####


  ###
  Umat=apply(Umat,2,normalize)


  newS=t(sapply(1:L,function(l) Ltrans(Umat[,l]%*%t(Umat[,l]))))

  S_cross_S=tcrossprod(newS)
  # Update A
  #newA = Y%*%t(newS) %*% solve(newS%*%t(newS))

  # Update A
  Sig_AconY=mapply(Sig_AconY_func,sigma2=sigma2,phi2=phi2,MoreArgs=list(S_cross_S=S_cross_S,p=p,L=L),SIMPLIFY = FALSE)

  newA_Grp =  mapply(mu_AconY_func,Sig_AconY= Sig_AconY,xFmat=xFmat,sigma2=sigma2,Y=YGrp,phi2=phi2,MoreArgs=list(SL=t(newS),p=p),SIMPLIFY = FALSE)

  newA=do.call(rbind, newA_Grp)


  #gamma_new=mapply(gamma_func,newA_Grp,MoreArgs = list(ni=ni), SIMPLIFY = FALSE)
  Fmat_new=t(newA)%*%x%*%solve(t(x)%*%x)
  xFmat_new=mapply(xFmat_func,xGrp,MoreArgs = list(Fmat=Fmat_new), SIMPLIFY = FALSE)
  #sigma2_new=mapply(sigma2_func,YGrp,MoreArgs=list(SL=newS),SIMPLIFY = FALSE)
  sigma2_new=mapply(sigma2_func,Sig_AconY=Sig_AconY,A=newA_Grp,xFmat=xFmat_new,H=HGrp,SIMPLIFY = FALSE)
  #phi2_new=mapply(phi2_func,Y=YGrp,A=newA_Grp,MoreArgs = list(S=newS), SIMPLIFY = FALSE)
  phi2_new=mapply(phi2_func,Y=YGrp,A=newA_Grp,Sig_AconY= Sig_AconY,MoreArgs = list(S=newS,S_cross_S=S_cross_S,p=p), SIMPLIFY = FALSE)



  return(list(A=newA,U=Umat,S=newS,S_cross_S=S_cross_S,Fmat=Fmat_new,xFmat=xFmat_new,sigma2=sigma2_new,phi2=phi2_new,Sig_AconY=Sig_AconY))

}

SCH_update_TLP_full <- function(Y,A,U,grp,xFmat,x=x,sigma2,phi2,Sig_AconY,HGrp, penalt = NULL,lambda_ch = 0.5,tau=0.3 * sqrt(log(L*V)/n),silent = FALSE,MaxIteration,espli,MaxIteration_U)
{

  if(is.null(penalt)){
    if(!silent)
      cat("SCH without penalty.")
  }else{
    if(!silent)
      cat(paste("SCH with", penalt,"penalty."))
  }



  YGrp <- split.data.frame(Y,grp)
  xGrp=split.data.frame(x,grp)


  p = dim(Y)[2]         # Number of edges
  V = (sqrt(1+8*p)-1)/2
  N = dim(Y)[1]
  L = ncol(U)

  rmat = matrix(rep(1:V,V),ncol=V)
  cmat = t(rmat)
  Umat<-U
  Lcoorsym = matrix(0,ncol=2,nrow=V*(V+1)/2)
  Lcoorsym[,1] = Ltrans(rmat,T)
  Lcoorsym[,2] = Ltrans(cmat,T)


  newS = array(dim=c(L,p))

  Iter_U=1
  Umat_new=matrix(ncol=L,nrow=V)




  while (Iter_U <= MaxIteration_U){

    b1<-sapply(1:nrow(A),function(i) Umat%*%diag(A[i,],ncol=L))
    b2<-t(b1)
    b3<-matrix(b2,ncol=L)
    b11<-sapply(1:nrow(A),function(i) (Umat^2)%*%diag(diag(Sig_AconY[[grp[i]]]),ncol=length(diag(Sig_AconY[[grp[i]]]))))
    b21<-t(b11)
    b31<-matrix(b21,ncol=L)
    phi2_vec=rep(unlist(phi2)[grp],V)
    yvpen=sapply(1:V,function(v)c(Y[,which(Lcoorsym[,1] == v | Lcoorsym[, 2] == v)])             # dat
    )

    b4 <- b3 %*% t(Umat)
    b5 <- lapply(1:length(Sig_AconY),function(x) Umat * (Umat %*% Sig_AconY[[x]] - Umat %*% diag(diag(Sig_AconY[[x]]),ncol=length(diag(Sig_AconY[[x]])))))
    b6 <-do.call('rbind',b5[grp])

    rho=matrix(0,nrow=V,ncol=L)
    z=colMeans((b3^2+b31)/phi2_vec)


    lambda<-lambda_ch/tau*ifelse(abs(Umat)<=tau,1,0)


    for (l in 1:L) {
      b4 <- b3 %*% t(Umat)
      rho[,l] <- colMeans((b3[,l]*(yvpen-(b4 - outer( b3[, l], t(Umat)[l, ])))-outer(b6[,l],t(Umat)[l, ]))/phi2_vec)
      Umat_new[,l]=sapply(1:V,function(v) soft_thre(rho[v,l],z[l],lambda[v,l]))
    }

    Umat_new=apply(Umat_new,2,normalize)

    err_U=norm(as.matrix(Umat_new-Umat))/norm(as.matrix(Umat))


    Umat=Umat_new



    if (err_U< espli){
      break
    }
    Iter_U=Iter_U+1
  }


  Umat=apply(Umat,2,normalize)


  newS=t(sapply(1:L,function(l) Ltrans(Umat[,l]%*%t(Umat[,l]))))
  S_cross_S=tcrossprod(newS)

  # Update A
  Sig_AconY=mapply(Sig_AconY_func,sigma2=sigma2,phi2=phi2,MoreArgs=list(S_cross_S=S_cross_S,p=p,L=L),SIMPLIFY = FALSE)

  newA_Grp =  mapply(mu_AconY_func,Sig_AconY= Sig_AconY,xFmat=xFmat,sigma2=sigma2,Y=YGrp,phi2=phi2,MoreArgs=list(SL=t(newS),p=p),SIMPLIFY = FALSE)

  newA=do.call(rbind, newA_Grp)


   Fmat_new=t(newA)%*%x%*%solve(t(x)%*%x)
  xFmat_new=mapply(xFmat_func,xGrp,MoreArgs = list(Fmat=Fmat_new), SIMPLIFY = FALSE)
   sigma2_new=mapply(sigma2_func,Sig_AconY=Sig_AconY,A=newA_Grp,xFmat=xFmat_new,H=HGrp,SIMPLIFY = FALSE)
   phi2_new=mapply(phi2_func,Y=YGrp,A=newA_Grp,Sig_AconY= Sig_AconY,MoreArgs = list(S=newS,S_cross_S=S_cross_S,p=p), SIMPLIFY = FALSE)


  return(list(A=newA,U=Umat,S=newS,S_cross_S=S_cross_S,Fmat=Fmat_new,xFmat=xFmat_new,sigma2=sigma2_new,phi2=phi2_new,Sig_AconY=Sig_AconY))
  #return(list(A=newA,U=Umat,S=newS))
}


### reparametrize mean effect

alpha_gamma_func<-function(Fmat,grp,n.batch,L,ni){

  q=ncol(Fmat)-n.batch

  if ((n.batch>1)&(q>0)){

    Zmat=Fmat[,c(1,(q+2):(q+n.batch)),drop=FALSE]
    beta=Fmat[,2:(q+1)]
    gamma=matrix(ncol=n.batch,nrow=L)

    gamma[,1]=rowSums(((rep(1,L)%*%t(-ni[2:n.batch]))*Zmat[,2:n.batch,drop=FALSE]))/(sum(ni))
    alpha=Zmat[,1]-gamma[,1]
    gamma[,2:n.batch]=Zmat[,2:n.batch]+gamma[,1]
  }else if((n.batch==1)&(q>0)){
    alpha=Fmat[,1]
    beta=Fmat[,2:(q+1)]
    gamma=NULL
  }else{
    Zmat=Fmat[,c(1,(q+2):(q+n.batch)),drop=FALSE]
    beta=NULL
    gamma=matrix(ncol=n.batch,nrow=L)

    gamma[,1]=rowSums(((rep(1,L)%*%t(-ni[2:n.batch]))*Zmat[,2:n.batch,drop=FALSE]))/(sum(ni))
    alpha=Zmat[,1]-gamma[,1]
    gamma[,2:n.batch]=Zmat[,2:n.batch]+gamma[,1]
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
xFmat=estimates$xFmat

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

xFmat_cb=matrix(0, nrow = n, ncol = L)

  for(i in seq_along(order.temp.batch)){
  xFmat_cb[order.temp.batch[[i]], ] =xFmat_cb[order.temp.batch[[i]], ] + t(xFmat[[i]])
  }


  dat_h=(sigma_weight*(A-xFmat_cb)+alpha_mat+modbeta)%*%S+phi_weight*E

  if (standarize){
    dat_h= t(scale*t(dat_h))+rep(1,n)%*%t(center)
  }

  return(dat_h)
}

