#' @title LASERC harmonization
#' @description Main function to perform LASERC harmonization.
#' @param dat Numeric matrix with subjects as rows \code{(n)},
#'                and functional connectivity features as columns \code{(p)}.
#' @param batch Numeric or character vector specifying the batch/site
#'     variable needed for harmonization.
#' @param mod Optional model matrix for outcome of interest and
#'     other covariates besides batch/scanner.
#' @param L Numeric specifying the number of sparsity connectivity patterns that needs to be tuned (default = 3).
#' @param standardize  Logical indicating whether to standardize the columns of \code{dat} prior to
#'   model fitting. Defaults to \code{TRUE}.
#' @param center A numeric vector of length \code{p} giving column means for centering (if \code{standardize=TRUE}).
#'   If \code{NULL}, column means of \code{dat} are used.
#' @param scale A numeric vector of length \code{p} giving column standard deviations for scaling (if \code{standardize=TRUE}).
#'   If \code{NULL}, column standard deviations of \code{dat} are used.
#' @param MaxIteration An integer specifying the maximum number of iterations for each stage of the
#'   iterative procedure. Default is 200.
#' @param lambda Numeric scalar controlling the TLP penalty (default = \code{log(n)/n}). Used as the constant part in TLP.
#' @param tau Numeric scalar controlling a threshold parameter in TLP (default = \code{0.3 * sqrt(1/(0.5*V))}).
#' @param espli A numeric convergence threshold for the relative change in the loss function across
#'   iterations (default = 1e-4).
#' @param silent Logical indicating whether to suppress iteration progress messages. Defaults to \code{FALSE}.
#' @param initial A character string specifying the initial estimation method for
#'   \code{(A, S, U)}. Possible values are \code{"eigen"} or \code{"prespecify"}.
#'   Defaults to \code{"eigen"}.
#' @param Theta_ini A list containing user-specified initial values for \code{A}, \code{S},
#'   \code{U}, \code{Fmat}, \code{sigma2}, and \code{phi2} if \code{initial = "prespecify"}.
#'   Ignored otherwise.
#' @return A list with the following components:
#' \describe{
#'   \item{\code{dat.LASERC}}{The harmonized (corrected) data matrix after adjusting for batch
#'   effects.}
#'   \item{\code{estimates}}{A list of final parameter estimates, including \code{A}, \code{S},
#'   \code{U}, \code{xFmat}, \code{sigma2}, \code{phi2}, \code{Fmat}, \code{gamma}, \code{beta},
#'   and \code{alpha}.}
#'   \item{\code{dat.original}}{The original data matrix (n x p) passed to the function.}
#'   \item{\code{batch}}{The batch vector used for correction.}
#' }
#'
#' @examples
#' \dontrun{
#' data(sim_data)
#' dat=sim_data$dat
#' batch=sim_data$batch
#' mod=sim_data$mod
#' n=nrow(dat)
#' p=ncol(dat)
#' V=(sqrt(1+8*p)-1)/2
#' results <- LASERC(
#'   dat = dat,
#'   batch = batch,
#'   mod = mod,
#'   L = 5
#' )
#' }
#' @importFrom stats model.matrix sd var weighted.mean
#' @export

LASERC<-function(dat,batch, mod=NULL,L=3,standardize=TRUE,center=NULL,scale=NULL, MaxIteration=100,lambda = log(n*V)/(n*V), tau= 0.5*sqrt(1/(V)),
                 espli=1e-4, silent=FALSE,initial='eigen',Theta_ini=NULL,penalty){

  n=nrow(dat);p=ncol(dat);ni=table(batch)
  M=length(ni);batch.id=unique(batch);n.batch=length(ni)
  V = (sqrt(1+8*p)-1)/2



  if(is.null(mod)){

    q=0
  }else{
    q=ncol(mod)-1
  }

  if ((M>1)&(q>0)){
    x=model.matrix(~mod[,-1]+batch)
  }else if(M==1){
    x=mod
  }else if(q==0){
    x=model.matrix(~batch)
  }

  dat.original=dat

  order.temp.batch=lapply(1:M,function(b) which(batch==batch.id[b]))
  if (standardize){
    if (is.null(center)&is.null(scale)){
      center=colMeans(dat)
      scale=apply(dat,2,sd)
    }
    dat=scale(dat)
  }

  datGrp <- split.data.frame(dat,batch)
  xGrp=split.data.frame(x,batch)

  H=x%*%solve(t(x)%*%x)%*%t(x)
  HGrp=split.data.frame(as.matrix(diag(H)),batch)



  if (initial=='eigen'){
    Theta_ini =  SCH_eigen_initial(dat,L,V)

    A = Theta_ini$A;
    S = Theta_ini$S
    U = Theta_ini$U

    AGrp <- split.data.frame(A,batch)


    Fmat=t(A)%*%x%*%solve(t(x)%*%x)
    xFmat=mapply(xFmat_func,xGrp,MoreArgs = list(Fmat=Fmat), SIMPLIFY = FALSE)

    sigma2=lapply(AGrp, function(x) apply(x,2,var))

    phi2=lapply(1:M,function(i) sum((datGrp[[i]]- AGrp[[i]]%*%S)^2)/(ni[i]*p))

    S_cross_S=tcrossprod(S)
    Sig_AconY=mapply(Sig_AconY_func,sigma2=sigma2,phi2=phi2,MoreArgs=list(S_cross_S=S_cross_S,p=p,L=L),SIMPLIFY = FALSE)

  }else if (initial=='warm'){

    A = Theta_ini$A;
    S = Theta_ini$S
    U = Theta_ini$U
    Fmat=Theta_ini$Fmat
    xFmat=mapply(xFmat_func,xGrp,MoreArgs = list(Fmat=Fmat), SIMPLIFY = FALSE)
    sigma2=Theta_ini$sigma2
    phi2=Theta_ini$phi2


    S_cross_S=tcrossprod(S)


    AGrp <- split.data.frame(A,batch)
    Sig_AconY=mapply(Sig_AconY_func,sigma2=sigma2,phi2=phi2,MoreArgs=list(S_cross_S=S_cross_S,p=p,L=L),SIMPLIFY = FALSE)



  }else if (initial=='warmL'){
    U_ini = Theta_ini$U
    A_ini = Theta_ini$A;
    S_ini = Theta_ini$S
    sigma2_ini=Theta_ini$sigma2
    xFmat_ini=Theta_ini$xFmat
    Fmat_ini=Theta_ini$Fmat

    L_ini=ncol(U_ini)
    L_upd=L-L_ini
    if (L<=V){
      Theta_ini = SCH_warmL_initial(dat,L,V,U_ini,A_ini,S_ini)
    }else{
      Theta_ini = SCH_warmL_random_initial(dat,L,V,U_ini,A_ini,S_ini)
    }

    A_upd = Theta_ini$A;
    S_upd = Theta_ini$S
    U_upd = Theta_ini$U
    #D = Theta_ini$D

    AGrp_upd <- split.data.frame(A_upd,batch)


    Fmat_upd=t(A_upd)%*%x%*%solve(t(x)%*%x)
    xFmat_upd=mapply(xFmat_func,xGrp,MoreArgs = list(Fmat=Fmat_upd), SIMPLIFY = FALSE)

    sigma2_upd=lapply(AGrp_upd, function(x) apply(x,2,var))

    Res=dat-A_ini%*%S_ini
    Res_Grp= split.data.frame(Res,batch)


    phi2=lapply(1:M,function(i) sum((Res_Grp[[i]]- AGrp_upd[[i]]%*%S_upd)^2)/(ni[i]*p))

    sigma2=lapply(1:M,function(i) c(sigma2_ini[[i]],sigma2_upd[[i]]))
    xFmat=lapply(1:M,function(i) rbind(xFmat_ini[[i]],xFmat_upd[[i]]))

    S=rbind(S_ini,S_upd)
    AGrp_upd <- split.data.frame(A_upd,batch)
    S_cross_S_upd=tcrossprod(S_upd)
    Sig_AconY_upd=mapply(Sig_AconY_func,sigma2=sigma2_upd,phi2=phi2,MoreArgs=list(S_cross_S=S_cross_S_upd,p=p,L=L_upd),SIMPLIFY = FALSE)

    S_cross_S=tcrossprod(S)

    A =cbind(A_ini,A_upd);
    U=cbind(U_ini,U_upd)
    Sig_AconY=mapply(Sig_AconY_func,sigma2=sigma2,phi2=phi2,MoreArgs=list(S_cross_S=S_cross_S,p=p,L=L),SIMPLIFY = FALSE)

  }



  loss=NULL

  S_cross_S=tcrossprod(S)
  Sig_AconY=mapply(Sig_AconY_func,sigma2=sigma2,phi2=phi2,MoreArgs=list(S_cross_S=S_cross_S,p=p,L=L),SIMPLIFY = FALSE)

  loss<-loss_func(sigma2=sigma2,phi2=phi2,S=S,S_cross_S=S_cross_S,xFmat=xFmat,datGrp,M=M,ni=ni)

  Iter = 1


  while(Iter <= MaxIteration){
    if(penalty=='OLS+warmL'){
      Theta_new = SCH_update_OLS(Res,A_upd,U_upd,batch,xFmat=xFmat_upd,x=x,sigma2=sigma2_upd,phi2=phi2,Sig_AconY=Sig_AconY_upd,HGrp=HGrp,penalt= penalty,lambda_ch = 0, tau=tau,silent = silent)

      S_upd = Theta_new$S;
      S_cross_S_upd=Theta_new$S_cross_S
      Sig_AconY_upd=Theta_new$Sig_AconY
      A_upd = Theta_new$A;
      U_upd=Theta_new$U
      #D_new=Theta_new$D
      xFmat_upd=Theta_new$xFmat
      Fmat_upd=Theta_new$Fmat
      sigma2_upd=Theta_new$sigma2
      phi2_upd=Theta_new$phi2


      S_new = rbind(S_ini,Theta_new$S);
      S_cross_S=tcrossprod(S_new)
      A_new =cbind(A_ini,Theta_new$A);
      U_new=cbind(U_ini,Theta_new$U)
      #D_new=Theta_new$D
      #xFmat_new=Theta_new$xFmat
      xFmat_new=lapply(1:M,function(i) rbind(xFmat_ini[[i]],Theta_new$xFmat[[i]]))

      Fmat_new=rbind(Fmat_ini, Theta_new$Fmat)
      sigma2_new=lapply(1:M,function(i) c(sigma2_ini[[i]],Theta_new$sigma2[[i]]))

      # sigma2_new=Theta_new$sigma2
      phi2_new=Theta_new$phi2


    }else if(penalty=='OLS'){
      Theta_new = SCH_update_OLS(dat,A,U,batch,xFmat=xFmat,x=x,sigma2=sigma2,phi2=phi2,Sig_AconY=Sig_AconY,HGrp=HGrp,penalt= penalty,lambda_ch = lambda, tau=tau,silent = silent,MaxIteration=MaxIteration,espli=espli)

      S_new = Theta_new$S;
      S_cross_S=Theta_new$S_cross_S
      Sig_AconY=Theta_new$Sig_AconY
      A_new = Theta_new$A;
      U_new=Theta_new$U
      #D_new=Theta_new$D
      xFmat_new=Theta_new$xFmat
      Fmat_new=Theta_new$Fmat
      sigma2_new=Theta_new$sigma2
      phi2_new=Theta_new$phi2


    }else if(penalty=='TLP'){
      Theta_new = SCH_update_TLP_full(dat,A,U,batch,xFmat=xFmat,x=x,sigma2=sigma2,phi2=phi2,Sig_AconY=Sig_AconY,HGrp=HGrp,penalt= penalty,lambda_ch =lambda, tau=tau,silent = silent,MaxIteration=MaxIteration,espli=espli,MaxIteration_U=25)

      S_new = Theta_new$S;
      S_cross_S=Theta_new$S_cross_S
      Sig_AconY=Theta_new$Sig_AconY
      A_new = Theta_new$A;
      U_new=Theta_new$U
      #D_new=Theta_new$D
      xFmat_new=Theta_new$xFmat
      Fmat_new=Theta_new$Fmat
      sigma2_new=Theta_new$sigma2
      phi2_new=Theta_new$phi2

    }



    loss_new=loss_func(sigma2=sigma2_new,phi2=phi2_new,S=S_new,S_cross_S=S_cross_S,xFmat=xFmat_new,datGrp,M=M,ni=ni)

    loss=c(loss,loss_new)

    err_loss=abs(loss_new- loss[length(loss)-1])/abs(loss[length(loss)-1])


    if(is.na(err_loss)>0){
      repara=alpha_gamma_func(Fmat,batch,n.batch,L,ni=ni)

      gamma=repara$gamma
      beta=repara$beta
      alpha=repara$alpha
      return(list(Conver=FALSE,A=A,S=S,U=U,xFmat=xFmat, sigma2=sigma2,phi2=phi2,gamma=gamma,beta=beta,alpha=alpha,Fmat=Fmat,loss=loss))}

    if(!silent){
      message(paste("Iter ",Iter,"; Percentage change on loss: " , format(round(err_loss,7), scientific = TRUE),".",sep=""))
    }
    A = A_new; S= S_new; U=U_new;

    xFmat=xFmat_new;
    sigma2=sigma2_new;
    Fmat=Fmat_new;
    phi2=phi2_new

    if(err_loss < espli)
    {

      break
    }
    Iter = Iter + 1


  }

  repara=alpha_gamma_func(Fmat,batch,n.batch,L,ni=ni)
  gamma=repara$gamma
  beta=repara$beta
  alpha=repara$alpha


  estimates=list(A=A,S=S,U=U,xFmat=xFmat,sigma2=sigma2,phi2=phi2,Fmat=Fmat,gamma=gamma,beta=beta,alpha=alpha)


  harmonized=harmonize_func(estimates,covariate=mod[,-1],ni,batch,dat,order.temp.batch,L,standardize,center,scale)

return(list(dat.LASERC=harmonized,estimates=estimates,dat.original=
              dat.original,batch=batch))

  }
