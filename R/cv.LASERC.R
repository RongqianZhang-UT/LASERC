#' @title Cross-validation for selecting L in LASERC
#' @description
#' \code{cv.LASERC} performs k-fold cross-validation over a grid of possible
#' ranks (\code{L_set}) for the \code{LASERC} function.
#' @param dat Numeric matrix with subjects as rows,
#'                and functional connectivity features as columns.
#' @param batch Numeric or character vector specifying the batch/site
#'     variable needed for harmonization.
#' @param mod Optional model matrix for outcome of interest and
#'     other covariates besides batch/scanner.
#' @param standardize  Logical indicating whether to standardize the columns of \code{dat} prior to
#'   model fitting. Defaults to \code{TRUE}.
#' @param center A numeric vector of length \code{p} giving column means for centering (if \code{standardize=TRUE}).
#'   If \code{NULL}, column means of \code{dat} are used.
#' @param scale A numeric vector of length \code{p} giving column standard deviations for scaling (if \code{standardize=TRUE}).
#'   If \code{NULL}, column standard deviations of \code{dat} are used.
#' @param \dots Additional parameters passed to \code{\link{LASERC}}, such as \code{lambda},
#'   \code{tau}, or \code{MaxIteration}.
#' @param L_set A numeric vector of candidate \code{L} to be evaluated via cross-validation.
#' @param K An integer specifying the number of folds for cross-validation. Default is 3.
#' @param parallel Logical indicating whether to parallelize the cross-validation loop
#'   using \pkg{doParallel} and \pkg{foreach}. Defaults to \code{TRUE}.
#' @param ncores An integer specifying the number of CPU cores to register when
#'   \code{parallel = TRUE}. Default is 3.
#' @param seed An integer seed used for reproducibility in creating folds. Default is 1234.
#' @return An object of class \code{"cv.LASERC"}, which is a list containing:
#' \describe{
#'   \item{\code{L_set}}{The vector of \code{L} values tested.}
#'   \item{\code{cv.mean}}{A numeric vector of mean loss across folds for each \code{L}.}
#'   \item{\code{cv.se}}{A numeric vector of standard deviations of the loss across folds
#'   for each \code{L}.}
#'   \item{\code{L_est}}{The chosen \code{L} (the value that minimizes the cross-validation loss).}
#' }
#' @examples
#' \dontrun{
#' data(sim_data)
#' dat=sim_data$dat
#' batch=sim_data$batch
#' mod=sim_data$mod
#' # Candidate L values
#' L_candidates <- 1:10
#'
#' # Perform 3-fold CV for each candidate L
#' cv_results <- cv.LASERC(
#'   dat = dat,
#'   batch = batch,
#'   mod = mod,
#'   L_set = L_candidates,
#'   K = 3,
#'   parallel = FALSE  # set to TRUE to test parallel mode
#' )
#' plot(cv_results)
#' }

#' @importFrom caret createFolds
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @export

cv.LASERC<-function(dat,batch, mod,standardize=TRUE,center=NULL,scale=NULL,...,L_set,K=3,parallel=TRUE,ncores=3,seed=1234){
  n=nrow(dat);p=ncol(dat);ni=table(batch)
  if(is.null(mod)){
    q=0
  }else{
    q=ncol(mod)-1
  }

  M=length(ni);batch.id=unique(batch)
  V = (sqrt(1+8*p)-1)/2
  set.seed(seed)
  folds <- createFolds(batch, k = K, list = TRUE)
  cv.res=NULL
  for (l in L_set){
    if (parallel){

      registerDoParallel(cores=ncores)
      cv.res[[l]]<- foreach (k = 1:K, .combine = 'c',.packages  = c('LASERC'))  %dopar% {

        test_indices <- folds[[k]]
        train_indices <- setdiff(seq_len(length(batch)), test_indices)

        # Split the data into training and testing sets
        train_mod <- mod[train_indices,,drop=FALSE ]
        test_mod <- mod[test_indices, ,drop=FALSE]

        train_dat <- dat[train_indices,,drop=FALSE ]
        test_dat <- dat[test_indices,,drop=FALSE ]

        if (standardize){
          center=colMeans(train_dat)
          scale=apply(train_dat,2,sd)

          test_dat=scale(test_dat,center=center,scale=scale)
        }

        train_batch<-batch[train_indices ]
        test_batch<-batch[test_indices ]
        train_ni=table(train_batch)
        test_ni=table(test_batch)


        RESULT_train<-LASERC(dat=train_dat,mod=train_mod,batch=train_batch,standardize = standardize,center=center,scale=scale, L=l, lambda = log(sum(train_ni))/(sum(train_ni)),tau=0.3 * sqrt(1/(0.5*V)),silent = TRUE,
        ...)

        if ((M>1)&(q>0)){
          test_X=model.matrix(~test_mod[,-1]+test_batch)
        }else if(M==1){
          test_X=test_mod
        }else{
          test_X=model.matrix(~test_batch)
        }

        test_XGrp=split.data.frame(test_X,test_batch)
        test_XFmat=mapply(XFmat_func,test_XGrp,MoreArgs = list(Fmat=RESULT_train$estimates$Fmat), SIMPLIFY = FALSE)


        test_datGrp=split.data.frame(test_dat,test_batch)
        S_cross_S=tcrossprod(RESULT_train$estimates$S)
        loss=loss_func(RESULT_train$estimates$sigma2,RESULT_train$estimates$phi2,RESULT_train$estimates$S,S_cross_S=S_cross_S,test_XFmat,datGrp=test_datGrp,M=M,ni=test_ni)

        loss
      }
      stopImplicitCluster()
      } else {
        loss_L=NULL
        for (k in 1:K){
          test_indices <- folds[[k]]
          train_indices <- setdiff(seq_len(length(batch)), test_indices)

          # Split the data into training and testing sets
          train_mod <- mod[train_indices,,drop=FALSE ]
          test_mod <- mod[test_indices, ,drop=FALSE]

          train_dat <- dat[train_indices,,drop=FALSE ]
          test_dat <- dat[test_indices,,drop=FALSE ]


          if (standardize){
            center=colMeans(train_dat)
            scale=apply(train_dat,2,sd)

             test_dat=scale(test_dat,center=center,scale=scale)
          }

          train_batch<-batch[train_indices ]
          test_batch<-batch[test_indices ]
          train_ni=table(train_batch)
          test_ni=table(test_batch)


          RESULT_train<-LASERC(dat=train_dat,mod=train_mod,standardize=standardize,center=center,scale=scale,batch=train_batch, L=l, lambda = log(sum(train_ni))/(sum(train_ni)),tau=0.3 * sqrt(1/(0.5*V)),silent = TRUE,
          ...)

          if ((M>1)&(q>0)){
            test_X=model.matrix(~test_mod[,-1]+test_batch)
          }else if(M==1){
            test_X=test_mod
          }else{
            test_X=model.matrix(~test_batch)
          }

          test_XGrp=split.data.frame(test_X,test_batch)
          test_XFmat=mapply(XFmat_func,test_XGrp,MoreArgs = list(Fmat=RESULT_train$estimates$Fmat), SIMPLIFY = FALSE)


          test_datGrp=split.data.frame(test_dat,test_batch)
          S_cross_S=tcrossprod(RESULT_train$estimates$S)
          loss=loss_func(RESULT_train$estimates$sigma2,RESULT_train$estimates$phi2,RESULT_train$estimates$S,S_cross_S=S_cross_S,test_XFmat,datGrp=test_datGrp,M=M,ni=test_ni)
          loss_L=c(loss_L,loss)

        }
        cv.res[[l]]=loss_L
      }
    }



  cv.res=do.call('cbind',cv.res)

  cv.mean=colMeans(cv.res)

  cv.se=apply(cv.res,2,sd)

  L_est=L_set[which.min(cv.mean)]


  out <- structure(list(
                         L_set=L_set,
                        cv.mean = cv.mean,
                        cv.se = cv.se,
                        L_est =  L_est

  ),
  class = "cv.LASERC")

  return(out)
}

