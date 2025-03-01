#' @title Plot method for the \code{cv.LASERC} function.
#' @description
#' \code{plot.cv.LASERC} provides a convenient way to visualize the cross-validation 
#' loss (with error bars) across different values of \code{L}. The plot also indicates 
#' the chosen \code{L_est} via a vertical dashed line.
#'
#' @param x An object of class \code{"cv.LASERC"} produced by \code{\link{cv.LASERC}}.
#' @param vertical.line Logical indicating whether to draw a vertical dashed line at 
#'   the chosen \code{L_est}. Defaults to \code{TRUE}.
#' @param \dots Additional graphical parameters passed to \code{\link[ggplot2]{ggplot}}.
#'
#' @details
#' If only one \code{L} value is tested, a warning is given and no plot is produced. 
#' Otherwise, the cross-validation mean loss and +/- 1 standard error are displayed.
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @method plot cv.LASERC
#' @export
#' 



plot.cv.LASERC<- function(x, vertical.line=TRUE, ...) {
  
  if (length(x$cv.mean) <= 1) {
    warning("No plot produced, since the length of cv sequence is <= 1.")
    return ()
  }
  df <- NULL
  df$index <- as.integer(x$L_set)
  df$cvm <- x$cv.mean
  df$cvl <- x$cv.mean - x$cv.se
  df$cvu <- x$cv.mean + x$cv.se
  df=as.data.frame(df)
  xlab='est L'
  ylab <- "Loss"
  title <- "Cross-validation Plot"
  g <- ggplot(df, aes(x = index, y = cvm)) + 
    xlab(xlab) + ylab(ylab) + ggtitle(title) + 
    geom_line(lwd=1, color = "red") + geom_point(pch=19, color = "red") + 
    theme_bw() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    geom_errorbar(aes(ymin = cvl, ymax = cvu), width=0.02, color = "gray")
  
  if (vertical.line) {
    g <- g + geom_vline(xintercept = x$L_est, linetype = "dashed")
  }
  
  g
}
