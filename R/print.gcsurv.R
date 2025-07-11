print.gcsurv <- function (x, digits=4, ...)
{
  cat("Call:", "\n", sep = "")
  dput(x$formula)
  cat("\n")
  
  tmp <- matrix(c(x$RMST$mean.deltaRMST, x$RMST$se.deltaRMST,
                  x$RMST$mean.deltaRMST/x$RMST$se.deltaRMST, x$RMST$p.value.deltaRMST),nrow=1)
  colnames(tmp) <- c("deltaRMST","se(deltaRMST)","z", "p")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, signif.stars = FALSE, ...)
  cat("\n")
  
  tmp <- matrix(c(x$AHR$mean.AHR, x$AHR$se.AHR,
                  x$AHR$mean.AHR/x$AHR$se.AHR, x$AHR$p.value.AHR),nrow=1)
  colnames(tmp) <- c("AHR","se(AHR)","z", "p")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, signif.stars = FALSE, ...)
  cat("\n")
  
  cat(paste0("n= ",x$n,", number of events= ",x$nevent))
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep="");cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep="");cat("\n") }
}