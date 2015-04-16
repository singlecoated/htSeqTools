rowLogRegLRT <- function(counts, exact=TRUE, p.adjust.method='none') {
  if ((!is.matrix(counts)) & (!is.data.frame(counts))) stop('Argument counts must be either a matrix or a data.frame')
  countsTotal <- rowSums(counts)
  n <- colSums(counts)
  p0 <- countsTotal/sum(countsTotal)
  l0 <- log(p0)*countsTotal + log(1-p0)*(sum(countsTotal)-countsTotal)
  p1 <- t(t(counts)/n)
  l1 <- log(p1)*counts + log(1-p1)*t(n-t(counts))
  l1 <- rowSums(l1,na.rm=TRUE)
  tstat <- 2*(l1-l0)
  pvals <- 1-pchisq(tstat,df=ncol(counts)-1)
  if (exact) {
    sel <- countsTotal<=5*ncol(counts)
    if (any(sel)) pvals[sel] <- apply(counts[sel,,drop=FALSE],1,function(z) chisq.test(rbind(z,n-z),simulate.p.value=TRUE)$p.value)
  }
  pvals <- p.adjust(pvals,method=p.adjust.method)
  data.frame(tstat=tstat,pvals=pvals)
}
