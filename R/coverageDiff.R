coverageDiff <- function(sample1, sample2, chrLength) {
  sample1 <- sample1[names(sample1) %in% names(chrLength)]
  sample2 <- sample2[names(sample2) %in% names(chrLength)]
  chrLengthSel <- as.list(chrLength[names(sample1)]); chrLengthSel <- chrLengthSel[!is.na(chrLengthSel)]
  cover1 <- coverage(sample1, width=chrLengthSel)
  chrLengthSel <- as.list(chrLength[names(sample2)]); chrLengthSel <- chrLengthSel[!is.na(chrLengthSel)]
  cover2 <- coverage(sample2, width=chrLengthSel)
  #
  sel <- names(cover1) %in% names(cover2)
  n <- names(cover1)[sel]
  coverdiff <- cover1[n]-cover2[n]
  n1 <- names(cover1)[!sel]
  coverdiff1 <- cover1[n1]
  n2 <- names(cover2)[!(names(cover2) %in% names(cover1))]
  if (length(n2)>0) coverdiff2 <- -1*cover2[n2] else coverdiff2 <- cover2[n2]
  #
  ans <- c(coverdiff,coverdiff1,coverdiff2)
  names(ans) <- c(n,n1,n2)
  return(ans)
}
