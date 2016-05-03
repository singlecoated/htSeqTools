
setMethod("listOverlap", signature(list1='character', list2='character', list3='missing', univ='character'),
function(list1, list2, univ, ...) {
  list1 <- unique(list1[!is.na(list1)])
  list2 <- unique(list2[!is.na(list2)])
  univ <- unique(univ[!is.na(univ)])   
  nboth <- sum(list1 %in% list2)
  nonly1 <- length(list1)-nboth
  nonly2 <- length(list2)-nboth
  nnone <- sum((!(univ %in% list1)) & (!(univ %in% list2)))
  xtab <- matrix(c(nnone,nonly1,nonly2,nboth),nrow=2,byrow=TRUE)
  rownames(xtab) <- c('notList1','List1')
  colnames(xtab) <- c('notList2','List2')
  ans <- chisq.test(xtab, ...)
  return(ans)
}
)

# linear models version
# setMethod("listOverlap", signature(list1='character', list2='character', list3='missing', univ='character'),
# function(list1, list2, univ, ...) {
#   list1 <- unique(list1[!is.na(list1)])
#   list2 <- unique(list2[!is.na(list2)])
#   univ <- unique(univ[!is.na(univ)])
#   xtab <- data.frame(Univ=univ, List1=as.integer(univ %in% list1), List2=as.integer(univ %in% list2))
#   ftable <- as.data.frame(table(xtab[,c('List1','List2')]))
#   glm1 <- glm(Freq ~ List1 * List2, family='poisson',data=ftable)
#   glm2 <- glm(Freq ~ List1 + List2, family='poisson',data=ftable)
#   anovatest <- anova(glm1,glm2,test='Chisq')
#   pvalue <- anovatest[['P(>|Chi|)']][2]
#   ans <- list(xtab,ftable,glm1,glm2,pvalue)
#   names(ans) <- c('xtab','ftable','glm1','glm2','pvalue')
#   return(ans)
# }
# )

setMethod("listOverlap", signature(list1='character', list2='character', list3='character', univ='character'),
function(list1, list2, list3, univ, ...) {
  list1 <- unique(list1[!is.na(list1)])
  list2 <- unique(list2[!is.na(list2)])
  list3 <- unique(list3[!is.na(list3)])
  univ <- unique(univ[!is.na(univ)])
  xtab <- data.frame(Univ=univ, List1=as.integer(univ %in% list1), List2=as.integer(univ %in% list2), List3=as.integer(univ %in% list3))
  ftable <- as.data.frame(table(xtab[,c('List1','List2','List3')]))
  glm1 <- glm(Freq ~ List1 * List2 * List3, family='poisson',data=ftable)
  glm2 <- glm(Freq ~ List1*List2 + List1*List3 + List2*List3, family='poisson',data=ftable)
  anovatest <- anova(glm1,glm2,test='Chisq')
  pvalue <- anovatest[['P(>|Chi|)']][2]
  ans <- list(xtab,ftable,glm1,glm2,pvalue)
  names(ans) <- c('xtab','ftable','glm1','glm2','pvalue')
  return(ans)
}
)

setMethod("listOverlap", signature(list1='factor', list2='character', list3='missing', univ='missing'),
function(list1, list2, univ, ...) {
  if (is.null(names(list1))) stop("list1 must have names for signature('numeric','numeric','missing')")
  inlist2 <- names(list1) %in% list2
  xtab <- table(list1, inlist2)
  colnames(xtab) <- c('notInList2','List2')
  ans <- chisq.test(xtab, ...)
  return(ans)
}
)
