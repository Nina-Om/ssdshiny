genFits <- function(x, type = "norm", method="mge", gof="CvM"){

  ini.par = mledist(x[x >= 0], "weibull")
  
  fit <- switch(type,
         "norm" = fitdist(x, "norm", method = method, gof = gof, start = list(mean=mean(x), sd=sd(x))),
         "logis" =  fitdist(x, "logis", method = method, gof = gof),
         "extreme" = fextreme(x, 'extreme', method=method, gof=gof, start=list(mean=mean(x), sd=sd(x))),
         "gumbel" = fgumbel(x, 'gumbel', method=method, gof=gof, start=list(mean=mode(x), sd=0.78*sd(x))),
         "weibull" = fitdist(x, "weibull", method = method, gof = gof, start=list(shape=ini.par$estimate[[1]], scale=ini.par$estimate[[2]])),
  return(fit)
)}

dgumbel <- function(x, mean, sd) 1/sd*exp((mean-x)/sd)*exp(-exp((mean-x)/sd)) #f
pgumbel <- function(q, mean, sd) exp(-exp((mean-q)/sd)) #F
qgumbel <- function(p, mean, sd) mean-sd*log(-log(p))
rgumble <- function(n) log(-log(runif(n)))

fgumbel <- function(x, type = "gumbel" , method="mge", gof="CvM", start=list(mean=mode(x), sd=0.78*sd(x))) {
  fitdist(x, 'gumbel', method=method, gof=gof, start=list(mean=mode(x), sd=0.78*sd(x)))
}

dextreme <- function(x, mean, sd) 1/sd*exp((-mean+x)/sd)*exp(-exp((-mean+x)/sd))
pextreme <- function(q, mean, sd) 1-exp(-exp((-mean+q)/sd))
qextreme <- function(p, mean, sd) mean+sd*log(-log(1-p))
rextreme <- function(n) log(-log(1-runif(n)))

fextreme <- function(x, type = 'extreme', method="mge", gof="CvM", start=list(mean=mean(x), sd=sd(x))) {
  fitdist(x, 'extreme', method=method, gof=gof, start=list(mean=mean(x), sd=sd(x)))
}

 # getParams <- function(fit, bootmethod = "nonparam", niter = 5000){
 #   b4 <-bootdist(fit, bootmethod = bootmethod, niter=niter)
 #   sumb4 <- summary(b4)
 #   par41 <- sumb4$fitpart$estimate[1]
 #   par42 <- sumb4$fitpart$estimate[2]
 #   return(list(par1 = par41, par2 = par42))
 # }

 getParams <- function(fit.nor){
   sumb4 <- summary(fit.nor)
   par41 <- sumb4$estimate[1]
   par42 <- sumb4$estimate[2]
   return(list(par1 = par41, par2 = par42))
 }

mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
} 

SSD <- function(x, p.sim) {
  inv.ecdf <- function(x){return(ecdf(x)(x))}
  p.data <- inv.ecdf(x) 
  ssd <- sum((p.sim - p.data)^2)
  return(ssd)
}

MSE <- function(x, p.sim) {
  inv.ecdf <- function(x){return(ecdf(x)(x))}
  p.data <- inv.ecdf(x)
  ssd <- sum((p.sim - p.data)^2)
  mse <- ssd/(length(x)-2)
return(mse)
}

hc <- function(fit, bootfit, p) {
  hc <- quantile(fit, probs = p)
  hcq <- quantile(bootfit, probs = p)
  return(list(hc=hc, hcq=hcq))
}

HC <- function(model.name, scale, hc.list, p) {
  if(scale == 2) {
    HCp <- as.data.frame(sapply(hc.list, function(i) c(10^i$hc$quantile$p, 10^i$hcq$quantCI$p)))
  }else{
    HCp <- as.data.frame(sapply(hc.list, function(i) c(i$hc$quantile$p, i$hcq$quantCI$p)))  
  }
  colnames(HCp) <- model.name
  rownames(HCp) <- list(paste0("HC", p*100), paste0("HC", p*100, "L"), paste0("HC", p*100, "H"))
  HCp <- format(HCp, digits=3, nsmall = 3, justify="left") 
  return(HCp)
}

fx <- function(scale, x, y) ifelse(scale == 2, return(x), return(y))

CI.boot <- function(scale, b1, CI.level) {
  qboot <- quantile(b1, probs = seq(0.001,0.999,0.001), CI.type = "two.sided", CI.level = CI.level)
  cis <- qboot$quantCI
  cis <- fx(scale, 10^cis, cis)
  Pr <- qboot$probs
  rownames(cis) <- c('lwr' ,'upr')
  CI.Pr <- data.frame(cbind(Pr, t(cis)))
  return(CI.Pr)
}