library(fitdistrplus)
library(MASS)
library(fBasics)
library(EnvStats)
library(readxl)
library(gridExtra)
library(broom)
library(latticeExtra)
library(knitr)
library(docxtools)

Curve.fitting <- function(Data, unit, CI.level, scale, exposure, samplsize=500) {
  source("genFunctions.R")
  
  # Data preparation (columns name in the uploaded data file are user defined)
  # Data <- read_excel("SSDdata.xlsx"); scale <- 2; unit <- "ug/l"; CI.level <- 95; exposure=2.5
  CI.level <- CI.level * 0.01
  
  ## rank=rank(LConc)
  ## cumdens.Hazen=(rank-0.5)/20
  ## cumdens.Filibens=(rank-0.3175)-(20+0.365)
  
  colnames(Data)<- c("X to include","Taxa Grouping","Species", "Concentration")
  Data <- Data[order(Data$Concentration),]
  Conc <- Data[[4]]
  LConc = log10(Conc)
  #fx <- function(scale, x, y) {ifelse(scale == "Log", return(x), return(y))}
  x <- fx(scale, LConc, Conc)
  summary(x)
  
  # Histogram and Kernel density of Observed data 
  par("mar") 
  par(mar=c(4,4,4,4))
  par(mfrow=c(1,2))
  hist(Conc, prob=T, cex.main=0.7, cex.axis=0.6, cex.lab=0.6, main=paste("Histogram and Kernel density \n(Arithmetic scale)"), xlim = c(min(Conc)-0.2*min(Conc),max(Conc)+0.2*max(Conc)), xlab = paste0("Concentration ", "(", unit , ")"), ylab="Density")
  lines(density(Conc), col="blue", lwd=2)
  hist(LConc, prob=T, cex.main=0.7, cex.axis=0.6, cex.lab=0.6, main="Histogram and Kernel density \n(Log scale)", xlim = c(min(LConc)+0.8*min(LConc),max(LConc)+0.8*max(LConc)), xlab = paste0("Log10(Concentration) ", "(" , unit , ")"), ylab="Density")
  lines(density(LConc), col="blue", lwd=2)
  
  p1 <- ggplot(Data, aes(Conc)) + 
    stat_ecdf(geom = "point", pad = FALSE) + 
    ggtitle("Scatterplot of Data in Arithmetic Space") + 
    theme(plot.title = element_text(size=8, hjust=0), axis.text=element_text(size=10),
          axis.title=element_text(size=8,face="bold"))+
    labs(x = paste0("Concentration ", "(", unit , ")"), y ="Proportion of Species")
  p2 <- ggplot(Data, aes(LConc)) + 
    stat_ecdf(geom = "point", pad = FALSE) + 
    ggtitle("Scatterplot of Data in Log Space") + 
    theme(plot.title = element_text(size=8, hjust=0), axis.text=element_text(size=10),
          axis.title=element_text(size=8,face="bold"))+
    labs(x = paste0("Log10(Concentration) ", "(", unit , ")"), y ="Proportion of Species")
  grid.arrange(p1, p2, nrow = 1) 

# Parameters’ estimate
  model.name <- list("Normal", "Logistic", "Extreme.Value", "Gumbel", "Weibull")
  fit.nor <- genFits(x, type = "norm", method="mge", gof="CvM")
  fit.log <- genFits(x, type = "logis", method="mge", gof="CvM")
  fit.ext <- genFits(x, type = "extreme", method="mge", gof="CvM")
  fit.gum <- genFits(x, type = "gumbel", method="mge", gof="CvM")
  fit.wbl <- genFits(x, type = "weibull", method="mge", gof="CvM")
  par1 <- getParams(fit.nor)
  par2 <- getParams(fit.log)
  par3 <- getParams(fit.ext)
  par4 <- getParams(fit.gum)
  par5 <- getParams(fit.wbl)
# Bootstrap estimations (65 seconds)
  fit.list <- list(fit.nor, fit.log, fit.ext, fit.gum, fit.wbl)
  b1 <- bootdist(fit.nor, bootmethod = "nonparam", niter=500)
  b2 <- bootdist(fit.log, bootmethod = "nonparam", niter=500)
  b3 <- bootdist(fit.ext, bootmethod = "nonparam", niter=500)
  b4 <- bootdist(fit.gum, bootmethod = "nonparam", niter=500)
  b5 <- bootdist(fit.wbl, bootmethod = "nonparam", niter=500)
  fit.bootlist <- list(b1,b2,b3,b4,b5)

# get Confidence Intervals from bootstraps 
  ci1 <- CI.boot(scale, b1, CI.level)
  ci2 <- CI.boot(scale, b2, CI.level)
  ci3 <- CI.boot(scale, b3, CI.level)
  ci4 <- CI.boot(scale, b4, CI.level)
  ci5 <- CI.boot(scale, b4, CI.level)
  CI <- list(ci1, ci2, ci3, ci4, ci5)
  names(CI) <- model.name
  CI <- as.data.frame(CI)
  
# New fits based on bootstrap estimations 
  fit.nor <- fitdist(x, "norm", method="mge", gof="CvM", start = list(mean=b1$fitpart$estimate[[1]], sd=b1$fitpart$estimate[[2]]))
  fit.log <- fitdist(x, "logis", method="mge", gof="CvM", start = list(location=b2$fitpart$estimate[[1]], scale=b2$fitpart$estimate[[2]]))
  fit.ext <- fitdist(x, "extreme", method="mge", gof="CvM", start=list(mean=b3$fitpart$estimate[[1]], sd=b3$fitpart$estimate[[2]]))
  fit.gum <- fitdist(x, "gumbel", method="mge", gof="CvM", start=list(mean=b4$fitpart$estimate[[1]], sd=b4$fitpart$estimate[[2]]))
  fit.wbl <- fitdist(x, "weibull", method="mge", gof="CvM", start=list(shape=b5$fitpart$estimate[[1]], scale=b5$fitpart$estimate[[2]]))
  par1 <- getParams(fit.nor)
  par2 <- getParams(fit.log)
  par3 <- getParams(fit.ext)
  par4 <- getParams(fit.gum)
  par5 <- getParams(fit.wbl)
  # Table of fitted distributions parameters
  par.fit <- as.data.frame(sapply(fit.list, function(i) c(i$estimate[[1]], i$estimate[[2]], i$loglik)))
  colnames(par.fit) <- model.name
  rownames(par.fit) <- list("location", "scale", "loglik")
  par.fit <- format(par.fit, digits=3, nsmall = 3, justify="left") 
  par.fit
  
  # Goodness-of-fit statistics
  gof.stat <- gofstat(fit.list, fitnames=model.name)
  par.gof <- as.data.frame(cbind(gof.stat$ad[1:length(model.name)], rank(gof.stat$ad[1:length(model.name)]), gof.stat$adtest[1:length(model.name)], 
                                 gof.stat$ks[1:length(model.name)],  rank(gof.stat$ks[1:length(model.name)]), gof.stat$kstest[1:length(model.name)], 
                                 gof.stat$chisq[1:length(model.name)],  rank(gof.stat$chisq[1:length(model.name)]), gof.stat$chisqpvalue[1:length(model.name)]))
  colnames(par.gof) <- c("Anderson-Darling", "Rank", "Reject", "Kolmogrov-Smirnov", "Rank", 
                         "Reject", "Chi-Squared", "Rank", "Chisq p-value")
  
  # Plots of Empirical and theoretical CDFs
  pro = seq(0.001,0.999,0.001)

   x.simul1 <- (qnorm(pro, par1$par1[[1]], par1$par2[[1]]))
   x.simul2 <- (qlogis(pro, par2$par1[[1]], par2$par2[[1]]))
   x.simul3 <- (qextreme(pro, par3$par1[[1]], par3$par2[[1]]))
   x.simul4 <- (qgumbel(pro, par4$par1[[1]], par4$par2[[1]]))
   x.simul5 <- (qweibull(pro, par5$par1[[1]], par5$par2[[1]]))
   
   x.simul1 <- fx(scale, 10^x.simul1, x.simul1)
   x.simul2 <- fx(scale, 10^x.simul2, x.simul2)
   x.simul3 <- fx(scale, 10^x.simul3, x.simul3)
   x.simul4 <- fx(scale, 10^x.simul4, x.simul4)
   x.simul5 <- fx(scale, 10^x.simul5, x.simul5)
   
   Table.sim <- data.frame(cbind(pro, x.simul1, x.simul2, x.simul3, x.simul4, x.simul5, CI))
   colnames(Table.sim)[1:(length(model.name)+1)] <- c("pro", "Normal", "Logistic", "Extreme.Value", "Gumbel", "Weibull")
   Table.sim$Normal
   
   # dev.off()
   # xx <- fx(scale, 10^x, x)
   # par(mar=c(4,4,4,4), mgp=c(2,1,0))
   # plot(ecdf(xx), axes = F, main=NA, cex.lab=0.6, col="black", xlab="Concentrations (units)",
   #      ylab= "Proportion of Taxa Affected", 
   #      log="x", xlim=c(0.00001,1000*max(xx)))
   # title(main= "Empirical and theoretical CDFs", line=1, cex.main=0.8)
   # axis(1, cex.axis = 0.6)
   # axis(2, cex.axis = 0.6)
   # lines(Table.sim[,2],pro, col="red")
   # lines(Table.sim[,3],pro, col="orange")
   # lines(Table.sim[,4],pro, col="purple")
   # lines(Table.sim[,5],pro, col="green")
   # lines(Table.sim[,6],pro, col="black")
   # lines(CI$Gumbel.lwr, pro, col="blue")
   # legend("bottomright", legend=c("Normal", "Logistic", "Extreme Value", "Gumbel", "Weibull", "lwr"),
   #        col=c("red", "orange", "purple", "green", "black", "blue"), lty=1:5, cex=0.8)
  
   
   # dev.off()
   # xx <- fx(scale, 10^x, x)
   # par(mar=c(4,4,4,4), mgp=c(2,1,0))
   # plot(ecdf(xx), axes = F, main=NA, cex.lab=0.6, col="black", xlab="Concentrations (units)",
   #      ylab= "Proportion of Taxa Affected", 
   #      xlim=c(-0.5*max(xx),max(xx)))
   # title(main= "Empirical and theoretical CDFs", line=1, cex.main=0.8)
   # axis(1, cex.axis = 0.6)
   # axis(2, cex.axis = 0.6)
   # lines((Table.sim[,2]),pro, col="red")
   # lines(Table.sim[,3],pro, col="orange")
   # lines(Table.sim[,4],pro, col="purple")
   # lines(Table.sim[,5],pro, col="green")
   # lines(CI$Gumbel.lwr, pro, col="blue")
   # legend("bottomright", legend=c("Normal", "Logistic", "Extreme Value", "Gumbel", "lwr"),
   #        col=c("red", "orange", "purple", "green", "blue"), lty=1:5, cex=0.8)

  # SSE and MSE calculations
  p.sim1 <- pnorm(x, par1$par1[[1]], par1$par2[[1]])
  p.sim2 <- plogis(x, par2$par1[[1]], par2$par2[[1]])
  p.sim3 <- pextreme(x, par3$par1[[1]], par3$par2[[1]])
  p.sim4 <- pgumbel(x, par4$par1[[1]], par4$par2[[1]])
  p.sim5 <- pweibull(x, par5$par1[[1]], par5$par2[[1]])
  
  ssd1 <- SSD(x, p.sim1)
  mse1 <- MSE(x, p.sim1)
  ssd2 <- SSD(x, p.sim2)
  mse2 <- MSE(x, p.sim2)
  ssd3 <- SSD(x, p.sim3)
  mse3 <- MSE(x, p.sim3)
  ssd4 <- SSD(x, p.sim4)
  mse4 <- MSE(x, p.sim4)
  ssd5 <- SSD(x, p.sim5)
  mse5 <- MSE(x, p.sim5)
  
  gof.measurs <- as.data.frame(rbind(cbind(ssd1, ssd2, ssd3, ssd4, ssd5), cbind(mse1, mse2, mse3, mse4, mse5))) 
  colnames(gof.measurs) <- model.name
  rownames(gof.measurs) <- list("SSD", "MSE")                          
  gof.measurs <-  format(gof.measurs, digits=3, nsmall = 3, justify="left")
  gof.measurs
  
  # hazardous concentration
  p=0.05
  hc1 <- hc(fit.list[[1]], fit.bootlist[[1]], p)
  hc2 <- hc(fit.list[[2]], fit.bootlist[[2]], p)
  hc3 <- hc(fit.list[[3]], fit.bootlist[[3]], p)
  hc4 <- hc(fit.list[[4]], fit.bootlist[[4]], p)
  hc5 <- hc(fit.list[[5]], fit.bootlist[[5]], p)
  hc.list <- list(hc1, hc2, hc3, hc4, hc5)
  HC5 <- HC(model.name, scale, hc.list, p)
  
  p=0.50
  hc1 <- hc(fit.list[[1]], fit.bootlist[[1]], p)
  hc2 <- hc(fit.list[[2]], fit.bootlist[[2]], p)
  hc3 <- hc(fit.list[[3]], fit.bootlist[[3]], p)
  hc4 <- hc(fit.list[[4]], fit.bootlist[[4]], p)
  hc5 <- hc(fit.list[[5]], fit.bootlist[[5]], p)
  hc.list <- list(hc1, hc2, hc3, hc4, hc5)
  HC50 <- HC(model.name, scale, hc.list, p) 
  
  hcTable <- rbind(HC5, HC50[1,])
  
  # Fraction Affected
  fa1 <- pnorm(log10(exposure), par1$par1[[1]], par1$par2[[1]])
  fa2 <- plogis(log10(exposure), par2$par1[[1]], par2$par2[[1]])
  fa3 <- pextreme(log10(exposure), par3$par1[[1]], par3$par2[[1]])
  fa4 <- pgumbel(log10(exposure), par4$par1[[1]], par4$par2[[1]])
  fa5 <- pgumbel(log10(exposure), par5$par1[[1]], par5$par2[[1]])
  fa.Table <- data.frame(cbind(fa1, fa2, fa3, fa4, fa5))
  colnames(fa.Table) <- model.name
  rownames(fa.Table) <- c("Fraction Affected")
  fa.Table <-  format(fa.Table, digits=3, nsmall = 3, justify="left")

  # Additional plots
  # CDF, pdf
  # dev.off()
  # c <- cdfcomp(fit.list, plotstyle = "ggplot", ylim= c(0,1), xlab="Concentrations (units)", ylab= "Proportion of Taxa Affected", legendtext=model.name)
  # denscomp(fit.list, plotstyle = "ggplot", xlab="Concentrations (units)", legendtext=model.name)
  # qqcomp(fit.list, xlab="Concentrations (units)", ylab= "Proportion of Taxa Affected", legendtext=model.name)
  # ppcomp(fit.list, legendtext=model.name)

   ## Bootstap dist. plots
    # par("mar") 
    # par(mar=c(2,2,2,2))
    # par(mfrow=c(2,2))
    # CIcdfplot(b1, CI.output = "quantile", CI.level = CI.level,  xlogscale = F, xlab="Concentrations (units)", ylab= "Proportion of Taxa Affected", main= "Normal distribution")
    # CIcdfplot(b2, CI.output = "quantile", CI.level = CI.level, xlogscale = F, xlab="Concentrations (units)", ylab= "Proportion of Taxa Affected", main= "Logistic distribution")
    # CIcdfplot(b3, CI.output = "quantile", CI.level = CI.level, xlogscale = F, xlab="Concentrations (units)", ylab= "Proportion of Taxa Affected", main= "Extreme Value distribution")
    # CIcdfplot(b4, CI.output = "quantile", CI.level = CI.level, xlogscale = F, xlab="Concentrations (units)", ylab= "Proportion of Taxa Affected", main= "Gumbel distribution")
    # CIcdfplot(b5, CI.output = "quantile", CI.level = CI.level, xlogscale = F, xlab="Concentrations (units)", ylab= "Proportion of Taxa Affected", main= "Weibull distribution")
    # 
  # plot(fit.nor)
  # plot(fit.log)
  # plot(fit.ext)
  # plot(fit.gum)
  
  #Table.sim <- format_engr(Table.sim)
  #Table.sim <- kable(Table.sim)
  saveRDS(Table.sim, "SSDfit.rds")
  return(list(SSD.table=as.data.frame(Table.sim), fit.parameters = par.fit, gof.test = par.gof, 
              sse = gof.measurs, df.hc=as.data.frame(hcTable), df.fa = as.data.frame(fa.Table)))
}
