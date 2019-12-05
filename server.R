###################################################################
# SSD Tool ########################################################
###################################################################
library(shiny)
library(shinycssloaders)
library(shinyjqui)
library(shinyBS)
library(fitdistrplus)
library(fBasics)
library(EnvStats)
library(docxtools)
library(magrittr)
library(DT)
gc()
source("genFunctions.R")
source("Curve-fitting-func.R")

shinyServer(function(input, output) {

  SSDdata <- reactive({
    req(input$file1)
    inFile <- input$file1
    Data <- read_excel(inFile$datapath)
    colnames(Data)<- c("X to include","Taxa Grouping","Species", "Concentration")
    Data <- Data[order(Data$Concentration),]
    Data$fraq <- ppoints(Data[[4]], 0.5)
    Data$fraq <- as.numeric(format(round(Data$fraq, 2), nsmall = 2))
    print(head(Data))
    return(as.data.frame(Data))
  })

  filteredData3 <- reactive({
    req(input$file2)
    inFile <- input$file2
    Data <- read.csv(inFile$datapath)
    return(as.data.frame(Data))
  })

  output$downTemplate <- downloadHandler(
    filename = function(){
      paste("Toxicity_Template-Download","xlsx",sep=".")
    },
    content = function(con){
      file.copy("SSDdata.xlsx", con)
    })

  output$downTemplate2 <- downloadHandler(
    filename = function(){
      paste("Exposure_Template-Download","csv",sep=".")
    },
    content = function(con){
      file.copy("exposuredata.csv", con)
    })

  output$SSDtable <- DT::renderDataTable(
    SSDdata())

  # Generate a summary of the dataset
  output$summary.tox <- renderPrint({
    toxdata <- cbind(SSDdata()[4], SSDdata()[5])
    summary(toxdata)
  })

  # SSD Options
  output$unitselssd <- renderUI({
    selectizeInput("unitssd", "Units", c("", "ng/l", "ug/l", "mg/l", "lbs/A"), multiple = F)
  })

  output$plotsel <- renderUI({
    selectizeInput("plotpos", "Plotting Position ", c("", "Hazen", "Filiben's"), multiple = F)
  })

  # Unit Conversion
  unitfactor <- reactive({
    req(input$unitssd)
    as.character(input$unitssd)
    unifact <- 1
    if(input$unitssd %in% c("mg/l")) {
      unifact <- 0.001}
    else{
      if(input$unitssd %in% c("ng/l")) {
        unifact <- 1000}
    }
    print(unifact)
    return(unifact)
  })

  output$samplsizesel <- renderUI({
    selectizeInput("samplsize", "Bootstrap Iterations", c("", 500, 1000, 5000, 10000), multiple = F)
  })

  # SSD Calculation
  fit.fun <- eventReactive(input$dist, {
    fitfun <- Curve.fitting(SSDdata(), input$unitssd, input$conf, input$scale, as.numeric(input$exposure), as.numeric(input$samplsize))
    fitfun
  })

  SSDfit <- reactive({
    req(fit.fun())
    aa <- fit.fun()$SSD.table
    print(aa[1:6])
    return(aa)
  })

  output$FitTable <- DT::renderDataTable({
    SSDfit() %>% format(digits=2)},
    options = list(
      dom = 'T<"clear">lfrtip',
      autoWidth = TRUE,
      columnDefs = list(list(width = '20%', targets = list(2,3,4))),
      deferRender=TRUE,
      scrollX=TRUE,scrollY=400,
      pageLength = 5,
      scrollCollapse=TRUE),
    caption="Species Sensitivity Distributions with upper and lower confidence interval")

  fit.pars <- eventReactive(input$act2, {
    fit.fun()$fit.parameters
  })

  # Models parameters
  output$fitpars.table <- DT::renderDataTable(
    fit.pars(), rownames = TRUE)

  # GoF tests
  gof.pars <- eventReactive(input$act3, {
    fit.fun()$gof.test
  })

  output$gofPars <- DT::renderDataTable(
    datatable(gof.pars() , rownames = TRUE) %>%
      formatRound(columns=c("Anderson-Darling","Kolmogrov-Smirnov","Chi-Squared", "Chisq p-value"), digits=3)
  )

  # SSE and MSE
  dev.pars <- eventReactive(input$act4, {
    fit.fun()$sse
  })

  output$sse.mse <- DT::renderDataTable(
    dev.pars(),
    rownames = TRUE)

  # HC% and HC50
  hc2 <- eventReactive(input$act5, {
    fit.fun()$df.hc
  })

  output$hcTabl <- DT::renderDataTable(
    hc2() %>% format(digits=3))

  # Fraction Affected
  fa <- eventReactive(input$dist, {
    fit.fun()$df.fa
  })

  #### Plot SSD
  output$modelsel <- renderUI({
    selectizeInput("model", "Select Distribution Name",
                   c("", "Normal", "Logistic", "Extreme.Value", "Gumbel", "Weibull"), multiple = F)
  })

  # "hc" for plot
  hc <- reactive({
    fit.fun()$df.hc
  })

  # Data for SSD plots
  distdata <- reactive({
    colname1 <- input$model
    colname3 <- paste0(input$model,".lwr")
    colname4 <- paste0(input$model,".upr")
    data <- cbind(SSDfit()[,colname1], SSDfit()$pro, SSDfit()[,colname3], SSDfit()[,colname4])
    colnames(data) <- c("dat", "pro", "lwr", "upr")
    return(as.data.frame(data))
  })

  output$taxasel <- renderUI({
    taxa.group <- c(sort(unique(as.character(SSDdata()$"Taxa Grouping"))))
    pickerInput("taxa", "Add label to plot", choices=taxa.group, options = list(`actions-box` = TRUE),multiple = T)
  })

  ## Plots #########################################################################################################
  output$plot_USGS <- renderPlot({
    req(filteredData3())
    ggplot(filteredData3(), aes(ResultMeasureValue)) +
      stat_ecdf(geom = "point", pad = FALSE, col="dodgerblue") +
      scale_x_log10(name = paste0("Concentration (ug/l)")) +
      scale_y_continuous(name = "Likelihood of exceeding exposure treshold")+
      ggtitle(paste0("Distribution of prometryn detects in surface water samples from the USGS database")) +
      theme(plot.title = element_text(size=14, hjust=0), legend.position="none", axis.text=element_text(size=12),
            plot.background = element_rect(fill = "lavender"), axis.title=element_text(size=12,face="bold"))
  }, height = 400, width = 600)

  plotl <- reactive({
    req(input$model, unitfactor())
    #savedhc <<- hc()
    n <- 6 # count of dist. models + 1
    xint <- as.numeric(c(hc()[[input$model]][[1]], hc()[[input$model]][[4]]))
    x.lim <- c(0.01*min(SSDfit()[2:n]), max(SSDfit()[2:n]))
    print(xint, x.lim)

    ggplot() +
      stat_ecdf(data=filteredData3(), aes(ResultMeasureValue*unitfactor()), geom = "point", pad = FALSE) +
      stat_ecdf(data = SSDdata(), aes(Concentration) ,geom = "point", pad = FALSE, col="dodgerblue") +
      geom_line(data=distdata(),  aes(x= dat, y=pro), colour="magenta3") +
      geom_line(data=distdata(), aes(x= lwr, y=pro), colour="magenta3", lty="dotted") +
      geom_line(data=distdata(), aes(x= upr, y=pro), colour="magenta3", lty="dotted") +
      geom_vline(aes(xintercept = xint), linetype = "dashed", colour = "lightslateblue") +
      geom_vline(aes(xintercept = as.numeric(input$exposure)), colour = "lightslateblue") +
      geom_text(aes(x=xint[1], y= 0.75*max(SSDdata()$fraq), label=paste0("HC5=", xint[1])), colour="black", angle=90, text=element_text(size=11)) +
      geom_text(aes(x=xint[2], y= 0.5*max(SSDdata()$fraq), label=paste0("HC50=", xint[2])), colour="black", angle=90, text=element_text(size=11)) +
      geom_text(aes(x=as.numeric(input$exposure), y= 0.25*max(SSDdata()$fraq), label=paste0("FA=", as.numeric(fa()[[input$model]])*100,"%")),
                colour="black", angle=90, text=element_text(size=11)) +
      geom_text(data = SSDdata()[SSDdata()$"Taxa Grouping" %in% input$taxa,],
                aes(x = Concentration, y = fraq, label = Species), hjust = 1, size = 4) +
      scale_x_log10(limits=x.lim) +
      #ggtitle(paste0("Water quality concentration (", input$unitssd,"), Toxicity concentration (", input$unitssd, ")")) +
      theme(plot.title = element_text(size=11, hjust=0), legend.position="none", axis.text=element_text(size=10),
            axis.title=element_text(size=10))+
      labs(x = paste0("Concentration"), y ="Fraction of species affected")
  })

  plota <- reactive({
    req(input$model)
    xint <- as.numeric(c(hc()[[input$model]][[1]], hc()[[input$model]][[4]]))
    print(xint)
    ggplot() +
      stat_ecdf(data=filteredData3(), aes(ResultMeasureValue*unitfactor()), geom = "point", pad = FALSE) +
      stat_ecdf(data = SSDdata(), aes(Concentration) ,geom = "point", pad = FALSE, col="dodgerblue") +
      geom_line(data=distdata(),  aes(x= dat, y=pro), colour="magenta3") +
      geom_line(data=distdata(), aes(x= lwr, y=pro), colour="magenta3", lty="dotted") +
      geom_line(data=distdata(), aes(x= upr, y=pro), colour="magenta3", lty="dotted") +
      geom_vline(aes(xintercept = xint), linetype = "dashed", colour = "lightslateblue") +
      geom_vline(aes(xintercept = as.numeric(input$exposure)), colour = "lightslateblue") +
      geom_text(aes(x=xint[1], y= 0.75*max(SSDdata()$fraq), label=paste0("HC5=", xint[1])), colour="black", angle=90, text=element_text(size=11)) +
      geom_text(aes(x=xint[2], y= 0.5*max(SSDdata()$fraq), label=paste0("HC50=", xint[2])), colour="black", angle=90, text=element_text(size=11)) +
      geom_text(aes(x=as.numeric(input$exposure), y= 0.25*max(SSDdata()$fraq), label=paste0("FA=", as.numeric(fa()[[input$model]])*100,"%")),
                colour="black", angle=90, text=element_text(size=11)) +
      scale_x_discrete(limits=c(-1*max(SSDdata()$Concentration), max(SSDdata()$Concentration))) +
      geom_text(data = SSDdata()[SSDdata()$"Taxa Grouping" %in% input$taxa,],
                aes(x = Concentration, y = fraq, label = Species), hjust = 1, size = 4) +
      #ggtitle(paste0("Water quality concentration (", input$unitssd,"), Toxicity concentration (", input$unitssd, ")")) +
      theme(plot.title = element_text(size=11, hjust=0), legend.position="none", axis.text=element_text(size=10),
            axis.title=element_text(size=10))+
      labs(x = paste0("Concentration"), y ="Fraction of species affected")
  })

  output$plot <- renderPlot({
    if(input$scale == 2){
      plotl() + ggtitle(paste0("Water quality concentration (", input$unitssd,"), Toxicity concentration (", input$unitssd, ")"))
    }else{
      plota() + ggtitle(paste0("Water quality concentration (", input$unitssd,"), Toxicity concentration (", input$unitssd, ")"))
    }
  })
  ###########################################################################################################
  #### Generating R Markdown report #########################################################################
  ###########################################################################################################
  wrapper <- function(x, ...)
  {
    paste(strwrap(x, ...), collapse = "\n")
  }

  plot_USGS2 <- reactive({
    req(filteredData3())
    p <- ggplot(filteredData3(), aes(ResultMeasureValue)) +
      stat_ecdf(geom = "point", pad = FALSE, col="dodgerblue") +
      scale_x_log10(name = paste0("Concentration (", input$unitssd,")")) +
      scale_y_continuous(name = "Likelihood of exceeding exposure treshold")+
      #labs(caption = "Source: the Lahman baseball database")+
      ggtitle(wrapper("Figure 6. Distribution of prometryn detects in surface water samples from the
                      \nUSGS database", width = 100)) +
      theme(plot.title = element_text(size=11, hjust=0), legend.position="none", axis.text=element_text(size=10),
            axis.title=element_text(size=10))
    return(p)
  })

  ssd.plotl <- reactive({
    x.lim <- c(min(SSDfit()[2:6]), max(SSDfit()[2:6]))
    p2 <- ggplot() +
      stat_ecdf(data = SSDdata(), aes(Concentration) ,geom = "point", pad = FALSE, col="dodgerblue") +
      geom_text(data = as.data.frame(SSDdata()),
                aes(x = Concentration, y = fraq, label = Species), hjust = 1, size = 2) +
      geom_line(data=SSDfit(),  aes(x= Normal, y=pro, colour="Normal")) +
      geom_line(data=SSDfit(),  aes(x= Logistic, y=pro, colour="Logistic")) +
      geom_line(data=SSDfit(),  aes(x= Extreme.Value, y=pro, colour="Extreme.Value")) +
      geom_line(data=SSDfit(),  aes(x= Gumbel, y=pro, colour="Gumbel")) +
      geom_line(data=SSDfit(),  aes(x= Weibull, y=pro, colour="Weibull")  ) +
      scale_x_log10(limits=x.lim) +
      scale_color_manual(name = NULL, values = c("Normal"="red", "Logistic"="yellow", "Extreme.Value"="green",
                                                 "Gumbel"="magenta3", "Weibull"="black")) +
      theme(plot.title = element_text(size=11, hjust=0), axis.text=element_text(size=10),
            axis.title=element_text(size=10), legend.position="right") +
      labs(x = paste0("Concentration (", input$unitssd,")"), y ="Fraction of species affected")
    return(p2)
  })

  ssd.plota <- reactive({
    x.lim <- c(min(SSDfit()[2:6]), max(SSDfit()[2:6]))
    p2 <- ggplot() +
      stat_ecdf(data = SSDdata(), aes(Concentration) ,geom = "point", pad = FALSE, col="dodgerblue") +
      geom_text(data = as.data.frame(SSDdata()),
                aes(x = Concentration, y = fraq, label = Species), hjust = 1, size = 2) +
      geom_line(data=SSDfit(),  aes(x= Normal, y=pro, colour="Normal")) +
      geom_line(data=SSDfit(),  aes(x= Logistic, y=pro, colour="Logistic")) +
      geom_line(data=SSDfit(),  aes(x= Extreme.Value, y=pro, colour="Extreme.Value")) +
      geom_line(data=SSDfit(),  aes(x= Gumbel, y=pro, colour="Gumbel")) +
      geom_line(data=SSDfit(),  aes(x= Weibull, y=pro, colour="Weibull")) +
      scale_x_discrete(limits=c(-1*max(SSDdata()$Concentration), max(SSDdata()$Concentration))) +
      scale_color_manual(name = NULL, values = c("Normal"="red", "Logistic"="yellow", "Extreme.Value"="green",
                                                 "Gumbel"="magenta3", "Weibull"="black"))+
      theme(plot.title = element_text(size=11, hjust=0), axis.text=element_text(size=10),
            axis.title=element_text(size=10), legend.position="right")+
      labs(x = paste0("Concentration (", input$unitssd,")"), y ="Fraction of species affected")
    return(p2)
  })

  ssd.plot1 <- reactive({
    if(input$scale == 2){
      ssd.plotl() + ggtitle(wrapper("Figure 7. Model fits visualization of the aquatic plant species sensitivity
                                    distributions \ndeveloped for prometryn. Dots on the graph represent endpoints
                                    for individual species, while the lines represent different model fits of the
                                    data distributions.", width = 100))
    }else{
      ssd.plota() + ggtitle(wrapper("Figure 7. Model fits visualization of the aquatic plant species sensitivity
                                    distributions \ndeveloped for prometryn. Dots on the graph represent endpoints
                                    for individual species, while the lines represent different model fits of the
                                    data distributions.", width = 100))
    }
    })

  ssd.plot2 <- reactive({
    if(input$scale == 2){
      plotl() + ggtitle(wrapper("Figure 8. Distribution of prometryn detects in surface water samples
                                from the \nUSGS database.", width = 100))
    }else{
      plota() + ggtitle(wrapper("Figure 8. Distribution of prometryn detects in surface water samples
                                from the \nUSGS database.", width = 100))
    }
  })

  mrkdwn.gof <- reactive({
    t(fit.fun()$gof.test)
  })

  ssd.summary <- reactive({
    rbind(fit.fun()$sse, fit.fun()$df.hc)
  })

  output$report <- downloadHandler(
    filename = function() {
      paste('report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'))
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "report2.Rmd")
      file.copy("report2.Rmd", tempReport, overwrite = TRUE)

      out <- rmarkdown::render(tempReport, switch(
        input$format,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    })

})
