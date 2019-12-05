library(shiny)
library(ggplot2)
library(shinyBS)
library(DT)
library(shinyWidgets)

source("respin.R", local = TRUE)

# Define UI for application that draws a histogram
shinyUI(navbarPage("Species Sensitivity Distributions Tool",
                   theme = "bootstrap.css",
                   inverse = TRUE,
          tabsetPanel(
                   tabPanel("SSD Tool", icon = icon("bar-chart-o"),
                   sidebarLayout(
                     sidebarPanel(width=3,
                                  SpinBox(),
                                  h4("1. Please upload input data files:"),
                                  br(),
                    fileInput("file1", label="Toxicity Data"),
                    fileInput("file2", label="Exposure Data"),
                    hr(),
                    h4("Example input data files:"),
                    br(),
                    downloadButton("downTemplate","Download Toxicity Data Template"),
                    br(),
                    br(),
                    downloadButton("downTemplate2","Download Exposure Data Template")
                    ),
         mainPanel(

           fluidRow(
             SpinBox(),
             h4("2. Specify Options for SSD"),
             br(),
             column(3,
                    fluidRow(
                      uiOutput("unitselssd"),
                      bsTooltip("unitselssd", "Units of measure for toxicity data", placement = "top",
                                trigger = "hover", options = NULL),
                      uiOutput("plotsel"),
                      bsTooltip("plotsel", "Only Hazen plotting position is active",
                                placement = "top", trigger = "hover", options = NULL),

                      radioButtons("scale", "Data Scale",
                                   list("Arthmetic"=1, "Logarithmic"=2), selected = 1),
                      uiOutput("samplsizesel")
                    )),
             column(6,
                    textInput("chem", ("Chemical Name"), value = ""),
                    bsTooltip("chem", "Active ingredient name (optional)", placement = "top",
                              trigger = "hover", options = NULL),
                    textInput("exposure", "Exposure Concentration", value ="0.0"),
                    bsTooltip("exposure", "Is required to calculate Fraction Affected (FA), unit of EC value must be identical to the unit of the toxicity values. (optional)",
                              placement = "top", trigger = "hover", options = NULL),
                    numericInput("conf", "Desired Confidence %", value = 95)),
             bsTooltip("conf", "Confidence interval for median by bootstrap",
                       placement = "top", trigger = "hover", options = NULL)
           ),

           fluidRow(
             SpinBox(),
             hr(),
             h4("3. Species Sensitivity Distributions"),
             br(),
             actionButton("dist","Calculate"),
             bsTooltip("dist", "The program calculates various distributions through toxicity values for the species",
                       placement = "top", trigger = "hover", options = NULL),
             br(),
             br(),
             DT::dataTableOutput("FitTable"),
             br(),
             br(),
             hr(),
             h4("4. Outputs"),
             br(),
             actionButton("act2","Distribution Parameters"),
             bsTooltip("act2", "Table of distributions parameters", placement = "top", trigger = "hover",
                       options = NULL),
             br(),
             br(),
             DT::dataTableOutput("fitpars.table"),
             br(),
             actionButton("act3","Goodness of Fit Tests"),
             bsTooltip("act3", "Goodness-of-fit tests to decide whether your data comes from the specific population distribution", placement = "top", trigger = "hover",
                       options = NULL),
             br(),
             br(),
             DT::dataTableOutput("gofPars"),
             br(),
             actionButton("act4","SSE and MSE"),
             bsTooltip("act4", "Sum of Squared Errors and Mean Squared Error to determine the distribution with minimum deviation", placement = "top", trigger = "hover",
                       options = NULL),
             br(),
             br(),
             DT::dataTableOutput("sse.mse"),
             br(),
             actionButton("act5","HC5 and HC50"),
             bsTooltip("act5", "Hazardous concentration affecting 5% and 50% of the species", placement = "top", trigger = "hover",
                       options = NULL),
             br(),
             br(),
             DT::dataTableOutput("hcTabl")
           ),

           fluidRow(
             hr(),
             h4("4. Plot Species Sensitivity Distribution"),
             br(),
             uiOutput("modelsel"),
             column(9,plotOutput("plot")),
             column(3, uiOutput("taxasel"))
           )
  ))),

   tabPanel("Table and Summary", icon = icon("table"),
  fluidRow(
    h3("Summary of Toxicity Data"),
    column(6, verbatimTextOutput("summary.tox")),
    br(),
    DT::dataTableOutput("SSDtable")
  )),
  # ),

  tabPanel("Report", icon=icon("pencil",lib="glyphicon"),
           SpinBox(),
           radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                        inline = TRUE),
           downloadButton("report", "Generate report")
  ),
  tabPanel("Help", icon=icon("pencil",lib="glyphicon"),
           h4("Input exposure Data"),
           "Exposure data file must have column 'ResultMeasureValue' for concentrations.
            All other columns in the input data file are redundant.
            Example input data file includes 'Prometon' concentrations from WQP website."),
  hr(),
  tags$p("Contact: Nina Omani (ninaomani@gmail.com)")
))) # Closes tabset



