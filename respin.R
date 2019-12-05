library(shinycssloaders)
library(shinyjqui)

# ui.R Shared resources and functions ----


#' Creates a resizable element with a css spinner
#'
#' @param x Shiny Item to resize and spin
#' @param type Int from 1-9 from https://projects.lukehaas.me/css-loaders/
#' @param color hex color for foreground
#' @param color.background Applicable only to types 2-3, set to be app background color
#'
#' @return resizable element with attached spinner
#' @export
#'
#' @examples
#'
respin <- function(x, type = 1, color = "#5f7800", color.background = '#fcfcfc'){
  jqui_resizable(x) %>%
    withSpinner(type = type,
                color = color,
                color.background = color.background)
}

#' Creates a spinning busy box when Shiny is refreshing
#'
#' @return ui element within Shiny
#' @export
#'
#' @examples
#' SpinBox()
#'
SpinBox <- function(){
  tagList(
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     absolutePanel(top="25%", left ="25%",
                                   wellPanel(style="background-color:green; color:white",
                                             draggable = TRUE,icon("refresh fa-spin fa-2x", lib="font-awesome"),
                                             "Loading..."
                                   )
                     )
    )
  )
}


