library(taxogeno)
library(shiny)
library(shinyTree)
library(ggplot2)
## https://christophergandrud.github.io/networkD3/
library(networkD3)
library(ggiraph)
library(stringr)
library(parallel)
library("DBI")
library("jsonlite")
library(taxogeno)

runApp(shinyApp(taxogenoWebAppUi, taxogenoWebAppServer), port=10081)
