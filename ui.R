# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(Rcpp)
library(RcppArmadillo)
library(shiny)
library(shinythemes)
library(shinydashboard)
library(visNetwork)

###############################################################################################
header<-dashboardHeader(title="Clustering Networks through ERGMs",titleWidth = 400)

sidebar<-dashboardSidebar(
  width = 300,
  sidebarMenu(id = "main",
              menuItem(text = "Clustering", 
                       tabName = "clust_net", 
                       icon = icon("dashboard"),
                       selected = TRUE,
                       startExpanded = TRUE
              ),
              menuSubItem(icon = NULL,tabName = "sub_start_ID",selectInput("nclust", "Number of clusters:",
                                                                           c("2" = "2",
                                                                             "3" = "3",
                                                                             "4" = "4")))
  )
)

###############################################################################################
demo_box<-fluidRow(
  column(width= 12, align="center",
         radioButtons("demo_indicator", "Do you want to see demo?",
                      c("Yes" = "yes",
                        "No" = "no"),selected = "yes"),
         conditionalPanel(condition = "(input.demo_indicator == 'no')",
                          fileInput("adjacency", "Choose CSV File with adjacency matrix")
         ))
)

update_button<-actionButton(inputId = 'go',label = h4("Analyze and Plot"),width = "175px",style="color: #fff; background-color: #337ab7; border-color: #2e6da4; margin: 0 auto")

update_button_row<-fluidRow(
  column(width= 12,icon("refresh"),update_button,align="center",style='padding:10px;')
)

tab_clust_net<-tabItem(tabName = "clust_net", h2("Clustering Networks through ERGMs"), demo_box, visNetworkOutput("NetPlot"),update_button_row)

###############################################################################################
body<- dashboardBody(tab_clust_net,
                     tags$head(tags$style(HTML('
                                               /* logo (dark grey) */
                                               .skin-blue .main-header .logo {
                                               background-color: #236EAF;
                                               font-family: "Arial";
                                               color: #FFFFFF;
                                               }
                                               
                                               /* logo when hovered */
                                               .skin-blue .main-header .logo:hover {
                                               background-color: #000000;
                                               }
                                               
                                               /* navbar (rest of the header) (dark grey) */
                                               .skin-blue .main-header .navbar {
                                               background-color: #236EAF;
                                               }
                                               
                                               /* main sidebar (dark grey)*/
                                               .skin-blue .main-sidebar {
                                               background-color: #EEEEEE;
                                               color: #EEEEEE;
                                               }
                                               
                                               /* active selected tab in the sidebarmenu */
                                               .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                                               background-color: #EEEEEE;
                                               color: #000000;
                                               font-weight: bold;
                                               }
                                               
                                               /* other links in the sidebarmenu 
                                               .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                               background-color: #EEEEEE;
                                               color: #EEEEEE;
                                               font-weight: bold;
                                               }*/
                                               
                                               /* other links in the sidebarmenu (background is light grey) */
                                               .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                               background-color: #EEEEEE;
                                               font-family: "Arial";
                                               color: #000000;
                                               }
                                               
                                               /* other links in the sidebarmenu when hovered */
                                               .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                                               background-color: #000000;
                                               color: #FFFFFF;
                                               font-weight: bold;
                                               }
                                               
                                               body, label, input, button, select {
                                               font-family: "Arial";
                                               }
                                               
                                               .content-wrapper,.right-side {
                                               background-color: #ffffff;
                                               }
                                               
                                               /* toggle button  */
                                               .skin-blue .main-header .navbar .sidebar-toggle{
                                               color: #FFFFFF;
                                               }
                                               
                                               /* toggle button when hovered  */
                                               .skin-blue .main-header .navbar .sidebar-toggle:hover{
                                               background-color: #000000;
                                               }
                                               
                                               .nav-tabs {
                                               background-color: #EEEEEE;
                                               border-color:#000000;
                                               }
                                               
                                               .nav-tabs-custom>.nav-tabs>li>a:hover {
                                               background-color: #000000;
                                               color: #FFFFFF;
                                               border-radius: 0;
                                               }
                                               
                                               .nav-tabs-custom>.tab-content {
                                               background: #FFFFFF;
                                               padding: 3px;
                                               border-bottom-right-radius: 1px;
                                               border-bottom-left-radius: 1px;
                                               border-color: #000000;
                                               }
                                               
                                               .nav-tabs-custom .nav-tabs li.active a {
                                               background-color: #F7800A;
                                               border-color: #236EAF;
                                               }
                                               
                                               .nav-tabs-custom .nav-tabs li.active:hover a{
                                               background-color: #000000;
                                               color: #FFFFFF;
                                               border-color: #236EAF;
                                               }
                                               
                                               .nav-tabs-custom .nav-tabs li.active {
                                               border-top-color: #236EAF;
                                               }
                                               
                                               .leaflet .legend i{
                                               border-radius: 50%;
                                               width: 10px;
                                               height: 10px;
                                               margin-top: 4px;
                                               }
                                               .test_class{color:#000000; text-align: center; font-weight: bold;}
                                               
                                               .box.box-solid.box-primary>.box-header {
                                               color:#000000;
                                               background:#EEEEEE;
                                               font-weight: bold;
                                               }
                                               
                                               .box.box-solid.box-primary>.box-header:hover {
                                               color:#FFFFFF;
                                               background:#000000;
                                               font-weight: bold;
                                               }
                                               
                                               .box.box-solid.box-primary{
                                               border-bottom-color:#666666;
                                               border-left-color:#666666;
                                               border-right-color:#666666;
                                               border-top-color:#666666;
                                               background-color:#FFFFFF;
                                               }
                                               
                                               .box.box-solid.box-primary>.box-header .btn, .box.box-solid.box-primary                               >.box-header a {
                                               color: #000000;
                                               }
                                               .box.box-solid.box-primary>.box-header .btn:hover {
                                               color: #FFFFFF;
                                               }
                                               
                                               .box.box-solid.box-primary>.box-header:hover a {
                                               color: #FFFFFF;
                                               }
                                               
                                               .skin-blue .wrapper {
                                               background-color: #EEEEEE;
                                               }
                                               
                                               ')))
                     )

shinyUI(dashboardPage(header,sidebar,body))
