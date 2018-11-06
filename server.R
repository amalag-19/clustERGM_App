
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
require(visNetwork)

shinyServer(function(input, output) {
  #Requiring packages to plot
  source("shinyfunc.R")
  require(grid)
  # Generating Plots
  output$NetPlot <- renderVisNetwork({
    input$go
    ## Plotting scree plot for spatial clustering by K-means
    isolate({
      if(input$demo_indicator=="yes"){
        ## Loading International trade network 1990
        net_result_trade<-wrapper_ERGM_stat_undir_Dens(sim.net=trade50[,,10],nclust = as.numeric(input$nclust),thres=10^(-6),theta_init=jitter(rep(0,as.numeric(input$nclust))),sim_indicator=0)
        plot_out<-net_plotter_visNetwork(input_demo = "yes",adjacency_mat = trade50[,,10],net_result = net_result_trade,nclust = as.numeric(input$nclust))
        #print("flag")
      }else if(input$demo_indicator=="no"){
        inFile <- input$adjacency
        if (is.null(inFile))
          return(NULL)
        df_adjacency<-read.csv(inFile$datapath)
        labels<-names(df_adjacency)
        adjacency_mat<-unname(data.matrix(df_adjacency))
        net_result<-wrapper_ERGM_stat_undir_Dens(sim.net=adjacency_mat,nclust = as.numeric(input$nclust),thres=10^(-6),theta_init=jitter(rep(0,as.numeric(input$nclust))),sim_indicator=0)
        plot_out<-net_plotter_visNetwork(input_demo = "no",adjacency_mat = adjacency_mat,net_result = net_result,nclust = as.numeric(input$nclust),node_labels=labels)
      }
      print(plot_out)
    })
  })
})

