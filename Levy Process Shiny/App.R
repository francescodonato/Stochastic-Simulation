library(shiny)
library(shinythemes)
library(shinyWidgets)
library(tidyverse)
library(ggrepel)
library(reshape2)
library(gghighlight)
library(ggforce)
library(latex2exp)

#######################################################################################
###########                       Auxiliary  functions                      ###########
#######################################################################################
hyperexp<-function(n,lambda1,lambda2,p1){
  times<-rep(0,n)
  exp_matr<-matrix(c(rexp(n,lambda1),rexp(n,lambda2)),ncol = 2)
  temp<-replicate(n,sample(1:2,size = 1,prob = c(p1,1-p1)))
  for(i in 1:n){
    times[i]<-exp_matr[i,temp[i]]
  }
  return(times)
}



## common number pareto
pareto<-function(n,jump_size,k){
  if(k<=1){stop("Error Inputs, k must be > 1")}
  # beta computed to give a mean = jump_size
  beta<-((k-1)*jump_size)/k
  times<-beta*(1-runif(n))^(-1/k)
  return(times)
}


theta_v<-function(mu,sigma,lambda){
  d= -(mu/sigma^2)+sqrt((mu^2/sigma^4)+(2*lambda/sigma^2))
  return(d)
}

theta_w<-function(mu,sigma,lambda){
  d=(mu/sigma^2)+sqrt((mu^2/sigma^4)+(2*lambda/sigma^2))
  return(d)
}

#######################################################################################
###########                          Main function                          ###########
#######################################################################################

simulation<-function(n,mu,sigma,lambda, Y="Exponential",return_val="final",
                     jump_size=1,#check if jump size is lambda poisson
                     rate1=NULL,rate2=NULL,p=NULL,sd=NULL,k=NULL){
  
  ############# init
  unif<-runif(n)
  A<-rep(0,n+1)
  P<-rep(0,n+1)
  M<-rep(0,n+1)
  rate_v<-theta_v(mu,sigma,lambda)
  rate_w<-theta_w(mu,sigma,lambda)
  
  ############# choiche F common numbers
  if(Y=="Exponential"){
    Y_vect<-rexp(n,1/jump_size)
  }
  else if(Y=="Hyperexponential"){#still work to do here
    Y_vect<-hyperexp(n,rate1,rate2,p)
  }
  else if(Y=="Normal"){
    Y_vect<-rnorm(n,sd)
  }
  else if(Y=="Pareto"){
    Y_vect<-pareto(n,jump_size,k)
  }
  
  ############# Variables
  # Y_vect<-rexp(n,jump_size)#built-in
  V<-rexp(n,rate_v)#built-in
  W<-rexp(n,rate_w)#built-in
  
  ## common number
  # V<-simple_exp(unif,rate_v)
  # W<-simple_exp(unif,rate_w)
  
  ############# simulation
  
  for(i in 2:(n+1)){
    P[i]<-A[i-1]+(V[i-1]-W[i-1])
    A[i]<-P[i]+Y_vect[i-1]
    M[i]<-max(M[i-1], P[i]-W[i-1], A[i])
  }
  ############# return
  if(return_val=="final"){
    return(c(P[n+1],A[n+1],M[n+1]))
  }
  else if(return_val=="path"){
    # return(list("P"=P,"A"=A,"M"=M,"Time"=c(0,cumsum(Y_vect))))
    return(data.frame("P"=P,"A"=A,"M"=M,"Iteration"=seq(0,n)))
  }
}

############################################

ui <- navbarPage("Lèvy Process", theme = shinytheme("cyborg"),
                 tabPanel(
                   "Description",
                   withMathJax(includeMarkdown("descr.md"))
                 ),
                 tabPanel("Simple process",fluid = TRUE,
                          sidebarLayout(
                            sidebarPanel(
                              sliderInput("n", label = "Number of interaction", min = 1000, 
                                          max = 10000, value = 1000,step = 100),
                              
                              helpText("Experiment with different distributions and parameters"),
                              
                              fluidRow(
                                column(6,
                                       numericInput("mu",
                                                    label = HTML("&mu;"),
                                                    value = 1)),
                                column(6,
                                       numericInput("sigma",
                                                    label = HTML("&sigma;"),
                                                    value = 1,
                                                    min = 0.001)),
                              ),
                              fluidRow(
                                column(6,
                                       numericInput("lambda_poi",
                                                    label = HTML("&lambda;<sub>poisson</sub>"),
                                                    value = 1,
                                                    min = 0.001)),
                                column(6,
                                       numericInput("jump_size",
                                                    label= "Mean of Y",
                                                    value = 1,
                                                    min = 0.001)),
                              ),
                              
                              selectInput("var",
                                          label = paste("Choose a distribution for Y", "(the mean is fixed)", sep = "\n"),
                                          choices = list("Exponential",
                                                         "Hyperexponential",
                                                         "Pareto",
                                                         "Normal"),
                                          selected = "Exponential"
                              ),
                              
                              
                              conditionalPanel(
                                condition = "input.var == 'Hyperexponential'",
                                fluidRow(
                                  column(4,
                                         numericInput("lambda1",
                                                      label = HTML("&lambda;<sub>1</sub>"),
                                                      value = 1,
                                                      min = 0.001)
                                  ),
                                  column(4,
                                         numericInput("lambda2",
                                                      label = HTML("&lambda;<sub>2</sub>"),
                                                      value = 1,
                                                      min = 0.001)
                                  ),
                                  column(4,
                                         numericInput("p1",
                                                      label = HTML("p<sub>1</sub>"),
                                                      value = 0.5,
                                                      min = 0,max = 1)
                                  )
                                )
                              ),
                              
                              conditionalPanel(
                                condition = "input.var == 'Pareto'",
                                fluidRow(
                                  column(6,
                                         numericInput("k",
                                                      label = HTML("k"),
                                                      value = 1, min = 1)
                                  )
                                )
                              ),
                              
                              conditionalPanel(
                                condition = "input.var == 'Normal'",
                                fluidRow(
                                  column(6,
                                         numericInput("sd",
                                                      label= "Standard dev.",
                                                      value = 1)
                                  )
                                )
                              ),
                              
                              actionButton("action", "Simulate"),
                              
                            ),

                            mainPanel(
                              
                              fluidRow(column(12,
                                              plotOutput("plot")
                              )
                              ),
                              fluidRow(
                                column(4,
                                       radioButtons("P_A", label = "Variable to display",
                                                    choices = list("P" = 1, "A" = 2),
                                                    selected = 1),
                                )
                              )
                            )
                              )
                           ),
                 
                 tabPanel("Multiple simulations",fluid = TRUE,
                          sidebarLayout(
                            sidebarPanel(
                              sliderInput("n_m", label = "Number of interaction", min = 1000, 
                                          max = 10000, value = 1000,step = 100),
                              
                              helpText("Experiment with different distributions and parameters"),
                              
                              fluidRow(
                                column(6,
                                       numericInput("mu_m",
                                                    label = HTML("&mu;"),
                                                    value = 1)),
                                column(6,
                                       numericInput("sigma_m",
                                                    label = HTML("&sigma;"),
                                                    value = 1,
                                                    min = 0.001)),
                              ),
                              fluidRow(
                                column(6,
                                       numericInput("lambda_poi_m",
                                                    label = HTML("&lambda;<sub>poisson</sub>"),
                                                    value = 1,
                                                    min = 0.001)),
                                column(6,
                                       numericInput("jump_size_m",
                                                    label= "Mean of Y",
                                                    value = 1,
                                                    min = 0.001)),
                              ),
                              
                              selectInput("var_m",
                                          label = paste("Choose a distribution for Y", "(the mean is fixed)", sep = "\n"),
                                          choices = list("Exponential",
                                                         "Hyperexponential",
                                                         "Pareto",
                                                         "Normal"),
                                          selected = "Exponential"
                              ),
                              
                              
                              conditionalPanel(
                                condition = "input.var_m == 'Hyperexponential'",
                                fluidRow(
                                  column(4,
                                         numericInput("lambda1_m",
                                                      label = HTML("&lambda;<sub>1</sub>"),
                                                      value = 1,
                                                      min = 0.001)
                                  ),
                                  column(4,
                                         numericInput("lambda2_m",
                                                      label = HTML("&lambda;<sub>2</sub>"),
                                                      value = 1,
                                                      min = 0.001)
                                  ),
                                  column(4,
                                         numericInput("p1_m",
                                                      label = HTML("p<sub>1</sub>"),
                                                      value = 0.5,
                                                      min = 0,max = 1)
                                  )
                                )
                              ),
                              
                              conditionalPanel(
                                condition = "input.var_m == 'Pareto'",
                                fluidRow(
                                  column(6,
                                         numericInput("k_m",
                                                      label = HTML("k"),
                                                      value = 1, min = 1)
                                  )
                                )
                              ),
                              
                              conditionalPanel(
                                condition = "input.var_m == 'Normal'",
                                fluidRow(
                                  column(6,
                                         numericInput("sd_m",
                                                      label= "Standard dev.",
                                                      value = 1)
                                  )
                                )
                              ),
                              
                              actionButton("action_m", "Simulate"),
                              
                            ),
                            
                            mainPanel(
                              
                              fluidRow(column(12,
                                              plotOutput("plot_m")
                              )
                              ),

                                fluidRow(column(4,
                                                sliderInput("Histogram_bins",
                                                            label = "Number of bins:",
                                                            min = 0, max = 50, value =20)
                                )
                              )
                            )

                            )
                          )
                 )



# Define server logic ----
server <- function(input, output) {
# first tab
    v <- reactiveValues(data = NULL)
  
  observeEvent(input$action, {
    v$data <- simulation(input$n,input$mu,input$sigma,input$lambda_poi, jump_size=input$jump_size,
                         Y = input$var_m,return_val = "path",
                         rate1 = input$lambda1_m, rate2 = input$lambda2_m, p = input$p1_m,
                         sd = input$sd_m, k=input$k_m)
  })
  
  output$plot <- renderPlot({
    if (is.null(v$data)) return()
    if(input$P_A==1){
      v$data%>%
        ggplot(aes(y=P,x=Iteration))+ #can choose P or A
        geom_step()+  #or use geom_path
        labs(title = "Evolution of the Lèvy process",
             subtitle = "Timestamp at each jump",
             y="Values of P")
    }
    else if(input$P_A==2){
      v$data%>%
        ggplot(aes(y=A,x=Iteration))+ #can choose P or A
        geom_step()+  #or use geom_path
        labs(title = "Evolution of the Lèvy process",
             subtitle = "Timestamp at each jump",
             y="Values of A")
      
    }
    
  })
  
  # second tab
  v_m <- reactiveValues(data_m = NULL)
  
  observeEvent(input$action_m, {
    v_m$data_m <- as.data.frame(t(replicate(100,simulation(input$n_m,input$mu_m,input$sigma_m,input$lambda_poi_m, jump_size=input$jump_size_m,
                         Y = input$var,rate1 = input$lambda1, rate2 = input$lambda2, p = input$p1,
                         sd = input$sd, k=input$k))))
  })
  
  output$plot_m <- renderPlot({
    if (is.null(v_m$data_m)) return()
    v_m$data_m%>%
      rename("P"="V1","A"="V2","M"="V3")%>%#renaming
      melt()%>%
      ggplot()+#plotting
      geom_histogram(aes(x=value,group=variable),bins = input$Histogram_bins,fill="white",color="black")+
      labs(title = "Histogram for Exponential process",
           x="Values", y="Count")+
      facet_wrap(vars(variable),scales = "free")
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)