#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("How long will it take to recover from a pandemic?", windowTitle = "Epidemic Look"),
  
  sidebarLayout(
    
    # Sidebar with a slider input
    sidebarPanel(
      sliderInput("beta",
                  "Social Distancing Factor (Infection Rate):",
                  min = 0,
                  max = 5,
                  value = 2.16,
                  step = .01
                  ),
      
      sliderInput("gamma",
                   "Rate of Recovery (1/days):",
                   min = 0,
                   max = 1,
                   value = .5,
                   step = .01
                  ),
      sliderInput("birth",
                  "Population Constant:",
                  min = 0,
                  max = .1,
                  value = .01,
                  step = .000001),
      sliderInput("epsilon",
                  "Wanning Rate:",
                  min = 0,
                  max = .1,
                  value = .001,
                  step = .00001),
      
      numericInput("time",
                   "Timeline:",
                   min = 0,
                   max = 1000,
                   value = 250),
      
      selectInput(inputId = "model", "Choose a Model",
                  list("SIR", "SIRS")
                  ,width = 1000),
      
      checkboxInput("vital", "Vital Dynamics",
                    TRUE)

      
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("epidemic"),
      strong("Instructions on how to use are as follows..."),
      p("The social distancing factor slider ranges from 0 to 5 with 0 being complete isolation as in the infection rate is effectivly zero and with 5 being 
        an almost uncontrollable infection. If the social distancing factor is low that means we are flattening the curve. If the rate is high then we aren't 
        effectivley doing that. The rate of recovery is, much like it sounds, how long does it take for someone to recover and thus be added to the recovered pool
        of people. The timeline input is how far out you would like to look based on the other initial conditions."),
      p("The initial conditions when you open the app are that the infected multiply by 2.16 times a day and it takes 2 days to recover from said infection."),
      p("The population constant is the birth and death rate as mu and nu are both equal to each other in most assumptions."),
      p("The wanning variable indicates a wanning immunity since you won't carry a lifelong immunity to the disease e.g. the seasonal flu"),
      strong("Play with the sliders and see how different birth/death rates, recovery rates, social distancing rates, and also a immunity impact epidemics. The SIR
             model is mainly used for (one and done) type diseases and the SIRS model is used for the other variety.")
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$epidemic <- renderPlot({
    library(deSolve)
    library(ggplot2)
    t <- input$time
    b <- input$beta
    g <- input$gamma
    model <- input$model
    vital <- input$vital
    m <- input$birth
    th <- m
    e<- input$epsilon
    
    ##ode solution
    if (model=="SIR" & isFALSE(vital)){
      init <- c(S=1-1e-6,I=1e-6,R=0)
      parameters <- c(bet=b,gamm=g)
      time <- seq(0,t,by=t/(2*length(1:t)))
      eqn <- function(time,state,parameters){
        with(as.list(c(state,parameters)),{
          dS <- -bet*S*I
          dI <- bet*S*I-gamm*I
          dR <- gamm*I
          return(list(c(dS,dI,dR)))})}
      out<-ode(y=init,times=time,eqn,parms=parameters)
      out.df<-as.data.frame(out)

    } else if (model=="SIR" & !isFALSE(vital)){
      init <- c(S=1-1e-6,I=1e-6,R=0)
      parameters <- c(bet=b,gamm=g, mu=m, theta=th)
      time <- seq(0,t,by=t/(2*length(1:t)))
      eqn <- function(time,state,parameters){
        with(as.list(c(state,parameters)),{
          dS <- m-bet*S*I-th*S
          dI <- bet*S*I-gamm*I-th*I
          dR <- gamm*I-m*R
          return(list(c(dS,dI,dR)))})}
      out<-ode(y=init,times=time,eqn,parms=parameters)
      out.df<-as.data.frame(out)

    } else if (model=="SIRS" & isFALSE(vital)){
      init <- c(S=1-1e-6,I=1e-6,R=0)
      parameters <- c(bet=b,gamm=g, epsi=e)
      time <- seq(0,t,by=t/(2*length(1:t)))
      eqn <- function(time,state,parameters){
        with(as.list(c(state,parameters)),{
          dS <- -bet*S*I+epsi*R
          dI <- bet*S*I-gamm*I
          dR <- gamm*I-epsi*R
          return(list(c(dS,dI,dR)))})}
      out<-ode(y=init,times=time,eqn,parms=parameters)
      out.df<-as.data.frame(out)
      out.df$W <- 1-out.df$S-out.df$I-out.df$R

    } else if (model=="SIRS" & !isFALSE(vital)){
      init <- c(S=1-1e-6,I=1e-6,R=0)
      parameters <- c(bet=b,gamm=g,  mu=m, theta=th, epsi=e)
      time <- seq(0,t,by=t/(2*length(1:t)))
      eqn <- function(time,state,parameters){
        with(as.list(c(state,parameters)),{
          dS <- mu-bet*S*I+epsi*R-theta*S
          dI <- bet*S*I-gamm*I-theta*I
          dR <- gamm*I-epsi*R-mu*R
          return(list(c(dS,dI,dR)))})}
      out<-ode(y=init,times=time,eqn,parms=parameters)
      out.df<-as.data.frame(out)
      out.df$W <- 1-out.df$S-out.df$I-out.df$R
      } 
    
    ## take to long to run for a nice clean app
    # else if (model=="SEIR" & !isFALSE(vital)){
    #   init <- c(S = 0.1, E = 1e-04, I = 1e-04, R = 1 - 0.1 - 1e-4 - 1e-4)
    #   parameters <- c(bet=b, gamm=g, sig=s)
    #   time <- seq(0,t,by=t/(2*length(1:t)))
    #   eqn <- function(time,state,parameters){
    #     with(as.list(c(state,parameters)),{
    #       dS <- -bet*S*I
    #       dE <- bet*S*I-sig*E
    #       dI <- sig*E-gamm*I
    #       dR <- gamm*I
    #       return(list(c(dS,dE,dI,dR)))})}
    #   out<-ode(y=init,times=time,eqn,parms=parameters)
    #   out.df<-as.data.frame(out)
    # 
    # } else if (model=="SEIR" & isFALSE(vital)){
    #   init <- c(S = 0.1, E = 1e-04, I = 1e-04, R = 1 - 0.1 - 1e-4 - 1e-4)
    #   parameters <- c(bet=b,gamm=g, sig=s, mu=m, theta=th)
    #   time <- seq(0,t,by=t/(2*length(1:t)))
    #   eqn <- function(time,state,parameters){
    #     with(as.list(c(state,parameters)),{
    #       dS <- mu-theta*S-bet*S*I
    #       dE <- bet*S*I-sig*E-theta*E
    #       dI <- sig*E-gamm*I-theta*I
    #       dR <- gamm*I-theta*R
    #       return(list(c(dS,dE,dI,dR)))})}
    #   out<-ode(y=init,times=time,eqn,parms=parameters)
    #   out.df<-as.data.frame(out)
    # 
    # } else if (model=="SEIRS" & !isFALSE(vital)){
    #   init <- c(S = 0.1, E = 1e-04, I = 1e-04, R = 1 - 0.1 - 1e-4 - 1e-4)
    #   parameters <- c(bet=b,gamm=g, sig=s, epsi=e)
    #   time <- seq(0,t,by=t/(2*length(1:t)))
    #   eqn <- function(time,state,parameters){
    #     with(as.list(c(state,parameters)),{
    #       dS <- -bet*S*I+epsi*E
    #       dE <- bet*S*I-sig*E
    #       dI <- sig*E-gamm*I
    #       dR <- gamm*I-epsi*R
    #       return(list(c(dS,dE,dI,dR)))})}
    #   out<-ode(y=init,times=time,eqn,parms=parameters)
    #   out.df<-as.data.frame(out)
    # 
    # } else if (model=="SEIRS" & isFALSE(vital)){
    #   init <- c(S = 0.1, E = 1e-04, I = 1e-04, R = 1 - 0.1 - 1e-4 - 1e-4)
    #   parameters <- c(bet=b,gamm=g, sig=s, epsi=e, mu=m, theta=th)
    #   time <- seq(0,t,by=t/(2*length(1:t)))
    #   eqn <- function(time,state,parameters){
    #     with(as.list(c(state,parameters)),{
    #       dS <- mu-bet*S*I+epsi*R-theta*S
    #       dE <- bet*S*I-sig*E-mu*E
    #       dI <- sig*E-gamm*I-theta*I
    #       dR <- gamm*I-epsi*R-theta*R
    #       return(list(c(dS,dE,dI,dR)))})}
    #   out<-ode(y=init,times=time,eqn,parms=parameters)
    #   out.df<-as.data.frame(out)
    # 
    # }
    
    
    
    mytheme4 <- theme_bw() +
      theme(text=element_text(colour="black")) +
      theme(panel.grid = element_line(colour = "white")) +
      theme(panel.background = element_rect(fill = "#B2B2B2"))
    theme_set(mytheme4)
    
    
    if(model=="SIR"& isFALSE(vital)){
      ## ggplot2 data
      title <- bquote("SIR Model w/o Vital Dynamics")
      subtit <- bquote(list(beta==.(parameters[1]),~gamma==.(parameters[2])))
      res<-ggplot(out.df,aes(x=time))+
        ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit))))))+
        geom_line(aes(y=S,colour="Susceptible"))+
        geom_line(aes(y=I,colour="Infected"))+
        geom_line(aes(y=R,colour="Recovered"))+
        ylab(label="Proportion")+
        xlab(label="Time (days)")+
        theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
        theme(legend.title=element_text(size=12,face="bold"),
              legend.background = element_rect(fill='#FFFFFF',
                                               size=0.5,linetype="solid"),
              legend.text=element_text(size=10),
              legend.key=element_rect(colour="#FFFFFF",
                                      fill='#C2C2C2',
                                      size=0.25,
                                      linetype="solid"))+
        scale_colour_manual("Compartments",
                            breaks=c("Susceptible","Infected","Recovered"),
                            values=c("blue","red","darkgreen"))
      
      print(res)
    } else if (model=="SIR" & !isFALSE(vital)){
      ## ggplot2 data
      title <- bquote("SIR Model w/ Vital Dynamics")
      subtit <- bquote(list(beta==.(parameters[1]),~gamma==.(parameters[2]), ~mu==~nu==.(parameters[3])))
      res<-ggplot(out.df,aes(x=time))+
        ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit))))))+
        geom_line(aes(y=S,colour="Susceptible"))+
        geom_line(aes(y=I,colour="Infected"))+
        geom_line(aes(y=R,colour="Recovered"))+
        ylab(label="Proportion")+
        xlab(label="Time (days)")+
        theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
        theme(legend.title=element_text(size=12,face="bold"),
              legend.background = element_rect(fill='#FFFFFF',
                                               size=0.5,linetype="solid"),
              legend.text=element_text(size=10),
              legend.key=element_rect(colour="#FFFFFF",
                                      fill='#C2C2C2',
                                      size=0.25,
                                      linetype="solid"))+
        scale_colour_manual("Compartments",
                            breaks=c("Susceptible","Infected","Recovered"),
                            values=c("blue","red","darkgreen"))
      
      print(res)
    } else if (model=="SIRS"& isFALSE(vital)){
      ## ggplot2 data
      title <- bquote("SIRS Model w/o Vital Dynamics")
      subtit <- bquote(list(beta==.(parameters[1]),~gamma==.(parameters[2])))
      res<-ggplot(out.df,aes(x=time))+
        ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit))))))+
        geom_line(aes(y=S,colour="Susceptible"))+
        geom_line(aes(y=I,colour="Infected"))+
        geom_line(aes(y=R,colour="Recovered"))+
        ylab(label="Proportion")+
        xlab(label="Time (days)")+
        theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
        theme(legend.title=element_text(size=12,face="bold"),
              legend.background = element_rect(fill='#FFFFFF',
                                               size=0.5,linetype="solid"),
              legend.text=element_text(size=10),
              legend.key=element_rect(colour="#FFFFFF",
                                      fill='#C2C2C2',
                                      size=0.25,
                                      linetype="solid"))+
        scale_colour_manual("Compartments",
                            breaks=c("Susceptible","Infected","Recovered"),
                            values=c("blue","red","darkgreen"))
      
      print(res)
    } else if (model=="SIRS"& !isFALSE(vital)){
      ## ggplot2 data
      title <- bquote("SIRS Model w/ Vital Dynamics")
      subtit <- bquote(list(beta==.(parameters[1]),~gamma==.(parameters[2]), ~mu==~nu==.(parameters[3]), ~epsilon==.(parameters[5])))
      res<-ggplot(out.df,aes(x=time))+
        ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit))))))+
        geom_line(aes(y=S,colour="Susceptible"))+
        geom_line(aes(y=I,colour="Infected"))+
        geom_line(aes(y=R,colour="Recovered"))+
        ylab(label="Proportion")+
        xlab(label="Time (days)")+
        theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
        theme(legend.title=element_text(size=12,face="bold"),
              legend.background = element_rect(fill='#FFFFFF',
                                               size=0.5,linetype="solid"),
              legend.text=element_text(size=10),
              legend.key=element_rect(colour="#FFFFFF",
                                      fill='#C2C2C2',
                                      size=0.25,
                                      linetype="solid"))+
        scale_colour_manual("Compartments",
                            breaks=c("Susceptible","Infected","Recovered"),
                            values=c("blue","red","darkgreen"))
      print(res)
      } 
    # else if (model=="SIER"& isFALSE(vital)){
    #   ## ggplot2 data
    #   title <- bquote("SIER Model w/ Vital Dynamics")
    #   subtit <- bquote(list(beta==.(parameters[1]),~gamma==.(parameters[2])))
    #   res<-ggplot(out.df,aes(x=time))+
    #     ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit))))))+
    #     geom_line(aes(y=S,colour="Susceptible"))+
    #     geom_line(aes(y=I,colour="Infected"))+
    #     geom_line(aes(y=R,colour="Recovered"))+
    #     geom_line(aes(y=E,colour="Exposed"))+
    #     ylab(label="Proportion")+
    #     xlab(label="Time (days)")+
    #     theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
    #     theme(legend.title=element_text(size=12,face="bold"),
    #           legend.background = element_rect(fill='#FFFFFF',
    #                                            size=0.5,linetype="solid"),
    #           legend.text=element_text(size=10),
    #           legend.key=element_rect(colour="#FFFFFF",
    #                                   fill='#C2C2C2',
    #                                   size=0.25,
    #                                   linetype="solid"))+
    #     scale_colour_manual("Compartments",
    #                         breaks=c("Susceptible","Infected","Recovered","Exposed"),
    #                         values=c("blue","red","darkgreen","yellow"))
    #   
    #   print(res)
    # } else if (model=="SIER"& !isFALSE(vital)){
    #   ## ggplot2 data
    #   title <- bquote("SIER Model w/o Vital Dynamics")
    #   subtit <- bquote(list(beta==.(parameters[1]),~gamma==.(parameters[2])))
    #   res<-ggplot(out.df,aes(x=time))+
    #     ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit))))))+
    #     geom_line(aes(y=S,colour="Susceptible"))+
    #     geom_line(aes(y=I,colour="Infected"))+
    #     geom_line(aes(y=R,colour="Recovered"))+
    #     geom_line(aes(y=E,colour="Exposed"))+
    #     ylab(label="Proportion")+
    #     xlab(label="Time (days)")+
    #     theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
    #     theme(legend.title=element_text(size=12,face="bold"),
    #           legend.background = element_rect(fill='#FFFFFF',
    #                                            size=0.5,linetype="solid"),
    #           legend.text=element_text(size=10),
    #           legend.key=element_rect(colour="#FFFFFF",
    #                                   fill='#C2C2C2',
    #                                   size=0.25,
    #                                   linetype="solid"))+
    #     scale_colour_manual("Compartments",
    #                         breaks=c("Susceptible","Infected","Recovered","Exposed"),
    #                         values=c("blue","red","darkgreen","yellow"))
    #   
    #   print(res)
    # }
     
    ## ggplot2 data
    # title <- bquote("SIR Model")
    # subtit <- bquote(list(beta==.(parameters[1]),~gamma==.(parameters[2])))
    # res<-ggplot(out.df,aes(x=time))+
    #   ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit))))))+
    #   geom_line(aes(y=S,colour="Susceptible"))+
    #   geom_line(aes(y=I,colour="Infected"))+
    #   geom_line(aes(y=R,colour="Recovered"))+
    #   ylab(label="Proportion")+
    #   xlab(label="Time (days)")+
    #   theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
    #   theme(legend.title=element_text(size=12,face="bold"),
    #         legend.background = element_rect(fill='#FFFFFF',
    #                                          size=0.5,linetype="solid"),
    #         legend.text=element_text(size=10),
    #         legend.key=element_rect(colour="#FFFFFF",
    #                                 fill='#C2C2C2',
    #                                 size=0.25,
    #                                 linetype="solid"))+
    #   scale_colour_manual("Compartments",
    #                       breaks=c("Susceptible","Infected","Recovered"),
    #                       values=c("blue","red","darkgreen"))
    # 
    # print(res)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

