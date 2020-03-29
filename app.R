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
      
      numericInput("time",
                   "Timeline:",
                   min = 0,
                   max = 1000,
                   value = 100)
      
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("epidemic"),
      strong("Instructions on how to use are as follows..."),
      p("The social distancing factor slider ranges from 0 to 5 with 0 being complete isolation as in the infection rate is effectivly zero and with 5 being 
        an almost uncontrollable infection. If the social distancing factor is low that means we are flattening the curve. If the rate is high then we aren't 
        effectivley doing that. The rate of recovery is, much like it sounds, how long does it take for someone to recover and thus be added to the recovered pool
        of people. The timeline input is how far out you would like to look based on the other initial conditions."),
      p("The initial conditions when you open the app are that the infected multiply by 2.16 times a day and it takes 2 days to recover from said infection.")
      
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
    
    ##ode solution
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
    
    mytheme4 <- theme_bw() +
      theme(text=element_text(colour="black")) +
      theme(panel.grid = element_line(colour = "white")) +
      theme(panel.background = element_rect(fill = "#B2B2B2"))
    theme_set(mytheme4)
    
    ## ggplot2 data
    title <- bquote("SIR Model")
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
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

