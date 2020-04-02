SIR.model <- function(t, b, g){
  library(deSolve)
  library(ggplot2)
  
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
  ggsave(plot=res,
         filename=paste0("SIRplot_","time",t,"beta",b,"gamma",g,".png"),
         width=8,height=6,dpi=180)
  }

SIR.model(70,1.4,0.3)



parameters <- c(mu = 1 / (70 * 365), beta = 520 / 365,
                sigma = 1 / 14, gamma = 1 / 7)
initials <- c(S = 0.1, E = 1e-04, I = 1e-04, R = 1 - 0.1 - 1e-4 - 1e-4)

# Uncomment the following lines (running it takes more than a few seconds):
seir <- SEIR(pars = parameters, init = initials, time = 0:(60 * 365))
PlotMods(seir)


# Parameters and initial conditions.
parameters <- c(beta = 1.4247, gamma = 0.14286)
initials <- c(S = 1 - 1e-06, I = 1e-06, R = 1 - (1 - 1e-06 - 1e-06))


# Solve and plot.
sir <- SIR(pars = parameters, init = initials, time = 0:70)
PlotMods(sir)

test <-seir[["results"]]
test
