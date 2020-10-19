# load shiny 
library(shiny)


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Simple Shiny!"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "Kpmu",
                  label = "Kpmu:",
                  min = 0.1,
                  max = 1,
                  value = 0.5,
                  step = 0.1)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "PBPKPlot")
      
    )
  )
)


##########################################################


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  output$PBPKPlot <- renderPlot({
    
    #load observed IV infusion data
    obs <- read.csv("../data/source/Adult_IV.csv")
    
    #set simulation conditions
    bw   <- 73
    amt  <- 4*bw
    rate <- 4*bw
    cmt  <- "VEN"
    ii   <- 12
    addl <- 13
    ss   <- 1
    
    modA <- mread("../model/voriPBPK")
    
    #run simulation
    sim <- 
      modA %>% 
      param(Kpmu = input$Kpmu) %>%
      ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss) %>% 
      mrgsim(delta = 0.1, end = 12) %>% 
      filter(row_number() != 1)  
    
    #plot prediction and compare to observed data
    gp <- ggplot() + 
      geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
      geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
      geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
      scale_colour_manual(name='', 
                          values=c('sim'='black', 'observed'='black'), 
                          breaks=c("observed","sim"),
                          labels=c("observed","predicted")) +
      guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
      labs(title="Adult 4 mg/kg IV", x="time (h)", y="Plasma concentration (mg/L)") +
      theme_bw()
    gp
    
  })
  
}


##########################################################


shinyApp(ui, server)

