#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

load("./pbmc_tutorial.rds")
exprs <- as.data.frame(pbmc@scale.data)
target <- as.factor((pbmc@ident))
rm(list = 'pbmc')

plot_binarization <- function(BLA, threshold, GOI){
  THR <- 1 / threshold
  intercept <- (max(BLA[GOI]) - min(BLA[GOI])) / THR + min(BLA[GOI])
  BLA$highexprs <-
    as.factor(BLA[GOI] > intercept)
  p <- ggplot(BLA, aes(x = get(GOI), fill = highexprs)) +
    scale_fill_manual(values = c("blue", "red")) +
    geom_histogram(binwidth = 0.01) +
    geom_vline(
      xintercept = intercept,
      colour = 'red',
      linetype = 'dashed'
    ) +
    theme(legend.position="none")+
    ggtitle(paste0('Binarization for threshold = ',threshold)) +
    ylab('Single Cell Counts') +
    xlab(paste(GOI, 'normalized expression'))
  
  p
  
  
}

BLA <-
  as.data.frame(t(exprs[,]))
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Binarization"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("bins",
                     "Threshold",
                     min = 0,
                     max = 1,
                     step = 0.025,
                     value = 0.5)
      ,
      selectizeInput(inputId = "gene",
                  choices = rownames(exprs),
                  label = "Gene of Interest"
                  
                  )
      
      ),
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
     plot_binarization(BLA, threshold = input$bins, GOI = input$gene)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

