# Analysis of Motif Distributions in Endocytic Proteins 
#
#
#

library(shiny)
library(tidyverse)
library(ggpubr)

data <- read_csv("human_proteome_motifs_across_domains.csv")
data %<>%
  group_by(sequence_id,domain_type,motif_type) %>%
  summarize(count = n()) %>%
  spread(domain_type, count) %>%
  replace_na(list(D = 0, O = 0)) %>% 
  arrange(desc(D))

ui <- fluidPage(
  
  titlePanel("Analysis of Motif Distributions in Endocytic Proteins"),
  
  sidebarPanel(
    
    helpText("This shiny app shows the distributions of selected motifs in disordered versus ordered regions of proteins across the human proteome."),
    ## -> input$motif
    selectInput('motif', 'Which motif to show?', choices=(c("NPXY","NPF","DX[FW]","YXX[LIMFV]","[ST]XXXX[LI]","X[DE]XXXL[LI]")), selected="NPF"),
    
    selectInput('x_max', 'What should the higher limit on the x axis be?', choices=(c(1:10)), selected=5),
    selectInput('y_max', 'What should the higher limit on the y axis be?', choices=(c(1:30)), selected=10)
    
  ),
  
  mainPanel(
    plotOutput('plot')
  )
)



server <- function(input, output) {
  
  dataset <- reactive({
    data %>% filter(motif_type == input$motif) %>%
      group_by(D,O) %>%
      tally() %>%
      mutate(xmin = O-0.5, ymin = D-0.5, xmax = xmin + 1, ymax = ymin+1)
      
  }) 
  
  output$plot <- renderPlot({
    
    p <- ggplot(dataset()) + 
      geom_rect(color = "white", aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = n)) + 
      geom_text(color = "white", aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label = n)) + 
      scale_x_continuous("Ordered Motif Count", breaks = 0:input$x_max) +
      scale_y_continuous("Disordered Motif Count", breaks = 0:input$y_max) +
      grids(linetype = "dashed") +
      scale_fill_continuous(name = "Protein Count")
    
    print(p)
    
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

