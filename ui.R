#-----------------------
# 
# Shiny app for visualizing the data in Figures 1A-1C of [Peck and Lauring 2018](http://jvi.asm.org/content/early/2018/04/26/JVI.01031-17.short)
#
# Developed by Kayla Peck, 18.05.23
#
#-----------------------

library(shiny)
shinyUI(fluidPage(
  #Sidebar
  sidebarLayout(position="left",
    sidebarPanel(
      selectInput("plot", "Choose figure to view:", 
                  choices=c("1A: Evolution vs. mutation rate (Baltimore classes)",
                            "1B: Evolution vs. mutation rate (individual viruses)",
                            "1C: Mutation rate vs. genome size")),
      #selectInput("print", "Choose a value to print:",
      #            choices=c("None","Evolution rate", "Mutation rate", "Virus class")),
      #br(),
      #hr()
      checkboxInput("all", "Include data collected after 2018 publication")
    ),
    #Spot for the plot
    mainPanel(
      h1("The complexities of viral mutation rates"),
      h3("Peck and Lauring 2018"),
      br(),
      h4("Abstract: Many viruses evolve rapidly. This is due, in part, to their high mutation rates. Mutation rate estimates for over 
        25 viruses are currently available. Here, we review the population genetics of virus mutation rates. We specifically 
        cover the topics of mutation rate estimation, the forces that drive the evolution of mutation rates, and how the optimal mutation rate can be context-dependent."),
      br(),
      h4("Choose a plot you would like to view using the left drop-down menu. Draw a box around points to explore the raw data.", style="color:blue"),
      br(),
      plotOutput("virusPlot",
                 height = "500px",
                 width = "600px",
                 dblclick = "virusPlot_dblclick",
                 brush = brushOpts(
                   id = "virusPlot_brush"
                 )),
      uiOutput("legend.ui"),
       fluidRow(
         tags$head(tags$style(type="text/css",
                              "#brush_info {font-size: 8px}")),
         br(),
         verbatimTextOutput("brush_info2"),
        fluidRow(
          tags$head(tags$style(type="text/css",
                               "#brush_info {font-size: 12px}")),
          dataTableOutput("brush_info"))
        )
      )
    )
  )
)

