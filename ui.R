#-----------------------
# 
# Shiny app for visualizing the data in Figures 1A-1C of [Peck and Lauring 2018](http://jvi.asm.org/content/early/2018/04/26/JVI.01031-17.short)
#
# Developed by Kayla Peck, 18.05.23
#
#-----------------------

library(shiny)
shinyUI(fluidPage(
  #Main title
  titlePanel("The complexities of viral mutation rates."),
  #Sidebar
  sidebarLayout(position="right",
    sidebarPanel(
      selectInput("plot", "Choose figure to view:", 
                  choices=c("1A: Evolution vs. mutation rate (Baltimore classes)",
                            "1B: Evolution vs. mutation rate (individual viruses)",
                            "1C: Mutation rate vs. genome size")),
      #selectInput("print", "Choose a value to print:",
      #            choices=c("None","Evolution rate", "Mutation rate", "Virus class")),
      br(),
      #h5("Virus key"),
      #p("TMV = tobacco mosaic virus"),
      #p("polio1 = Poliovirus type 1"),
      #p("HCV = hepatitis C virus"),
      #p("fluva = Influenza A virus"),
      #p("HIV = Human immunodeficiency virus"),
      #p("HSV1 = herpes simplex virus 1"),
      #p("AHBV = avian hepatitis B virus"),
      hr()
    ),
    #Spot for the plot
    mainPanel(
      h5("Peck and Lauring 2018."),
      p("Many viruses evolve rapidly. This is due, in part, to their high mutation rates. Mutation rate estimates for over 
        25 viruses are currently available. Here, we review the population genetics of virus mutation rates. We specifically 
        cover the topics of mutation rate estimation, the forces that drive the evolution of mutation rates, and how the optimal mutation rate can be context-dependent."),
      p("Choose a plot you would like to view using the right drop-down menu. Hover over points for more information.", style="color:blue"),
      plotOutput("virusPlot")
    )
  )))

