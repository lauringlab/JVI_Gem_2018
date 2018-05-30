#-----------------------
# 
# Shiny app for visualizing the data in Figures 1A-1C of [Peck and Lauring 2018](http://jvi.asm.org/content/early/2018/04/26/JVI.01031-17.short)
#
# Developed by Kayla Peck, 18.05.23
#
#-----------------------

library(shiny)

#Read in the data and prepare the data frames needed for each plot
dat <- read.csv("Figure_1_mu_and_K_data.csv")
dat <- dat[,c(1:6,8,10)]

datL <- read.csv("Lynch_2016_mu_data.csv")
names(datL)[1] <- "group"
datL$Reference <- "Lynch et al. 2016"

dat$U <- dat$G*1000*dat$mu 
datL$U <- datL$G*1000000*datL$mu 

K <- mu <- G <- 0
for(i in 1:length(unique(dat$group))){
  sub <- subset(dat, group==unique(dat$group)[i])
  K[i] <- 10^mean(log10(sub$K), na.rm=T)
  mu[i] <- 10^mean(log10(sub$mu), na.rm=T)
  G[i] <- 10^mean(log10(sub$G), na.rm=T)
}

#Set up color palettes
cols <- c("gold","forestgreen","orangered","firebrick","dodgerblue", "darkorchid")
colsL <- c("black", "white","grey70")
dat$colors <- "dodgerblue"
for(i in 1:length(dat$group)){
  if(dat$group[i]=="ss_neg_RNA"){
    dat$colors[i] <- "firebrick"
  }
  if(dat$group[i]=="dsRNA"){
    dat$colors[i] <- "forestgreen"
  }
  if(dat$group[i]=="retro"){
    dat$colors[i] <- "orangered"
  }
  if(dat$group[i]=="ssDNA"){
    dat$colors[i] <- "darkorchid"
  }
  if(dat$group[i]=="dsDNA"){
    dat$colors[i] <- "gold"
  }
}
datL$colors <- "black"
for(i in 1:length(datL$group)){
  if(datL$group[i]=="eubacteria"){
    datL$colors[i] <- "grey70"
  }
  if(datL$group[i]=="unicellular"){
    datL$colors[i] <- "white"
  }
}

#create organized data frames for shiny data retrieval
#Figure 1 data
fig1dat <- data.frame(K=K, mu=mu, group=unique(dat$group))
fig1dat$colors <- cols

#Figure 2 data
fig2dat <- na.omit(dat)

#Figure 3 data
fig3dat <- data.frame(group=c(as.character(dat$group),as.character(datL$group)), 
                           species=c(as.character(dat$virus),as.character(datL$species)), G=c(dat$G*1000, datL$G*1e6),
                           mu = c(dat$mu, datL$mu), U=c(dat$U,datL$U), 
                      Reference=c(as.character(dat$Full.mu.reference),as.character(datL$Reference)),
                      colors=c(as.character(dat$colors),as.character(datL$colors)))
fig3dat <- na.omit(fig3dat)

#-----------------
# Start shiny plot!
#-----------------

shinyServer(function(input,output,session)
{
  par(mar=c(5.1,4.1,0,0))
  #renderPlot indicates that the function is "reactive" - it should automatically
  #re-execute when the input changes
  current_data <- reactive({
    if(input$plot == "1A: Evolution vs. mutation rate (Baltimore classes)"){
      return(fig1dat)
    }
    if(input$plot == "1B: Evolution vs. mutation rate (individual viruses)"){
      return(fig2dat)
    }
    if(input$plot == "1C: Mutation rate vs. genome size"){
      return(fig3dat)
    }
  })
  current_x <- reactive({
    if(input$plot == "1A: Evolution vs. mutation rate (Baltimore classes)" | input$plot == "1B: Evolution vs. mutation rate (individual viruses)"){
      return("mu")
    }
    if(input$plot == "1C: Mutation rate vs. genome size"){
      return("G")
    }
  })
  current_y <- reactive({
    if(input$plot == "1A: Evolution vs. mutation rate (Baltimore classes)" | input$plot == "1B: Evolution vs. mutation rate (individual viruses)"){
      return("K")
    }
    if(input$plot == "1C: Mutation rate vs. genome size"){
      return("mu")
    }
  })
  
  output$virusPlot <- renderPlot({
    
    #render the plot
    if(input$plot == "1A: Evolution vs. mutation rate (Baltimore classes)"){
      par(mar=c(5,6,1,1))
      plot(K~mu, dat=fig1dat, ylim=c(1e-5,1e-2), xlim=c(1e-7,1e-4), log="xy", pch=21, bg=fig1dat$colors, cex=2,
           ylab=expression(paste("Evolutionary rate (s/n/y)")), xlab=expression(paste("Mutation rate (s/n/c)")), 
           xaxt='n', yaxt='n', cex.lab=2)
      axis(1, cex.axis=1.25, at=c(1e-7,1e-6,1e-5,1e-4), labels=c(expression(paste(10^{-7})), expression(paste(10^{-6})), expression(paste(10^{-5})), expression(paste(10^{-4}))))
      axis(2, cex.axis=1.25, las=2, at=c(1e-5,1e-4,1e-3,1e-2), labels=c(expression(paste(10^{-5})), expression(paste(10^{-4})), expression(paste(10^{-3})), expression(paste(10^{-2}))))
      legend("topleft", ncol=2, c("dsDNA","dsRNA","retro","(-)ssRNA", "(+)ssRNA", "ssDNA"), 
             pch=21, pt.bg=cols, cex=1.5, bty='n')
    }
    if(input$plot == "1B: Evolution vs. mutation rate (individual viruses)"){
      par(mar=c(5,6,1,1))
      plot(K~mu, dat=fig2dat, ylim=c(1e-5,1e-1), xlim=c(2e-8,2e-4), log="xy", pch=c(21,21,23,21,22,21,22,23,24,25),
           bg=fig2dat$colors,
           ylab=expression(paste("Evolutionary rate (s/n/y)")), xlab=expression(paste("Mutation rate (s/n/c)")), 
           xaxt='n', yaxt='n',cex=2,cex.lab=2)
      axis(1, cex.axis=1.25, at=c(1e-7,1e-6,1e-5,1e-4), labels=c(expression(paste(10^{-7})), expression(paste(10^{-6})), expression(paste(10^{-5})), expression(paste(10^{-4}))))
      axis(2, cex.axis=1.25, las=2, at=c(1e-5,1e-4,1e-3,1e-2,1e-1), labels=c(expression(paste(10^{-5})), expression(paste(10^{-4})), expression(paste(10^{-3})), expression(paste(10^{-2})), expression(paste(10^{-1}))))
      legend("topleft", ncol=2, c("Tobacco mosaic virus", "Human rhinovirus", "Poliovirus 1",
                                  "Human norovirus", "Hepatitis C virus", "Influenza A virus",
                                  "Measles virus", "Avian HBV", "HIV-1", "Herpes simplex virus"), 
             pt.bg=c(rep("dodgerblue",5),rep("firebrick",2),rep("orangered",2),"gold"), pch=c(21,22,23,24,25,21,22,21,23,21),
             cex=1.2, bty='n')
    }
    if(input$plot == "1C: Mutation rate vs. genome size"){
      par(mar=c(5,6,1,1))
      plot(mu~G, dat=fig3dat, ylim=c(1e-11,1e-1), xlim=c(1e3,1e10), log="xy",
           ylab=expression(paste("Mutation rate (s/n/c or s/n/g)")), xlab=expression(paste("G (bp)")), 
           xaxt='n', yaxt='n', bg=as.character(fig3dat$colors), pch=21, cex=1.5, cex.lab=2)
      axis(1, cex.axis=1.25, at=c(1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10), labels=c(expression(paste(10^{3})),expression(paste(10^{4})),expression(paste(10^{5})),expression(paste(10^{6})),expression(paste(10^{7})),expression(paste(10^{8})),expression(paste(10^{9})),expression(paste(10^{10}))))
      axis(2, cex.axis=1.25, las=2, at=c(1e-11, 1e-9,1e-7,1e-5,1e-3,1e-1), labels=c(expression(paste(10^{-11})), expression(paste(10^{-9})), expression(paste(10^{-7})), expression(paste(10^{-5})), expression(paste(10^{-3})), expression(paste(10^{-1}))))
      legend("topright", ncol=2, c("multicellular","unicellular", "eubacteria","", "", "", "dsDNA", "dsRNA", "retro", "(-)ssRNA", "(+)ssRNA", "ssDNA"), 
             pch=21, pt.bg=c(colsL,rep("white",3),cols), col=c(rep(1,3),rep("white",3),rep(1,6)),cex=1.5, bty='n')
    }
    
    #Make the plot brushable
    output$brush_info <- renderDataTable({
      if(is.null(input$virusPlot_brush))
        return()
      if(input$plot == "1A: Evolution vs. mutation rate (Baltimore classes)"){
        selected_points <- brushedPoints(current_data(), input$virusPlot_brush, current_x(), current_y())
        temp <- dat[dat$group %in% selected_points$group,]
        temp
      }
      else{
        brushedPoints(current_data(), input$virusPlot_brush, current_x(), current_y())
      }
    })
    
    #Plot summary box for selected points
     output$brush_info2 <- renderPrint({
       if(is.null(input$virusPlot_brush))
         return()
       if(input$plot == "1A: Evolution vs. mutation rate (Baltimore classes)"){
         brushedPoints(current_data(), input$virusPlot_brush, current_x(), current_y())
       }
     })
    
  })
})

#--------------------------------
# Old code for showing data when the mouse hovers over a data point
#--------------------------------
# output$info <- renderText({
#   #Function to obtain the closest coordinate point
#   hoverValue <- function(hover,x,y,other=NULL,tolerance=0.05){
#     if(!is.null(hover)){
#       x0 <- hover$x #x coordinate in user space
#       y0 <- hover$y #y coordinate in user space
#       xrange <- hover$domain$right - hover$domain$left
#       yrange <- hover$domain$top - hover$domain$bottom
#       #find the observation closest to the user coordinates
#       dist <- abs(x0-x)/xrange + abs(y0-y)/yrange
#       i <- which.min(dist)
#       #return corresponding index if close enough
#       if(dist[i] < tolerance){
#         #cat(captions[i])
#         convx <- formatC(10^as.numeric(x[i]),format="e", digits=2)
#         convy <- formatC(10^as.numeric(y[i]), format="e", digits=2)
#         paste0("Group = ", other[i], "\n",
#                "mu = ", convx, "\n",
#                "K = ", convy)
#       }
#     }
#   }
#   
#   hoverValue(input$plot_hover, muKdat$mu, muKdat$K, muKdat$group)