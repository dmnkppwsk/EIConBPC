

#Packages
library(shiny)
library(shinysense)
library(ggplot2)
library(DT)
library(plotly)
library(xcms)

#XChromatograms to data frame functions
.ChromatogramLongFormat <- function(x){
  data.frame(rtime = x@rtime, 
             intensity = x@intensity,
             spectrum_id = names(x@rtime),
             fromFile = x@fromFile, 
             mz_interval = paste0(x@filterMz[1]," - ",paste0(x@filterMz[2])), 
             row.names = NULL, 
             stringsAsFactors = FALSE
  )
}

.ChromatogramsLongFormat <- function(x, pdata = NULL){
  table <- lapply(x@.Data, .ChromatogramLongFormat)
  table <- do.call(rbind.data.frame,table)
  
  table$mz_interval <- as.factor(table$mz_interval)
  
  # add filename
  filename <- cbind.data.frame(fromFile =  1:attributes(x)$dim[2], filename = attributes(x)$dimnames[[2]], stringsAsFactors = FALSE)
  table <- merge(table, filename, all.y = TRUE,  by = "fromFile")
  
  # if we supplied pdata add the corresponding annotation to each row
  if(!is.null(pdata)){
    
    # we merge the pdata columns to be able to supply a text column that can be shown in plotly tooltips
    text <- lapply(names(pdata@data), function(x) paste0(x,": " , as.character(as.matrix(pdata@data[x])),"<br>"))
    text <- apply(as.data.frame(text),1,paste, collapse="")
    
    pdata <- cbind.data.frame(fromFile =  1:nrow(pdata), pdata@data, text = text, stringsAsFactors = FALSE)
    
    table <- merge(table, pdata, all.y = TRUE,  by = "fromFile")
  }
  
  return(table)
}


#Set working directory
setwd("D:/Natalia/Dataset")

#Load raw data (after peak filling)
load("xdata.Rdata")
#Load metadata
load("short_Experimental.Rdata")
load("short_Experimental2.Rdata")
load("short_Experimental3.Rdata")
#Load BPC data frames
for (i in short_Experimental3$df_filename){load(i)}
#Load Tinderesting results (with isotopes removed)
load("AllResults_zonderisotopen.RData")
#Generate choices for selecting the sample
time<-unique(short_Experimental$`short_Experimental$time`)
Replicate<-unique(short_Experimental$`short_Experimental$Replicate`)


# UI
ui <- fluidPage(

    # Application title
    titlePanel("Metabolite candidate revision on BPC"),

    # Inputs
    sidebarLayout(
        sidebarPanel(
            numericInput("mz","Target m/z:", value = 623.10),
            sliderInput("pm_mz", "+/- m/z:", min=0, max=0.1, value=0.05, step=0.001),
            selectInput("time", "Timepoint (h):", selected=min(time), choices=time),
            #selectInput("repilcate", "Replicate", selected=Replicate[1], choices=Replicate),
            numericInput("rowname", "Selected sample (row_name)", value = 10),
            actionButton("plot_button", "Generate EIC on BPC")),
    #Output tabs        
    mainPanel(tabsetPanel(
          tabPanel("Choose m/z",DT::DTOutput('table_metabolites')),
          tabPanel("Choose sample", DT::DTOutput('table_exp'), textOutput("selected")),
          tabPanel("Chromatogram",plotly::plotlyOutput('plot'))
        ))
    ))


# Define server logic
server <- function(input, output) {
  #Metadata table to choose sample number (rowname)
    table_exp <- function(){
      short_Experimental3[,c(9,7,4,1)] %>% 
      filter(time == input$time)
    }
    output$table_exp <- DT::renderDT({table_exp()})
    #Tinderesting results table
    table_metabolites <- AllResults_zonderisotopen[,c(17,1,15,14,13,11,12,10,16)]
    table_metabolites <- table_metabolites%>%
      arrange(tind_score_order)
    output$table_metabolites <- DT::renderDT({table_metabolites})
    
    
   #Plot
    
    observeEvent(input$plot_button, {
      sn<-input$rowname
      mzr<-input$mz+c(-(input$pm_mz), input$pm_mz)
      chr_sn <- chromatogram(filterFile(xdata, sn), mz = mzr, aggregationFun = "max")
      df_EIC <- .ChromatogramsLongFormat(chr_sn)
      dm2<-df_EIC
      dm2[is.na(dm2<-df_EIC)] <- 0
      df_EIC<-dm2
      nam_df_bpis <- paste("df_bpis", sn, sep = "_")
      target_df_bpis<-get(nam_df_bpis)
      target_df_bpis$intensityEIC <- df_EIC$intensity
      
      
      plotEoB<-function () {
        plot<-ggplotly(
        ggplot(data = target_df_bpis, aes(rtime/60)) +
        geom_line(aes(y=intensity, colour="BPC"))+
        geom_line(aes(y=intensityEIC, colour="EIC"))+
        theme_classic() +
        labs(x = "Retention time", y = "Intensity"))}
      # }

    
    output$plot <-plotly::renderPlotly({plotEoB()})
    })
      
      
  
}

# Run the application 
shinyApp(ui = ui, server = server)

