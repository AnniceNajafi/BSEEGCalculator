#' Authors@R: person("Annice", "Najafi", email = "annice-najafi@uiowa.edu")
#' The University of Iowa, Carver College of Medicine, Department of Psychiatry
#' Spring 2020
#' 
#'Discription: This is a web-based application designed to detect delirium in mice. For the application to work, 
#'bispectral EEG files in the European Data Format have to be uploaded to the application.

#Load related libraries
library(shiny)
library(shinycssloaders)
library(ggplot2)
library(dplyr)
library(future)
library(remotes)
library(devtools)
library(shinythemes)
library(eegUtils)
library(edfReader)
library(gridExtra)
library(DT)
library(shinymanager)
library(data.table)
library(tidyverse)
library(shinyWidgets)
library(shinyjs)
#library(edfReader)
#library(suncalc)
#library(lubridate)
#increase the maximum size of the file that can be received by the app to 300 MB
options(shiny.maxRequestSize=2000*1024^2)
#' Function to process EEG files for a day for the detection of delirium
#'
#' @param mydata EEG file in the form of .edf stored in a variable
#' @return A dataframe with calculated BSEEG scores
#' 

slReadEDF <- function(mydata, filldis, filldisrupt) {
  if(is.null(colnames(mydata$signals))){
    #MAKE SURE YOU UNCOMMENT THE COMMAND BELOW
    #eeg_example <- import_raw("Matt-task-spatcue.bdf")
  }
  #get the number of hours
  srate <- mydata$srate
  totalSeconds <- nrow(mydata$timings)
  numHrs <- totalSeconds/(srate*3600)
  #parse the data into 24 hours and 10 minutes intervals
  numHrs <- floor(numHrs)
  thresh <- floor(totalSeconds/(numHrs*6))
  #define columns for storing ratios and means
  myTable <- data.frame(c(1:thresh))
  #initialize columns to hold the average for the first EEG electrode
  #variable to hold low frequency
  myAvg3EEG1 <- c(NA)
  #variable to hold high frequency
  myAvg10EEG1 <- c(NA)
  #initialize columns to hold the average for the second EEG electrode
  #variable to hold low frequency
  myAvg3EEG2 <- c(NA)
  #variable to hold high frequency
  myAvg10EEG2 <- c(NA)
  #initialize two variables for holding mean values
  finalAvg1 <- c(NA)
  finalAvg2 <- c(NA)
  
  #mydata <- round(mydata$signals,5) 
  #go through all the data for every 10 minute window / 
  thresh2 <- (numHrs*6) - 1
  for(i in 0:thresh2){
    if(is.null(colnames(mydata$signals))){
      storeHere <- eeg_example
      ###split data 
      storeHere$signals <- tibble(file$signal$EEG1$data[c(((i*thresh)+1):(thresh*(i+1)))], file$signal$EEG2$data[c(((i*thresh)+1):(thresh*(i+1)))])
      storeHere$timings = tibble(file$signal$EEG1$t[c(((i*thresh)+1):(thresh*(i+1)))], file$signal$EEG2$t[c(((i*thresh)+1):(thresh*(i+1)))])
    }else{
      
      storeHere <- mydata
      ###split data 
      storeHere$signals <- mydata$signals[c(((i*thresh)+1):(thresh*(i+1))), ]
      storeHere$timings <- mydata$timings[c(((i*thresh)+1):(thresh*(i+1))), ]
    }
    if(filldisrupt){
      h3 <- round(storeHere$signals,5) 
      colnames(h3)[1]<-"EEG1"
      colnames(h3)[2]<-"EEG2"
      which(h3$EEG1 == 0.00045)-> idx
      h3[!(h3$EEG2 == 0.00045),] -> h4
      #h3[(h3$`EEG EEG1A-B` == 0.00045)] <- NA
      storeHere$signals <- h4
      storeHere$timings[idx,]
    }
    if(length(storeHere$signals[1]) != 0){
      ###compute power spectral density 
      holdPSD<- compute_psd(storeHere)
      colnames(holdPSD)[1]<-"EEG1"
      colnames(holdPSD)[2]<-"EEG2"
      ###find range
      hold3Hz <- holdPSD %>% filter(frequency>=3 & frequency<=3.5)
      hold10Hz <- holdPSD %>% filter(frequency>=9.5 & frequency<=10)
      ###find the mean of hold3Hz and hold8Hz for the first electrode
      myAvg3EEG1<- mean(hold3Hz$EEG1)
      myAvg10EEG1<- mean(hold10Hz$EEG1)
      ###find the mean of hold3Hz and hold8Hz for the second electrode
      myAvg3EEG2 <- mean(hold3Hz$EEG2)
      myAvg10EEG2 <- mean(hold10Hz$EEG2)
      ###find the ratio of high to low frequency for the two electrodes 
      finalAvg1[i] = myAvg3EEG1/myAvg10EEG1
      finalAvg2[i] = myAvg3EEG2/myAvg10EEG2
    }else{
      finalAvg1[i] = NA
      finalAvg2[i] = NA
    }
    
  }
  #Initialize variables for holding results of computation
  avgavg1 <- c(NA)
  avgavg2 <- c(NA)
  normavg1 <- c(NA)
  normavg2 <- c(NA)
  #find the average of ratios for every hour
  numhm <- numHrs - 1
  for(k in 0:numhm){
    #find the average for every hour 
    #***Note: We found the ratio of every 10 minutes, there are a total
    #of 6 10-min intervals in every hour, find the avg for those ratios
    valVal1 <- finalAvg1[c((6*(k)+1):(6*(k+1)+1))]
    valVal2 <- finalAvg2[c((6*(k)+1):(6*(k+1)+1))]
    avgavg1[k] = mean(valVal1)
    avgavg2[k] = mean(valVal2)
  }
  #avgavg1 <- avgavg1 %>% drop_na()
  #avgavg2 <- avgavg2 %>% drop_na()
  #Specify a 24 interval for the BSEEG plot
  #myTime <- c(1:numhm)
  myTime <- c(1:23)
  x <- myTime
  #combine ratios and time for the BSEEG plot
  coolEDFTable <- data.frame(time=c(1:24))
  #coolEDFTable <- data.frame(time=c(1:numhm))
  coolEDFTable["EEG1"]<- c(NA)
  coolEDFTable["EEG2"]<- c(NA)
  if(numhm<=24){
    coolEDFTable$"EEG1"[1:numhm]<- avgavg1
    coolEDFTable$"EEG2"[1:numhm]<- avgavg2
  }else{
    coolEDFTable$"EEG1"[1:24]<- avgavg1
    coolEDFTable$"EEG2"[1:24]<- avgavg2
  }
  hold.na.dropped <- coolEDFTable %>% drop_na()
  hold.na.dropped <- subset(hold.na.dropped, !is.nan(hold.na.dropped$EEG1))
  length(hold.na.dropped) -> implen
  coolEDFTableF <- data.frame(time=c(1:12))
  coolEDFTableS <- data.frame(time=c(13:24))
  coolEDFTableF["EEG1"]<- c(NA)
  coolEDFTableF["EEG2"]<- c(NA)
  coolEDFTableS["EEG1"]<- c(NA)
  coolEDFTableS["EEG2"]<- c(NA)
  coolEDFTableF$EEG1 <- coolEDFTable$EEG1[1:12]
  coolEDFTableS$EEG1 <- coolEDFTable$EEG1[13:24]
  coolEDFTableF$EEG2 <- coolEDFTable$EEG2[1:12]
  coolEDFTableS$EEG2 <- coolEDFTable$EEG2[13:24]
  hold.na.dropped.F <- coolEDFTableF
  hold.na.dropped.S <- coolEDFTableS
  hold.na.dropped.F.count <- coolEDFTableF %>% drop_na()
  hold.na.dropped.F.count <- subset(hold.na.dropped.F.count, !is.nan(hold.na.dropped.F.count$EEG1))
  hold.na.dropped.S.count <- coolEDFTableS %>% drop_na()
  hold.na.dropped.S.count <- subset(hold.na.dropped.S.count, !is.nan(hold.na.dropped.S.count$EEG1))
  if(filldis){
    
    if(length(hold.na.dropped.F.count$EEG1) >= 3){
      hold.na.dropped.F$EEG1 <- as.numeric(hold.na.dropped.F$EEG1)
      hold.na.dropped.F$EEG2 <- as.numeric(hold.na.dropped.F$EEG2)
      hold.na.dropped.F$EEG1[is.nan(hold.na.dropped.F$EEG1)] <- NA
      hold.na.dropped.F$EEG2[is.nan(hold.na.dropped.F$EEG2)] <- NA
      hold.na.dropped.F$EEG1[is.na(hold.na.dropped.F$EEG1)] <- mean(hold.na.dropped.F$EEG1, na.rm = TRUE)
      hold.na.dropped.F$EEG2[is.na(hold.na.dropped.F$EEG2)] <- mean(hold.na.dropped.F$EEG2, na.rm = TRUE)
      print(hold.na.dropped.F)
    }
    if(length(hold.na.dropped.S.count$EEG1) >= 3){
      hold.na.dropped.S$EEG1 <- as.numeric(hold.na.dropped.S$EEG1)
      hold.na.dropped.S$EEG2 <- as.numeric(hold.na.dropped.S$EEG2)
      hold.na.dropped.S$EEG1[is.nan(hold.na.dropped.S$EEG1)] <- NA
      hold.na.dropped.S$EEG2[is.nan(hold.na.dropped.S$EEG2)] <- NA
      hold.na.dropped.S$EEG1[is.na(hold.na.dropped.S$EEG1)] <- mean(hold.na.dropped.S$EEG1, na.rm = TRUE)
      hold.na.dropped.S$EEG2[is.na(hold.na.dropped.S$EEG2)] <- mean(hold.na.dropped.S$EEG2, na.rm = TRUE)
      print(hold.na.dropped.S)
    }
    coolEDFTable <- rbind(hold.na.dropped.F, hold.na.dropped.S)
  }
  #if(!is.na(coolEDFTable$time))
  print(coolEDFTable)
  #return the dataframe as result
  return(c(coolEDFTable, implen))
}

#' The shiny app consists of two parts
#' A user interface part for the user to upload files and see the results
#' server section which processes the edf files received from the user
###########################################USER INTERFACE SIDE######################################################
credentials <- data.frame(
  user = c("admin", "research", "iowa", "iowa", "guest"),
  password = c("pass", "1234", "1234", "iowa", "1234"),
  # comment = c("alsace", "auvergne", "bretagne"), %>% 
  stringsAsFactors = FALSE
)
ui <- secure_app(
  theme = shinythemes::shinytheme("flatly"),
  head_auth = tags$script()
  ,
  tags_bottom = tags$div(
    tags$p(
      "To request for access please contact ",
      tags$a(
        href = "mailto:annicenajafi@tamu.edu?Subject=BSEEG Access Request",
        target="_top", "administrator at annicenajafi@tamu.edu"
      )
    )
  ),
  fluidPage(
    theme=shinytheme("yeti"),
    navbarPage(
      "BSEEG Calculator", id = "del",
      tabPanel("Instructions",
               sidebarLayout(
                 sidebarPanel(
                   tags$img(src='University-of-Iowa-Top-50-Most-Affordable-Executive-MBA-Online-Programs-2019.png', height=100, width=150),
                   br(),
                   tags$style("

             .progress-bar {
             background-color: #CFB53B;
             }

             "),          div(
               tags$head(
                 tags$style(HTML('#refresh1{background-color: #CFB53B; border-color: #CFB53B}'))
               ),
               shinyjs::useShinyjs(),
               id = "filefile1",
               circleButton(
                 inputId = "refresh1",
                 icon = icon("refresh"),
                 status = "primary"
               ),
               checkboxInput("filldis", "Fill incomplete recordings based on recorded hours", value = FALSE),
               checkboxInput("filldisrupt", "Delete disruptions in data", value = FALSE),
               fileInput("fileEdf", "Please upload day 1 file in .edf format"),br(),
               fileInput("fileEdf1", "Please upload day 2 file in .edf format"),br(),
               fileInput("fileEdf2", "Please upload day 3 file in .edf format"),br(),
               fileInput("fileEdf3", "Please upload day 4 file in .edf format"),br(),
               fileInput("fileEdf4", "Please upload day 5 file in .edf format"),br(),
               fileInput("fileEdf5", "Please upload day 6 file in .edf format"),br(),
               fileInput("fileEdf6", "Please upload day 7 file in .edf format"),
               actionButton('submit',   strong(em('Submit')), style="color: #000000; background-color: #CFB53B; border-color: #CFB53B")
             )
             
                 ), 
             mainPanel(
               h4(strong("Instructions:")),
               p("1. Please upload EEG recordings of mouse in order using the buttons on the left.", strong("If you keep getting disconnected from the server try waiting for one file to fully upload before uploading the next one.")),
               p("2. Click on the submit button and wait for the results to appear on page.You can upload as many files as you want.", em("You can switch to the 'Standardized BSEEG Plot' tab without clicking the submit button.")),
               p("3. To download a sample .edf file or if you have encountered an error navigate to the 'more' tab."),
               tags$img(src='Screen Shot 2020-05-31 at 4.27.55 PM.png', height=100, width=300),
               p("Please make sure you are uploading files in the European Data Format (edf)."),
               #p("To know more about the algorithm used in the application please click here to access our paper."
               #, href ="https://www.sciencedirect.com/science/article/pii/S0022395620311444")
               p("4. The program is able to fill an incomplete recording for daytime or nighttime if at least 3 hours of recording is present for the daytime or the nighttime. Check the box on the top to try."),
               p("5. Disruptions in the recordings can be deleted when taking the PSD for every 10 minute window by checking the box at the top of the page."),
               p("6. The refresh button at the top of the page can be used to clear the file input slots without having to log out."),
               br(),
               tags$a(href ="https://www.sciencedirect.com/science/article/pii/S0022395620311444", "For more information regarding the application click here to see our publication", tags$style(HTML("a {color: #CFB53B}")))
             )
               )
      ),
      tabPanel("BSEEG Plot", title = "Standardized BSEEG Plot", value = "BSEEG Plot", 
               br(),
               p(strong(em("Please wait while the app is processing. It may take a few minutes to show the results after all files are uploaded."))),
               br(),
               p("You can move the sliders to change the range of Y axis. Make sure you move the sliders if part of the plot is missing."),
               #tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: yellow}")),
               #tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: yellow}")),
               tags$head(tags$style(HTML('.js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {
                                                  background: #CFB53B;
                                                  border-top: 1px solid #CFB53B ;
                                                  border-bottom: 1px solid #CFB53B ;}

                            /* changes the colour of the number tags */
                           .irs-from, .irs-to, .irs-single { background: #CFB53B }'
               ))
               ),
               h4(strong("BSEEG Plots for Every Day")),
               plotOutput("gridDiurnal")  %>% withSpinner(color="#CFB53B"),
               sliderInput("ylimBase", "Y axis scale:",
                           min = 0, max = 10, value = c(0,7), step = 0.5
               ), 
               #plotOutput("edfPlot") %>% withSpinner(color="#CFB53B"),
               #plotOutput("edfFullPlot"),
               #once all files received plot the final normalized BSEEG plot    
               tags$head(tags$style(HTML('.js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {
                                                  background: #CFB53B;
                                                  border-top: 1px solid #CFB53B ;
                                                  border-bottom: 1px solid #CFB53B ;}

                            /* changes the colour of the number tags */
                           .irs-from, .irs-to, .irs-single { background: #CFB53B }'
               ))
               ),
               h4(strong("Table of BSEEG Scores for Every Day")), 
               tabsetPanel(
                 id="myGrid",
                 tabPanel("day 1", DT::dataTableOutput("cool0") %>% withSpinner(color="#CFB53B"), downloadButton("downloadData", "Download as excel")
                 ),
                 tabPanel("day 2", DT::dataTableOutput("cool1") %>% withSpinner(color="#CFB53B")
                 ),
                 tabPanel("day 3", DT::dataTableOutput("cool2") %>% withSpinner(color="#CFB53B")
                 ),
                 tabPanel("day 4", DT::dataTableOutput("cool3") %>% withSpinner(color="#CFB53B")
                 ),
                 tabPanel("day 5", DT::dataTableOutput("cool4") %>% withSpinner(color="#CFB53B")
                 ),
                 tabPanel("day 6", DT::dataTableOutput("cool5") %>% withSpinner(color="#CFB53B")
                 ),
                 tabPanel("day 7", DT::dataTableOutput("cool6") %>% withSpinner(color="#CFB53B")
                 )
                 
               ),
               h4(strong("Standardized BSEEG Plot")),
               plotOutput("normalizedPlot") %>% withSpinner(color="#CFB53B"),
               sliderInput("ylim", "Y axis scale:",
                           min = -1, max = 3, step = 0.5, value = c(-1, 2)),
               checkboxInput("EEG1", "EEG1", value = TRUE),
               checkboxInput("EEG2", "EEG2", value = TRUE),
               #selectInput("Ticker", "select x-axis scaling", choices=c("12-hours", "Days"), selected="12-hours"),
               h4(strong("Table of Standardized BSEEG Scores for Every 12 Hours")),
               tableOutput("normalizedTable")  %>% withSpinner(color="#CFB53B")
      ),
      
      navbarMenu("more",
                 tabPanel("Sample .edf Files",
                          #p("Please wait while the app is processing. It may take as long as 10 minutes to show the results after all files are uploaded."),
                          strong('Click on the links below to download sample edf files.'),
                          br(),
                          tags$a(href="https://drive.google.com/file/d/184F0WqwsIqifVmFc-6w36KSO5C7B46La/view?usp=sharing", strong("Day 1 .edf file.")),
                          br(),
                          tags$a(href="https://drive.google.com/file/d/1eSALPgDlhA0-H_pgo-bvaPupDWAdg92s/view?usp=sharing", strong("Day 1 .edf file.")),
                          br(),
                          tags$a(href="https://drive.google.com/file/d/1gm_yNfPDRVW-VauN6tIPc4w6i5ju0DPp/view?usp=sharing", strong("Day 3 LPS injected .edf file.")),
                          br(),
                          tags$a(href="https://drive.google.com/file/d/165jq8QG3WAmsHsXzzSA-fETiOFdigKxL/view?usp=sharing", strong("Day 4 .edf file.")),
                          br(),
                          tags$a(href="https://drive.google.com/file/d/1B31Z9RxVp3ArWvrq2RsnHseVYIGyYCh-/view?usp=sharing", strong("Day 5 .edf file.")),
                          br(),
                          tags$a(href="https://drive.google.com/file/d/1mIMzqlkiHTxAYNja6qnfiSLYEPnMb15q/view?usp=sharing", strong("Day 6 .edf file.")),
                          br(),
                          h6(strong('EEG Recordings Less Than 24 Hours')),
                          tags$a(href="https://drive.google.com/file/d/1oLWgiXKuFQfejw8tgPYhriyXUrmZCEQm/view?usp=sharing", strong("21 Hour 18 minute EEG Recording File")),
                          br(),
                          tags$a(href="https://drive.google.com/file/d/1WLqY4IclK-MZ61HUnQc6rGHlw6R46U-I/view?usp=sharing", strong("Less Than 24 Hr EEG Recording File")),
                          br(),
                          #plotOutput("gridDiurnal")  %>% withSpinner(color="#CFB53B")
                 ),
                 #tabPanel("Excel Input", sidebarPanel(
                 #fileInput("file", "Please upload Day 0 file in .csv format"),
                 #fileInput("file1", "Please upload Day 1 file in .csv format"),
                 #fileInput("file2", "Please upload Day 2 file in .csv format"),
                 #fileInput("file3", "Please upload Day 3 file in .csv format"),
                 #fileInput("file4", "Please upload Day 4 file in .csv format"),
                 #fileInput("file5", "Please upload Day 5 file in .csv format"),
                 #fileInput("file6", "Please upload Day 6 file in .csv format")
                 #),
                 #mainPanel(plotOutput("fullPlot"))), 
                 
                 tabPanel("Contact Us", 
                          p("1. I keep getting disconnected from the server while uploading the files. What should I do?"),
                          p("Try uploading files one by one if the problem continues, contact us."),
                          p("2. I get the error \"Unsupported filetype\" once I submit the files."),
                          p("The program is currently only able to process .edf files therefore other file formats are not accepted."),
                          p("For other questions or errors regarding the app please contact"), tags$a("annicenajafi@tamu.edu", href =("mailto:annice-najafi@tamu.edu?Subject=Questions%20Regarding The BSEEG app"))
                          # ,tags$a(href="https://github.com/AnniceNajafi/DelDetector", strong("Click here to view app code on Github.")))
                 )
      )
    )))
###############################################END of USERINTERFACE##################################################
#'yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
#'yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
#'yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
#'yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
#'yyyyyyyyyyyyyyyyyyyssssssyyyyyyyyyyyyyyyyyyyyyyyyy
#'yyyyyyyyyyyso+//::::::::::::://+osyyyyyyyyyyyyyyyy
#'yyyyyyys+/:::::::::::::::::::::::::/+syyyyyyyyyyyy
#'yyyyy+:::::::::::::::::::::::::::::::::+syyyyyyyyy
#'yyy+:::::::::::::::::::::::::::::::::::::oyyyyyyyy
#'ys/::::::::::::::::::::/+::::::::o/:::::::syyyyyyy
#'y:::::::::::::::::::/osyo::::::::oy+::::::ysoyyyyy
#'+:::::::::::::::::+o/::sy+::::::+yyyyo++++::::+yyy
#'so::::::::::::::oyo+/:::+sysoooyyyyyyys/::::::::sy
#'y:::::::::::::oyyyyyyy/::::////:yyyyo/:::::::::::s
#'yoy:::::::::oyyyyyyys+::::::::::syy/:::::::::::::+
#'yyy:/:::::oyyyyyys+:::::::::::::yyossyyyyso+:::::/
#'yyyyo:::+yyyyyyy/::::::::::::::sys+oosyyyyyyy/:::o
#'yyyyo:/syyyyyyy/:::::::::::::/syyy+:::/yyyyyyy::+y
#'yyyyyoyyyyyyyyy/:::::::::::/oyyyyyys+:::/oyyyyosyy
#'yyyyyyyyyyyyyyys+:::::::/oyyyyyyyyyyyyyyyyyyyyyyyy
#'yyyyyyyyyyyyyyyyyyssssyyyyyyyyyyyyyyyyyyyyyyyyyyyy
#'yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
##############################################SERVER SIDE#############################################################
#'Function server receives an input from the user, processes and shows results
#'@param input the edf files received from the user
#'@param output the plots shown on the window
server <- function(input, output, session){
  
  
  result_auth <- secure_server(check_credentials = check_credentials(credentials))
  values <- reactiveValues(
    upload_state = NULL
  )
  values1 <- reactiveValues(
    upload_state = NULL
  )
  values2 <- reactiveValues(
    upload_state = NULL
  )
  values3 <- reactiveValues(
    upload_state = NULL
  )
  values4 <- reactiveValues(
    upload_state = NULL
  )
  values5 <- reactiveValues(
    upload_state = NULL
  )
  values6 <- reactiveValues(
    upload_state = NULL
  )
  
  output$res_auth <- renderPrint({
    reactiveValuesToList(result_auth)
  })
  observeEvent(input$fileEdf, {
    values$upload_state <- 'uploaded'
  })
  observeEvent(input$fileEdf1, {
    values1$upload_state <- 'uploaded'
  })
  observeEvent(input$fileEdf2, {
    values2$upload_state <- 'uploaded'
  })
  observeEvent(input$fileEdf3, {
    values3$upload_state <- 'uploaded'
  })
  observeEvent(input$fileEdf4, {
    values4$upload_state <- 'uploaded'
  })
  observeEvent(input$fileEdf5, {
    values5$upload_state <- 'uploaded'
  })
  observeEvent(input$fileEdf6, {
    values6$upload_state <- 'uploaded'
  })
  
  observeEvent(input$refresh1, {
    values$upload_state <- 'reset'
    values1$upload_state <- 'reset'
    values2$upload_state <- 'reset'
    values3$upload_state <- 'reset'
    values4$upload_state <- 'reset'
    values5$upload_state <- 'reset'
    values6$upload_state <- 'reset'
    #reset('fileEdf')
    #reset('fileEdf1')
    reset("filefile1")
    #attr(input$fileEdf1, "readonly") <- FALSE
    #input$fileEDF <- NULL
  })
  
  
  
  #The submit button is made to be reactive so once the user clicks on it they are automatically directed to the results page
  observeEvent(input$submit, {
    updateTabsetPanel(session, "del",
                      selected = "BSEEG Plot")
  })
  
  
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste("eeg_data_downloaded", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(normBack()[[12]], file)
    }
    
    #filename = "downloaded_data.csv",
  )
  datasetInput <- reactive({
    switch(normBack()[[12]],
           "time" = time,
           "EEG1" = EEG1,
           "EEG2" = EEG2)
  })
  #*********These plots are shown as tabs of a data table
  #TO change ====>Table of BSEEG Scores for Every Day edit the lines below this comment
  #optimized to received the calculated results from the function.
  
  output$cool0<- DT::renderDataTable({
    
    coolcool0 <- normBack()[[3]]
    #return if no file is uploaded./ No error message
    if(is.null(coolcool0$EEG1)){
      return()
    }
    colnames(coolcool0)<- c('Time(Hr)', 'EEG1', 'EEG2')
    coolcool0
  })
  output$cool1<- DT::renderDataTable({
    
    coolcool1 <- normBack()[[4]]
    if(is.null(coolcool1$EEG1)){
      return()
    }
    colnames(coolcool1)<- c('Time(Hr)', 'EEG1', 'EEG2')
    coolcool1
  })
  output$cool2<- DT::renderDataTable({
    
    coolcool2 <- normBack()[[5]]
    if(is.null(coolcool2$EEG1)){
      return()
    }
    colnames(coolcool2)<- c('Time(Hr)', 'EEG1', 'EEG2')
    coolcool2
  })
  output$cool3<- DT::renderDataTable({
    
    coolcool3 <- normBack()[[6]]
    if(is.null(coolcool3$EEG1)){
      return()
    }
    colnames(coolcool3)<- c('Time(Hr)', 'EEG1', 'EEG2')
    coolcool3
  })
  output$cool4<- DT::renderDataTable({
    
    coolcool4 <- normBack()[[7]]
    if(is.null(coolcool4$EEG1)){
      return()
    }
    colnames(coolcool4)<- c('Time(Hr)', 'EEG1', 'EEG2')
    coolcool4
  })
  output$cool5<- DT::renderDataTable({
    
    coolcool5 <- normBack()[[8]]
    if(is.null(coolcool5$EEG1)){
      return()
    }
    colnames(coolcool5)<- c('Time(Hr)', 'EEG1', 'EEG2')
    coolcool5
  })
  output$cool6<- DT::renderDataTable({
    
    coolcool6 <- normBack()[[9]]
    if(is.null(coolcool6$EEG1)){
      return()
    }
    colnames(coolcool6)<- c('Time(Hr)', 'EEG1', 'EEG2')
    coolcool6
  })
  #This plot shows the standardized BSEEG score over time
  output$normalizedTable<- renderTable({
    midt<- normBack()[[1]]
    #return in case no file is uploaded/ avoid error message
    if(is.null(midt)){
      return()
    }
    colnames(midt)<- c('12 hour period', 'EEG1', 'EEG2')
    midt
  })
  #This is the algorithm for the 7 plots shown in the grid. 
  output$gridTable <- renderTable({
    coolcool0 <- normBack()[[3]]
    coolcool1 <- normBack()[[4]]
    coolcool2 <- normBack()[[5]]
    coolcool3 <- normBack()[[6]]
    coolcool4 <- normBack()[[7]]
    coolcool5 <- normBack()[[8]]
    coolcool6 <- normBack()[[9]]
    #Define variables to store how many files are uploaded
    count1=0;
    count2=0;
    count3=0;
    count4=0;
    count5=0;
    count6=0;
    if(is.null(coolcool0$EEG1)){
      return()
    }
    if(is.null(coolcool1$EEG1)){
      count1=1
    }
    if(is.null(coolcool2$EEG1)){
      count2=1
    }
    if(is.null(coolcool3$EEG1)){
      count3=1
    }
    if(is.null(coolcool4$EEG1)){
      conut4=1
    }
    if(is.null(coolcool5$EEG1)){
      count5=1
    }
    if(is.null(coolcool6$EEG1)){
      count6=1
    }
    if(count1==1){
      colnames(coolcool0)<- c('Time(Hr)', 'EEG1', 'EEG2')
      grid.table(coolcool0, nrow = 1)
    }
    else if(count2==1){
      colnames(coolcool0)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool1)<- c('Time(Hr)', 'EEG1', 'EEG2')
      grid.table(coolcool0, coolcool1, nrow = 1)
    }
    else if(count3==1){
      colnames(coolcool0)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool1)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool2)<- c('Time(Hr)', 'EEG1', 'EEG2')
      grid.table(coolcool0, coolcool1, coolcool2, nrow = 1)
    }
    else if(count4==1){
      colnames(coolcool0)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool1)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool2)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool3)<- c('Time(Hr)', 'EEG1', 'EEG2')
      grid.table(coolcool0, coolcool1, coolcool2, coolcool3, nrow = 1)
    }
    else if(count5==1){
      colnames(coolcool0)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool1)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool2)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool3)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool4)<- c('Time(Hr)', 'EEG1', 'EEG2')
      grid.table(coolcool0, coolcool1, coolcool2, coolcool3, coolcool4, nrow = 1)
    }
    else if(count6==1){
      colnames(coolcool0)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool1)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool2)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool3)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool4)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool5)<- c('Time(Hr)', 'EEG1', 'EEG2')
      grid.table(coolcool0, coolcool1, coolcool2, coolcool3, coolcool4, coolcool5, nrow = 1)
    }
    else{
      colnames(coolcool0)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool1)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool2)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool3)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool4)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool5)<- c('Time(Hr)', 'EEG1', 'EEG2')
      colnames(coolcool6)<- c('Time(Hr)', 'EEG1', 'EEG2')
      grid.table(coolcool0, coolcool1, coolcool2,coolcool3, coolcool4, coolcool5, coolcool6, nrow = 2)
    }
  })
  
  output$normalizedPlot <- renderPlot({
    df <- normBack()[[1]]
    if(is.null(df)){
      return()
    }
    if(input$EEG1 && input$EEG2){
      df %>% ggplot(aes(x=x, y=y1), xaxt='n') + geom_line(aes(color = 'EEG1'), size=1.5)+ geom_line(aes(y=y2, color ='EEG2'), size=1.5) +
        scale_y_continuous(name="Standardized BSEEG Score", limits=input$ylim) + ggtitle("Standardized BSEEG") + scale_x_continuous(name="Time (12 hr)",limits=c(1, normBack()[[2]]), breaks = normBack()[[11]],labels=normBack()[[10]]
                                                                                                                                    
        )  +
        theme_bw() + theme(panel.border = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           text = element_text(size = 17))
      #xis(line=1,col="red",col.ticks="red",col.axis="red", labels=c('Day1'))
    }
    else if(input$EEG1){
      df %>% ggplot(aes(x = x, y=y1)) + geom_line(color = 'red', size=1.5) +
        scale_y_continuous(name="Standardized BSEEG Score", limits=input$ylim) + ggtitle("Standardized BSEEG") + scale_x_continuous(name="Time (12 hr)", limits=c(1, normBack()[[2]]),  breaks = normBack()[[11]],labels=normBack()[[10]]
                                                                                                                                    #, labels=normBack()[[10]]
        ) +
        theme_bw() + theme(panel.border = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           text = element_text(size = 17))
      
    }
    else if(input$EEG2){
      df %>% ggplot(aes(x = x, y=y2)) + geom_line(color = 'blue', size=1.5) +
        scale_y_continuous(name="Standardized BSEEG Score", limits=input$ylim) + ggtitle("Standardized BSEEG") + scale_x_continuous(name="Time (12 hr)", limits=c(1, normBack()[[2]]),  breaks = normBack()[[11]],labels=normBack()[[10]]
                                                                                                                                    #, labels=normBack()[[10]]
        ) +
        theme_bw() + theme(panel.border = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           text = element_text(size = 17))
    }
    
  })
  #This plot in the second tab outputs the BSEEG score over time for every single day
  output$gridDiurnal<- renderPlot({
    #hold receieved files from user in dataframes
    #day 0
    # holdedf = input$fileEdf
    # #day 1
    # holdedf1 = input$fileEdf1
    # #day 2
    # holdedf2 = input$fileEdf2
    # #day 3
    # holdedf3 = input$fileEdf3
    # #day 4
    # holdedf4 = input$fileEdf4
    # #day 5
    # holdedf5 = input$fileEdf5
    # #day 6
    # holdedf6 = input$fileEdf6
    holdedf <- holdFiles()[8:13]
    #put if statements to not continue if no file was received
    #day 0
    count1=0;
    count2=0;
    count3=0;
    count4=0;
    count5=0;
    count6=0;
    if(is.null(input$fileEdf)){
      return()
    }
    #day 1
    if(holdedf[1] == 1){
      count1=1;
    }
    #day 2
    if(holdedf[2] == 1){
      count2=1;
    }
    #day 3
    if(holdedf[3] == 1){
      count3=1;
    }
    #day 4
    if(holdedf[4] == 1){
      count4=1;
    }
    #day 5
    if(holdedf[5] == 1){
      count5=1;
    }
    #day 6
    if(holdedf[6] == 1){
      count6=1;
    }
    #use import_raw to import recordings for day 0
    #dataday0 <- import_raw(holdedf$datapath)
    #use the written function slReadEDF to find the BSEEG score
    coolTable0 <- normBack()[[3]]
    print(coolTable0)
    #rm(dataday0)
    #store time in a new vector 
    #lenlen <- length(coolTable0$time)
    
    myTime <- c(1:24)
    #plot for day 0, bind values for two electrodes
    x <- myTime
    y10 <- coolTable0$EEG1
    y20 <- coolTable0$EEG2
    df0 <- tbl_df(data.frame(x, y10, y20))
    myPlot0 <- df0 %>% ggplot(aes(x = x, y=y10)) + geom_line(aes(color = 'EEG1'), size=1.5)+ geom_line(aes(y=y20, color ='EEG2'), size=1.5)+ scale_x_continuous(name="Time (Hr)") +
      scale_y_continuous(name="BSEEG Score", limits=input$ylimBase) + ggtitle("BSEEG - Baseline Day") + theme_bw() + theme(panel.border = element_blank(),
                                                                                                                           panel.grid.major = element_blank(),
                                                                                                                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                           text = element_text(size = 17))
    #+ expand_limits(x = c(1:23))
    #+ coord_cartesian(xlim = c(1:23))
    coolTable1 <- c(NA)
    coolTable1$EEG1 <- c(NA)
    coolTable1$EEG2 <- c(NA)
    if(count1==0){
      #use import_raw to import recordings for day 1
      #dataday1 <- import_raw(holdedf1$datapath)
      #use the written function slReadEDF to find the BSEEG score
      coolTable1 <- normBack()[[4]]
      #rm(dataday1)
    }
    
    #plot for day 2, bind values for two electrodes
    #lenlen <- length(coolTable1$time)
    
    
    #plot for day 0, bind values for two electrodes
    x <- myTime
    
    #plot for day 1, bind values for two electrodes
    y11 <- coolTable1$EEG1
    y21 <- coolTable1$EEG2
    df1 <- tbl_df(data.frame(x, y11, y21))
    myPlot1 <- df1 %>% ggplot(aes(x = x, y=y11)) + geom_line(aes(color = 'EEG1'), size=1.5)+ geom_line(aes(y=y21, color ='EEG2'), size=1.5)+ scale_x_continuous(name="Time (Hr)") +
      scale_y_continuous(name="BSEEG Score", limits=input$ylimBase) + ggtitle("BSEEG - Day 2")+ theme_bw() + theme(panel.border = element_blank(),
                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                   text = element_text(size = 17))
    coolTable2 <- c(NA)
    coolTable2$EEG1 <- c(NA)
    coolTable2$EEG2 <- c(NA)
    if(count2==0){
      #use import_raw to import recordings for day 2
      #dataday2 <- import_raw(holdedf2$datapath)
      #use the written function slReadEDF to find the BSEEG score
      coolTable2 <- normBack()[[5]]
      #rm(dataday2)
    }
    #plot for day 2, bind values for two electrodes
    #lenlen <- length(coolTable2$time
    
    #myTime <- c(1:23)
    #plot for day 0, bind values for two electrodes
    x <- myTime
    y12 <- coolTable2$EEG1
    y22 <- coolTable2$EEG2
    df2 <- tbl_df(data.frame(x, y12, y22))
    myPlot2 <- df2 %>% ggplot(aes(x = x, y=y12)) + geom_line(aes(color = 'EEG1'), size=1.5)+ geom_line(aes(y=y22, color ='EEG2'), size=1.5)+ scale_x_continuous(name="Time (Hr)") +
      scale_y_continuous(name="BSEEG Score", limits=input$ylimBase) + ggtitle("BSEEG - Day 3")+ theme_bw() + theme(panel.border = element_blank(),
                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                   text = element_text(size = 17))
    coolTable3 <- c(NA)
    coolTable3$EEG1 <- c(NA)
    coolTable3$EEG2 <- c(NA)
    if(count3==0){
      #use import_raw to import recordings for day 3
      #dataday3 <-import_raw(holdedf3$datapath)
      #use the written function slReadEDF to find the BSEEG score
      coolTable3 <- normBack()[[6]]
      #rm(dataday3)
    }
    
    #lenlen <- length(coolTable3$time)
    #myTime <- c(1:lenlen)
    #x <- myTime
    #plot for day 3, bind values for two electrodes
    y13 <- coolTable3$EEG1
    y23 <- coolTable3$EEG2
    df3 <- tbl_df(data.frame(x, y13, y23))
    myPlot3 <- df3 %>% ggplot(aes(x = x, y=y13)) + geom_line(aes(color = 'EEG1'), size=1.5)+ geom_line(aes(y=y23, color ='EEG2'), size=1.5)+ scale_x_continuous(name="Time (Hr)") +
      scale_y_continuous(name="BSEEG Score", limits=input$ylimBase) + ggtitle("BSEEG - Day 4")+ theme_bw() + theme(panel.border = element_blank(),
                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                   text = element_text(size = 17))
    coolTable4 <- c(NA)
    coolTable4$EEG1 <- c(NA)
    coolTable4$EEG2 <- c(NA)
    if(count4==0){
      #use import_raw to import recordings for day 1
      #dataday4 <- import_raw(holdedf4$datapath)
      #use the written function slReadEDF to find the BSEEG score
      coolTable4 <- normBack()[[7]]
      #rm(dataday4)
    }
    
    #lenlen <- length(coolTable4$time)
    #myTime <- c(1:lenlen)
    #x <- myTime
    #plot for day 4, bind values for two electrodes
    y14 <- coolTable4$EEG1
    y24 <- coolTable4$EEG2
    df4 <- tbl_df(data.frame(x, y14, y24))
    myPlot4 <- df4 %>% ggplot(aes(x = x, y=y14)) + geom_line(aes(color = 'EEG1'), size=1.5)+ geom_line(aes(y=y24, color ='EEG2'), size=1.5)+ scale_x_continuous(name="Time (Hr)") +
      scale_y_continuous(name="BSEEG Score", limits=input$ylimBase) + ggtitle("BSEEG - Day 5")+ theme_bw() + theme(panel.border = element_blank(),
                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                   text = element_text(size = 17))
    coolTable5 <- c(NA)
    coolTable5$EEG1 <- c(NA)
    coolTable5$EEG2 <- c(NA)
    if(count5==0){
      #use import_raw to import recordings for day 5
      #dataday5 <- import_raw(holdedf5$datapath)
      #use the written function slReadEDF to find the BSEEG score
      coolTable5 <- normBack()[[8]]
      #rm(dataday5)
    }
    #lenlen <- length(coolTable5$time)
    #myTime <- c(1:lenlen)
    #x <- myTime
    #plot for day 5, bind values for two electrodes
    y15 <- coolTable5$EEG1
    y25 <- coolTable5$EEG2
    df5 <- tbl_df(data.frame(x, y15, y25))
    myPlot5 <- df5 %>% ggplot(aes(x = x, y=y15)) + geom_line(aes(color = 'EEG1'), size=1.5)+ geom_line(aes(y=y25, color ='EEG2'), size=1.5)+ scale_x_continuous(name="Time (Hr)") +
      scale_y_continuous(name="BSEEG Score", limits=input$ylimBase) + ggtitle("BSEEG - Day 6")+ theme_bw() + theme(panel.border = element_blank(),
                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                   text = element_text(size = 17))
    coolTable6 <- c(NA)
    coolTable6$EEG1 <- c(NA)
    coolTable6$EEG2 <- c(NA)
    if(count6==0){
      #use import_raw to import recordings for day 6
      #dataday6 <- import_raw(holdedf6$datapath)
      #use the written function slReadEDF to find the BSEEG score
      coolTable6 <- normBack()[[9]]
      #rm(dataday6)
    }
    #lenlen <- length(coolTable6$time)
    #myTime <- c(1:lenlen)
    #x <- myTime
    #plot for day 6, bind values for two electrodes
    y16 <- coolTable6$EEG1
    y26 <- coolTable6$EEG2
    df6 <- tbl_df(data.frame(x, y16, y26))
    myPlot6 <- df6 %>% ggplot(aes(x = x, y=y16)) + geom_line(aes(color = 'EEG1'), size=1.5)+ geom_line(aes(y=y26, color ='EEG2'), size=1.5)+ scale_x_continuous(name="Time (Hr)") +
      scale_y_continuous(name="BSEEG Score", limits=input$ylimBase) + ggtitle("BSEEG - Day 7")+ theme_bw() + theme(panel.border = element_blank(),
                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                   text = element_text(size = 17))
    #arrange plots for every single day of the week in a grid
    if(count1==1){
      grid.arrange(myPlot0, nrow = 1)
    }
    else if(count2==1){
      grid.arrange(myPlot0, myPlot1, nrow = 1)
    }
    else if(count3==1){
      grid.arrange(myPlot0, myPlot1, myPlot2, nrow = 1)
    }
    else if(count4==1){
      grid.arrange(myPlot0, myPlot1, myPlot2, myPlot3, nrow = 1)
    }
    else if(count5==1){
      grid.arrange(myPlot0, myPlot1, myPlot2, myPlot3, myPlot4, nrow = 1)
    }
    else if(count6==1){
      grid.arrange(myPlot0, myPlot1, myPlot2, myPlot3, myPlot4, myPlot5, nrow = 1)
    }
    else{
      grid.arrange(myPlot0, myPlot1, myPlot2, myPlot3, myPlot4, myPlot5, myPlot6, nrow = 2)
    }
  })
  # output$coolPlot<- renderPlot({
  #   hold = input$file
  #   if(is.null(hold)){
  #     return()
  #   }
  #  
  #   data <- read.csv2(hold$datapath, sep=",", header = FALSE)
  #   data <- data[-c(1:5),]
  #   colnames(data)<- as.matrix(data[1,])
  #   data = data[-1, ]
  #   data[]<- lapply(data, function(x) type.convert(as.character(x)))
  #   data <- data[-1,]
  #   data["Ratio1"] <- NA
  #   data["Ratio2"] <- NA
  #   names(data)[1] <- "Date"
  #   names(data)[2] <- "Time"
  #   names(data)[3] <- "TimeStamp"
  #   names(data)[4] <- "TimefromStart"
  #   names(data)[5] <- "EEG1_gs_10hz"
  #   names(data)[6] <- "EEG1_gs_3hz"
  #   names(data)[7] <- "EEG2_gs_10hz"
  #   names(data)[8] <- "EEG2_gs_3hz"
  #   
  #   data <- transform(data, EEG1_gs_10hz = as.numeric(as.character(EEG1_gs_10hz)), 
  #                     EEG1_gs_3hz = as.numeric(as.character(EEG1_gs_3hz)), EEG2_gs_10hz = as.numeric(as.character(EEG2_gs_10hz)), 
  #                     EEG2_gs_3hz = as.numeric(as.character(EEG2_gs_3hz)))
  #   data <- mutate(data, Ratio1 = EEG1_gs_3hz/EEG1_gs_10hz)
  #   data <- mutate(data, Ratio2 = EEG2_gs_3hz/EEG2_gs_10hz)
  #   #Create new table
  #   coolTable <- data.frame(time=c(-18:-1, 0, 1, 2, 3, 4, 5))
  #   newTable <-  data.frame(time=c(-18:-1, 0, 1, 2, 3, 4, 5))
  #   coolTable["EEG1"]<- NA
  #   coolTable["EEG2"]<- NA
  #   coolAverage1 <- c(NA)
  #   coolAverage2 <- c(NA)
  #   #data$Time <- as.numeric(strptime(data$Time, format="%H:%M:%S") - as.POSIXct(format(Sys.Date())), units="secs")
  #   i=1
  #   for(j in 0:23){
  #     
  #     holdVal1 <- data$Ratio1[c((360*(j)+1):(360*(j+1)+1))]
  #     holdVal2 <- data$Ratio2[c((360*(j)+1):(360*(j+1)+1))]
  #     coolAverage1[i]<- mean(holdVal1)
  #     coolAverage2[i]<- mean(holdVal2)
  #     i=i+1
  #   }
  #   
  #   coolTable["EEG1"]<- coolAverage1
  #   coolTable["EEG2"]<- coolAverage2
  #   x <- coolTable$time
  #   y1 <- coolTable$EEG1
  #   y2 <- coolTable$EEG2
  #   df <- tbl_df(data.frame(x, y1, y2))
  #   coolPlot <- df %>% ggplot(aes(x = x, y=y2)) + geom_line(color = 'red', size=2)+geom_line(y=y1, color ='blue', size=2) + scale_x_continuous(name="Time") +
  #     scale_y_continuous(name="BSEEG Score", limits=c(0, 20)) + ggtitle("BSEEG")
  #   coolPlot
  # })
  
  # output$fullPlot<- renderPlot({
  #   
  #   holdday0 = input$file
  #   holdday1 = input$file1
  #   holdday2 = input$file2
  #   holdday3 = input$file3
  #   holdday4 = input$file4
  #   holdday5 = input$file5
  #   holdday6 = input$file6
  #   
  #   if(is.null(holdday0)){
  #     return()
  #   }
  #   if(is.null(holdday1)){
  #     return()
  #   }
  #   if(is.null(holdday2)){
  #     return()
  #   }
  #   if(is.null(holdday3)){
  #     return()
  #   }
  #   if(is.null(holdday4)){
  #     return()
  #   }
  #   if(is.null(holdday5)){
  #     return()
  #   }
  #   if(is.null(holdday6)){
  #     return()
  #   }
  #   
  #   dataday0 <- read.csv2(holdday0$datapath, sep=",", header = FALSE)
  #   dataday1 <- read.csv2(holdday1$datapath, sep=",", header = FALSE)
  #   dataday2 <- read.csv2(holdday2$datapath, sep=",", header = FALSE)
  #   dataday3 <- read.csv2(holdday3$datapath, sep=",", header = FALSE)
  #   dataday4 <- read.csv2(holdday4$datapath, sep=",", header = FALSE)
  #   dataday5 <- read.csv2(holdday5$datapath, sep=",", header = FALSE)
  #   dataday6 <- read.csv2(holdday6$datapath, sep=",", header = FALSE)
  #   
  #   coolTable0 <- slRead(dataday0)
  #   coolTable1 <- slRead(dataday1)
  #   coolTable2 <- slRead(dataday2)
  #   coolTable3 <- slRead(dataday3)
  #   coolTable4 <- slRead(dataday4)
  #   coolTable5 <- slRead(dataday5)
  #   coolTable6 <- slRead(dataday6)
  #   
  #   coolTable0$time <- c(-18:5)
  #   coolTable1$time <- c(6:29)
  #   coolTable2$time <- c(30:53)
  #   coolTable3$time <- c(54:77)
  #   coolTable4$time <- c(78:101)
  #   coolTable5$time <- c(102:125)
  #   coolTable6$time <- c(126:149)
  #   
  #   
  #   
  #   
  #   fullTable <- rbind(coolTable0, coolTable1, coolTable2, coolTable3, coolTable4, coolTable5, coolTable6)
  #   
  #   
  #   x <- fullTable$time
  #   y1 <- fullTable$EEG1
  #   y2 <- fullTable$EEG2
  #   df <- tbl_df(data.frame(x, y1, y2))
  #   fullPlot <- df %>% ggplot(aes(x = x, y=y2)) + geom_line(color = 'red', size=2)+geom_line(y=y1, color ='blue', size=2) + scale_x_continuous(name="Time") +
  #     scale_y_continuous(name="EEG Ratio", limits=c(0, 30)) + ggtitle("Standardized BSEEG Combined Average Score")
  #   print(fullTable)
  #   fullPlot
  #   
  # })
  
  infile <- eventReactive({}, {
    data<-as.data.frame(normalizedPlot)       
    updateSelectInput(session, inputId = "ylim", choices = names(data), selected = names(data)[2])
    return(data)
  }, ignoreNULL = FALSE)
  
  infile <- eventReactive({}, {
    data<-as.data.frame(gridDiurnal[1])       
    updateSelectInput(session, inputId = "filldis", choices = names(data), selected = names(data)[2])
    return(data)
  }, ignoreNULL = FALSE)
  
  infile <- eventReactive({}, {
    data<-as.data.frame(gridDiurnal[1])       
    updateSelectInput(session, inputId = "filldisrupt", choices = names(data), selected = names(data)[2])
    return(data)
  }, ignoreNULL = FALSE)
  
  infile <- eventReactive({}, {
    data<-as.data.frame(edfPlot)    
    updateSelectInput(session, inputId = "ylimBase", choices = names(data), selected = names(data)[2])
    return(data)
  }, ignoreNULL = FALSE)
  
  
  #holdFiles <- eventReactive(input$submit,{
  holdFiles <- function(){
    
    if (is.null(values$upload_state)) {
      holdedf = NULL
    } else if (values$upload_state == 'uploaded') {
      holdedf = input$fileEdf
    } else if (values$upload_state == 'reset') {
      holdedf = NULL
    }
    
    if (is.null(values1$upload_state)) {
      holdedf1 = NULL
    } else if (values1$upload_state == 'uploaded') {
      holdedf1 = input$fileEdf1
    } else if (values1$upload_state == 'reset') {
      holdedf1 = NULL
    }
    
    if (is.null(values2$upload_state)) {
      holdedf2 = NULL
    } else if (values$upload_state == 'uploaded') {
      holdedf2 = input$fileEdf2
    } else if (values$upload_state == 'reset') {
      holdedf2 = NULL
    }
    
    if (is.null(values3$upload_state)) {
      holdedf3 = NULL
    } else if (values3$upload_state == 'uploaded') {
      holdedf3 = input$fileEdf3
    } else if (values3$upload_state == 'reset') {
      holdedf3 = NULL
    }
    
    if (is.null(values4$upload_state)) {
      holdedf4 = NULL
    } else if (values4$upload_state == 'uploaded') {
      holdedf4 = input$fileEdf4
    } else if (values4$upload_state == 'reset') {
      holdedf4 = NULL
    }
    
    if (is.null(values5$upload_state)) {
      holdedf5 = NULL
    } else if (values5$upload_state == 'uploaded') {
      holdedf5 = input$fileEdf5
    } else if (values5$upload_state == 'reset') {
      holdedf5 = NULL
    }
    
    if (is.null(values6$upload_state)) {
      holdedf6 = NULL
    } else if (values$upload_state == 'uploaded') {
      holdedf6 = input$fileEdf6
    } else if (values$upload_state == 'reset') {
      holdedf6 = NULL
    }
    
    
    
    # #hold receieved files from user in dataframes
    # #day 0
    # holdedf = input$fileEdf
    # #day 1
    # holdedf1 = input$fileEdf1
    # #day 2
    # holdedf2 = input$fileEdf2
    # #day 3
    # holdedf3 = input$fileEdf3
    # #day 4
    # holdedf4 = input$fileEdf4
    # #day 5
    # holdedf5 = input$fileEdf5
    # #day 6
    # holdedf6 = input$fileEdf6
    #put if statements to not continue if no file was received
    #day 0
    count1=0;normalEdf1=0;
    count2=0;normalEdf2=0;
    count3=0;normalEdf3=0;
    count4=0;normalEdf4=0;
    count5=0;normalEdf5=0;
    count6=0;normalEdf6=0;
    if(is.null(holdedf)){
      return()
    }
    #day 1
    if(is.null(holdedf1)){
      count1=1;
    }
    #day 2
    if(is.null(holdedf2)){
      count2=1;
    }
    #day 3
    if(is.null(holdedf3)){
      count3=1;
    }
    #day 4
    if(is.null(holdedf4)){
      count4=1;
    }
    #day 5
    if(is.null(holdedf5)){
      count5=1;
    }
    #day 6
    if(is.null(holdedf6)){
      count6=1;
    }
    #day 0
    holdedf = input$fileEdf
    
    #use import_raw to import recordings for day 0
    #dataday0 <- import_raw(holdedf$datapath)
    
    dataday0 <- tryCatch(
      {
        import_raw(holdedf$datapath)
      }, 
      error= function(cond){
        edf::read.edf(holdedf$datapath)
      }
    )
    #day 1
    if(!is.null(holdedf1)){
      holdedf1 = input$fileEdf1
      #dataday1 <- import_raw(holdedf1$datapath)
      dataday1 <- tryCatch(
        {
          import_raw(holdedf1$datapath)
        }, 
        error= function(cond){
          edf::read.edf(holdedf1$datapath)
        }
      )
    }else{
      dataday1 <- NULL
    }
    if(!is.null(holdedf2)){
      holdedf2 = input$fileEdf2
      #dataday2 <- import_raw(holdedf2$datapath)
      dataday2 <- tryCatch(
        {
          import_raw(holdedf2$datapath)
        }, 
        error= function(cond){
          edf::read.edf(holdedf2$datapath)
        }
      )
    }else{
      dataday2 <- NULL
    }
    if(!is.null(holdedf3)){
      holdedf3 = input$fileEdf3
      #dataday3 <- import_raw(holdedf3$datapath)
      dataday3 <- tryCatch(
        {
          import_raw(holdedf3$datapath)
        }, 
        error= function(cond){
          edf::read.edf(holdedf3$datapath)
        }
      )
    }else{
      dataday3 <- NULL
    }
    if(!is.null(holdedf4)){
      holdedf4 = input$fileEdf4
      #dataday4 <- import_raw(holdedf4$datapath)
      dataday4 <- tryCatch(
        {
          import_raw(holdedf4$datapath)
        }, 
        error= function(cond){
          edf::read.edf(holdedf4$datapath)
        }
      )
    }else{
      dataday4 <- NULL
    }
    if(!is.null(holdedf5)){
      holdedf5 = input$fileEdf5
      #dataday5 <- import_raw(holdedf5$datapath)
      dataday5 <- tryCatch(
        {
          import_raw(holdedf5$datapath)
        }, 
        error= function(cond){
          edf::read.edf(holdedf5$datapath)
        }
      )
    }else{
      dataday5 <- NULL
    }
    if(!is.null(holdedf6)){
      holdedf6 = input$fileEdf6
      #dataday6 <- import_raw(holdedf6$datapath)
      dataday6 <- tryCatch(
        {
          import_raw(holdedf6$datapath)
        }, 
        error= function(cond){
          edf::read.edf(holdedf6$datapath)
        }
      )
    }else{
      dataday6 <- NULL
    }
    return(list(dataday0, dataday1, dataday2, dataday3, dataday4, dataday5, dataday6,  count1, count2, count3, count4, count5, count6))
  }
  showTest<- renderTable({
    
  })
  normBack <- eventReactive(input$submit,{
    # #hold receieved files from user in dataframes
    # #day 0
    # holdedf = input$fileEdf
    holdedf <- holdFiles()[8:13]
    # #day 1
    #  holdedf1 = holdFiles()[[8]]
    # # #day 2
    #  holdedf2 = holdFiles()[[9]]
    # # #day 3
    #  holdedf3 = holdFiles()[[10]]
    # # #day 4
    #  holdedf4 = holdFiles()[[11]]
    # # #day 5
    #  holdedf5 = holdFiles()[[12]]
    # # #day 6
    #  holdedf6 = holdFiles()[[13]]
    #put if statements to not continue if no file was received
    count1=0;
    count2=0;
    count3=0;
    count4=0;
    count5=0;
    count6=0;
    if(is.null(input$fileEdf)){
      return()
    }
    #day 1
    if(holdedf[1] == 1){
      count1=1;
    }
    #day 2
    if(holdedf[2] == 1){
      count2=1;
    }
    #day 3
    if(holdedf[3] == 1){
      count3=1;
    }
    #day 4
    if(holdedf[4] == 1){
      count4=1;
    }
    #day 5
    if(holdedf[5] == 1){
      count5=1;
    }
    #day 6
    if(holdedf[6] == 1){
      count6=1;
    }
    breakHold <- c(1,2)
    labelHold <- c('Day 1 - hour 12','Day1 - hour 24')
    #hold receieved files from user in dataframes
    #day 0
    ####holdedf = input$fileEdf
    #use import_raw to import recordings for day 0
    #dataday0 <- import_raw(holdedf$datapath)
    dataday0 <- holdFiles()[[1]]
    print(dataday0)
    #use the written function slReadEDF to find the BSEEG score
    coolTable0 <- as.data.table(slReadEDF(dataday0, input$filldis, input$filldisrupt)[1:3])
    
    
    srate <- dataday0$srate
    totalSeconds <- nrow(dataday0$timings)
    numHrs <- totalSeconds/(srate*3600)
    #parse the data into 24 hours and 10 minutes intervals
    numHrs <- floor(numHrs)
    numHrs <- numHrs - 2
    
    rm(dataday0)
    #add a column to coolTable0 including time
    coolTable0$time <- c(-18:5)
    #At this point we have the BSEEG scores for every single day of the week
    #We normalize the plot according to the baseline day
    #first we have to filter data for every 12 hours 
    #initialize vector for first 12 hours of day 0 EEG1
    normAvg0firstHalf1 <- c(NA)
    #initialize vector for second 12 hours of day 0 EEG1
    normAvg0secondHalf1 <- c(NA)
    normAvg0firstHalf2 <- c(NA)
    normAvg0secondHalf2 <- c(NA)
    #coolTable0 <- coolTable0 %>% drop_na()
    
    # lenT <- length(coolTable0)/2
    # normAvg0firstHalf1 <- coolTable0$EEG1[c(1:lenT)]
    # normAvg0secondHalf1 <- coolTable0$EEG1[c(lenT+1:2*lenT)]
    # normAvg0firstHalf2 <- coolTable0$EEG2[c(1:lenT)]
    # normAvg0secondHalf2 <- coolTable0$EEG2[c(lenT+1:2*lenT)]
    
    normAvg0firstHalf1 <- coolTable0$EEG1[c(1:11)]
    normAvg0secondHalf1 <- coolTable0$EEG1[c(12:22)]
    normAvg0firstHalf2 <- coolTable0$EEG2[c(1:11)]
    normAvg0secondHalf2 <- coolTable0$EEG2[c(12:22)]
    
    Baseline1 <- mean(normAvg0firstHalf1)
    holdNorm1 <- c(NA)
    holdNorm1[1]<- 0
    holdNorm1[2]<- (mean(normAvg0secondHalf1) - Baseline1)/Baseline1
    Baseline2 <- mean(normAvg0firstHalf2)
    holdNorm2 <- c(NA)
    holdNorm2[1]<- 0
    holdNorm2[2]<- (mean(normAvg0secondHalf2) - Baseline2)/Baseline2
    ####holdedf1 = input$fileEdf1
    #use import_raw to import recordings for day 1
    res <- coolTable0
    coolTable1 <- c(NA)
    coolTable1$time <- c(6, 29)
    if(count1==0){
      #dataday1 <- import_raw(holdedf1$datapath)
      dataday1 <- holdFiles()[[2]]
      # 
      # srate <- dataday1$srate
      # totalSeconds <- nrow(dataday1$timings)
      # numHrs2 <- totalSeconds/(srate*3600)
      # #parse the data into 24 hours and 10 minutes intervals
      # numHrs2 <- floor(numHrs2)
      # numHrs2 <- numHrs2 - 2
      # 
      # #use the written function slReadEDF to find the BSEEG score
      coolTable1 <- as.data.table(slReadEDF(dataday1, input$filldis, input$filldisrupt)[1:3])
      # coolTable1$time <- c(((-18+numHrs)+1):(-17+numHrs+numHrs2))
      rm(dataday1)
      breakHold <- c(1, 2, 3, 4)
      labelHold <- c('Day 1 - hour 12','Day1 - hour 24','Day2 - hour 12', 'Day2 - hour 24')
      res <- rbind(coolTable0, coolTable1)
    }
    # else{
    
    #}
    rm(holdedf)
    #add a column to coolTable1 including time
    normAvg1firstHalf1 <- c(NA)
    normAvg1secondHalf1 <- c(NA)
    normAvg1firstHalf2 <- c(NA)
    normAvg1secondHalf2 <- c(NA)
    
    # lenT <- length(coolTable1)/2
    # normAvg1firstHalf1 <- coolTable1$EEG1[c(1:lenT)]
    # normAvg1secondHalf1 <- coolTable1$EEG1[c(lenT+1:2*lenT)]
    # normAvg1firstHalf2 <- coolTable1$EEG2[c(1:lenT)]
    # normAvg1secondHalf2 <- coolTable1$EEG2[c(lenT+1:2*lenT)]
    
    normAvg1firstHalf1 <- coolTable1$EEG1[c(1:11)]
    normAvg1secondHalf1 <- coolTable1$EEG1[c(12:22)]
    normAvg1firstHalf2 <- coolTable1$EEG2[c(1:11)]
    normAvg1secondHalf2 <- coolTable1$EEG2[c(12:22)]
    
    holdNorm1[3]<- (mean(normAvg1firstHalf1) - Baseline1)/Baseline1
    holdNorm1[4]<- (mean(normAvg1secondHalf1) - Baseline1)/Baseline1
    holdNorm2[3]<- (mean(normAvg1firstHalf2) - Baseline2)/Baseline2
    holdNorm2[4]<- (mean(normAvg1secondHalf2) - Baseline2)/Baseline2
    
    coolTable2 <- c(NA)
    numHrs3 <- 1
    coolTable2$time <- c(30, 53)
    if(count2==0){
      #   #use import_raw to import recordings for day 2
      dataday2 <- holdFiles()[[3]]
      #   #use the written function slReadEDF to find the BSEEG score
      #   
      #   srate <- dataday2$srate
      #   totalSeconds <- nrow(dataday2$timings)
      #   numHrs3 <- totalSeconds/(srate*3600)
      #   #parse the data into 24 hours and 10 minutes intervals
      #   numHrs3 <- floor(numHrs3)
      #   numHrs3 <- numHrs3 - 2
      coolTable2 <- as.data.table(slReadEDF(dataday2, input$filldis, input$filldisrupt)[1:3])
      #   coolTable2$time <- c(((-18+numHrs2)+1):(-17+numHrs3+numHrs2))
      rm(dataday2)
      breakHold <- c(1, 2, 3, 4, 5, 6)
      labelHold <- c('Day 1 - hour 12','Day1 - hour 24','Day2 - hour 12', 'Day2 - hour 24', 'Day3 - hour 12', 'Day3 - hour 24')
      res <- rbind(coolTable0, coolTable1, coolTable2)
    }
    # else{
    
    #}
    #add a column to coolTable2 including time
    normAvg2firstHalf1 <- c(NA)
    normAvg2secondHalf1 <- c(NA)
    normAvg2firstHalf2 <- c(NA)
    normAvg2secondHalf2 <- c(NA)
    lenT <- length(coolTable2)/2
    # normAvg2firstHalf1 <- coolTable2$EEG1[c(1:lenT)]
    # normAvg2secondHalf1 <- coolTable2$EEG1[c(lenT+1:2*lenT)]
    # normAvg2firstHalf2 <- coolTable2$EEG2[c(1:lenT)]
    # normAvg2secondHalf2 <- coolTable2$EEG2[c(lenT+1:2*lenT)]
    normAvg2firstHalf1 <- coolTable2$EEG1[c(1:11)]
    normAvg2secondHalf1 <- coolTable2$EEG1[c(12:22)]
    normAvg2firstHalf2 <- coolTable2$EEG2[c(1:11)]
    normAvg2secondHalf2 <- coolTable2$EEG2[c(12:22)]
    
    holdNorm1[5]<- (mean(normAvg2firstHalf1) - Baseline1)/Baseline1
    holdNorm1[6]<- (mean(normAvg2secondHalf1) - Baseline1)/Baseline1
    holdNorm2[5]<- (mean(normAvg2firstHalf2)  - Baseline2)/Baseline2
    holdNorm2[6]<- (mean(normAvg2secondHalf2)  - Baseline2)/Baseline2
    
    
    coolTable3 <- c(NA)
    coolTable3$time <- c(54:77)
    if(count3==0){
      #   #use import_raw to import recordings for day 3
      dataday3 <- holdFiles()[[4]]
      #   #use the written function slReadEDF to find the BSEEG score
      #   srate <- dataday3$srate
      #   totalSeconds <- nrow(dataday3$timings)
      #   numHrs4 <- totalSeconds/(srate*3600)
      #   #parse the data into 24 hours and 10 minutes intervals
      #   numHrs4 <- floor(numHrs4)
      #   numHrs4 <- numHrs4 - 2
      coolTable3 <- as.data.table(slReadEDF(dataday3, input$filldis, input$filldisrupt)[1:3])
      #   coolTable3$time <- c(((-18+numHrs3)+1):(-17+numHrs3+numHrs4))
      rm(dataday3)
      breakHold <- c(1, 2, 3, 4, 5, 6, 7, 8)
      labelHold <- c('Day 1 - hour 12','Day1 - hour 24','Day2 - hour 12', 'Day2 - hour 24', 'Day3 - hour 12', 'Day3 - hour 24', 'Day4 - hour 12', 'Day4 - hour 24')
      res <- rbind(coolTable0, coolTable1, coolTable2, coolTable3)
    }
    # else{
    #}
    #add a column to coolTable3 including time
    normAvg3firstHalf1 <- c(NA)
    normAvg3secondHalf1 <- c(NA)
    normAvg3firstHalf2 <- c(NA)
    normAvg3secondHalf2 <- c(NA)
    lenT <- length(coolTable3)/2
    normAvg3firstHalf1 <- coolTable3$EEG1[c(1:11)]
    normAvg3secondHalf1 <- coolTable3$EEG1[c(12:22)]
    normAvg3firstHalf2 <- coolTable3$EEG2[c(1:11)]
    normAvg3secondHalf2 <- coolTable3$EEG2[c(12:22)]
    holdNorm1[7]<- (mean(normAvg3firstHalf1) - Baseline1)/Baseline1
    holdNorm1[8]<- (mean(normAvg3secondHalf1) - Baseline1)/Baseline1
    holdNorm2[7]<- (mean(normAvg3firstHalf2)  - Baseline2)/Baseline2
    holdNorm2[8]<- (mean(normAvg3secondHalf2)  - Baseline2)/Baseline2
    
    coolTable4 <- c(NA)
    coolTable4$time <- c(78:101)
    if(count4==0){
      #use import_raw to import recordings for day 4
      dataday4 <- holdFiles()[[5]]
      #use the written function slReadEDF to find the BSEEG score
      #srate <- dataday4$srate
      #totalSeconds <- nrow(dataday4$timings)
      #numHrs5 <- totalSeconds/(srate*3600)
      #parse the data into 24 hours and 10 minutes intervals
      #numHrs5 <- floor(numHrs5)
      #numHrs5 <- numHrs5 - 2
      coolTable4 <- as.data.table(slReadEDF(dataday4, input$filldis, input$filldisrupt)[1:3])
      #coolTable4$time <- c(((-18+numHrs4)+1):(-17+numHrs5+numHrs4))
      
      rm(dataday4)
      breakHold <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
      labelHold <- c('Day 1 - hour 12','Day1 - hour 24','Day2 - hour 12', 'Day2 - hour 24', 'Day3 - hour 12', 'Day3 - hour 24', 'Day4 - hour 12', 'Day4 - hour 24', 'Day5 - hour 12', 'Day5 - hour 24')
      res <- rbind(coolTable0, coolTable1, coolTable2, coolTable3, coolTable4)
    }
    #else{
    #}
    #add a column to coolTable4 including time
    
    normAvg4firstHalf1 <- c(NA)
    normAvg4secondHalf1 <- c(NA)
    normAvg4firstHalf2 <- c(NA)
    normAvg4secondHalf2 <- c(NA)
    lenT <- length(coolTable4)/2
    normAvg4firstHalf1 <- coolTable4$EEG1[c(1:11)]
    normAvg4secondHalf1 <- coolTable4$EEG1[c(12:22)]
    normAvg4firstHalf2 <- coolTable4$EEG2[c(1:11)]
    normAvg4secondHalf2 <- coolTable4$EEG2[c(12:22)]
    holdNorm1[9]<- (mean(normAvg4firstHalf1) - Baseline1)/Baseline1
    holdNorm1[10]<- (mean(normAvg4secondHalf1) - Baseline1)/Baseline1
    holdNorm2[9]<- (mean(normAvg4firstHalf2)  - Baseline2)/Baseline2
    holdNorm2[10]<- (mean(normAvg4secondHalf2)  - Baseline2)/Baseline2
    coolTable5 <- c(NA)
    coolTable5$time <- c(102:125)
    if(count5==0){
      #use import_raw to import recordings for day 5
      dataday5 <- holdFiles()[[6]]
      #use the written function slReadEDF to find the BSEEG score
      coolTable5 <- as.data.table(slReadEDF(dataday5, input$filldis, input$filldisrupt)[1:3])
      
      # srate <- dataday5$srate
      # totalSeconds <- nrow(dataday5$timings)
      # numHrs6 <- totalSeconds/(srate*3600)
      #parse the data into 24 hours and 10 minutes intervals
      # numHrs6 <- floor(numHrs6)
      # numHrs6 <- numHrs6 - 2
      # coolTable5$time <- c(((-18+numHrs5)+1):(-17+numHrs5+numHrs6))
      
      rm(dataday5)
      breakHold <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
      labelHold <- c('Day 1 - hour 12','Day1 - hour 24','Day2 - hour 12', 'Day2 - hour 24', 'Day3 - hour 12', 'Day3 - hour 24', 'Day4 - hour 12', 'Day4 - hour 24', 'Day5 - hour 12', 'Day5 - hour 24', 'Day6 - hour 12', 'Day6 - hour 24')
      res <- rbind(coolTable0, coolTable1, coolTable2, coolTable3, coolTable4, coolTable5)
    }
    #else{
    
    #}
    #add a column to coolTable5 including time
    
    normAvg5firstHalf1 <- c(NA)
    normAvg5secondHalf1 <- c(NA)
    normAvg5firstHalf2 <- c(NA)
    normAvg5secondHalf2 <- c(NA)
    lenT <- length(coolTable5)/2
    normAvg5firstHalf1 <- coolTable5$EEG1[c(1:11)]
    normAvg5secondHalf1 <- coolTable5$EEG1[c(12:22)]
    normAvg5firstHalf2 <- coolTable5$EEG2[c(1:11)]
    normAvg5secondHalf2 <- coolTable5$EEG2[c(12:22)]
    holdNorm1[11]<- (mean(normAvg5firstHalf1) - Baseline1)/Baseline1
    holdNorm1[12]<- (mean(normAvg5secondHalf1) - Baseline1)/Baseline1
    holdNorm2[11]<- (mean(normAvg5firstHalf2)  - Baseline2)/Baseline2
    holdNorm2[12]<- (mean(normAvg5secondHalf2)  - Baseline2)/Baseline2
    coolTable6 <- c(NA)
    coolTable6$time <- c(126:149)
    if(count6==0){
      #use import_raw to import recordings for day 6
      dataday6 <- holdFiles()[[7]]
      #use the written function slReadEDF to find the BSEEG score
      coolTable6 <- as.data.table(slReadEDF(dataday6, input$filldis, input$filldisrupt)[1:3])
      
      #srate <- dataday6$srate
      #totalSeconds <- nrow(dataday6$timings)
      #numHrs7 <- totalSeconds/(srate*3600)
      #parse the data into 24 hours and 10 minutes intervals
      #umHrs7 <- floor(numHrs7)
      #numHrs7 <- numHrs7 - 2
      #coolTable6$time <- c(((-18+numHrs6)+1):(-17+numHrs7+numHrs6))
      
      rm(dataday6)
      breakHold <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
      labelHold <- c('Day 1 - hour 12','Day1 - hour 24','Day2 - hour 12', 'Day2 - hour 24', 'Day3 - hour 12', 'Day3 - hour 24', 'Day4 - hour 12', 'Day4 - hour 24', 'Day5 - hour 12', 'Day5 - hour 24', 'Day6 - hour 12', 'Day6 - hour 24', 'Day7 - hour 12', 'Day 7 - hour 24')
      res <- rbind(coolTable0, coolTable1, coolTable2, coolTable3, coolTable4, coolTable5, coolTable6)
    }
    #else{
    
    #}
    #add a column to coolTable6 including time
    normAvg6firstHalf1 <- c(NA)
    normAvg6secondHalf1 <- c(NA)
    normAvg6firstHalf2 <- c(NA)
    normAvg6secondHalf2 <- c(NA)
    lenT <- length(coolTable6)/2
    normAvg6firstHalf1 <- coolTable6$EEG1[c(1:11)]
    normAvg6secondHalf1 <- coolTable6$EEG1[c(12:22)]
    normAvg6firstHalf2 <- coolTable6$EEG2[c(1:11)]
    normAvg6secondHalf2 <- coolTable6$EEG2[c(12:22)]
    holdNorm1[13]<- (mean(normAvg6firstHalf1) - Baseline1)/Baseline1
    holdNorm1[14]<- (mean(normAvg6secondHalf1) - Baseline1)/Baseline1
    holdNorm2[13]<- (mean(normAvg6firstHalf2)  - Baseline2)/Baseline2
    holdNorm2[14]<- (mean(normAvg6secondHalf2)  - Baseline2)/Baseline2
    #plot the standardized BSEEG score over time
    normX <- c(1:14)
    x <- normX
    #ratios for EEG1 and EEG2 on the same plot
    y1 <- holdNorm1*2
    y2 <- holdNorm2*2
    df <- tbl_df(data.frame(x, y1, y2))
    normalizedTable <- df
    lims <- 1
    limlim <- 14
    if(count6==1){
      limlim <- 12
      normalizedTable <- normalizedTable[c(1:12),]
    }
    if(count5==1){
      limlim <- 10
      normalizedTable <- normalizedTable[c(1:10),]
    }
    if(count4==1){
      limlim <- 8
      normalizedTable <- normalizedTable[c(1:8),]
    }
    if(count3==1){
      limlim <- 6
      normalizedTable <- normalizedTable[c(1:6),]
    }
    if(count2==1){
      limlim <- 4
      normalizedTable <- normalizedTable[c(1:4),]
    }
    if(count1==1){
      limlim <- 2
      normalizedTable <- normalizedTable[c(1:2),]
    }
    normalizedPlot<- df %>% ggplot(aes(x = x, y=y1)) + geom_line(aes(color = 'EEG1'), size=1)+ geom_line(aes(y=y2, color ='EEG2'), size=1) +
      scale_y_continuous(name="BSEEG Score", limits=input$ylim) + ggtitle("Standardized BSEEG") + scale_x_continuous(name="Time (12 hr)", limits=c(lims, limlim)) 
    #return the final plot
    # Make dataset_1
    # Make dataset_2
    return(list(normalizedTable, limlim, coolTable0, coolTable1, coolTable2, coolTable3, coolTable4, coolTable5, coolTable6, labelHold, breakHold, res))
  },ignoreNULL = FALSE)
  
}
#run the shiny app
shinyApp(ui = ui, server=server)





