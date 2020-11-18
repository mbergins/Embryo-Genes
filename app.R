library(shiny)
library(tidyverse)
library(here)
library(ggplot2)
library(reactable)

parnell_data = read_rds(here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_forShinyApp.rds'))
plot_range = c(floor(min(parnell_data$Mean-parnell_data$SE)), ceiling(max(parnell_data$Mean+parnell_data$SE)))

parnell_data_list = read_rds(here('Parnell_data_split.rds'))
parnell_data_full_list = read_rds(here('Parnell_data_full_split.rds'))

gene_list = parnell_data %>% 
    group_by(Gene) %>% 
    summarise(gene_mean_exp = mean(Mean)) %>% 
    arrange(desc(gene_mean_exp)) %>%
    pull(Gene)

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("FAS Data Browser"),
    
    sidebarLayout(
        sidebarPanel(
            # selectizeInput(inputId = "gene",
            #             label = h2("Select Gene of Interest"),
            #             choices = gene_list,
            #             options = list(maxOptions = length(gene_list))),
            selectizeInput(inputId = "gene",
                           label = h2("Select Gene of Interest"),
                           choices = gene_list,
                           multiple = F,
                           options = list(maxOptions = 50,
                                          placeholder = "Please Select a Gene",
                                          maxItems = 1,
                                          onInitialize = I('function() { this.setValue(""); }'))),            
            checkboxInput("include_PAE",
                          label = 'Include prenatal alcohol exposure (PAE) data',
                          value = F),
            checkboxGroupInput("mouse_strains", label = h3("Mouse Strains"), 
                               choices = list("C57BL6J" = "C57BL6J", "C57BL6N" = "C57BL6N"),
                               selected = c("C57BL6J","C57BL6N")),
        ),
        
        mainPanel(
            plotOutput("expressionPlot")
            # downloadButton('download_plot',"Download Plot")
        )
    ),
    hr(),
    fluidRow(
        column(12,
               reactableOutput("data_summary")
        )
    ),
    hr(),
    fluidRow(id = "download-button",
             column(2,
                    downloadButton("download_data_summary", "Download Selected Data")
             ),
             column(10,
                    downloadButton("download_full_data", "Download Full Data Set")
             )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    selected_data <- reactive({
        if (input$gene == "") {
            return(data.frame())
        } else {
            parnell_data_list[[input$gene]] %>%
                filter(Strain %in% input$mouse_strains) %>%
                #slightly strange if else here, turns out inline if-else like this
                #is cool in a tidyverse pipe
                filter(Treatment %in% if(input$include_PAE) c("Vehicle","PAE") else c("Vehicle")) %>%
                group_by(Timepoint,Strain,Treatment) %>%
                identity()
        }
    })
    
    selected_data_full <- reactive({
        if (input$gene == "") {
            return(data.frame())
        } else {
        parnell_data_full_list[[input$gene]] %>%
            filter(Strain %in% input$mouse_strains) %>%
            #slightly strange if else here, turns out inline if-else like this
            #is cool in a tidyverse pipe
            filter(Treatment %in% if(input$include_PAE) c("Vehicle","PAE") else c("Vehicle")) %>%
            group_by(Timepoint,Strain,Treatment) %>%
            identity()
        }
    })
    
    output$expressionPlot <- renderPlot({
        if (dim(selected_data())[1] > 0) {
            ggplot(selected_data(), aes(x=Strain,y=Mean,fill=Treatment)) +
                geom_bar(stat="identity",position="dodge",alpha=0.5) +
                geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=0.2,position=position_dodge(width=0.9)) +
                geom_point(data=selected_data_full(), 
                           mapping=aes(x=Strain,y=values,color=Treatment), 
                           position=position_jitterdodge(jitter.width=0.4)) +
                facet_wrap(~Timepoint) +
                ggtitle(input$gene) +
                ylim(plot_range)
        } else if (input$gene == "") {
            ggplot(data.frame(text = "Please Select a Gene"),
                   aes(x=0,y=0,label=text)) +
                geom_text() +
                theme_void()
        } else {
            ggplot(data.frame(text = "Please Make a Different Selection\nNo Data Found"),
                   aes(x=0,y=0,label=text)) +
                geom_text() +
                theme_void()
        }
    })
    
    output$data_summary <- renderReactable({
        if (dim(selected_data())[1] > 0) {
            reactable(selected_data(), 
                      filterable = TRUE)
        } else {
            
        }
    })
    
    output$download_data_summary <- downloadHandler(
        filename = "FAS_data_summary.csv",
        content = function(file) {
            write_csv(selected_data(), file)
        }
    )
    
    output$download_full_data <- downloadHandler(
        filename = "FAS_full_data_set.csv.gz",
        content = function(file) {
            write_csv(read_rds(here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_table.rds')), file)
        }
    )
    
    # output$foo = downloadHandler(
    #     filename = 'FAS_plot.png',
    #     content = function(file) {
    #         device <- function(..., width, height) {
    #             grDevices::png(..., width = width, height = height,
    #                            res = 300, units = "in")
    #         }
    #         ggsave(file, plot = plotInput(), device = device)
    #     })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
