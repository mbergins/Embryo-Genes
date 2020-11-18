library(shiny)
library(tidyverse)
library(here)
library(ggplot2)
library(reactable)

parnell_data = read_rds(here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_forShinyApp.rds'))
plot_range = c(floor(min(parnell_data$Mean-parnell_data$SE)), ceiling(max(parnell_data$Mean+parnell_data$SE)))

gene_list = parnell_data %>% 
    group_by(Gene) %>% 
    summarise(gene_mean_exp = mean(Mean)) %>% 
    arrange(desc(gene_mean_exp))

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("FAS Data Browser"),
    
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "gene",
                        label = h2("Select Gene of Interest"),
                        choices = gene_list),
            checkboxInput("include_EtOH",
                          label = 'Include EtOH Treatment',
                          value = F),
            checkboxGroupInput("mouse_strains", label = h3("Mouse Strains"), 
                               choices = list("C57BL6J" = "C57BL6J", "C57BL6N" = "C57BL6N"),
                               selected = c("C57BL6J","C57BL6N")),
        ),
        # Show a plot of the generated distribution
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
             column(12,
                    downloadButton("download_data_summary", "Download Data")
             )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    selected_data <- reactive({
        parnell_data %>%
            filter(Gene == input$gene) %>%
            filter(Strain %in% input$mouse_strains) %>%
            #slightly strange if else here, turns out inline if-else like this
            #is cool in a tidyverse pipe
            filter(Treatment %in% if(input$include_EtOH) c("Control","EtOH") else c("Control")) %>%
            identity()
    })
    
    output$expressionPlot <- renderPlot({
        if (dim(selected_data())[1] > 0) {
            ggplot(selected_data(), aes(x=Strain,y=Mean,fill=Treatment)) +
                geom_bar(stat="identity",position="dodge") +
                geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=0.2,position=dodge) +
                facet_wrap(~Timepoint) +
                ggtitle(input$gene) +
                ylim(plot_range)
        } else {
            ggplot(data.frame(text = "Please Make a Different Selection\nNo Data Found"),
                   aes(x=0,y=0,label=text)) +
                geom_text() +
                theme_void()
        }
    })
    
    output$data_summary <- renderReactable({
        reactable(selected_data(), 
                  filterable = TRUE)
    })
    
    output$download_data_summary <- downloadHandler(
        filename = "FAS_data_summary.csv",
        content = function(file) {
            write_csv(selected_data(), file)
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
