library(shiny)
library(tidyverse)
library(here)
library(ggplot2)
library(reactable)
library(dqshiny)
library(shinylogs)

parnell_data_list = read_rds(here('Parnell_data_split.rds'))
parnell_data_full_list = read_rds(here('Parnell_data_full_split.rds'))

plot_range = lapply(parnell_data_list, 
                           function(x) {
                               return(c(min(x$Mean - x$SE), max(x$Mean + x$SE)))
                           }) %>%
    unlist() %>%
    range()
plot_range[1] = floor(plot_range[1])
plot_range[2] = ceiling(plot_range[2])

lower_case_match = read_rds(here('lower_case_match.rds'))

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Gastrulation-stage mouse embryo transcriptome browser"),
    
    tags$head(tags$style(".rightAlign{float:right;}")),
    
    tags$p("This visualization tool contains data supporting our manuscript CITATION LINK HERE."),

    tags$p("RNA-seq was performed on normally developing and alcohol-exposed mouse embryos from either the C57BL/6J or C57BL/6NHsd strains. Samples were collected during gastrulation at either embryonic day (E) 7, E7.25, or E7.5. Alcohol-exposed embryos were exposed to alcohol beginning at E7. Our study provides a view of how expression of specific genes changes across gastrulation during normal development and the impact of strain and/or alcohol on key developmental genes. Data are normalized based on the mean VST-normalized values from each strain and time point to allow comparison across age, strain, and treatment. Negative values indicate very low/no expression of a particular gene, while values ~0 indicates low expression, and values above 0 indicates moderate to high levels of expression. VST-normalized values for each replicate can be found in the supplemental files of the linked publication."),
    sidebarLayout(
        sidebarPanel(
            autocomplete_input("gene", 
                               h2("Select a Gene of Interest"), 
                               c(sort(names(parnell_data_list)), sort(tolower(names(parnell_data_list)))),
                               placeholder = "Start Typing to Find a Gene",
                               max_options = 100),
            tags$p("*Searches are case-sensitive, format gene name with capital first letter"),
            # selectizeInput(inputId = "gene",
            #                label = h2("Select Gene of Interest"),
            #                choices = gene_list,
            #                multiple = F,
            #                options = list(maxOptions = 10,
            #                               placeholder = "Please Select a Gene",
            #                               items = "",
            #                               maxItems = 1,
            #                               onInitialize = I('function() { this.setValue(""); }'),
            #                               createFilter = I('function(input) { return input.length >= 1;}'),
            #                               openOnFocus = F,
            #                               closeAfterSelect = T
            #                               )
            #                ),            
            checkboxGroupInput("mouse_strains", label = h3("Mouse Strains"), 
                               choices = list("C57BL/6J" = "C57BL6J", "C57BL/6N" = "C57BL6N"),
                               selected = c("C57BL6J","C57BL6N")),
            hr(),
            checkboxInput("include_PAE",
                          label = 'Include prenatal alcohol exposure (PAE) data',
                          value = F)
            
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
             column(9,
                    downloadButton("download_full_data", "Download Full Data Set")
             ),
             column(1,
                    tags$a(href="https://github.com/mbergins/Embryo-Genes", icon("github", class="rightAlign fa-2x")))
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    #Blocking tracking until after publication
    # track_usage(storage_mode = store_rds(path = here('logs')))
    
    selected_data <- reactive({
        selected_gene = input$gene
        if (is.null(parnell_data_list[[input$gene]])) {
            if (is.null(lower_case_match[[selected_gene]])) {
                return(data.frame())
            } else {
                selected_gene = lower_case_match[[selected_gene]]
            }
        }
        
        if (input$gene == "") {
            return(data.frame())
        } else {
            parnell_data_list[[selected_gene]] %>%
                mutate(Treatment = fct_relevel(Treatment,c("PAE","Vehicle"))) %>%
                mutate(str_treat = paste0(Strain,"\n",Treatment)) %>%
                mutate(str_treat = fct_relevel(str_treat, c("C57BL6J\nVehicle","C57BL6N\nVehicle"))) %>%
                filter(Strain %in% input$mouse_strains) %>%
                #slightly strange if else here, turns out inline if-else like this
                #is cool in a tidyverse pipe
                filter(Treatment %in% if(input$include_PAE) c("Vehicle","PAE") else c("Vehicle")) %>%
                group_by(Timepoint,Strain,Treatment) %>%
                identity()
        }
    })
    
    selected_data_full <- reactive({
        selected_gene = input$gene
        if (is.null(parnell_data_list[[input$gene]])) {
            if (is.null(lower_case_match[[selected_gene]])) {
                return(data.frame())
            } else {
                selected_gene = lower_case_match[[selected_gene]]
            }
        }
        
        if (input$gene == "") {
            return(data.frame())
        } else {
            parnell_data_full_list[[selected_gene]] %>%
                mutate(Treatment = fct_relevel(Treatment,c("PAE","Vehicle"))) %>%
                mutate(str_treat = paste0(Strain,"\n",Treatment)) %>%
                mutate(str_treat = fct_relevel(str_treat, c("C57BL6J\nVehicle","C57BL6N\nVehicle"))) %>%
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
            cols <- c("C57BL6J\nPAE" = "grey30", 
                      "C57BL6J\nVehicle" = "darkgreen", 
                      "C57BL6N\nPAE" = "grey60", 
                      "C57BL6N\nVehicle" = "darkblue")

            ggplot(selected_data(),aes(x=Strain,y=Mean,color=str_treat, fill=str_treat)) +
                geom_hline(yintercept = 0, color = "black") +
                geom_bar(stat="identity",position="dodge",alpha=0.5) +
                geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=0.2,position=position_dodge(width=0.9)) +
                geom_point(data=selected_data_full(),
                           mapping=aes(x=Strain,y=values,color=str_treat, fill = str_treat),
                           position=position_jitterdodge(jitter.width=0.25)) +
                labs(y="Mean-centered expression", 
                     color = "Strain\nTreatment",
                     fill = "Strain\nTreatment") +
                facet_wrap(~Timepoint) +
                ggtitle(selected_data()$Gene[1]) +
                theme(text = element_text(size=20)) + 
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(), 
                      axis.line = element_line(colour = "black")) +
                scale_colour_manual(values = cols) +
                scale_fill_manual(values = cols) +
                ylim(plot_range)

        } else if (input$gene == "") {
            ggplot(data.frame(text = "Please Select a Gene"),
                   aes(x=0,y=0,label=text)) +
                geom_text(size=20) +
                theme_void()
        } else {
            ggplot(data.frame(text = "Please Make a Different Selection\nNo Data Found"),
                   aes(x=0,y=0,label=text)) +
                geom_text(size=20) +
                theme_void()
        }
    })
    
    output$data_summary <- renderReactable({
        if (dim(selected_data())[1] > 0) {
            reactable(selected_data() %>% 
                          mutate(Mean = signif(Mean,3),
                                 SE = signif(SE, 3)) %>%
                          select(-str_treat), 
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
        filename = "FAS_full_data_set.csv",
        content = function(file) {
            write_csv(read_rds(here('Parnell_data_reformated.rds')), file)
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
