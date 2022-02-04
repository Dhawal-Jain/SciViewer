

tabBox(title = "",width = NULL,
       tabPanel(title="FAQs",
                icon = icon("info"),
                fluidRow(column(includeMarkdown("documents/faq.md"), width = 10, offset = 0))
       ),
       tabPanel(title="Glossary",
                icon = icon("info"),
                fluidRow(column(includeMarkdown("documents/glossary.md"), width = 10, offset = 0))
       )
)
