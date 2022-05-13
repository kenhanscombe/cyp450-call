library(rmarkdown)

output_dir <- "workflow/reports"
render("workflow/scripts/pgxpop_out.Rmd",
    output_dir = output_dir,
    params = list(output_dir = output_dir)
)