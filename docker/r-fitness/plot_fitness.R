#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 3) {
  stop("Usage: Rscript your_script_name.R <style_value> <tsv_file> <r_image>")
}

style_value <- args[1]
tsv_file <- args[2]
r_image <- args[3]

# Different styles group by different columns
if (style_value == "top-proteins") {
  group = "expDesc"
} else if (style_value == "chosen-winners"){
  group = "group"
} else {
  stop("The style value is not valid. Only 'top-proteins' and 'chosen-winners' are allowed.")
}

load(r_image)

generate_plot <- function(group, locus_tag) {
  cleaned_group_name <- gsub("\\s+", "-", trimws(group))
  filename <- paste0(locus_tag, ".svg")
  directory <- file.path(style_value, cleaned_group_name)
  output_path <- file.path(directory, filename)

  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  sprintf("Saving plot to %s", output_path)

  svg(output_path, width = 7, height = 7) # Adjust width and height in inches
  par(mar = c(4, 4.5, 2.5, 1), mgp = c(3, 1, 0), cex.main = 2, cex.axis = 1, cex.lab = 1.5)
  show(group, locus = locus_tag, background = "Ec", ymax = 20,
      col.bg = "darkgrey", main = group)
}

data <- read.table(tsv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Iff data has a column named plasmid_ids, filter out any rows that equal '[]'
# because if it tries to plot and there's no reads in that region it breaks.
if ("plasmid_ids" %in% colnames(data)) {
  data <- data[data$plasmid_ids != "[]", ]
}

mapply(generate_plot, data[, group], data$locus_tag)
dev.off()
