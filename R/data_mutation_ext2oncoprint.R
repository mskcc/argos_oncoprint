library("argparse")
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

# Takes data_mutation_extended.txt and creates the oncoprint

default_output_dir <- normalizePath(getwd())
default_output_file <- file.path(default_output_dir, "Oncoprint.pdf")

parser <- ArgumentParser()
parser$add_argument("-m","--mutations", dest = 'mut_file', default="data_mutations_extended.txt", help="Mutations file")
parser$add_argument("-o", "--output_filename", default=default_output_file, help="Output filename")


args <- parser$parse_args()

input_filename=args$mut_file
output_filename=args$output_filename

cols_needed <- c("Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification")
data_mutations_extened = read.csv(input_filename, header = TRUE, stringsAsFactors = FALSE, comment.char = '#', sep = "\t")[,cols_needed]

data_mutations_extened[data_mutations_extened=="5'UTR"] <- "UTR5prime"
data_mutations_extened[data_mutations_extened=="5'Flank"] <- "Flank5prime"

data <- aggregate(data_mutations_extened[3],data_mutations_extened[-3],FUN = function(X) paste(X, collapse=",")) %>% 
      spread(Hugo_Symbol,Variant_Classification) %>% 
      replace(is.na(.), "")

rownames(data) = data[, 1]
data = data[, -1]
data <- t(data)

col=c("Missense_Mutation"= "darkgreen",
      "Nonsense_Mutation" = "black",
      "Frame_Shift_Del" = "darkorchid1",
      "Frame_Shift_Ins" = "darkorchid3",
      "In_Frame_Del"= "darkred",
      "In_Frame_Ins"="firebrick",
      "Splice_Site"= "orange",
      "Splice_Region"="darkorange",
      "Nonstop_Mutation" = "pink",
      "RNA" = "cyan",
      "Flank5prime" = "coral1",
      "UTR5prime" ="chocolate1",
      "Translation_Start_Site" ="lawngreen"
      )      

heatmap_legend_param = list(
  title = "Alterations",
  at= c("Missense_Mutation","Nonsense_Mutation",
        "Frame_Shift_Del","Frame_Shift_Ins",
        "In_Frame_Del", "In_Frame_Ins",
        "Splice_Site", "Splice_Region",
        "Nonstop_Mutation",
        "RNA" ,
        "Flank5prime" , "UTR5prime" ,
        "Translation_Start_Site"),
  labels= c("Missense Mutation","Nonsense Mutation",
        "Frameshift Deletion","Frameshift Insertion",
        "In-frame Deletion", "In-frame Insertion",
        "Splice Site", "Splice Region",
        "Nonstop Mutation",
        "RNA" ,
        "5'Flank" , "5'UTR" ,
        "Translation Start Site"),
  title_position ="topcenter"
)

alter_fun = list(
  Missense_Mutation         = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Missense_Mutation"],col = NA))},
  Nonsense_Mutation         = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Nonsense_Mutation"],col = NA))},
  Frame_Shift_Del           = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Frame_Shift_Del"],col = NA))},
  Frame_Shift_Ins           = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Frame_Shift_Ins"],col = NA))},
  In_Frame_Del              = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["In_Frame_Del"],col = NA))},
  In_Frame_Ins              = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["In_Frame_Ins"],col = NA))},
  Splice_Site               = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Splice_Site"],col = NA))},
  Splice_Region             = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Splice_Region"],col = NA))},
  Nonstop_Mutation          = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Nonstop_Mutation"],col = NA))},
  RNA                       = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["RNA"],col = NA))},
  Flank5prime               = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Flank5prime"],col = NA))},
  UTR5prime                 = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["UTR5prime"],col = NA))},
  Translation_Start_Site    = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Nonstop_Mutation"],col = NA))},
  background                = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), h- unit(0.5,"mm"),  gp = gpar(fill = "#CCCCCC", col = NA))}
)


pdf(output_filename)

oncoPrint(data,
          top_annotation=NULL,
          width = ncol(data)*unit(8, "pt"),
          alter_fun = alter_fun, 
          col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          show_pct=TRUE, pct_side = "right", pct_gp = gpar(fontsize = 9),
          right_annotation=NULL,
          column_title = "Paired Clinical - PDX Samples", column_title_gp = gpar(fontsize = 9, fontface='bold'),
          heatmap_legend_param = heatmap_legend_param,
          row_names_side = "left", row_names_gp = gpar(fontsize = 9),
          show_column_names=TRUE,column_names_gp = gpar(fontsize = 8)
          )

dev.off()




