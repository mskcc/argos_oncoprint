setwd("~/Work/oncoprint_Rudin_lab_BIC/")

library(ComplexHeatmap)
library(tidyr)
library(dplyr)
library(readr)



cols_needed <- c("Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification")
dmp_df = read.csv("./dmp_samples_data_mutation_extended_reduced.oncokb_annotated_filtered.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")[,cols_needed]
pdx_df = read.csv("./pdx_mutations_reduced_cols.oncokb_annotated_only_Oncogenic_reduced_columns.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")[,cols_needed]

# they didn't wanted these samples perhaps due to NO IMPACT dmp match
pdx_df = pdx_df[!(pdx_df$Tumor_Sample_Barcode %in% c("s_C_E00XWX_X001_d","s_CR_pdx_239_X","s_CR_pdx_151_X","s_C_37RVCF_X001_d01","s_C_659M3V_X002_d02")),]


combined_df <- dmp_df %>% bind_rows(pdx_df)

#mat<-dmp_df
mat <- combined_df


x <- aggregate(mat[3],mat[-3],FUN = function(X) paste(X, collapse=",")) %>% 
      spread(Hugo_Symbol,Variant_Classification) %>% 
      replace(is.na(.), "")

rownames(x) = x[, 1]
x = x[, -1]
x <- t(x)





#convert sample names to pi names
x<- as.data.frame(x)
x <- x %>% rename(
  "MSK_LX_537"  = "s_CR_pdx_537_X",
  "MSK_LX_740"  ="s_C_92WDH3_X001_d",
  "MSK_LX_871"  ="s_C_90UNXC_X001_d",
  "MSK_LX_726"  ="s_C_U5WXWP_X001_d",
  "MSK_LX_1405a"="s_C_XJ0RYH_X001_d01",
  "MSK_LX_68"   ="s_CR_pdx_068_X",
  "MSK_LX_867"  = "s_C_JKM4C1_X001_d",
  "MSK_LX_1042" = "s_C_60DMC0_X001_d",
  "MSK_LX_1305B" = "s_C_P4RW3W_X001_d01",
  "MSK_LX_148b" = "s_CR_pdx_342_X"
)

ordered_columns <- c(
  "MSK_LX_148b","P-0014813-T01-IM6",
  "MSK_LX_740","P-0014003-T02-IM6",
  "MSK_LX_68", "P-0010571-T01-IM5",
  "MSK_LX_537","P-0016779-T01-IM6",
  "MSK_LX_726","P-0022354-T01-IM6",
  "MSK_LX_1042","P-0023310-T02-IM6",
  "MSK_LX_867","P-0028903-T01-IM6",
  "MSK_LX_871","P-0028176-T01-IM6",
  "MSK_LX_1405a","P-0040730-T02-IM6",
  "MSK_LX_1305B","P-0052717-T02-IM7"
)
x[1,1]<-"no_mutation"

col=c("Missense_Mutation"= "darkgreen",
      "Nonsense_Mutation" = "black",
      "Frame_Shift_Indel" = "purple",
      "Frame_Shift_Del" = "darkorchid1",
      "Frame_Shift_Ins" = "darkorchid3",
      "In_Frame_Indel" = "brown",
      "In_Frame_Del"= "darkred",
      "Splice_Site"= "orange",
      "Nonstop_Mutation" = "pink",
      "no_mutation"= "#CCCCCC"
      )      

heatmap_legend_param = list(
  title = "Alternations",
  at= c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation","no_mutation"),
  labels= c("Missense Mutation","Nonsense Mutation","Frameshift Del","Frameshift Ins","In-frame Del","Splice Site","Nonstop Mutation","No Mutation")
)

alter_fun = list(
  Missense_Mutation   = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Missense_Mutation"],col = NA))},
  Nonsense_Mutation = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Nonsense_Mutation"],col = NA))},
  Frame_Shift_Del    = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Frame_Shift_Del"],col = NA))},
  Frame_Shift_Ins    = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Frame_Shift_Ins"],col = NA))},
  In_Frame_Del    = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["In_Frame_Del"],col = NA))},
  Splice_Site    = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Splice_Site"],col = NA))},
  Nonstop_Mutation    = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["Nonstop_Mutation"],col = NA))},
  no_mutation    = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), (h- unit(0.5,"mm"))*0.33,  gp = gpar(fill = col["no_mutation"],col = NA))},
  background          = function(x, y, w, h) { grid.rect(x, y, w- unit(0.9,"mm"), h- unit(0.5,"mm"),  gp = gpar(fill = "#CCCCCC", col = NA))}
)

pdf("clical-pdx_paired-our_var_classes.pdf")

oncoPrint(x,
          width = ncol(x)*unit(8, "pt"),
          alter_fun = alter_fun, 
          col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          top_annotation=NULL,
          row_names_side = "left", row_names_gp = gpar(fontsize = 9),
          show_pct=TRUE, pct_side = "right", pct_gp = gpar(fontsize = 9),
          right_annotation=NULL,
          column_title = "Paired Clinical - PDX Samples", column_title_gp = gpar(fontsize = 9, fontface='bold'),
          heatmap_legend_param = heatmap_legend_param,
          show_column_names=TRUE,column_names_gp = gpar(fontsize = 9),column_order =ordered_columns,
          )

dev.off()

