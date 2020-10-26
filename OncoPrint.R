# ONCOPRINT
# BY: LIAM FLINN SPURR
# Refer to this link for help with changing any options in the OncoPrint function:
#   https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html

onc_only <- TRUE # logical to specify whether to include just oncogenic mutations

library(tidyverse)
library(data.table)
library(ComplexHeatmap)

clinical <- fread("CLINICAL ANNOTATION/SPECIMEN FILE") %>% # sample metadata/clinical annotations/TMB
  mutate(SAMPLE_ACCESSION_NBR = gsub("-", "_", SAMPLE_ACCESSION_NBR)) # removes the dashes in sample names
samples_list <- clinical$SAMPLE_ACCESSION_NBR

# a list of the genes convered on all panel versions in the data set, you may remove the filters for any versions that are not included in your data
genes_conc <- (fread("profile_versions.txt") %>% filter(V1 == "T", V2 == "T", V3 == "T", Type == "E"))$HUGO_SYMBOL

if(onc_only) {
  mutations <- fread("ONCOGENIC MUTATIONS ONCDRS FILE") %>%
    mutate(SAMPLE_ACCESSION_NBR = gsub("-", "_", SAMPLE_ACCESSION_NBR)) # removes the dashes in sample names
} else {
  mutations <- fread("GERMLINE FILTERED MUTATIONS ONCDRS FILE") %>%
    mutate(SAMPLE_ACCESSION_NBR = gsub("-", "_", SAMPLE_ACCESSION_NBR))
}

# categorize the variants and format into just sample, gene, and event columns
mutations <- mutations %>% filter(SAMPLE_ACCESSION_NBR %in% samples_list, CANONICAL_GENE %in% genes_conc) %>%  # filter to only include the correct genes/samples
  mutate(CANONICAL_VARIANT_CLASS = ifelse(CANONICAL_VARIANT_CLASS %in% c("In_Frame_Del", "Inframe_Del", "Inframe_Ins", "In_Frame_Ins", "protein_altering"), "In-frame mutation",
                                          ifelse(CANONICAL_VARIANT_CLASS %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Frameshift", "Initiator_Codon", 
                                                                                "Splice_Acceptor", "Splice_Donor", "Splice_Site", "Nonsense", 
                                                                                "Nonsense_Mutation", "Stop_Lost"), "Truncating mutation",
                                                 ifelse(CANONICAL_VARIANT_CLASS %in% c("Missense_Mutation", "Missense"), "Missense mutation", 
                                                        ifelse(CANONICAL_VARIANT_CLASS == "Rearrangement", "Rearrangement", "Other mutation"))))) %>%
  dplyr::select(SAMPLE_ACCESSION_NBR, CANONICAL_GENE, CANONICAL_VARIANT_CLASS) %>%
  distinct()
names(mutations) <- c("sample", "gene", "event")

# NOTE: If you have rearrangements to include, just add them to the "mutations" object in the same sample, gene, event (= Rearrangement) format

cna <- fread("COPY NUMBER VARIANTS FILE") %>%
  mutate(SAMPLE_ACCESSION_NBR = gsub("-", "_", SAMPLE_ACCESSION_NBR)) %>%
  filter(CNV_TYPE_CD %in% c("2DEL", "HA"), # keep only 2 copy deletions and high amps
         GENE %in% genes_conc, SAMPLE_ACCESSION_NBR %in% samples_list) %>%  # filter to only include the correct genes/samples
  mutate(Event = ifelse(CNV_TYPE_CD == "2DEL", "Deep deletion", "High amplification")) %>%
  select(SAMPLE_ACCESSION_NBR, GENE, Event) # format into just sample, gene, and event columns
names(cna) <- c("sample", "gene", "event")

# combine the mutations, CNVs, (and SVs if included)
df <- bind_rows(mutations, cna)
df <- df %>% group_by(sample, gene) %>% 
  summarize(event = paste(event, collapse = ';')) %>% 
  ungroup()

# turn the 3 column file into a matrix preserving column and row nuames
mat <- df %>% spread(sample, event) 
r <- unlist(mat[,1])
c <- colnames(mat)
mat <- data.frame(mat[,-1])
rownames(mat) <- gsub("-", "_", r) # remove any dashes in gene names

# transpose matrix to facilitate clinical annotation
mat.t <- t(mat)
mat.t <- data.frame(mat.t)
mat.t$sample <- c[-1]

clin_colors <- clinical %>% select(SAMPLE_ACCESSION_NBR, "INSERT CLINICAL TRACKS (REMOVE QUOTES)") # ADD WHICHEVER CLINICAL TRACKS YOU WANT INTO THE SELECT CALL
names(clin_colors) <- c("sample", "CLINICAL1", "CLINICAL2") # give them names that are easier to type (optional)
mat.t <- left_join(clin_colors, mat.t) # add on the clinical annotation

# order the genes in the OncoPrint based on frequency
gene.order <- as.character((data.frame(gene = rownames(mat), freq = (rowSums(!is.na(mat)) / ncol(mat))) %>%
                              filter(freq >= 0.01) %>% # set a lower limit of gene alteration frequency for inclusion in the OncoPrint
                              arrange(desc(freq)))$gene)
gene.order.names <- sapply(gene.order, as.name)

mat.t <- mat.t %>% arrange("INSERT CLINICAL TRACKS TO SORT ON HERE (REMOVE QUOTES)", !!!gene.order.names) # only include clinical tracks here that you want to use in the sorting of your sample order

sample.order <- unlist(mat.t$sample)

# add back in any samples that dropped off because they had no alterations
missing <- mat.t$sample[!(mat.t$sample %in% colnames(mat))]
mat[,missing] <- NA

sample.order <- sample.order[sample.order %in% colnames(mat)]

# enforce final sample and gene order and convert to matrix data type
mat <- mat[gene.order,sample.order]
mat <- as.matrix(mat)

# give your clinical annotation the same sample order as what will be used in the OncoPrint
s <- data.frame(sample = colnames(mat)) 
clin_colors <- left_join(s, clin_colors) 

get_type_fun <- function(x) strsplit(x, ";")[[1]] # Don't touch this, it splits the mutliple alterations internally in the OncoPrint

# Change the colors to what you want, you can add or remove alterations as needed
col <- c(`Missense mutation` = "#1DA048", 
         `In-frame mutation` = "#F26430",
         `Other mutation` = "#C3C3C3",
         `Rearrangement` = "#8E489B",
         `Truncating mutation` = "#000000",
         `High amplification` = "#D00000", 
         `Deep deletion` = "#226CE0")

# These define the shapes that are plotted on the OncoPrint (ex full cell for deletions vs small boxes for mutations), you can add or remove alterations as needed
alter_fun <- list(background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                                              gp = gpar(fill = "gray97", col = "white")), 
                  `High amplification` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                        gp = gpar(fill = col["High amplification"], col = NA)),
                  `Deep deletion` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                   gp = gpar(fill = col["Deep deletion"], col = NA)),
                  `Rearrangement` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                   gp = gpar(fill = col["Rearrangement"], col = NA)),
                  `Other mutation` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                    gp = gpar(fill = col["Other mutation"], col = NA)),
                  `Missense mutation` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                       gp = gpar(fill = col["Missense mutation"], col = NA)),
                  `In-frame mutation` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                       gp = gpar(fill = col["In-frame mutation"], col = NA)),
                  `Truncating mutation` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                         gp = gpar(fill = col["Truncating mutation"], col = NA)
                  )
)

# get the TMB annotation for the top of the plot
tmb <- left_join(s, clinical %>% select(SAMPLE_ACCESSION_NBR, TMB), by = c("sample" = "SAMPLE_ACCESSION_NBR")) %>% mutate(tmb = as.numeric(TMB))

# choose a name for the file
if(onc_only) name <- "OncoPrint_OncogenicOnly" else name <- "OncoPrint"

# plot the oncoprint
pdf(paste0(name, ".pdf"), width = 15, height = 15, onefile = F) # you'll want to play around with the dimensions based on how many samples/genes you have
o <- oncoPrint(mat, alter_fun = alter_fun, col = col, get_type = get_type_fun, # don't change this line
               column_order = colnames(mat), row_order = rownames(mat), # this uses the sample and gene orders we determined earlier
               top_annotation = HeatmapAnnotation(height = unit(6, "lines"), 
                                                  TMB = anno_barplot(tmb$TMB, which = "column", height = unit(3, "cm"), # plot TMB bars on the top of the plot
                                                                     gp = gpar(fill = "black", fontsize = 24), border = F), 
                                                  `Annotation 1` = clin_colors$ann1, # add the annotations we defined earlier in the clin_colors object, you can name them whatever you want
                                                  `Annotation 2` = clin_colors$ann2, # add in each annotation you want individually
                                                  col = list(`Annotation 1` = c("Category1" = "#1A535C",  # define colors for each category in each annotation you defined (you would have to bin continuous data)
                                                                                "Category2" = "#46ACC2"),
                                                             `Annotation 2` = c("Category1" = "#C73E1D", # make sure the names here match whatever you wrote above
                                                                                "Category2" = "#F26419"))),
               show_pct = TRUE, # show alteration percentages
               show_column_names = FALSE, # hide sample names
               row_names_gp = gpar(fontsize = 20), # change the font size of the gene names
               pct_gp = gpar(fontsize = 20), # change the font size of the percent
               column_split = clin_colors$split_variable, # if you want to split up the oncoprint by a clinical annotation (ex. different cohorts) do it here
               column_gap = unit(0.5, "lines") # choose size of gap between groups
)
o
dev.off()
