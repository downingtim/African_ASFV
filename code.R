library(ape)         # For reading FASTA files
library(seqinr)      # For parsing GenBank files
library(Biostrings)  # For working with DNA sequences in alignment
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggExtra)
library(data.table)
library(stringr) 
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

# Step 1: Read the GenBank file
gb_file <- "TAN_KEN_SUM.gb"
gb_data <- readLines(gb_file)

# Initialize empty lists to store information
samples <- list()
cds_entries <- list()

current_sample <- NULL
current_gene <- NULL  # To keep track of the current gene

# Step 2: Parse the GenBank file to extract CDS and Gene information
for (line in gb_data) {
  # Extract the sample name from the LOCUS entry
  if (grepl("^LOCUS", line)) {
    current_sample <- str_trim(strsplit(line, "\\s+")[[1]][2])  # Sample name
    samples[[length(samples) + 1]] <- current_sample  # Add sample to the list
  }
  
  # Extract Gene information
  if (grepl("^\\s+/gene=", line)) {
    current_gene <- gsub("/gene=\"|\"", "", line)  # Update the current gene
  }
  
  # Extract CDS information
  if (grepl("^\\s+CDS\\s+", line)) {
    cds_line <- gsub("CDS\\s+", "", line)
    if (grepl("complement", cds_line)) {
      coords <- gsub("[^0-9\\.]", "", cds_line)
      coord_split <- unlist(strsplit(coords, "\\.\\."))
      current_cds <- list(
        Start = as.integer(coord_split[2]),
        End = as.integer(coord_split[1]),
        Complement = TRUE,
        Gene_Name = current_gene,
        Sample = current_sample
      )
    } else {
      coords <- gsub("[^0-9\\.]", "", cds_line)
      coord_split <- unlist(strsplit(coords, "\\.\\."))
      current_cds <- list(
        Start = as.integer(coord_split[1]),
        End = as.integer(coord_split[2]),
        Complement = FALSE,
        Gene_Name = current_gene,
        Sample = current_sample
      )
    }
    
    # Add the CDS entry to the list
    cds_entries[[length(cds_entries) + 1]] <- current_cds
  }
}

set1 = trimws(unlist(sapply(cds_entries, function(x) x$Gene_Name)))
set1[516]="pMGF360-21R"

# Step 3: Convert the list of CDS entries into a dataframe
cds_df <- data.frame(
  Start = sapply(cds_entries, function(x) x$Start),
  End = sapply(cds_entries, function(x) x$End),
  Complement = sapply(cds_entries, function(x) x$Complement),
  Gene_Name = set1,
  Sample = sapply(cds_entries, function(x) x$Sample) )

# Display the first few rows of the dataframe to verify correctness
#print(head(cds_df))

alignment_file <- "TAN_SUM_KEN.aln"
aligned_sequences <- readDNAStringSet(alignment_file)
  mapped_positions <- data.frame(
    Start = integer(), End = integer(), Gene_Name = character(),
    Sample = character(), stringsAsFactors = FALSE )

sample_names=c("ASFV_TAN_1987_1", "ASFV_ZMB_SUM14_11", "Kenya1950")

for (sample_name in sample_names){
  cds_subset <- cds_df[cds_df$Sample == sample_name, ] 
  aligned_seq <- as.character(aligned_sequences[[sample_name]]) 
 
for (i in 1: nrow(cds_subset)) {
    cds <- cds_subset[i, ] 
    original_start <- cds$Start
    original_end <- cds$End
    gene_name <- cds$Gene_Name 
    aligned_start <- original_start
    aligned_end <- original_end 
    current_pos <- 1
    for (j in 1:nchar(aligned_seq)) {
      if (substr(aligned_seq, j, j) != "-"){ current_pos <- current_pos + 1}
      if (current_pos == original_start) { aligned_start <- j   }
      if (current_pos == original_end) {
        aligned_end <- j   } }
    
    # Store the results in the dataframe
    mapped_positions <- rbind(mapped_positions, data.frame(
    Gene_Start = original_start,
    Gene_End = original_end,
      Align_Start = aligned_start,
      Align_End = aligned_end,
      Gene_Name = gene_name,
      Sample = sample_name   ))
} }

head(mapped_positions)
write.csv(mapped_positions, "mapped_positions.csv")

# Read the alignment file
aln <- read.alignment("TAN_SUM_KEN.aln", format = "fasta") 
dna <- as.DNAbin(aln) 
window_size <- 1000
step_size <- 50 
results <- data.frame(position = integer(), dist_AB = numeric(),
    dist_AC = numeric())

# Calculate distances in sliding windows
for (i in seq(1000, nchar(aln$seq[[1]]) - window_size + 1 - 1000, by = step_size)) {
  window <- dna[, i:(i + window_size - 1)]
  dist_matrix <- dist.dna(window, model = "K80")
  results <- rbind(results, data.frame(
    position = i,
    dist_AB = dist_matrix[2],
    dist_AC = dist_matrix[3] ))  }
 
plot_data <- results %>%
  pivot_longer(cols = c(dist_AB, dist_AC), names_to = "comparison",
     values_to = "distance") %>%
  mutate(comparison = ifelse(comparison == "dist_AB",
      "TAN_1987_1 vs B", "TAN_1987_1 vs C"))
highlight_bases <- scan("bases.txt", what = integer(), sep = " ")
highlight_regions <- data.frame(xmin=highlight_bases,
    xmax=highlight_bases+step_size-1)
 
highlighted_data <- plot_data %>% filter(position %in% highlight_bases)

##### add MGF360 plot

highlight_bases <- scan("bases.txt", what = integer(), sep = " ")
highlight_regions <- data.frame(xmin = highlight_bases,
         xmax = highlight_bases + step_size - 1)
highlight_regions_dt <- data.table(highlight_regions)
plot_data_wide_dt <- as.data.table(plot_data_wide)
plot_data_wide_dt_old <- plot_data_wide_dt
plot_data_wide_dt <- read.csv("plot_data_wide_dt.csv")
plot_data_wide_dt <- plot_data_wide_dt[,-c(1)]

## subplots ##

plotit <- function(a1, b11){ 
  plot_data2 <- subset(plot_data, position >= a1 & position <= b11)  
  t1 <- subset(mapped_positions, ((Align_End >= a1 & Align_End <= b11) |
    (Align_Start >= a1 & Align_Start <= b11)) & Sample == "ASFV_TAN_1987_1")
  
  if (nrow(t1) == 0) { warning("No data")
      return(NULL) }
  
  # Create the data frame for geom_rect
  b1 <- data.frame(xmin = t1$Align_Start, xmax = t1$Align_End)
  size2 = 4
  angle1=0
  size <-  log2(6*(b11 - a1 + 1))
  if(size > 11){ size2 = 2.5 } 
    if(size > 13){ size2 = 2 } 
    if(size > 15){ size2 = 1.5 } 
    if(size > 17){ size2 =1
                 angle1=90}  
  p2 <- ggplot(plot_data2, aes(x = position, y = distance, color = comparison)) + 
    geom_rect(data = b1, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) + 
    geom_line() +  geom_point(size = 0.001, alpha = 0.01) + 
    geom_text(data = t1, aes(x = (Align_Start + Align_End) / 2,
              y = 0, label = Gene_Name),
              vjust = -0.5, hjust = 0.5, size =size2, angle =angle1, color="black")+ 
    labs(x = "Alignment Position", y = "Genetic Distance", color = " ") +
    theme_minimal() +
    scale_color_manual(values=c("TAN_1987_1 vs B"="blue", "TAN_1987_1 vs C"="red"))+
    theme(legend.position = c(0.35, 0.88), legend.justification = c("right", "top"),
          legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) + 
    scale_y_continuous(limits = c(0, NA)) 
  size <-    log2(0.4*(b11 - a1 + 1))
  pdf(paste0("d_subplot_", a1, "_", b11, ".pdf"), width = size, height = 4)
  print(p2)
  dev.off() }

plotit(4100,6450)  
plotit(14477,21800) 
plotit(42700,190510) 
plotit(122000,132400)
plotit(190550,194633) 

#   Filter for relevant gene entries containing "MGF360" and sample "ASFV_TAN_1987_1"
mgf360_positions <- mapped_positions %>%
  filter(grepl("MGF360", Gene_Name) & Sample == "ASFV_TAN_1987_1")
 
  t1 <- mgf360_positions  
  t1$Gene_Name[5] = "MGF360-5L"
  b1 <- data.frame(
  xmin = t1$Align_Start, 
  xmax = t1$Align_End, 
  Gene_Name = t1$Gene_Name,
  xmid = (t1$Align_Start + t1$Align_End) / 2  )
  
  p <- ggplot(plot_data, aes(x = position, y = distance, color = comparison)) +
  geom_rect(data=b1, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill= Gene_Name),
            alpha = 0.2, inherit.aes = F, show.legend = F) +
  geom_line(linewidth = 0.3) + 
  labs(x = "Position", y = "Genetic distance", color = "") + 
  theme_minimal() +
  theme(legend.position = c(0.7, 0.8),  legend.box.just = "right",
        legend.margin = margin(1,1,1,1),
        legend.background=element_rect(fill="white", linewidth=0,
        linetype="solid", colour="white")) +
  geom_point(size = 0.001, alpha = 0.01) + 
  scale_x_continuous(breaks = pretty(plot_data$position, n = 10)) +  
  scale_y_continuous(limits = c(0, 0.38)) +  
  geom_text(data = b1, aes(x = xmid, y = 0.35, label = Gene_Name), 
            angle = 89, size = 2.5, hjust = 0.5, vjust = 0.5, inherit.aes=F) + 
    scale_color_manual(values = c("TAN_1987_1 vs B" = "blue", 
                                  "TAN_1987_1 vs C" = "red"),
                       labels = c("TAN 1987/1 vs B (ZMB/SUM14/11)", 
                                  "TAN 1987/1 vs C (Kenya 1950)")) 
  pdf(paste0("mgf360_subplot.pdf"), width = 12, height = 4 )
  print(p)
  dev.off()

plotit2 <- function(a1, b11){ 
  plot_data2 <- subset(plot_data, position >= a1 & position <= b11)  
  replace(plot_data2$comparison, plot_data2$comparison=="TAN 1987/1 vs B",
                                             "TAN 1987/1 vs B (ZMB/SUM14/11)") 
  t1 <- subset(mapped_positions, ((Align_End >= a1 & Align_End <= b11) |
    (Align_Start >= a1 & Align_Start <= b11)) & Sample == "ASFV_TAN_1987_1") 
  if (nrow(t1) == 0) { warning("No data in subset t1 for the specified range.")
    return(NULL) } 
  b1 <- data.frame(xmin = t1$Align_Start, xmax = t1$Align_End)
  size2 = 5
  angle1=0  
  p2 <-  ggplot(plot_data2, aes(x = position, y = distance, color = comparison)) + 
    geom_rect(data = b1, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) + 
    geom_line() +  geom_point(size = 0.001, alpha = 0.01) + 
    geom_text(data = t1, aes(x = (Align_Start + Align_End) / 2, y = 0, label = Gene_Name),
              vjust = 0.5, hjust = 0.5, size = size2, angle = angle1, color = "black") + 
    labs(x = "Alignment Position", y = "Genetic Distance", color = "") +
    theme_minimal() +   
    scale_y_continuous(limits = c(0, 0.25)) + 
    scale_color_manual(values = c("TAN_1987_1 vs B" = "blue", 
                                  "TAN_1987_1 vs C" = "red"),
                       labels = c("TAN 1987/1 vs B (ZMB/SUM14/11)", 
                                  "TAN 1987/1 vs C (Kenya 1950)")) + 
    theme(legend.position = "none")
  size <-    log2(0.1*(b11 - a1 + 1))
  return(p2) }
  
# plto each of the MGF360 genes  
plots1 <- c()
for (e in 1:(length(b1$xmin)-1)){
  if(b1$xmax[e] < b1$xmin[e]){ start1 = b1$xmax[e]
                                end1 =   b1$xmin[e] }
  if(b1$xmax[e] > b1$xmin[e]){ start1 = b1$xmin[e]
                                end1 =   b1$xmax[e] }
  plots1[[e]] <- plotit2(start1, end1) }
#str(plots1)

plots1 <- Filter(Negate(is.null), plots1) 
combined_plot <- wrap_plots(plots1, ncol = 6)  # You can adjust ncol for different layouts
ggsave("MGF360_combined.pdf", combined_plot, width = 16, height = 6.4)

# combine plots     

final_plot <- plot_grid(
  p, combined_plot, 
  labels = c("(A)", "(B)"), # Add labels directly using plot_grid
  label_size = 15,          # Adjust label size
  ncol = 1,                 # Stack vertically
  rel_heights = c(0.37, 0.63) # Assign relative heights
)
ggsave("combined_plots.pdf", plot = final_plot, width = 17, height =10)

