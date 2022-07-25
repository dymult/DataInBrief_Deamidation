library(readxl)
library(stringr)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(RColorBrewer)

##NOTE: Currently only works with GPM peptide outputs
#Standard GPM outputs are .csv but this script is written for .xlsx

#Create data frames for compiling deamidation ratios
deamtable <- data.frame(stringsAsFactors = FALSE, 
                       sample_name = character(),
                       "ratioN" = numeric(),
                       "ratioQ" = numeric())

ker <- data.frame(stringsAsFactors = FALSE, sample_name = character(),
                  "N_k" = numeric(),
                  "Q_k" = numeric())

nker <- data.frame(stringsAsFactors = FALSE, sample_name = character(),
                  "N_nk" = numeric(),
                  "Q_nk" = numeric())

#Reading in all .xlsx files in working directory
files <- lapply(list.files(pattern = "^[^~]*.xlsx"), read_excel) 

#Generate sample names based on .xlsx name
names <- gsub(".xlsx", "", list.files(pattern = "^[^~]*.xlsx"))

# start of processing loop for each file
for (i in 1:length(files)){ 
#Read in peptide data and filter cRAP contaminants 
  input <- files[[i]][,c("sequence", "modifications", "protein", "description")]

#Additional contaminants can be added to this list as needed
  contam <- c("BOVIN", "TRYP", "SHEEP", "reversed", "CHICK")
  filt <- input[!grepl(paste0(contam, collapse = "|"), input$protein),]
  
#Unlist the GPM peptide modification column and extract all observed deamidation events
  mod <- unlist(strsplit(filt$modifications, ";"))
  deam <- mod[grepl("Deamidated", mod)]
  
#Count number of observed deamidated Qs and Ns respectively
  eventQ <- length(deam[grepl("Q", deam)])
  eventN <- length(deam[grepl("N", deam)])
  
#Count total number of Qs and Ns in the peptide data
  numQ <- sum(str_count(filt$sequence,"Q"))
  numN <- sum(str_count(filt$sequence,"N"))
  
#Calculate ratio of deamidation events for Q and N
  ratioQ <- eventQ/numQ
  ratioN <- eventN/numN
 
#Compile deamidation data into the deamtable data frame created earlier 
  deamtable[i,] <- c(names[i], ratioN, ratioQ)
  
  # filtering by keratin
  keratin <- filt[grepl("keratin", filt$description),]
  notker <- filt[!grepl("keratin", filt$description),]
  
  dataframes <- list(keratin, notker)

#This loop calculates deamidation ratios for keratin and non-keratin protein groups
    for (j in 1:2){
  data <- dataframes[[j]]
  mod <- unlist(strsplit(data$modifications, ";"))
  deam <- mod[grepl("Deamidated", mod)]
  
  eventQ <- length(deam[grepl("Q", deam)])
  eventN <- length(deam[grepl("N", deam)])
  
  numN <- sum(str_count(data$sequence,"N"))
  numQ <- sum(str_count(data$sequence,"Q"))
  
  ratioN <- eventN/numN 
  ratioQ <- eventQ/numQ 

#This compiles data for keratin protein groups
  if (j == 1){
    ker[i,] <- c(names[i], ratioN, ratioQ)
  }
#This compiles data for non-keratin protein groups
  if (j == 2){
    nker[i,] <- c(names[i], ratioN, ratioQ)
  }
    }
#Merge keratin and non-keratin deamidation ratios in a single data frame by sample
  pgdeam <- merge(ker, nker, by = "sample_name")
}

deammelt <- reshape2::melt(deamtable, id.vars = "sample_name", value.name = "ratio") 
deammelt$ratio <- as.numeric(deammelt$ratio) 
colnames(deammelt)[2] <- "amino_acid" 
deammelt$age <- ifelse(grepl("Modern", x = deammelt$sample_name) == "TRUE", "Modern", "Ancient")
deammelt$residue <- ifelse(grepl("ratioN", x = deammelt$amino_acid) == "TRUE", "Asparagine", "Glutamine")
deammelt$sampleID <- gsub(pattern = "_.*", replacement = "", deammelt$sample_name)

##PLOTTING
#Plot of the bulk deamidation ratios across all samples faceted by sample age
png(filename = "bulkdeam.png", res = 300, height = 1800, width = 2500)
ggplot(data = deammelt, aes(x = sample_name, y = 1-ratio, fill = residue)) + 
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), color = "black", size = 0.2) + 
  scale_y_continuous(limits = c(-0.03,1), expand = c(0,0)) + 
  labs(x = "", y = "Non-Deamidated N and Q (%)", title = "Bulk Deamidation Ratios Across Samples") + 
  facet_grid(~age, scales = "free") + 
  scale_fill_manual(values = brewer.pal(8, "Pastel2")[c(1,4)]) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = "white"), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 16), 
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.spacing.y = unit(1, "lines"),
        plot.margin = margin(1, 1, 1, 2, "lines"), 
        axis.title.y = element_text(vjust = 2), 
        axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size = 10, hjust = 1, angle = 45), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16))
dev.off()


pgdeammelt <- reshape2::melt(pgdeam, id.vars = "sample_name", value.name = "ratio") 
pgdeammelt$ratio <- as.numeric(pgdeammelt$ratio)
colnames(pgdeammelt)[2] <- "amino_acid"
pgdeammelt$age <- ifelse(grepl("Modern", x = pgdeammelt$sample_name) == "TRUE", "Modern", "Ancient")
pgdeammelt$protein <- ifelse(grepl("nk", x = pgdeammelt$amino_acid) == "TRUE", "Non-Keratins", "Keratins")
pgdeammelt$residue <- ifelse(grepl("N_", x = pgdeammelt$amino_acid) == "TRUE", "Asparagine", "Glutamine")
pgdeammelt$sampleID <- gsub(pattern = "_.*", replacement = "", pgdeammelt$sample_name)

#Plotting deamidation ratios for all samples faceted by sample age and protein group
png(filename = "pgdeam.png", res = 300, height = 2000, width = 2800)
ggplot(data = pgdeammelt, aes(x = sample_name, y = 1-ratio, fill = residue)) + 
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), color = "black", size = 0.2) + 
  scale_y_continuous(limits = c(-0.03,1), expand = c(0,0)) + 
  labs(x = "", y = "Non-Deamidated N and Q (%)", title = "Deamidation Ratios of Protein Groups Across\n Modern and Ancient Samples") + 
  facet_grid(protein~age, scales = "free") + 
  scale_fill_manual(values = brewer.pal(8, "Pastel2")[c(1,4)]) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = "white"), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 16), 
        plot.title = element_text(size = 23, hjust = 0.5),
        panel.spacing.y = unit(1, "lines"),
        plot.margin = margin(1, 1, 1, 2, "lines"),
        axis.title.y = element_text(vjust = 2), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 10, hjust = 1, angle = 45), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16)) 
dev.off()

#Create a summary function to calculate SD and mean values
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  
#Subset the protein group data frame to only include non-keratin proteins
  nkdeammelt <- subset(pgdeammelt, pgdeammelt$protein != "Keratins")
#Apply the summary function to subsetted nkdeammelt data frame 
  avdeammelt <- ddply(nkdeammelt, c("sampleID", "residue"), .fun=summary_func,
                  "ratio")

#Plotting average values for merged interior and exterior datasets for each sample
png(filename = "avdeamnonkeratin.png", res = 300, width = 2000, height = 1500)
  ggplot(data = avdeammelt, aes(x = sampleID, y = 1-mean, fill = residue)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.2) + 
  geom_errorbar(aes(ymin=1-mean-sd, ymax=1-mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(limits = c(-0.03,1), expand = c(0, 0)) + 
  labs(x = "", y = "Non-Deamidated N and Q (%)", title = "Average Deamidation Ratio Across\n Samples of Non-Keratin Proteins") + 
  scale_fill_manual(values = brewer.pal(8, "Pastel2")[c(1,4)]) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white"), 
        axis.title = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5),
        strip.text = element_text(size = 16), 
        plot.margin = margin(1, 1, 1, 2, "lines"),
        axis.title.y = element_text(vjust = 2), 
        axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size = 14, hjust = 1, angle = 45), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16))
  dev.off()
