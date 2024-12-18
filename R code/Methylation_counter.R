# The data loaded in these files corresponds to the extraction of information 
# from a .bam file, using samtools mpileup. In the original files only the 
# 1st, 2nd and 5th columns were kept. So the 'Bases' column corresponds to:
# · if the base matches the forward strand reference
# , if the base matches the reverse strand reference
# C for methylated cytocines
# T for unmetylated cytocines
# A,G or other bases in case of mismatch
library(RColorBrewer)

diseases <- c('Down', 'Alzheimer', 'Lewy', 'Parkinson')
Methylation_percentage <- c(0,0,0,0)
start_prot <- 68511780
end_prot <- 68551319
Methylation_prot <- c(0,0,0,0)

# Read the data and compute the Methylation_percentage
for (n in 1:4) {
  path <- paste0(diseases[n], '_base_counts.txt')   #Creates the file name
  data <- read.table(path, header = FALSE)  #Reads each file
  colnames(data) <- c('Chromosome', 'Position', 'Bases') #Asigns column names
  
  #Clean dataframe for possible NAs
  data$Bases <- sapply(data$Bases, function(x) ifelse(is.na(x), "0", x))
  
  #Select only the protein of interest
  prot_data <- data[data$Position >= start_prot & data$Position <= end_prot, ]
  
  # Count methylated (C) and unmethylated (T) reads for all Chr17
  data$Methylated <- sum(sapply(data$Bases, function(x) sum(tolower(strsplit(x, '')[[1]]) == 'c')))
  data$Unmethylated <- sum(sapply(data$Bases, function(x) sum(tolower(strsplit(x, '')[[1]]) == 't')))
  
  Methylation_percentage[n] <- data$Methylated[1] / (data$Methylated[1] + data$Unmethylated[1]) * 100
  
  # Count methylated (C) and unmethylated (T) reads for region of interest
  prot_Methylated <- sum(sapply(prot_data$Bases, function(x) sum(tolower(strsplit(x, '')[[1]]) == 'c')))
  prot_Unmethylated <- sum(sapply(prot_data$Bases, function(x) sum(tolower(strsplit(x, '')[[1]]) == 't')))
  
  Methylation_prot[n] <- prot_Methylated[1] / (prot_Methylated[1] + prot_Unmethylated[1]) * 100
  
  print(paste('Methylation percentage complete for',diseases[n]))
}

print(diseases)
print(Methylation_percentage)
print(Methylation_prot)

#Plot of the Methylation percentage for all Chr17:
coul <- brewer.pal(5, "Set2") 
barplot(height=Methylation_percentage, names=diseases,col=coul, 
        main='Methylation percentage of Chr17', xlab='Diseases', 
        ylab='Methylation percentage (%)', ylim=c(0,12))

#Plot of the Methylation percentage for region of interest:
coul <- brewer.pal(5, "Set2") 
barplot(height=Methylation_prot, names=diseases,col=coul, 
        main='Methylation percentage of PRKAR1A region', xlab='Diseases', 
        ylab='Methylation percentage (%)', ylim=c(0,2.5))
