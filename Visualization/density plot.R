wd<- "D:/OneDrive/Lake_virome/checkV"
setwd(wd)

# library
library(ggplot2)
library(ggExtra)
library(dplyr)

data <- read.delim('contig_length.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,na.strings="na")

data %>%
  ggplot( aes(x=coverage)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)
