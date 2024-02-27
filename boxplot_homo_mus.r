library(ggplot2)
library(data.table)
library(openxlsx)
library(tidyverse)
library(fs)
library(dplyr)
library(stringr)
library(reshape2)
library(hrbrthemes)


homolist <- read.table("homolist_path.txt")
homolist_long <- melt(homolist)

plot_homo <- ggplot(homolist_long, aes(x = variable, y = value)) +          
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=20), alpha= 0.02, fill = "#FEC7C7") +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=20, ymax=28), alpha= 0.02, fill = "#F6DDCC") +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=28, ymax=Inf), alpha= 0.02, fill = "#C8F0D3")+
  geom_boxplot(fill="#7B7D7D") + 
  ggtitle("Homo Sapiens per base quality distribution boxplot")+ 
  theme_ipsum() + xlab("Position in read") 
plot_homo


muslist <- read.table("muslist_path.txt")
muslist_long <- melt(muslist)

plot_mus <- ggplot(muslist_long, aes(x = variable, y = value)) +          
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=20), alpha= 0.02, fill = "#FEC7C7") +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=20, ymax=28), alpha= 0.02, fill = "#F6DDCC") +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=28, ymax=Inf), alpha= 0.02, fill = "#C8F0D3")+
  geom_boxplot(fill="#7B7D7D") + 
  ggtitle("Mus musculus per base quality distribution boxplot")+ 
  theme_ipsum() + xlab("Position in read") 
plot_mus
