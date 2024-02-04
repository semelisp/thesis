library(openxlsx)
library(dbplyr)
library(ggplot2)
library(hrbrthemes)
library(ggExtra)
library(tidyverse)

table_graph <- read.xlsx("Total_Table_graphs.xlsx")
median_data <- read.xlsx("Total_Median_Quality.xlsx")

#aligned vs uliquely aligned reads
scatterplot1 <- ggplot(table_graph, aes(x=Aligned.Reads, y=Uniquely.Aligned.Reads, color=Organism)) + 
  geom_point(size=1.5) + theme_ipsum()+ ggtitle("Plot of aligned reads vs uniquely aligned reads") + 
  xlab("Aligned Reads") + ylab("Uniquely Aligned Reads")
scatterplot1

#aligned vs uliquely aligned reads Homo
tHomo <- table_graph[table_graph$Organism == "Homo Sapiens",]
scatterplotH <- ggplot(tHomo, aes(x=Aligned.Reads, y=Uniquely.Aligned.Reads)) + 
  geom_point(color="#7FB3D5",
             fill="#7FB3D5",
             size=1.5,
             alpha=0.9) + 
  theme_ipsum()+ 
  ggtitle("Plot of aligned reads vs uniquely aligned reads Homo Sapiens") + 
  xlab("Aligned Reads") + 
  ylab("Uniquely Aligned Reads")+
  xlim(0.4,1)+
  ylim(0.4,1)
scatterplotH


#aligned vs uliquely aligned reads Mus
tMus <- table_graph[table_graph$Organism == "Mus musculus",]
scatterplotM <- ggplot(tMus, aes(x=Aligned.Reads, y=Uniquely.Aligned.Reads)) + 
  geom_point(color="#AF7AC5",
             fill="#AF7AC5",
             size=1.6,
             alpha=0.9) + 
  theme_ipsum()+ 
  ggtitle("Plot of aligned reads vs uniquely aligned reads Mus musculus") + 
  xlab("Aligned Reads") + 
  ylab("Uniquely Aligned Reads") +
  xlim(0,1) +
  ylim(0.25,1)
scatterplotM


#Median Quality
m_q_graph <- ggplot(median_data, aes(x=Quality.Median)) +
  geom_density(fill="#AB82AB", color="#800080", alpha=0.8) +
  ggtitle("Median Quality Distribution per Dataset") +
  theme_ipsum()
m_q_graph

#Median quality  histogram
q_hist <- ggplot(median_data, aes(x=Quality.Median)) +
  geom_histogram( binwidth=0.3, fill="#AB82AB", color="#800080", alpha=0.8) +
  ggtitle("Per base median quality histogram") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
q_hist


#Median Quality Homo Sapiens
qHomo <- median_data[median_data$Organism == "Homo Sapiens",]
m_q_graph <- ggplot(qHomo, aes(x=Quality.Median)) +
  geom_density(fill="#AD9FAD", color="#800080", alpha=0.8) +
  ggtitle("Median Quality Distribution for Homo Sapiens") +
  theme_ipsum() +
  ylim(0,0.7)
m_q_graph

#Median Quality  Histogram Homo Sapiens
qHomo <- median_data[median_data$Organism == "Homo Sapiens",]
q_hist <- ggplot(qHomo, aes(x=Quality.Median)) +
  geom_histogram( binwidth=0.55, fill="#AD9FAD", color="#800080", alpha=0.8) +
  ggtitle("Per base median quality histogram for Homo Sapiens") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  ) +
  xlim(30,40) +
  ylim(0,10)
q_hist


#Median Quality Mus musculus
qMus <- median_data[median_data$Organism == "Mus musculus",]
m_q_graph <- ggplot(qMus, aes(x=Quality.Median)) +
  geom_density(fill="#D8BFD8", color="#800080", alpha=0.8) +
  ggtitle("Median Quality Distribution for Mus musculus") +
  theme_ipsum() +
  xlim(30,40)
m_q_graph

#Median Quality Mus musculus Hist
qMus <- median_data[median_data$Organism == "Mus musculus",]
q_hist <- ggplot(qMus, aes(x=Quality.Median)) +
  geom_histogram( binwidth=0.55, fill="#D8BFD8", color="#800080", alpha=0.8) +
  ggtitle("Per base median quality histogram for Mus musculus") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )+
  xlim(30,40) +
  ylim(0,10)
q_hist


#Quality distribution
q_graph <- ggplot(table_graph, aes(x=Quality)) + 
  geom_density(fill="#6E6FB8", color="#45448E", alpha=0.8)  +
  ggtitle("Per base Quality distribution") +
  theme_ipsum()
q_graph


#quality distribution histogram
q_hist <- ggplot(table_graph, aes(x=Quality)) +
  geom_histogram( binwidth=0.5, fill="#6E6FB8", color="#45448E", alpha=0.8) +
  ggtitle("Per base quality histogram") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
q_hist


#Homo Sapiens
tHomo <- table_graph[table_graph$Organism == "Homo Sapiens",]
q_graph <- ggplot(tHomo, aes(x=Quality)) + 
  geom_density(fill="#A2A2D8", color="#45448E", alpha=0.8)  +
  ggtitle("Per base Quality distribution on Homo Sapiens") +
  theme_ipsum() +
  xlim(28,40) +
  ylim(0,0.4)
q_graph

#Homo Sapiens Hist
tHomo <- table_graph[table_graph$Organism == "Homo Sapiens",]
q_hist <- ggplot(tHomo, aes(x=Quality)) +
  geom_histogram( binwidth=0.55, fill="#A2A2D8", color="#45448E", alpha=0.8) +
  ggtitle("Per base quality histogram Homo Sapiens") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  ) +
  xlim(28,40) +
  ylim(0,160)
q_hist


#Mus musculus
tMus <- table_graph[table_graph$Organism == "Mus musculus",]
q_graph <- ggplot(tMus, aes(x=Quality)) + 
  geom_density(fill="#ADBBFF", color="#45448E", alpha=0.8)  +
  ggtitle("Per base Quality distribution on Mus musculus") +
  theme_ipsum() +
  xlim(28,40) +
  ylim(0,0.4)
q_graph

#Mus musculus Hist
tMus <- table_graph[table_graph$Organism == "Mus musculus",]
q_hist <- ggplot(tMus, aes(x=Quality)) +
  geom_histogram( binwidth=0.55, fill="#ADBBFF", color="#45448E", alpha=0.8) +
  ggtitle("Per base quality histogram Mus musculus") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  ) +
  xlim(28,40) +
  ylim(0,160)
q_hist

#Aligned Reads Distribution
ar_graph <- ggplot(table_graph, aes(x=Aligned.Reads)) + 
  geom_density(fill="#608054", color="#426135", alpha=0.8)  +
  ggtitle("Aligned Reads (%) distribution") +
  theme_ipsum()
ar_graph

#Aligned Reads histogram
q_hist <- ggplot(table_graph, aes(x=Aligned.Reads)) +
  geom_histogram( binwidth=0.03, fill="#608054", color="#426135", alpha=0.8) +
  ggtitle("Aligned Reads histogram") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
q_hist

#Aligned Reads Distribution Homo
tHomo <- table_graph[table_graph$Organism == "Homo Sapiens",]
ar_graph <- ggplot(tHomo, aes(x=Aligned.Reads)) + 
  geom_density(fill="#8FB082", color="#426135", alpha=0.8)  +
  ggtitle("Aligned Reads (%) distribution on Homo Sapiens") +
  theme_ipsum() +
  xlab("Aligned Reads") +
  ylim(0,40)
ar_graph

#Aligned Reads Histogram Homo
tHomo <- table_graph[table_graph$Organism == "Homo Sapiens",]
q_hist <- ggplot(tHomo, aes(x=Aligned.Reads)) +
  geom_histogram( binwidth=0.02, fill="#8FB082", color="#426135", alpha=0.8) +
  ggtitle("Aligned Reads histogram for Homo Sapiens") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
q_hist

#Aligned Reads Distribution Mus
tMus <- table_graph[table_graph$Organism == "Mus musculus",]
ar_graph <- ggplot(tMus, aes(x=Aligned.Reads)) + 
  geom_density(fill="#86967F", color="#426135", alpha=0.8)  +
  ggtitle("Aligned Reads (%) distribution on Mus musculus") +
  theme_ipsum()+
  xlab("Aligned Reads") +
  ylim(0,40)
ar_graph

#Aligned Reads Hist Mus
tMus <- table_graph[table_graph$Organism == "Mus musculus",]
q_hist <- ggplot(tMus, aes(x=Aligned.Reads)) +
  geom_histogram( binwidth=0.02, fill="#86967F", color="#426135", alpha=0.8) +
  ggtitle("Aligned Reads histogram for Mus musculus") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  ) + ylim(0,500)
q_hist


#Uniquely Aligned Reads Distribution
ar_graph <- ggplot(table_graph, aes(x=Uniquely.Aligned.Reads)) + 
  geom_density(fill="#C4734B", color="#8C3F17", alpha=0.8)  +
  ggtitle("Uniquely Aligned Reads (%) distribution") +
  theme_ipsum() +
  xlim(0,1)
ar_graph

#uniquely Aligned Reads Histogram
q_hist <- ggplot(table_graph, aes(x=Uniquely.Aligned.Reads)) +
  geom_histogram( binwidth=0.02, fill="#C4734B", color="#8C3F17", alpha=0.8) +
  ggtitle("Uniquely Aligned Reads (%) histogram") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
q_hist



#Uniquely Aligned Reads Distribution Homo
tHomo <- table_graph[table_graph$Organism == "Homo Sapiens",]
ar_graph <- ggplot(tHomo, aes(x=Uniquely.Aligned.Reads)) + 
  geom_density(fill="#E1B39B", color="#8C3F17", alpha=0.8)  +
  ggtitle("Uniquely Aligned Reads (%) distribution on Homo Sapiens") +
  theme_ipsum() +
  xlim(0,1) 
ar_graph

#uniquely Aligned Reads Histogram Homo
tHomo <- table_graph[table_graph$Organism == "Homo Sapiens",]
q_hist <- ggplot(tHomo, aes(x=Uniquely.Aligned.Reads)) +
  geom_histogram( binwidth=0.02, fill="#E1B39B", color="#8C3F17", alpha=0.8) +
  ggtitle("Uniquely Aligned Reads (%) histogram Homo Sapiens") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  ) +
  xlim(0,1) +
  ylim(0,120)
q_hist


#Uniquely Aligned Reads Distribution Mus
tMus <- table_graph[table_graph$Organism == "Mus musculus",]
ar_graph <- ggplot(tMus, aes(x=Uniquely.Aligned.Reads)) + 
  geom_density(fill="#F9C9A5", color="#8C3F17", alpha=0.8)  +
  ggtitle("Uniquely Aligned Reads (%) distribution on Mus musculus") +
  theme_ipsum() +
  ylim(0,8) +
  xlim(0,1)
ar_graph

#uniquely Aligned Reads Histogram Mus
tMus <- table_graph[table_graph$Organism == "Mus musculus",]
q_hist <- ggplot(tMus, aes(x=Uniquely.Aligned.Reads)) +
  geom_histogram( binwidth=0.02, fill="#F9C9A5", color="#8C3F17", alpha=0.8) +
  ggtitle("Uniquely Aligned Reads (%) histogram Mus musculus") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  ) +
  xlim(0,1) +
  ylim(0,120)
q_hist


#Aligned Reads Median
m_q_graph <- ggplot(median_data, aes(x=Aligned.Reads)) +
  geom_density(fill="#8F545E", color="#74303C", alpha=0.8) +
  ggtitle("Median Aligned Reads Distribution per Dataset") +
  theme_ipsum()
m_q_graph

q_hist <- ggplot(median_data, aes(x=Aligned.Reads)) +
  geom_histogram( binwidth=0.01, fill="#8F545E", color="#74303C", alpha=0.8) +
  ggtitle("Aligned Reads histogram") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
q_hist

#Median Aligned Reads Homo Sapiens
qHomo <- median_data[median_data$Organism == "Homo Sapiens",]
m_q_graph <- ggplot(qHomo, aes(x=Aligned.Reads)) +
  geom_density(fill="#B59096", color="#74303C", alpha=0.8) +
  ggtitle("Median Aligned Reads Distribution for Homo Sapiens") +
  theme_ipsum() + ylim(0,120)
m_q_graph

qHomo <- median_data[median_data$Organism == "Homo Sapiens",]
q_hist <- ggplot(qHomo, aes(x=Aligned.Reads)) +
  geom_histogram( binwidth=0.005, fill="#B59096", color="#74303C", alpha=0.8) +
  ggtitle("Median Aligned Reads histogram Homo Sapiens") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
q_hist

#Median Aligned Reads Mus musculus
qMus<- median_data[median_data$Organism == "Mus musculus",]
m_q_graph <- ggplot(qMus, aes(x=Aligned.Reads)) +
  geom_density(fill="#E9B6BF", color="#74303C", alpha=0.8) +
  ggtitle("Median Aligned Reads Distribution for Mus musculus") +
  theme_ipsum() + ylim(0,120)
m_q_graph

qMus<- median_data[median_data$Organism == "Mus musculus",]
q_hist <- ggplot(qMus, aes(x=Aligned.Reads)) +
  geom_histogram( binwidth=0.01, fill="#E9B6BF", color="#74303C", alpha=0.8) +
  ggtitle("Median Aligned Reads histogram Mus musculus") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  ) + ylim(0,17)
q_hist

#Uniquely Aligned Reads Median
m_q_graph <- ggplot(median_data, aes(x=Uniquely.Aligned.Reads)) +
  geom_density(fill="#5C8A76", color="#4F816E", alpha=0.8) +
  ggtitle("Median Uniquely Aligned Reads Distribution per Dataset") +
  theme_ipsum()
m_q_graph

q_hist <- ggplot(median_data, aes(x=Uniquely.Aligned.Reads)) +
  geom_histogram( binwidth=0.02, fill="#5C8A76", color="#4F816E", alpha=0.8) +
  ggtitle(" Median Uniquely Aligned Reads histogram per dataset") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
q_hist


#Median Uniquely Aligned Reads Homo Sapiens
qHomo <- median_data[median_data$Organism == "Homo Sapiens",]
m_q_graph <- ggplot(qHomo, aes(x=Uniquely.Aligned.Reads)) +
  geom_density(fill="#A7C6BA", color="#4F816E", alpha=0.8) +
  ggtitle("Median Uniquely Aligned Reads Distribution for Homo Sapiens") +
  theme_ipsum()
m_q_graph

qHomo <- median_data[median_data$Organism == "Homo Sapiens",]
q_hist <- ggplot(qHomo, aes(x=Uniquely.Aligned.Reads)) +
  geom_histogram( binwidth=0.015, fill="#A7C6BA", color="#4F816E", alpha=0.8) +
  ggtitle(" Median Uniquely Aligned Reads histogram per dataset Homo sapiens") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  ) + ylim(0,6)  + xlim(0.4,1)
q_hist

#Median Uniquely Aligned Reads Mus musculus
qMus<- median_data[median_data$Organism == "Mus musculus",]
m_q_graph <- ggplot(qMus, aes(x=Uniquely.Aligned.Reads)) +
  geom_density(fill="#B2EAD5", color="#4F816E", alpha=0.8) +
  ggtitle("Median Uniquely Aligned Reads Distribution for Mus musculus") +
  theme_ipsum() + ylim(0,8)
m_q_graph


qMus<- median_data[median_data$Organism == "Mus musculus",]
q_hist <- ggplot(qMus, aes(x=Uniquely.Aligned.Reads)) +
  geom_histogram( binwidth=0.015, fill="#B2EAD5", color="#4F816E", alpha=0.8) +
  ggtitle(" Median Uniquely Aligned Reads histogram per dataset Mus musculus") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15) 
  ) + ylim(0,6) + xlim(0.4,1)
q_hist
