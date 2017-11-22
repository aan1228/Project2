# Generating Figures using Rstudio

## Ocean Chemical Curves (Figure 1)
```Rstudio
install.packages("ggplot2")
library(ggplot2)

geochem_gradients <- read.csv(file = "geochem_gradients.csv", header = TRUE, sep = ",")
head(geochem_gradients)

ggplot(geochem_gradients, aes(x=PO4, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Phosphate Concentrations at Different Depths") +
  labs(x="[PO4] (uM)", y="Depth (m)")

ggplot(geochem_gradients, aes(x=SI, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Silicate Concentrations at Different Depths") +
  labs(x="[SI] (uM)", y="Depth (m)")

ggplot(geochem_gradients, aes(x=NO3, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Nitrate Concentrations at Different Depths") +
  labs(x="[NO3] (uM)", y="Depth (m)")

ggplot(geochem_gradients, aes(x=Mean_NH4, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Ammonium Concentrations at Different Depths") +
  labs(x="[NH4] (uM)", y="Depth (m)")

ggplot(geochem_gradients, aes(x=Mean_NO2, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Nitrite Concentrations at Different Depths") +
  labs(x="[NO2] (uM)", y="Depth (m)")

ggplot(geochem_gradients, aes(x=Mean_N2O, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Nitrous Oxide Concentrations at Different Depths") +
  labs(x="[N2O] (uM)", y="Depth (m)")

ggplot(geochem_gradients, aes(x=CTD_O2, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Dissolved Oxygen Concentrations at Different Depths") +
  labs(x="[O2] (uM)", y="Depth (m)")

ggplot(geochem_gradients, aes(x=Mean_H2S, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Sulfide Concentrations at Different Depths") +
  labs(x="[H2S] (uM)", y="Depth (m)")

T_salinity_historical <- read.csv(file = "T_Salinity_historical.csv", header = TRUE, sep=",")
head(T_salinity_historical)

ggplot(T_salinity_historical, aes(x=Salinity, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Water Salinity at Different Depths") +
  labs(x="Salinity (psu)", y="Depth (m)")

ggplot(T_salinity_historical, aes(x=Temperature, y=Depth)) +  
  geom_point() +
  scale_y_reverse(lim=c(200,10))+
  geom_path(colour="blue")+
  ggtitle("Water Temperature at Different Depths") +
  labs(x="T (degrees celcius)", y="Depth (m)")
```

## CheckM Output (Figure 2)

Generated using the checkM csv file
```
group4checkM <- 
  read.csv(file.choose())

library(ggplot2)

ggplot(group4checkM, aes(x=Completeness, y=Contamination)) +
  geom_point(aes(colour=Marker_lineage)) +
  geom_vline(xintercept = c(50, 70, 90), colour = "gray45", linetype=2) +
  geom_hline(yintercept = c(5, 10, 15),  colour = "gray45", linetype=2) + 
  xlab("Completeness (%)") +
  ylab("Contamination (%)") +
  xlim(0,100) +
  ylim(0,70)
```

## Annotations (Figure 3)

After annotating with MASH and SILVA, the number of times each phyla showed up in either MASH or SILVA was compiled into a table to generate Figure 3.

|Phylum	|Method	|Count|
|-------|--------|-------|
|Firmicutes	|SILVA	|3|
|Firmicutes	|MASH	|1|
|Proteobacteria	|SILVA|	4|
|Proteobacteria|	MASH	|9|
|Parcubacteria|	SILVA|	1|
|Parcubacteria|	MASH|	0|
|Actinobacteria|	SILVA|	2|
|Actinobacteria|	MASH|	0|
|Planctomycetes|	SILVA|	1|
|Planctomycetes|	MASH|	0|
|Marinimicrobia|	SILVA|	0|
|Marinimicrobia|	MASH	|1|
|Bacteroidetes|	SILVA|	1|
|Bacteroidetes|	MASH	|0|
|Thaumarchaeota|	SILVA|	0|
|Thaumarchaeota|	MASH	|1|
|Nematoda	|SILVA|	1|
|Nematoda	|MASH|	0|
|Tracheophyta	|SILVA|	1|
|Tracheophyta	|MASH|	0|


```
library(ggplot2)
phylum = read.table('Annotations.csv', header=TRUE, sep=',')

ggplot(phylum, aes(factor(Phylum), Count, fill = Method)) +
  geom_bar(stat="identity", position="dodge") +
  xlab("Phylum") +
  ylab("Counts") +
  scale_fill_brewer(palette="Set1") +
  scale_y_continuous(breaks=seq(0,10,2)) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
```
