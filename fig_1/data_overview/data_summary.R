library(tidyverse)
library(catmaply)


# Rearrange 
# Bulk - mixed
# Bulk methylation below atac
s1 <- read_csv('table.csv')
s1w <- pivot_wider(s1,names_from = c("Experiment"), values_from=c("Tissue"))

s1l <- pivot_longer(s1w, !Patient, names_to = "y")
s1l <- s1l[!is.na(s1l$value),]

# figfile <- 'pituitary_transcriptomic_data'

ylevs <- c("Sex","Age","USP8 Mutation","10x Multiome","Mass Spectroscopy","Methylation","rChIP", "PLA", "mIHC") 

ylevs <- rev(ylevs)
plevs <- unique(s1l$Patient)

# s1l$value <- gsub(pattern = "Tumor",replacement = "Core",x=s1l$value)
# s1l$value <- gsub(pattern = "Normal",replacement = "Margin",x=s1l$value)
# s1l$value <- gsub(pattern = "A",replacement = "Adult",x=s1l$value)
# s1l$value <- gsub(pattern = "P",replacement = "Pediatric",x=s1l$value)
# s1l$value <- gsub(pattern = "F",replacement = "Female",x=s1l$value)
# s1l$value <- gsub(pattern = "M",replacement = "Male",x=s1l$value)
# s1l$value <- gsub(pattern = "\\ \\+\\ ",replacement = "/",x=s1l$value)
# vlevs <- c("F","M","Adult","Pediatric","Core","Core + Margin","CD","sCD","GH","PRL","NFPA")
vlevs <- c("CD","Female","Male","Adult","Pediatric", "Mutated", "No Mutation", "Tumor","Tumor + Normal", "Normal")

s1l$Patient <- factor(s1l$Patient, levels=plevs)
s1l$y <- factor(s1l$y, ylevs)
s1l$value <- factor(s1l$value, vlevs)

library(scales)

# colors <- c("#FEE08B","turquoise", "limegreen", "#E6F598", "lightcoral","dodgerblue" )
colors <- c("turquoise","#FAA300", "limegreen", "#E6F598","#7209B7","#F3B391", "#FF1654", "#02A9EA")

pal <- scale_fill_manual(values = colors)

ggplot(data=s1l, aes(x=Patient, y=y, fill=value))  +
  geom_tile(colour = "grey25", show.legend=T) + 
  scale_x_discrete(position = "bottom", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  pal + 
  coord_equal(ratio=1)+ theme_bw() + labs(y = NULL) + theme(axis.title.x = element_text(size = 24, face='bold')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(vjust=0.4)) +
  theme(legend.position="right") +
  guides(fill = guide_legend(title = NULL, keywidth=8, keyheight=8, default.unit="mm")) +
  theme(axis.text = element_text(size = 16, face="bold")) + theme(legend.text = element_text(size = 16, face='bold'))

ggsave('data_summary.pdf', width = 14, height = 7, units = 'in')
