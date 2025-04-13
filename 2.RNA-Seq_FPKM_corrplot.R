setwd("~/linxing2/workspace/DMY_sci_data/results/DEGs")
getwd()
a <- read.csv("DMY_all_FPKM.csv")
rownames(a) <- a[,1]
a<-a[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28)]
c <- cor(a)
library(pheatmap)
pheatmap(c, cluster_row = TRUE,clustering_method = "ward.D2",display_numbers = FALSE,number_color = "black"
         ,cellwidth = 80, cellheight =40 )
