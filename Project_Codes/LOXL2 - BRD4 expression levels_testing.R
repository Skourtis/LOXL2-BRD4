if(!require(pacman)){install.packages("pacman");require(pacman)}
libraries <- c("openxlsx", "tidyverse")
p_load(libraries, character.only = T)
Breast_lines <- read.xlsx("./Project_Datasets/CCLE Table 1.xlsx" ,sheet = 2) %>%
    .[.$Tissue.of.Origin == "Breast", 2]

Proteins <- read.xlsx("./Project_Datasets/CCLE Table 2.xlsx" ,sheet = 2)
enzymes <- read.xlsx("./Project_Datasets/uniprot_hsa01100.xlsx")$Entry
#POI <- c("NAMPT", "LDHA", "CD274")
POI <- c("P12268")
Proteins <- Proteins[Proteins$Uniprot_Acc %in% POI, str_detect(colnames(Proteins), "_BREAST_|Gene_Symbol")]
CCLE_enzymes <- Proteins[Proteins$Uniprot_Acc %in% enzymes,]
CCLE_enzymes[is.na(CCLE_enzymes)] <- mean(CCLE_enzymes$CAL120_BREAST_TenPx02, na.rm = T)
Metaboverse_test <- data.frame(Enzymes = CCLE_enzymes$Uniprot_Acc, Fold_Change = CCLE_enzymes$AU565_BREAST_TenPx01/ CCLE_enzymes$PATU8988T_PANCREAS_TenPx03)
write_tsv(Metaboverse_test, "Metaboverse_CCLE_enzymes_test.tsv")
test <- t(Proteins)
colnames(Proteins) <-  str_remove_all(colnames(Proteins), "_BREAST_TenPx..")
#Proteins$Cell_line <- str_remove_all(Proteins$Cell_line, "_SKIN_TenPx..")
Proteins[is.na(Proteins)] <- 0

Proteins <- pivot_longer(Proteins, cols =  colnames(Proteins)[-1],
                              names_to = "Cell_line", values_to = "Abundance")
Proteins$Cell_line[Proteins$.copy == "2"] <- "CAL120_2"

TN_score <- Proteins %>% group_by(Cell_line) %>% summarise(TN_Pos_Score = mean(Abundance[Gene_Symbol %in% c("MCM5", "STMN1", "GLS", "RCL1", "SPOUT1", "ENO1")]),
              PGR_ERBB2 = mean(Abundance[Gene_Symbol %in% c("PGR", "ERBB2")]))

TN_score <- TN_score[order(-TN_score$PGR_ERBB2, TN_score$TN_Pos_Score),]

#Proteins$Cell_line <- factor(Proteins$Cell_line , levels = TN_score$Cell_line)
Proteins_test <- left_join(Proteins, CCLE_cell_lines[,c(1,4)]%>% mutate(CCLE_ID = str_remove_all(CCLE_cell_lines$CCLE_ID, "_BREAST")), 
          by = c("Cell_line" = "CCLE_ID"))

Laura_Ama_cell_lines <- c("MDAMB231","MDAMB468", "CAL51", "T47D", "MCF7")
Proteins[str_detect(Proteins$Cell_line, paste(Laura_Ama_cell_lines, collapse = "|")),] %>%
    group_by(Gene_Symbol) %>%
    ggplot(aes(x = factor(Cell_line, levels = c("T47D", "MCF7", "MDAMB231","MDAMB468", "CAL51")), y = Abundance, colour = Gene_Symbol, group = Gene_Symbol)) +
    geom_point(color = "blue", size = 3)+
    geom_line(aes(color = "blue"), color = "blue") + theme_classic(base_size = 15)+ 
    labs(title = "Chromatome-Associated Protein Abundance")+ xlab("Cell Lines") + ylab("Protein Abundance")+
    scale_color_discrete(name = "Accession", labels = "IMPDH2")
ggsave("P12268_CCLE_Expression.png",width= 300, height = 212 , units = "mm")



+ facet_grid(Pathology, scales= 'free_x', space = "free")

CCLE_cell_lines <- read_tsv("./Project_Datasets/Cell_lines_annotations_20181226.txt")
CCLE_skin <- CCLE_cell_lines[str_detect(CCLE_cell_lines$CCLE_ID,paste(Proteins$Cell_line, collapse = "|")),c(1,4)]
CCLE_skin$CCLE_ID <- str_remove_all(CCLE_skin$CCLE_ID,"_SKIN")

Proteins_skin <- left_join(Proteins, CCLE_skin,by = c("Cell_line"="CCLE_ID"))
# 
# level_order <- unique(unlist(Proteins_skin[Proteins_skin$Gene_Symbol == "NAMPT",] %>% .[order(.$Pathology,.$BRAF_mut, .$Abundance),2])) #this vector might be useful for other plots/analyses
# #Proteins_skin$Cell_line <- as.factor(Proteins_skin$Cell_line)
# 
# Proteins_skin$Group <- as.factor(paste(Proteins_skin$Gene_Symbol,Proteins_skin$Pathology,sep = "_"))
# test <- read_tsv("./Project_Datasets/gene_set_library_crisp.gmt", col_names = F)%>%.[.$X2 == 'skin',]
# test_test = test %>% unite("Mutations",colnames(test[-1]), sep = "_", remove = TRUE)#[,1:2]
# test_test$BRAF_mut <- if_else(str_detect(test_test$Mutations,"BRAF"),"TRUE","FALSE")
# Proteins_skin <- left_join(Proteins_skin,test_test,by = c("Cell_line"= "X1"))
# 
# Proteins_skin %>%
#     subset(Proteins_skin$Gene_Symbol == "NAMPT") %>%
#     ggplot(aes(x = factor(Cell_line, levels =  level_order), y = Abundance, colour = Gene_Symbol, group = Gene_Symbol))+#, group = interaction(Group, BRAF_mut))) +
#     geom_line()+ facet_grid(BRAF_mut ~ Pathology, scales= 'free_x', space = "free")
