####Cell cycle genes #####
library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
library(ComplexUpset)
pacman::p_load(tidyverse)
HUMAN_9606_idmapping <-
    readr::read_tsv(fs::path_rel(here::here("Project_Datasets","HUMAN_9606_idmapping.dat")),
                    col_names = FALSE) %>%
    setNames(c("Uniprot", "Type", "ID"))
HUMAN_9606_idmapping_ENS <- inner_join(HUMAN_9606_idmapping %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type),
                                       HUMAN_9606_idmapping %>% subset(Type == "Ensembl") %>% dplyr::select(-Type),
                                       by = "Uniprot") %>% 
    select(-Uniprot) %>% set_names(c("Gene_Name","Ensembl")) %>% 
    subset(!duplicated(Ensembl))
# 1: the list of all the dream target genes
# 2: the list of DREAM target genes affected by downregulation of LOXL2 and losing BRD4s occupancy at promoters
                                          
# I would suppose that list 1 is bigger than list 2, meaning that not all the DREAM target genes are equally affected by the disruption of LOXL2-BRD4 interaction. So the questions that we could try to solve bioinformatically are
# A: are genes in list 2 only early or late dream target genes?
#     B. are promoters of genes in list 2 bound, or not bound, by specific TFs?
#     C: are promoters of genes in list 2 enriched for loxl2 DNA binding regions? 
#     D: ?
#     
DREAM_targets <- read.delim(here::here("Project_Datasets","geneset.txt")) %>% unlist(use.names = F)#https://www.gsea-msigdb.org/gsea/msigdb/cards/FISCHER_DREAM_TARGETS
Cell_variation_HPA_RNA <- openxlsx::read.xlsx(here::here("Project_Datasets","41586_2021_3232_MOESM1_ESM.xlsx"), sheet = 2) %>% 
    pull("X3")
Cell_cycle_proteins <- openxlsx::read.xlsx(here::here("Project_Datasets","41586_2021_3232_MOESM1_ESM.xlsx"), sheet = 3) %>% 
    set_names(.,.[1,]) %>% 
    .[-1,]
Cell_cycle_genes <- openxlsx::read.xlsx(here::here("Project_Datasets","41586_2021_3232_MOESM1_ESM.xlsx"), sheet = 8) %>% 
    set_names(.,.[1,]) %>% 
    .[-1,]

LOXL2_lost_promoter_BRD4S <- read.delim(here::here("Project_Datasets","unique_CT_unique_AB.Homer.promotergenes.txt"),header = F) %>% 
    unlist(use.names = F)
testing2 <- RNA_seq_dream <-  read.delim(here::here("Project_Datasets","DESeq2ALLShrink.txt"))
RNA_seq_dream <-  read.delim(here::here("Project_Datasets","DESeq2ALLShrink.txt")) %>% 
    rownames_to_column("Ensembl") %>% 
    subset(!is.na(pvalue)) %>% 
    mutate(Ensembl = str_remove_all(Ensembl, '\\.[:graph:]*$'),
           padj = p.adjust(pvalue, "BH")) %>% 
    #left_join(HUMAN_9606_idmapping_ENS) %>% 
        left_join(
            getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                  values = HUMAN_9606_idmapping_ENS$Ensembl, mart = mart), by = c("Ensembl" ="ensembl_gene_id" ))  %>% 
    mutate(Lost_BRD4S_loxl2 = hgnc_symbol %in% LOXL2_lost_promoter_BRD4S)
    
UpSetR::fromList(list(All_dream_targets = DREAM_targets,
                      LOXL2_dep_Dream_BRD4S_prom  = LOXL2_lost_promoter_BRD4S %>% .[.%in%DREAM_targets])) %>% 
    UpSetR::upset()

testing <- RNA_seq_dream %>% 
    mutate(hgnc_symbol = if_else(is.na(hgnc_symbol), Ensembl,hgnc_symbol)) %>% 
    mutate(isSignificant = if_else(padj<0.05 & abs(log2FoldChange)>1,T,F),
           DREAM_index = hgnc_symbol %in% DREAM_targets) %>% 
    subset(!is.na(log2FoldChange)) #%>% 
    ggplot(testing,aes(x = log2FoldChange, y = -log10(pvalue), colour = Lost_BRD4S_loxl2,alpha = isSignificant, label = hgnc_symbol)) +
    geom_point(data = . %>% subset(Lost_BRD4S_loxl2 == F))+
    geom_point(data = . %>% subset(Lost_BRD4S_loxl2 == T))+
    ggrepel::geom_label_repel(data = . %>% subset(isSignificant == T & Lost_BRD4S_loxl2 == F))+
    ggrepel::geom_label_repel(data = . %>% subset(isSignificant == T & Lost_BRD4S_loxl2 == T))+
        ggtitle("LOXL2 KD RNA-seq with Lost_BRD4S_loxl2",
                "Genes which lose BRD4S on their promoter upon LOXL2KD are shown in blue")
LOXL2_lost_promoter_BRD4S_cell_cycle <- Cell_cycle_proteins %>% 
    left_join(HUMAN_9606_idmapping_ENS, by = c("ENSG"="Ensembl")) %>% 
    mutate(LOXL2_lost_promoterBRD4S = if_else(Gene_Name %in% LOXL2_lost_promoter_BRD4S,T,F),
           DREAM = if_else(Gene_Name %in% DREAM_targets,T,F))
LOXL2_lost_promoter_BRD4S_cell_cycle %>% 
    subset(ccd_reason !="No") %>% 
    mutate(peakExpressionPseudotime= as.numeric(peakExpressionPseudotime)) %>% 
    ggplot(aes(x = DREAM,  y = peakExpressionPseudotime, colour = LOXL2_lost_promoterBRD4S ))+
    geom_boxplot()+
    theme_bw()+
    geom_jitter(alpha = 0.5,width = 0.2,height = 0)+
    #geom_point(alpha = 0.5, position=position_dodge(width=0.75) )+
    facet_wrap(~LOXL2_lost_promoterBRD4S)+
    ggtitle("Proteins - unique_CT_unique_AB promotergenes with their Pseudo times \n out of 1200 genes in the dataset, 350 have pseudotimes and are considered CCD",
        subtitle = "Pseudo time Expression (low = G1, Middle = S-tr, High = S&G2)")
LOXL2_lost_promoter_BRD4S_cell_cycle_genes <- Cell_cycle_genes %>% 
     mutate(LOXL2_lost_promoterBRD4S = if_else(name %in% LOXL2_lost_promoter_BRD4S,T,F),
           DREAM = if_else(name %in% DREAM_targets,T,F))
LOXL2_lost_promoter_BRD4S_cell_cycle_genes %>% 
    #subset(ccd_reason !="No") %>% 
    mutate(peakExpressionPseudotime= as.numeric(peakExpressionPseudotime)) %>% 
    ggplot(aes(x = DREAM,  y = peakExpressionPseudotime, colour = LOXL2_lost_promoterBRD4S ))+
    geom_boxplot()+
    theme_bw()+
    geom_jitter(alpha = 0.5,width = 0.2,height = 0)+
    #geom_point(alpha = 0.5, position=position_dodge(width=0.75) )+
    facet_wrap(~LOXL2_lost_promoterBRD4S)+
    ggtitle("Genes - unique_CT_unique_AB promotergenes with their Pseudo times \n out of 1200 genes in the dataset, 350 have pseudotimes and are considered CCD",
            subtitle = "Pseudo time Expression (low = G1, Middle = S-tr, High = S&G2)")
DREAM_LOXL2_BRD4S <- data.frame(
   name = HUMAN_9606_idmapping %>% subset(Type=="Gene_Name") %>% pull(ID) %>% unique()) %>% 
    mutate(DREAM = if_else(name %in% DREAM_targets,T,F),
           LOXL2_hit = if_else(name %in% LOXL2_lost_promoter_BRD4S,T,F))
Convictions <- DREAM_LOXL2_BRD4S %>% group_by(DREAM,LOXL2_hit) %>% 
   add_count() %>% 
    dplyr::select(-name) %>% distinct() %>% 
    mutate(DREAM = if_else(DREAM == T,"DREAM","Not_DREAM"),
           LOXL2_hit = if_else(LOXL2_hit == T,"LOXL2","Not_LOXL2")) %>% 
    pivot_wider(names_from = LOXL2_hit,values_from = n) %>% column_to_rownames("DREAM") %>% as.matrix() #%>% 
    fisher.test(Convictions)
DREAM_LOXL2_BRD4S %>% 
    ggplot(aes(x = LOXL2_hit, fill = DREAM))+
    geom_bar(position = "fill")+
    ggtitle("Statistically Significant more Dream Genes in \nunique_CT_unique_AB.Homer.promotergenes than all genes possible\np-value < 2.2e-16")

# BRD4_DREAM <- data.frame(
#     name = HUMAN_9606_idmapping %>% subset(Type=="Gene_Name") %>% pull(ID) %>% unique()) %>% 
#     mutate(DREAM = if_else(name %in% DREAM_targets,T,F),
#            BRD4_promoter_peak = if_else(name %in% Annotated$CT_AB$gene_name,T,F))
# Convictions <- BRD4_DREAM %>% group_by(DREAM,BRD4_promoter_peak) %>% 
#     add_count() %>% 
#     dplyr::select(-name) %>% distinct() %>% 
#     mutate(DREAM = if_else(DREAM == T,"DREAM","Not_DREAM"),
#            BRD4_promoter_peak = if_else(BRD4_promoter_peak == T,"BRD4","Not_BRD4")) %>% 
#     pivot_wider(names_from = BRD4_promoter_peak,values_from = n) %>% column_to_rownames("DREAM") %>% as.matrix() #%>% 
# fisher.test(Convictions)
# BRD4_DREAM %>% 
#     ggplot(aes(x = BRD4_promoter_peak, fill = DREAM))+
#     geom_bar(position = "fill")+
#     ggtitle("Statistically Significant more Dream Genes in \nBRD4_AB promoters than all genes possible\np-value < 2.2e-16")
# ggsave(here::here("Project_Output","BRD4_AB_Dream.png"))
# over <- getEnrichedGO(Annotated, orgAnn="org.Hs.eg.db", condense=TRUE,multiAdjMethod = "BH",maxP = 1)
# testing_AB_enrichment <- over %>% purrr::reduce(rbind) %>% 
#     left_join(goterms %>% enframe("go.id","Go.name")) %>%
#     mutate(Cell_cycle = str_detect(Go.name,"cell cycle"),
#            Significant = if_else(BH.adjusted.p.value<0.05,T,F)) 
# testing_AB_enrichment %>% 
#     subset(Ontology == "BP") %>% 
#     ggplot(aes(x = go.id, y = -log10(pvalue), label =Go.name, colour =Cell_cycle, alpha  = Significant))+
#     geom_point()+
#     ggrepel::geom_label_repel(data = . %>% subset(Cell_cycle ==T& BH.adjusted.p.value<0.05))+
#     ggrepel::geom_label_repel(data = . %>% subset(Cell_cycle ==F& BH.adjusted.p.value<0.003))+
#     #facet_wrap("Ontology")
#     ggtitle("BP BRD4 AB antibody promoters")
# ggsave(here::here("Project_Output","BRD4_AB_BP.png"))
# 
# library(reactome.db)
# path <- getEnrichedPATH(Annotated, "org.Hs.eg.db", "reactome.db", maxP=1,multiAdjMethod = "BH")
ggsave(here::here("Project_Output","BRD4L_BRD4_PromoterPeaks.png"))
ol <- findOverlapsOfPeaks(Chip_peaks_combined$CT_AB, Chip_peaks_combined$LX_AB)
overlaps <- ol$peaklist$`Chip_peaks_combined.CT_AB///Chip_peaks_combined.LX_AB`
out <- genomicElementDistribution(overlaps, 
                                  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                  promoterRegion=c(upstream=2000, downstream=500),
                                  geneDownstream=c(upstream=0, downstream=2000),
                                  promoterLevel=list(
                                      # from 5' -> 3', fixed precedence 3' -> 5'
                                      breaks = c(-2000, -1000, -500, 0, 500),
                                      labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                 "upstream <500b", "TSS - 500b"),
                                      colors = c("#FFE5CC", "#FFCA99", 
                                                 "#FFAD65", "#FF8E32")))

# LX_specific_gene <- setdiff(c(Conditions$LX_BB,Conditions$LX_AB),c(Conditions$CT_BB,Conditions$CT_AB))
# CT_specific_gene <- setdiff(c(Conditions$CT_BB,Conditions$CT_AB),c(Conditions$LX_BB,Conditions$LX_AB))
# LOXL2_independent <- Conditions %>% unlist() %>% setdiff(.,c(LX_specific_gene,CT_specific_gene))
# BRD4_short_specific <-setdiff(c(Conditions$LX_AB,Conditions$CT_AB),c(Conditions$LX_BB,Conditions$CT_BB))
# BRD4_long_specific <- setdiff(c(Conditions$LX_BB,Conditions$CT_BB),c(Conditions$LX_AB,Conditions$CT_AB))
# BRD4_shared <- Conditions %>% unlist() %>% setdiff(.,c(BRD4_long_specific,BRD4_short_specific))
# Promoter_Genes <- data.frame(
#     Gene = Annotated %>% map(.x = .,~.x$gene_name) %>% unlist()  %>% unique()) %>% 
#     mutate(IsDreamtarget = if_else(Gene %in% DREAM_targets,"Dream_Target","Non_Dream")) %>% 
#     mutate(
#         Condition = case_when(
#             (Gene %in% LOXL2_independent) ~ "LOXL2-independent",
#             (Gene %in% CT_specific_gene) ~ "CT",
#             Gene %in% LX_specific_gene ~ "LX"),
#         BRD4_isoform = case_when(
#             
#             Gene %in% BRD4_long_specific ~ "Long",
#             Gene %in% BRD4_short_specific ~ "Short",
#             Gene %in% BRD4_shared ~ "Both isoforms",
#             
#         )) #%>% na.omit()
# 
# Promoter_Genes <- data.frame(
#     Gene = Annotated %>% map(.x = .,~.x$gene_name) %>% unlist()  %>% unique()) %>% 
#     mutate(IsDreamtarget = if_else(Gene %in% DREAM_targets,"Dream_Target","Non_Dream")) %>% 
#     mutate(
#         Condition = case_when(
#             (Gene %in% LOXL2_independent) ~ "LOXL2-independent",
#             (Gene %in% CT_specific_gene) ~ "CT",
#             Gene %in% LX_specific_gene ~ "LX"),
#         BRD4_isoform = case_when(
#        
#         Gene %in% BRD4_long_specific ~ "Long",
#         Gene %in% BRD4_short_specific ~ "Short",
#         Gene %in% BRD4_shared ~ "Both isoforms",
#         
#     )) #%>% na.omit()
# Alluvia <- Promoter_Genes %>% 
#     group_by(IsDreamtarget,Condition,BRD4_isoform) %>% 
#     summarise(Freq = n()) 
# ggplot(as.data.frame(Alluvia),
#        aes(y = Freq, axis1 = IsDreamtarget, axis2 = Condition, axis3 = BRD4_isoform)) +
#     geom_alluvium(aes(fill = IsDreamtarget), width = 1/12) +
#     geom_stratum(width = 1/12, fill = "black", color = "grey") +
#     geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#     #scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
#     scale_fill_brewer(type = "qual", palette = "Set1") +
#     ggtitle("UC Berkeley admissions and rejections, by sex and department")

