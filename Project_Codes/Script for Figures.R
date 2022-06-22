#figure 1#
library("ggpubr")
library("readr")
library("tidyverse")
library("openxlsx")
### CPTAC Metastasis
#https://cptac-data-portal.georgetown.edu/study-summary/S039

#panel 1F, BRD4 and LOXL2 between breast subtypes ##
S039_BRCA_propective_clinical_data <- read.xlsx("./Project_Datasets/S039_BRCA_propective_clinical_data_r1.xlsx") %>% 
.[,c("Participant_ID","Estrogen_Receptor_Status_by_IHC","Progesterone_Receptor_Status_by_IHC","HER2_ERBB2_Status_by_IHC")] %>% na.omit()

S039_BRCA_propective_clinical_data<- S039_BRCA_propective_clinical_data %>%
mutate(IHC = case_when(
    Estrogen_Receptor_Status_by_IHC == "Positive" & Progesterone_Receptor_Status_by_IHC == "Positive" & HER2_ERBB2_Status_by_IHC == "Positive" ~ "TP",
    Estrogen_Receptor_Status_by_IHC == "Positive" & Progesterone_Receptor_Status_by_IHC == "Negative" & HER2_ERBB2_Status_by_IHC == "Negative" ~ "HR+",
    Estrogen_Receptor_Status_by_IHC == "Positive" & Progesterone_Receptor_Status_by_IHC == "Positive" & HER2_ERBB2_Status_by_IHC == "Negative" ~ "HR+",
    Estrogen_Receptor_Status_by_IHC == "Negative" & Progesterone_Receptor_Status_by_IHC == "Positive" & HER2_ERBB2_Status_by_IHC == "Negative" ~ "HR+",
    Estrogen_Receptor_Status_by_IHC == "Negative" & Progesterone_Receptor_Status_by_IHC == "Negative" & HER2_ERBB2_Status_by_IHC == "Positive" ~ "HER2+",
    Estrogen_Receptor_Status_by_IHC == "Negative" & Progesterone_Receptor_Status_by_IHC == "Negative" & HER2_ERBB2_Status_by_IHC == "Negative" ~ "TN",
    TRUE ~ "Other")) %>% 
    subset(!(IHC == "Other")) 

CPTAC_S039 <- read_tsv("./Project_Datasets/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv") %>%
    dplyr::select(matches("Gene|Ratio")) %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance" ) %>%
    subset(!str_detect(Sample,"Unshared")) %>%
    mutate(Sample = str_remove(Sample," Log Ratio"))
sample_ids <- read.xlsx(here::here("Project_Datasets","S039_Breast_Cancer_Prospective_Collection_Specimens.xlsx"))[,c(1:4)]
sample_ids <- left_join(sample_ids, S039_BRCA_propective_clinical_data[,c("Participant_ID","IHC")], 
                        by = c(`Participant.Protocol.Identifier.:.Collection.Protocol.Registration` = "Participant_ID")) %>% na.omit() %>%
    mutate(IHC = if_else(Sample.Type == "Adjacent_Normal", "Adjacent", IHC))%>%
    mutate(IHC = factor(IHC, levels = c("Adjacent","HR+","HER2+", "TP","TN")))



CPTAC_S039 <- CPTAC_S039 %>% left_join(sample_ids[,c("Specimen.Label","IHC")] %>% set_names(c("Sample","IHC")))

Genes_of_interest <- c("BRD4" ,"LOXL2" , "EPCAM" ,"ERBB2" ,"ESR1" , "VIM","PGR")
CPTAC_S039_GOI<-CPTAC_S039 %>% subset(Gene %in% Genes_of_interest) %>% na.omit() %>% 
    mutate(Gene = factor(Gene,levels =Genes_of_interest ))

CPTAC_S039_GOI%>%
    ggplot((aes(x = IHC,y = Abundance, colour  =IHC)))+
    geom_boxplot()+geom_point(alpha = 0.5) + facet_wrap("Gene", ncol = 2)  +
    #      geom_rect(data = subset(tp,Gene %in% c('BRD4',"LOXL2")),fill = "#D54FA4",xmin = -Inf,xmax = Inf, 
    #             ymin = -Inf,ymax = Inf,alpha = 0.1) + 
    theme_bw() +
    scale_color_manual(values = c(Adjacent = "#F7A8D7",
                                  `HR+` = "#FF9EE2",
                                  `HER2+` = "#FF67B0",
                                  TP = "#F62681",
                                  TN = "#BD0055"))+
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    ggtitle("CPTAC proteomics data", subtitle =  "Breast cancer Patient samples")

kruskal.test(Abundance ~ IHC, data = CPTAC_S039_GOI %>% subset(Gene == "LOXL2"))


### panel 1g & S3D ###

##### Cell Sensitivity to BRD4 inh ####
sample_info <-  read.csv(fs::path_rel(here::here("Project_Datasets","sample_info.csv")), stringsAsFactors = FALSE)
HUMAN_9606_idmapping <-
    readr::read_tsv(fs::path_rel(here::here("Project_Datasets","HUMAN_9606_idmapping.dat")),
                    col_names = FALSE) %>%
    setNames(c("Uniprot", "Type", "ID"))
##panel 1g##
loading_CCLE_TPM <- function(x,y){
    RNA_seq <- read_tsv(x) %>% dplyr::select(-transcript_ids) %>% 
        mutate(gene_id = str_remove_all(gene_id,"\\.[:graph:]*$")) %>%
        column_to_rownames("gene_id") %>% 
        set_names(str_remove_all(colnames(.),"_[:graph:]*$")) 
    RNA_seq <- RNA_seq[rowSums(RNA_seq == 0) <= (ncol(RNA_seq)*0.75), ] %>%
        .[rowSums(.)>100,] %>% 
        rownames_to_column("gene_id")%>%
        left_join(inner_join(y %>% 
                                 subset(Type == "Ensembl"), 
                             y %>% 
                                 subset(Type == "Gene_Name"), 
                             by = "Uniprot")  %>% 
                      dplyr::select(contains("ID")) %>% 
                      distinct() %>% 
                      set_names(c("gene_id","Gene_Name")), by = "gene_id")
    top_id <- data.frame(Gene_Name = RNA_seq$Gene_Name,
                         gene_id = RNA_seq$gene_id,
                         Abundance = RNA_seq %>% dplyr::select(-c(gene_id,Gene_Name)) %>%  
                             rowSums())%>%  group_split(Gene_Name) %>% 
        map_df(.x = .,~.x %>% arrange(desc(Abundance)) %>% head(1)) %>% pull(gene_id)
    RNA_seq <- RNA_seq %>% subset((gene_id %in% top_id) & !duplicated(Gene_Name)) %>% 
        na.omit() %>% dplyr::select(-gene_id)%>% remove_rownames()%>%
        column_to_rownames("Gene_Name")  %>% `+`(.,0.01)%>%
        `/`(.,matrixStats::rowMedians(as.matrix(.),na.rm = T)) %>% log2() %>% as.matrix
    RNA_seq
    
}
CCLE_RNA_seq_file <- fs::path_rel(here::here("Project_Datasets","CCLE_RNAseq_rsem_genes_tpm_20180929.txt"))#Raw from https://portals.broadinstitute.org/ccle/data 3/18/2020
CCLE_RNA_seq<-loading_CCLE_TPM(CCLE_RNA_seq_file,HUMAN_9606_idmapping)
extract_top_btm_expressing_cell_lines<-function(Gene_Name){
    #Gene_Name = "LOXL2"
    ordered <- CCLE_RNA_seq[Gene_Name,]%>% unlist()%>%#.[names(.) %in% (sample_info %>% subset(lineage == "breast") %>% pull(stripped_cell_line_name))] %>% 
        sort(decreasing = T)
    list(
        top_cell_lines = head(ordered,round(length(ordered)/10)) %>% names(),
        btm_cell_lines = tail(ordered,round(length(ordered)/10)) %>% names()) 
}
extreme_cell_lines<- map(c("LOXL2", "BRD4"),extract_top_btm_expressing_cell_lines) %>%
    set_names(c("LOXL2", "BRD4"))
sec_screen_info <- read_csv(here::here("Project_Datasets","secondary-screen-replicate-collapsed-treatment-info.csv"))
BRD4_inh <- sec_screen_info %>% subset(str_detect(target, "BRD4")) %>% pull("broad_id") %>% paste(., collapse = "|")
sec_screen <- data.table::fread(here::here("Project_Datasets","secondary-screen-replicate-collapsed-logfold-change.csv")) %>% 
    inner_join(.,sample_info[,1:2], by = c("V1" = "DepMap_ID"))%>%
    #column_to_rownames(var = "stripped_cell_line_name")  %>% 
    dplyr::select(matches(BRD4_inh)|contains("stripped"))  %>% 
    pivot_longer(cols = -stripped_cell_line_name, names_to = "column_name", values_to = "Sensitivity") %>%
    left_join(sec_screen_info)
extreme_cell_lines_LOXL2<-rbind(data.frame(stripped_cell_line_name = extreme_cell_lines$LOXL2$top_cell_lines,
                                           Expression = "High"),
                                data.frame(stripped_cell_line_name = extreme_cell_lines$LOXL2$btm_cell_lines,
                                           Expression = "Low"))

sec_screen_BRDi <- sec_screen %>% inner_join(extreme_cell_lines_LOXL2) 
compounds <- sec_screen_BRDi$name %>% unique() %>% .[c(7,1,3,5,6)] %>% set_names(.,.)

#mutate(dose = factor(dose,unique(dose))) %>%
# sec_screen %>% subset(str_detect(name,"JQ") & dose == 10) %>%
#     group_by(stripped_cell_line_name) %>%
#     summarise(Mean_Sensitivity = median(Sensitivity)) %>%
#     inner_join(CCLE_Protein %>% t() %>% .[,"Q9Y4K0"] %>% as.data.frame()%>% set_names( c("LOXL2")) %>% rownames_to_column("stripped_cell_line_name") )%>%
#     left_join(sample_info[,c("stripped_cell_line_name","lineage", "lineage_sub_subtype","lineage_molecular_subtype")])%>%
#     ggplot(aes(x = Mean_Sensitivity,y = LOXL2))+
#     geom_point(alpha = 0.5)+
#     geom_point(data = . %>% subset(lineage == "breast"), aes( colour = lineage_molecular_subtype))+
#     facet_wrap("lineage")+
#     theme_bw()+
#     ggtitle("LOXL2 proteomic Abundance with JQ1 sensitivity")


sec_screen_BRDi %>%
    subset((name %in% compounds) & between(log2(dose), -4,4))%>% 
    mutate(name = factor(name, levels = compounds)) %>% 
    group_by(stripped_cell_line_name) %>%
    ggplot(aes(x = log2(dose), y = Sensitivity, colour = Expression, group = stripped_cell_line_name))+
    geom_line(alpha = 0.1, )+
    stat_summary(aes(group=Expression,colour = Expression), fun=median,  geom="line")+
    geom_point(alpha = 0.2) + 
    facet_grid("name")+
    scale_colour_manual(values = c(High = "black", Low = "#8C1DB8"))+
    ggtitle("BRD4i Sensitivity", 
            subtitle = "Between High and Low RNA expressing LOXL2 Cell Lines")+
    # labs(caption =paste(paste(compounds,p_val, sep = "_"), collapse = "     "))+
    theme_bw()

##panel S3D##
loading_CCLE_prot <- function(x){
    CCLE_proteins <- openxlsx::read.xlsx(x, sheet =2) %>%        .[,!str_detect(colnames(.),"Peptides|Column")]
    colnames(CCLE_proteins) <- str_remove_all(colnames(CCLE_proteins),"_TenPx..")
    CCLE_proteins$Uniprot_Acc <- str_remove_all(CCLE_proteins$Uniprot_Acc,"-.")
    CCLE_proteins <- CCLE_proteins %>% subset(!duplicated(Uniprot_Acc))
    rownames(CCLE_proteins) <- NULL
    CCLE_proteins <- CCLE_proteins %>%
        column_to_rownames(var = "Uniprot_Acc") %>%
        dplyr::select(-c(1:5)) 
    colnames(CCLE_proteins) <- str_match(colnames(CCLE_proteins),"^([:graph:]*?)_")[,2]
    CCLE_proteins <- CCLE_proteins[,which(!duplicated(colnames(CCLE_proteins)))] ##Removing cell_lines with the same name
    index_to_keep <- CCLE_proteins %>% is.na() %>%
        rowSums()<(ncol(CCLE_proteins)*0.75)
    CCLE_proteins %>% .[index_to_keep,] # removing lines with more than 75% NA]
}
CCLE_Protein<-loading_CCLE_prot(fs::path_rel(here::here("Project_Datasets","Table_S2_Protein_Quant_Normalized.xlsx")))

extract_top_btm_protein_expressing_cell_lines<-function(Gene_Name){
    #Gene_Name = "Q9Y4K0"
    ordered <- CCLE_Protein[Gene_Name,]%>%unlist() %>%  na.omit() %>%  #.[names(.) %in% (sample_info %>% subset(lineage == "breast") %>% pull(stripped_cell_line_name))] %>% 
        sort(decreasing = T)
    list(
        top_cell_lines = head(ordered,round(length(ordered)/10)) %>% names(),
        btm_cell_lines = tail(ordered,round(length(ordered)/10)) %>% names()) 
}

extreme_protein_cell_lines<- map(c("Q9Y4K0", "O60885"),extract_top_btm_protein_expressing_cell_lines) %>%
    set_names(c("LOXL2", "BRD4"))

extreme_prot_cell_lines_LOXL2<-rbind(data.frame(stripped_cell_line_name = extreme_protein_cell_lines$LOXL2$top_cell_lines,
                                                Expression = "High"),
                                     data.frame(stripped_cell_line_name = extreme_protein_cell_lines$LOXL2$btm_cell_lines,
                                                Expression = "Low"))
sec_screen_BRDi <- sec_screen %>% inner_join(extreme_prot_cell_lines_LOXL2) 
#mutate(dose = factor(dose,unique(dose))) %>%



sec_screen_BRDi %>%
    subset((name %in% compounds) & between(log2(dose), -4,4))%>% 
    mutate(name = factor(name, levels = compounds)) %>% 
    group_by(stripped_cell_line_name) %>%
    ggplot(aes(x = log2(dose), y = Sensitivity, colour = Expression, group = stripped_cell_line_name))+
    geom_line(alpha = 0.1, )+
    stat_summary(aes(group=Expression,colour = Expression), fun=median,  geom="line")+
    geom_point(alpha = 0.2) + 
    facet_grid("name")+
    scale_colour_manual(values = c(High = "black", Low = "#8C1DB8"))+
    theme_bw()+
    ggtitle("BRD4i Sensitivity", 
            subtitle = "Between High and Low Prot expressing LOXL2 Cell Lines")
###S3 A, B ###
selected_genes <- c("BRD4","LOXL2")
HUMAN_9606_idmapping_selected <- HUMAN_9606_idmapping %>% subset(ID %in% selected_genes) %>% dplyr::select(-Type) %>% distinct()
Breast <- sample_info %>% subset(lineage == "breast") %>% pull(stripped_cell_line_name)
scale_colours <- c(luminal = "#F8ABD9",
                   luminal_HER2_amp = "#FF9EE2",
                   HER2_amp = "#FF68B0",
                   basal_A = "#F62983",
                   basal_B = "#BD0055")
my.formula <- y ~ x

library(grid)
empty_panel <- grid.rect(gp=gpar(col="white"))
interesting_lineage <-  sample_info$lineage %>% unique() %>% .[.%in% c("bone","blood","breast","central_nervous_system","colorectal","esophagus","gastric","kidney","liver","lung","lymphocyte","ovary","pancreas","plasma_cell","prostate")]

p1 <- CCLE_Protein %>% rownames_to_column("Uniprot")%>% 
    subset(Uniprot %in% HUMAN_9606_idmapping_selected$Uniprot) %>% left_join(HUMAN_9606_idmapping_selected) %>% 
    dplyr::select(-Uniprot) %>% 
    column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
    rownames_to_column("stripped_cell_line_name") %>% 
    left_join(sample_info[,c("lineage","stripped_cell_line_name","lineage_molecular_subtype")]) %>% 
    subset(lineage %in% interesting_lineage) %>% 
    mutate(lineage= str_to_title(lineage), lineage_molecular_subtype = if_else(lineage == "Breast",lineage_molecular_subtype,NA_character_)) %>% 
    #pivot_longer(-ID, names_to = "Cell_line", values_to = "Protein_Abundance") %>% 
    add_count(lineage) %>% subset(n>4) %>% 
    ggplot(aes(x = BRD4, y = LOXL2, colour = lineage_molecular_subtype))+
    scale_colour_manual(values = scale_colours)+
    geom_point()+ggtitle("CCLE proteomics")+   facet_wrap("lineage", nrow = 3, ncol = 5 )+
    # strip.background parameter of theme
    # function is used to remove the facet wrap box
    # element_blank() makes the box background blank
    theme_bw() +
    theme(        strip.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA))
p2 <- CCLE_RNA_seq %>% t() %>% as.data.frame() %>% 
    rownames_to_column("stripped_cell_line_name") %>% 
    left_join(sample_info[,c("lineage","stripped_cell_line_name","lineage_molecular_subtype")]) %>% 
    subset(lineage %in% interesting_lineage) %>% 
    mutate(lineage= str_to_title(lineage), lineage_molecular_subtype = if_else(lineage == "Breast",lineage_molecular_subtype,NA_character_)) %>% 
    #pivot_longer(-ID, names_to = "Cell_line", values_to = "Protein_Abundance") %>% 
    add_count(lineage) %>% subset(n>4) %>% 
    ggplot(aes(x = BRD4, y = LOXL2, colour = lineage_molecular_subtype))+
    scale_colour_manual(values = scale_colours)+
    geom_point()+ggtitle("CCLE transcriptomics")+
    theme_bw() +
    theme(        strip.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA))+
    facet_wrap("lineage", nrow = 3, ncol = 5 )

###S3 C###
library(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)
# view and print all cancer studies
studies <- getCancerStudies(mycgds)

# select breast cancer study from query menu name
cancer_study <- "Breast Invasive Carcinoma (TCGA, Firehose Legacy)"
study_id <- studies[studies$name==cancer_study, "cancer_study_id"]
# view and print all patient cases
cases <- getCaseLists(mycgds, study_id)
# view and print all genetic profiles
profiles <- getGeneticProfiles(mycgds, study_id)
# select mass spec cases
case_id <- cases[cases$case_list_name=="Samples with protein data (Mass Spec)", "case_list_id"]
profile_id <- profiles[profiles$genetic_profile_name=="Protein levels (mass spectrometry by CPTAC)", "genetic_profile_id"]
# get genetic profile data from PAM50
# get genetic profile data from PAM50
genes <- c("UBE2T", "BIRC5", "NUF2", "CDC6", "CCNB1", "TYMS", "MYBL2", "CEP55", "MELK",
           "NDC80", "RRM2", "UBE2C", "CENPF", "PTTG1", "EXO1", "ORC6L", "ANLN", "CCNE1", "CDC20",
           "MKI67", "KIF2C", "ACTR3B", "MYC", "EGFR", "KRT5", "PHGDH", "CDH3", "MIA", "KRT17", "FOXC1",
           "SFRP1", "KRT14", "ESR1", "SLC39A6", "BAG1", "MAPT", "PGR", "CXXC5", "MLPH", "BCL2", "MDM2",
           "NAT1", "FOXA1", "BLVRA", "MMP11", "GPR160", "FGFR4", "GRB7", "TMEM45B", "ERBB2","LOXL2","BRD4")
GOI <- c("LOXL2", "BRD4","ESR1","ERBB2", "PGR",#"MYBL2","FOXM1"
         "VIM","EPCAM")
scale_colours <- c(LumA = "#FF9EE2",
                   LumB = "#FF68B0",
                   Her2 = "#F62983",
                   Basal = "#BD0055")
ms_data <- getProfileData(mycgds, GOI, profile_id, case_id)
# filter out uninformative genes
# fetch TCGA subtypes
# URL: https://tcgadata.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt
subtypes <- read.table(here::here("Project_Datasets", "BRCA.547.PAM50.SigClust.Subtypes.txt"), sep="\t", header=TRUE,
                       colClasses="character") %>% 
    mutate(Sample = Sample %>% 
               stringr::str_sub(start = 1, end=15) %>% 
               str_replace_all("-","."))
ms_data %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample") %>% 
    left_join(subtypes)  %>%
    ggplot(aes(x= BRD4,y = LOXL2, colour  = PAM50 ))+
    geom_point(size = 3)+
    theme_bw()+
    scale_colour_manual(values = scale_colours)+
    ggtitle("Patient Study MS proteomics")

### F4 ###
####Panel B ####

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
# install.packages('rlang')
get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}

Operreta_Objects_Population_Nuclei <- read.delim(here::here("Datasets","210812_H3P.txt"),skip = 9,header = T)%>% 
    janitor::clean_names()%>% dplyr::select(matches("ean|column|std_dev")) %>% 
    set_names(c("Sample","Mean_488","SD_488")) %>% 
    #subset(SD_488<(Mean_488/6)) %>% 
    mutate(Condition = case_when(
        Sample %in% 1:6 ~"Ctrl",
        Sample %in% 7:12~"LOXL2KD",
        TRUE ~ "Mistake"),
        Replicate = case_when(
            Sample %in% 1:2 ~"Cltr_1",
            Sample %in% 3:4~"Cltr_2",
            Sample %in% 5:6 ~"Cltr_3",
            Sample %in% 7:8~"LOXL2KD_1",
            Sample %in% 9:10 ~"LOXL2KD_2",
            Sample %in% 11:12~"LOXL2KD_3",
            TRUE ~ "Mistake"),
        Sample = Replicate,
        linear_Mean_488 = Mean_488,
        Mean_488 = log2(Mean_488))
lowest_replicate <- Operreta_Objects_Population_Nuclei$Sample %>% table() %>% min 
Operreta_Objects_Population_Nuclei <- Operreta_Objects_Population_Nuclei %>% 
    group_by(Sample) %>% 
    sample_n(lowest_replicate, replace = F ) %>% 
    ungroup() %>% 
    #pull(Sample) %>% table()
    #group_by(Replicate) %>% 
    #mutate(Mean_488 = mean(Mean_488)) %>% 
    #ungroup() %>% 
    group_by(Condition) %>% 
    mutate(median_Condition = median(Mean_488)) %>% 
    ungroup() %>% 
    group_by(Sample) %>% 
    mutate(median_Sample = median(Mean_488)) %>% 
    ungroup() %>% 
    mutate(Norm_median =Mean_488+(median_Condition-median_Sample),
           SD_total = sd(Norm_median),
           all_groups_median = median(Norm_median),
           Sample = factor(Sample, levels = unique(Sample)),
           Mitotic = if_else(Norm_median> (all_groups_median+3*SD_total),T,F)) %>% 
    group_by(Sample) %>% 
    mutate(Density = get_density(Norm_median,Norm_median, n = 100)) %>% 
    ungroup()

Operreta_Objects_Population_Nuclei %>% 
    #group_by(Sample) %>% 
    #head(10000) %>% 
    ggplot(aes( x = Condition, y = Norm_median, colour = Density, fill = Mitotic))+
    geom_jitter(width = 0.2,shape=21, alpha = 0.5)+ 
    scale_color_gradient2(low = "grey50", high = "white", mid = "grey20")+
    theme_bw()+
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    #scale_alpha_manual(values = c(0.1,1))+
    scale_fill_manual(values =  c("#CFD1D2","black"))+
    # facet_wrap(~Condition, scales = "free_x") +
    ggtitle("H3 Derived Mitotic cells in Ctrl and LOXL2-KD")
####Panel D####


####Panel F####

Achilles_file <- fs::path_rel(here::here("Project_Datasets","Achilles_gene_effect.csv")) #Raw from https://depmap.org/portal/download/ 2/11/2020
Achilles <- inner_join(sample_info[,1:2],read.csv(Achilles_file), by = "DepMap_ID")[,-1] %>%
    column_to_rownames(var = "stripped_cell_line_name") %>% t() %>%
    magrittr::set_rownames(str_match(rownames(.), "^([:graph:]*?)\\.")[,2])

plan(multisession, workers = 5)

check_diff_sensitivity <- function(Gene_Name){
    #Gene_Name = "LOXL2"
    Cell_lines<- extreme_cell_lines[[Gene_Name]]
    top_Achilles <- Achilles[,colnames(Achilles) %in% Cell_lines$top_cell_lines]
    btm_Achilles <- Achilles[,colnames(Achilles) %in% Cell_lines$btm_cell_lines]
    testing<-list(p_value = furrr::future_map_dbl(.x =rownames(Achilles),
                                                  ~t.test(top_Achilles[.x,],
                                                          btm_Achilles[.x,]) %>% pluck("p.value")) %>% 
                      set_names(rownames(Achilles)),
                  plot_df = furrr::future_map(.x =rownames(Achilles),
                                              ~rbind( top_Achilles[.x,]  %>% as.data.frame() %>% 
                                                          pivot_longer(cols = everything(), names_to = "Cell_line", values_to = "Sensitivity") %>% 
                                                          mutate(Expression = "High"),
                                                      btm_Achilles[.x,]  %>% as.data.frame() %>% 
                                                          pivot_longer(cols = everything(), names_to = "Cell_line", values_to = "Sensitivity") %>% 
                                                          mutate(Expression = "Low")), .progress = TRUE) %>% set_names(rownames(Achilles)))
}
Diff_sens<- map(c("LOXL2","BRD4"), check_diff_sensitivity) %>% set_names(c("LOXL2","BRD4"))

mean_diff<-Diff_sens %>% map(.x = .,~.x[["plot_df"]] %>% 
                                 imap_dfr(.x = .,~data.frame(Gene = .y,
                                                             Diff = diff(c(.x %>% subset(Expression == "High") %>% pull(Sensitivity) %>% mean(na.rm = T),
                                                                           .x %>% subset(Expression == "Low") %>% pull(Sensitivity) %>% mean(na.rm = T))),
                                                             Base_essentiality = .x$Sensitivity %>% median(na.rm = T)
                                 )))
plan(sequential) 
Diff_sens_LOXL_adj_pval <- Diff_sens$LOXL2$p_value %>% p.adjust(method="BH") %>% sort() 
interested_Genes <- c(Diff_sens$LOXL2$p_value %>% names() %>% str_subset("^MED[:digit:]{1,}"), "MYBL2", "MYC", "GMPS", "PAICS","B-Myb","BRD4","MDM2", "TP53")

Diff_sens_adj<-Diff_sens %>% map(.x = ., p.adjust,method = "BH" )
mean_diff$LOXL2 %>% inner_join(enframe(Diff_sens$LOXL2$p_value,"Gene","pval")) %>% 
    ggplot(aes(x = -Diff, y = -log10(pval), label = Gene))+
    geom_point()+theme_bw()+
    ggrepel::geom_text_repel(data = . %>% subset(Gene %in% interested_Genes))

### SF6###
####Panel F####

for( i in interested_Genes){
    Diff_sens$LOXL2$plot_df[[i]] %>%
        
        ggplot(aes(x = Expression, y= Sensitivity, colour = Expression))+
        geom_boxplot() + geom_point(alpha = 0.5)+ theme_bw()+
        scale_colour_manual(values = c(High = "#F62983", Low = "Black"))+
        ggtitle(glue::glue("Sensitivity to {i} KO"), subtitle =
                    glue::glue("Cell lines wtih Highest and Lowest LOXL2 transcript expression, \nadj_pvalue = {round(Diff_sens_LOXL_adj_pval[i],4)},"))
    ggsave(filename = here::here("Project_Output",glue::glue("{i} LOXL2 Diff_Essentiality.svg" )))
    
}

#### END ####
