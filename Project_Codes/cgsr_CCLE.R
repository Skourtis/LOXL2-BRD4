#####cgdsr####
pacman::p_load("tidyverse","cgdsr", "furrr","clusterProfiler","enrichplot","ggupset",magrittr,eulerr,data.table,ComplexUpset,ComplexHeatmap)
#library(cgdsr)
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
ms_data <- ms_data %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample") %>% 
    left_join(subtypes)  %>%
  pivot_longer(cols = c(BRD4,LOXL2, ERBB2, ESR1, PGR,# MYBL2,FOXM1
                        ),
                                         names_to = "Gene",
                                         values_to = "MS_Abundance") %>% 
  mutate(Gene = factor(Gene,levels = c("BRD4","LOXL2","ERBB2","ESR1","PGR")))
tp <- unique(ms_data[,c("Gene")])
tp$MS_Abundance <- tp$PAM50 <- 1 
ggplot(ms_data, aes(x = PAM50, y = MS_Abundance#, colour = PAM50
                    ))+
    geom_boxplot()+
    geom_point(alpha = 0.2)+
  # geom_rect(data = subset(tp,Gene %in% c('BRD4',"LOXL2")),fill = "#D54FA4",xmin = -Inf,xmax = Inf, 
  #           ymin = -Inf,ymax = Inf,alpha = 0.1) + 
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
    facet_wrap("Gene", ncol = 2)+
  ggtitle("Protein Abundances", subtitle = "74 TCGA Samples MS CPTAC")

Testing_kruskal_proteins <- function(Gene_name){
for_test<-ms_data %>% subset(Gene == Gene_name) %>%
    mutate(PAM50 = factor(PAM50,levels = unique(PAM50)))
kruskal.test(MS_Abundance ~ PAM50, data = for_test) %>% pluck("p.value")}
map_dbl(GOI,Testing_kruskal_proteins) %>% set_names(GOI)
###corr #####
BRD4_LOXL2 <- ms_data %>% subset(Gene %in% c("BRD4","LOXL2")) %>% 
    pivot_wider(names_from = "Gene",values_from = "MS_Abundance")

    ggplot(BRD4_LOXL2, aes(x = BRD4, y = LOXL2, colour = PAM50))+
    geom_point(size = 2) + 
        ggtitle("Protein Abundance",
                subtitle = "FireHose Legacy TCGA")
    BRD4_LOXL2 %$% cor.test(.$LOXL2,
                       .$BRD4, method = "spearman")
    map(unique(BRD4_LOXL2$PAM50),~BRD4_LOXL2 %>% subset(PAM50 == .x) %$% cor.test(.$LOXL2,
                            .$BRD4, method = "pearson")) %>% set_names()
#####essentiality###
    
    
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
    #####LOXL2- Protein Abudnaces### 
    CCLE_Protein %>% as.data.frame() %>% rownames_to_column("GeneName") %>% pivot_longer(-GeneName, names_to = "Cell_line", values_to = "Abundance") %>% 
        subset(GeneName == "Q9Y4K0") %>% mutate(Cell_line_interest = if_else(Cell_line %in% c("MDAMB468","MDAMB231","BT549"),T,F)) %>% 
        ggplot(aes(x = fct_reorder(Cell_line, desc(Abundance)), y= Abundance, fill = Cell_line_interest,  label = Cell_line ))+
        geom_col()+theme_bw()+ggrepel::geom_label_repel(data = . %>% subset(Cell_line_interest == T))+
        theme(legend.position = "none",axis.ticks.x = element_blank(), axis.text.x = element_blank())+
        ggtitle("LOXL2 Protein Abundance across lineage Cell lines" ,subtitle = "CCLE Proteomics Gygi (2019)")
    
     ######
    testing <- CCLE_Protein %>% as.data.frame() %>% rownames_to_column("GeneName") %>% 
      pivot_longer(-GeneName, names_to = "stripped_cell_line_name", values_to = "Abundance") %>% 
        left_join(HUMAN_9606_idmapping %>% subset(Type == "Gene_Name"), by =c("GeneName"= "Uniprot")) %>% 
        subset(ID %in% c("ESR1","VIM","BRD4","LOXL2","MYBL2","FOXM1")) %>% 
      mutate(Cell_line_interest = if_else(stripped_cell_line_name %in% c("MDAMB468","MDAMB231","BT549"),"T","F")) %>% 
        left_join(sample_info[,c("stripped_cell_line_name","lineage", "lineage_sub_subtype","lineage_molecular_subtype")])%>%
        mutate(Breast_cell_line =  if_else(lineage == "breast", T,F))%>%
      subset(lineage == "breast" )%>% 
      mutate(ID = factor(ID, c("BRD4","LOXL2","FOXM1","MYBL2","VIM")))
        ggplot(testing,aes(x = ID, y= Abundance,label = stripped_cell_line_name))+
        geom_boxplot(alpha= 0.1)+
   
          theme_bw()+
          theme(axis.line = element_line(color='black'),
                plot.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank())+
         # ggrepel::geom_label_repel(data = . %>% subset(Cell_line_interest == "T"))+
          geom_point(data = . %>% subset(Cell_line_interest == "T"),aes( colour =stripped_cell_line_name  ), size = 2)+ 
          # scale_colour_manual(values = c("T" = "#D64FA4",
          #                                "F" = "grey20"))+
        #theme(#legend.position = "none",
         #     axis.ticks.x = element_blank(), axis.text.x = element_blank())+
         
        ggtitle("LOXL2-BRD4 Protein Abundance across lineage Cell lines" ,subtitle = "CCLE Proteomics Gygi (2019)")
    
      ####
        
        CCLE_Protein %>% as.data.frame() %>% rownames_to_column("GeneName") %>% pivot_longer(-GeneName, names_to = "stripped_cell_line_name", values_to = "Abundance") %>% 
            left_join(HUMAN_9606_idmapping %>% subset(Type == "Gene_Name"), by =c("GeneName"= "Uniprot")) %>% 
        subset(ID %in% c("ESR1","VIM","BRD4","LOXL2","PGR","ERBB2","EPCAM","MYBL2","FOXM1")) %>% 
            mutate(Cell_line_interest = if_else(stripped_cell_line_name %in% c("MDAMB468","MDAMB231","BT549"),T,F)) %>% 
            left_join(sample_info[,c("stripped_cell_line_name","lineage", "lineage_sub_subtype","lineage_molecular_subtype")])%>%
            mutate(Breast_cell_line =  if_else(lineage == "breast", T,F))%>%
            subset(lineage == "breast")%>%
            ggplot(aes(x = lineage_molecular_subtype, y= Abundance,label = stripped_cell_line_name ))+
            geom_boxplot(alpha= 0.1)+
            geom_point(aes(alpha = Cell_line_interest))+theme_bw()+
            # ggrepel::geom_label_repel(data = . %>% subset(Cell_line_interest == T))+
            #theme(#legend.position = "none",
            #     axis.ticks.x = element_blank(), axis.text.x = element_blank())+
            facet_wrap("ID")+
            ggtitle("LOXL2-BRD4 Protein Abundance across lineage breast subtypes" ,
                    subtitle = "CCLE Proteomics Gygi (2019)\n          VIM         BRD4        LOXL2        ERBB2        EPCAM        MYBL2        FOXM1 \n
    pval = 0.0017       0.02580       0.00030       0.0047       0.0075       0.0648       0.3546")
      krustal_df <-  CCLE_Protein %>% as.data.frame() %>% rownames_to_column("GeneName") %>% pivot_longer(-GeneName, names_to = "stripped_cell_line_name", values_to = "Abundance") %>% 
        left_join(HUMAN_9606_idmapping %>% subset(Type == "Gene_Name"), by =c("GeneName"= "Uniprot")) %>% 
        subset(ID %in% c("ESR1","VIM","BRD4","LOXL2","PGR","ERBB2","EPCAM","MYBL2","FOXM1")) %>% 
        mutate(Cell_line_interest = if_else(stripped_cell_line_name %in% c("MDAMB468","MDAMB231","BT549"),T,F)) %>% 
        left_join(sample_info[,c("stripped_cell_line_name","lineage", "lineage_sub_subtype","lineage_molecular_subtype")])%>%
        mutate(Breast_cell_line =  if_else(lineage == "breast", T,F))%>%
        subset(lineage == "breast")
          Testing_kruskal_proteins <- function(Gene_name){
            #Gene_name = "VIM"
          for_test<-krustal_df %>% subset(ID == Gene_name) %>%
            mutate(lineage_molecular_subtype = factor(lineage_molecular_subtype,levels = unique(lineage_molecular_subtype)))
          kruskal.test(Abundance ~ lineage_molecular_subtype, data = for_test) %>% pluck("p.value")}
          krustal_df_pval <-   map_dbl(c("VIM","BRD4","LOXL2","ERBB2","EPCAM","MYBL2","FOXM1"),Testing_kruskal_proteins) %>% 
            set_names(c("VIM","BRD4","LOXL2","ERBB2","EPCAM","MYBL2","FOXM1"))
        ####
    
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
    HUMAN_9606_idmapping <-
    readr::read_tsv(fs::path_rel(here::here("Project_Datasets","HUMAN_9606_idmapping.dat")),
                    col_names = FALSE) %>%
        setNames(c("Uniprot", "Type", "ID"))


sample_info <-  read.csv(fs::path_rel(here::here("Project_Datasets","sample_info.csv")), stringsAsFactors = FALSE)
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
extract_top_btm_protein_expressing_cell_lines<-function(Gene_Name){
  #Gene_Name = "Q9Y4K0"
  ordered <- CCLE_Protein[Gene_Name,]%>%unlist() %>%  na.omit() %>%  #.[names(.) %in% (sample_info %>% subset(lineage == "breast") %>% pull(stripped_cell_line_name))] %>% 
    sort(decreasing = T)
  list(
    top_cell_lines = head(ordered,round(length(ordered)/10)) %>% names(),
    btm_cell_lines = tail(ordered,round(length(ordered)/10)) %>% names()) 
}

extreme_cell_lines<- map(c("LOXL2", "BRD4"),extract_top_btm_expressing_cell_lines) %>%
  set_names(c("LOXL2", "BRD4"))
extreme_protein_cell_lines<- map(c("Q9Y4K0", "O60885"),extract_top_btm_protein_expressing_cell_lines) %>%
  set_names(c("LOXL2", "BRD4"))
CCLE_RNA_seq %>% t()%>% as.data.frame() %>% 
  rownames_to_column("stripped_cell_line_name") %>% 
  left_join(sample_info[,c("stripped_cell_line_name","lineage")]) %>% 
  add_count(lineage) %>% 
  subset(n>15 & stripped_cell_line_name %in% c(extreme_cell_lines$LOXL2$top_cell_lines,extreme_cell_lines$LOXL2$btm_cell_lines)) %>% 
  ggplot(aes(x = LOXL2,y = BRD4, colour = lineage))+
  geom_point()+
  ylim(c(-4,4))+
  ggtitle("No significant Diff between LOXL2 high and LOXL2 Low RNA expressing cell lines p-value = 0.3264")
High_low_LOXL2 <- CCLE_RNA_seq %>% t()%>% as.data.frame() %>% 
  rownames_to_column("stripped_cell_line_name") %>% 
  left_join(sample_info[,c("stripped_cell_line_name","lineage")]) %>% 
  add_count(lineage) %>% 
  subset(n>15 & stripped_cell_line_name %in% c(extreme_cell_lines$LOXL2$top_cell_lines,extreme_cell_lines$LOXL2$btm_cell_lines)) %>% 
  dplyr::select(stripped_cell_line_name, BRD4) %>% 
  mutate(LOXL2_Expression = if_else(stripped_cell_line_name %in% extreme_cell_lines$LOXL2$top_cell_lines, "High","Low" ) )
  t.test(High_low_LOXL2 %>% subset(LOXL2_Expression == "High") %>% pull(BRD4),
         High_low_LOXL2 %>% subset(LOXL2_Expression == "Low") %>% pull(BRD4))
  #https://depmap.org/portal/download/all/?release=DepMap+Public+20Q3&file=Achilles_gene_effect.csv
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


Diff_sens$BRD4$p_value %>% p.adjust(method="BH") %>% sort() %>% head(10)
interested_Genes <- c(Diff_sens$LOXL2$p_value %>% names() %>% str_subset("^MED[:digit:]{1,}"), "MYBL2", "MYC", "GMPS", "PAICS","B-Myb","BRD4","MDM2", "TP53")

for( i in interested_Genes){
  Diff_sens$LOXL2$plot_df[[i]] %>%

    ggplot(aes(x = Expression, y= Sensitivity, colour = Expression))+
    geom_boxplot() + geom_point(alpha = 0.5)+ theme_bw()+
    scale_colour_manual(values = c(High = "#F62983", Low = "Black"))+
    ggtitle(glue::glue("Sensitivity to {i} KO"), subtitle =
              glue::glue("Cell lines wtih Highest and Lowest LOXL2 transcript expression, \nadj_pvalue = {round(Diff_sens_LOXL_adj_pval[i],4)},"))
  ggsave(filename = here::here("Project_Output",glue::glue("{i} LOXL2 Diff_Essentiality.svg" )))

}
Diff_sens$LOXL2$plot_df$BRD4 %>% 
    ggplot(aes(x = Expression, y= Sensitivity, colour = Expression))+
    geom_boxplot() + geom_point()+ theme_bw()+
    ggtitle("Sensitivity to BRD4 KO", subtitle = "Cell lines wtih Highest and Lowest LOXL2 transcript expression")
Diff_sens$LOXL2$plot_df$MYC %>% 
    ggplot(aes(x = Expression, y= Sensitivity, colour = Expression))+
    geom_boxplot() + geom_point()+ theme_bw()+
    ggtitle("Sensitivity to MYC KO", subtitle = "Cell lines wtih Highest and Lowest LOXL2 transcript expression")

Diff_sens$BRD4$plot_df$SCD %>% 
    ggplot(aes(x = Expression, y= Sensitivity, colour = Expression))+
    geom_boxplot() + geom_point()+
    ggtitle("Sensitivity to BRD4 KO", subtitle = "Cell lines with Highest and Lowest BRD4 transcript expression, pvalue = 0.033, adjusted = 0.45")
Diff_sens_adj<-Diff_sens %>% map(.x = ., p.adjust,method = "BH" )
plan(sequential)

mean_diff$LOXL2 %>% pull(Diff,Gene) %>% sort(decreasing = TRUE)
ego3 <- gseGO(geneList     = mean_diff$BRD4 %>% pull(Diff,Gene) %>% sort(decreasing = TRUE),
              OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
              ont          = "ALL",
              #nPerm        = 1000,
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              keyType = "SYMBOL",
              verbose      = FALSE)
p2 <- dotplot(ego3, showCategory=30) + ggtitle("dotplot for BRD4 co-essentaiality GSEA")


enrichplot::upsetplot(ego3,n = 10)

##### Cell Sensitivity to BRD4 inh ####
sec_screen_info <- read_csv(here::here("Project_Datasets","secondary-screen-replicate-collapsed-treatment-info.csv"))
BRD4_inh <- sec_screen_info %>% subset(str_detect(target, "BRD4")) %>% pull("broad_id") %>% paste(., collapse = "|")
sec_screen <- fread(here::here("Project_Datasets","secondary-screen-replicate-collapsed-logfold-change.csv")) %>% 
    inner_join(.,sample_info[,1:2], by = c("V1" = "DepMap_ID"))%>%
    #column_to_rownames(var = "stripped_cell_line_name")  %>% 
    dplyr::select(matches(BRD4_inh)|contains("stripped"))  %>% 
    pivot_longer(cols = -stripped_cell_line_name, names_to = "column_name", values_to = "Sensitivity") %>%
    left_join(sec_screen_info)
sec_screen %>% subset(str_detect(name,"JQ") & dose == 10) %>%
    group_by(stripped_cell_line_name) %>%
    summarise(Mean_Sensitivity = median(Sensitivity)) %>%
    inner_join(CCLE_Protein %>% t() %>% .[,"Q9Y4K0"] %>% as.data.frame()%>% set_names( c("LOXL2")) %>% rownames_to_column("stripped_cell_line_name") )%>%
    left_join(sample_info[,c("stripped_cell_line_name","lineage", "lineage_sub_subtype","lineage_molecular_subtype")])%>%
    ggplot(aes(x = Mean_Sensitivity,y = LOXL2))+
    geom_point(alpha = 0.5)+
    geom_point(data = . %>% subset(lineage == "breast"), aes( colour = lineage_molecular_subtype))+
    facet_wrap("lineage")+
    theme_bw()+
    ggtitle("LOXL2 proteomic Abundance with JQ1 sensitivity")

breast_JQ1_LOXL2<- sec_screen %>% subset(str_detect(name,"JQ") & dose == 10) %>%
    group_by(stripped_cell_line_name) %>%
    summarise(Mean_Sensitivity = median(Sensitivity)) %>%
    inner_join(CCLE_Protein %>% t() %>% .[,"Q9Y4K0"] %>% as.data.frame()%>% set_names( c("LOXL2")) %>% rownames_to_column("stripped_cell_line_name") )%>%
    left_join(sample_info[,c("stripped_cell_line_name","lineage", "lineage_sub_subtype","lineage_molecular_subtype")])%>%
    subset(lineage == "breast")

cor.test(breast_JQ1_LOXL2$Mean_Sensitivity,breast_JQ1_LOXL2$LOXL2)
compounds <- sec_screen_BRDi$name %>% unique() %>% .[c(7,1,3,5,6)] %>% set_names(.,.)
extreme_cell_lines_LOXL2<-rbind(data.frame(stripped_cell_line_name = extreme_cell_lines$LOXL2$top_cell_lines,
                                           Expression = "High"),
                                data.frame(stripped_cell_line_name = extreme_cell_lines$LOXL2$btm_cell_lines,
                                           Expression = "Low"))
sec_screen_BRDi <- sec_screen %>% inner_join(extreme_cell_lines_LOXL2) 
    #mutate(dose = factor(dose,unique(dose))) %>%



t.test_comparison <- function(compound_name){

  test_JQ1 <- sec_screen_BRDi %>% subset(name == compound_name)
  max_dose <- test_JQ1$dose %>% max()
  test_JQ1 <- test_JQ1 %>% subset(dose ==max_dose)
  t.test(test_JQ1 %>% subset(Expression == "High") %>% pull(Sensitivity),
         test_JQ1 %>% subset(Expression == "Low") %>% pull(Sensitivity)) %>% .$p.value
}
p_val <- map_dbl(compounds,t.test_comparison) %>% p.adjust(method = "fdr") %>% round(3)
sec_screen_BRDi %>%
  subset((name %in% compounds) & between(log2(dose), -4,4))%>% 
  mutate(name = factor(name, levels = compounds)) %>% 
  group_by(stripped_cell_line_name) %>%
  ggplot(aes(x = log2(dose), y = Sensitivity, colour = Expression, group = stripped_cell_line_name))+
  geom_line(alpha = 0.1, )+
  stat_summary(aes(group=Expression,colour = Expression), fun=median,  geom="line")+
  geom_point(alpha = 0.2) + 
  facet_grid("name")+
  scale_colour_manual(values = c(High = "#F62983", Low = "Black"))+
  ggtitle("BRD4i Sensitivity", 
          subtitle = "Between High and Low RNA expressing LOXL2 Cell Lines")+
  labs(caption =paste(paste(compounds,p_val, sep = "_"), collapse = "     "))
#PRoteins

extreme_prot_cell_lines_LOXL2<-rbind(data.frame(stripped_cell_line_name = extreme_protein_cell_lines$LOXL2$top_cell_lines,
                                           Expression = "High"),
                                data.frame(stripped_cell_line_name = extreme_protein_cell_lines$LOXL2$btm_cell_lines,
                                           Expression = "Low"))
sec_screen_BRDi <- sec_screen %>% inner_join(extreme_prot_cell_lines_LOXL2) 
#mutate(dose = factor(dose,unique(dose))) %>%



p_val <- map_dbl(compounds,t.test_comparison) %>% p.adjust(method = "fdr") %>% round(3)
sec_screen_BRDi %>%
  subset((name %in% compounds) & between(log2(dose), -4,4))%>% 
  mutate(name = factor(name, levels = compounds)) %>% 
  group_by(stripped_cell_line_name) %>%
  ggplot(aes(x = log2(dose), y = Sensitivity, colour = Expression, group = stripped_cell_line_name))+
  geom_line(alpha = 0.1, )+
  stat_summary(aes(group=Expression,colour = Expression), fun=median,  geom="line")+
  geom_point(alpha = 0.2) + 
  facet_grid("name")+
  scale_colour_manual(values = c(High = "#F62983", Low = "Black"))+
  ggtitle("BRD4i Sensitivity", 
          subtitle = "Between High and Low Prot expressing LOXL2 Cell Lines")+
  labs(caption =paste(paste(compounds,p_val, sep = "_"), collapse = "     "))
####DREAM###
REACTOME<-
    read_csv2(here::here("Project_Datasets","Reactome_pathways.csv")) %>% .[,-1]
DREAM_targets <- read.delim(here::here("Project_Datasets","geneset_DREAM.txt"))[-1,] %>% unlist() # https://www.gsea-msigdb.org/gsea/msigdb/cards/FISCHER_DREAM_TARGETS
E2F_targets_DREAM <- REACTOME %>% subset(str_detect(Reactome_Pathway,"DREAM")) %>% pull(Gene_name)
set.seed(1234)
library("scales")
tesing <- mean_diff$LOXL2 %>% inner_join(enframe(Diff_sens$LOXL2$p_value,"Gene","pval")) %>% 
                                 mutate(adj = p.adjust(pval,method="BH"),
                                        Significant_type = case_when(
                                          adj<0.05 & (-Diff)>0.15 ~"Less_relative_proliferation",
                                          adj<0.05 & (-Diff)<(-0.15 )~"More_relative_proliferation",
                                          TRUE~"No_Sign_Diff"),
                                        
                                        Significant = if_else(abs(Diff) > 0.4 & adj<0.05,T,F))%>% arrange(Diff) %>% 
   mutate(Gene = if_else(Gene == "MYBL2","B-Myb",Gene),
          thick2 = case_when(
            between(adj,0.2,1)~1,
            between(adj,0.05,0.2)~2,
            between(adj,0.005,0.05)~3,
            between(adj,0.0005,0.005)~4,
            TRUE~5,
          ))
    #                                         Dream = if_else(Gene %in% DREAM_targets, T,F)) %>%
    ggplot(tesing, aes(x = -Diff, y = -log10(pval), label = Gene,  colour = Significant_type))+
      #  scale_colour_gradient2(high = "#EA2A8A", low = "#FFCECE", mid = "#F39AC2")+
       # scale_colour_gradientn(colours = rev(c("#EA2A8A","#F39AC2","#FFCECE")), 
       #                      values = rescale(c(tesing$Diff %>% min(),0,tesing$Diff %>% max())),
                            # guide = "colorbar", limits=c(tesing$Diff %>% min(),tesing$Diff %>% max()))+
      geom_point(alpha = 0.3, aes(size =thick2))+
      scale_colour_manual(values = c(Less_relative_proliferation = "#000000",
                                     More_relative_proliferation = "#000000",
                                     No_Sign_Diff = "#D1D3D4"))+
      scale_size(breaks = c(1,2,3,4,5),labels = c('0-20000','20000-40000','40000+','5',"4"))+
    #geom_point(data = . %>% subset(Dream == TRUE),alpha = 1)+
        
    ggrepel::geom_text_repel(data = . %>% subset(Gene %in% interested_Genes & (Significant_type == "No_Sign_Diff")), fontface = 'bold',seed = 42, max.overlaps = 100, colour = "#B3B3B3")+
      ggrepel::geom_text_repel(data = . %>% subset(Gene %in% interested_Genes & (Significant_type != "No_Sign_Diff" )) , fontface = 'bold',seed = 42, max.overlaps = 100, colour = "#000000")+ 
      
      # ggrepel::geom_text_repel(data = . %>% subset(Significant == T & !(Gene %in% interested_Genes)),seed = 42, max.overlaps = 10)+ 
      
  theme_bw()+
  xlab("Essentiality Difference Between LOXL2 High and Low")+
    ggtitle("Genes with Diff Sensitivy LOXL2 High and Low",
            subtitle = "LOXL2 low cell lines proliferate less with gene KOs on the right, compared to LOXL2 High cell lines")
ggsave("Diff_Essentiality_LOXL2.pdf")

mean_diff$BRD4 %>% inner_join(enframe(Diff_sens$BRD4$p_value,"Gene","pval")) %>% 
  mutate(adj = p.adjust(pval,method="BH"),
         Significant = if_else(abs(Diff) > 0.1 & adj<0.05,T,F))%>% arrange(Diff) %>% 
  # mutate(Gene = factor(Gene,levels = unique(Gene)),
  #                                         Dream = if_else(Gene %in% DREAM_targets, T,F)) %>%
  ggplot(aes(x = -Diff, y = -log10(pval), label = Gene, colour = Base_essentiality))+
  scale_colour_gradient2(low = "red", high = "blue", mid = "#E6DAC1")+
  geom_point(aes(alpha = Significant))+
  #geom_point(data = . %>% subset(Dream == TRUE),alpha = 1)+
  ggrepel::geom_text_repel(data = . %>% subset(Significant ==T),seed = 42, max.overlaps = 100)+ 
  ggrepel::geom_text_repel(data = . %>% subset(Significant ==T),seed = 42, max.overlaps = 100)+ 
  ggrepel::geom_text_repel(data = . %>% subset(Significant ==T),seed = 42, max.overlaps = 100)+ 
  theme_bw()+
  ggtitle("Genes with Diff Sensitivy BRD4")

CCLE_RNA_seq %>% t() %>% 
    .[,c("LOXL2","BRD4")] %>% as.data.frame()%>% set_names( c("LOXL2","BRD4")) %>% 
    rownames_to_column("stripped_cell_line_name") %>%
    left_join(sample_info[,c("stripped_cell_line_name","lineage", "lineage_sub_subtype","lineage_molecular_subtype")])%>%
    mutate(is_breast = if_else(lineage=="breast",T,F))%>%
    na.omit()%>%
    subset(!str_detect(lineage,"neered")) %>%
    ggplot(aes(x = BRD4,y = LOXL2))+
    geom_point(alpha = 0.1)+
    geom_point(data = . %>% subset(lineage == "breast"), colour = "red")+
    scale_colour_manual(values= c("grey50","red"))+
    #facet_wrap("lineage")+
    theme_bw()+
    ggtitle("LOXL2 Transcriptomic BRD4 correlation", subtitle = "Breast correlation = 0.25 pvalue 0.08, non-breast = -0.12 pvalue 0.000183")

lists <- list.files(here::here("Project_Datasets","new_peaks")) %>% str_subset("\\.genes\\.") %>% 
  set_names(.,.) %>% 
  purrr::map(.x = .,~read.delim(here::here("Project_Datasets","new_peaks",.x)) %>% unlist(use.names = F))
s4 <- list(CT_AB = lists$CT_AB.Homer.genes.txt,
           CT_BB= lists$CT_BB.Homer.genes.txt)

plot(eulerr::euler(s4, shape = "ellipse"), quantities = TRUE)

s4 <- list(LX_AB = lists$LX_AB.Homer.genes.txt,
           LX_BB= lists$LX_BB.Homer.genes.txt)

plot(euler(s4, shape = "ellipse"), quantities = TRUE)
s4 <- list(LX_AB = lists$LX_AB.Homer.genes.txt,
           LX_BB= lists$LX_BB.Homer.genes.txt,
           CT_AB =lists$CT_AB.Homer.genes.txt, 
             CT_BB =lists$CT_AB.Homer.genes.txt )

plot(euler(s4, shape = "ellipse"), quantities = TRUE)

lists <- list.files(here::here("Project_Datasets","new_peaks")) %>% str_subset("\\.promotergenes\\.") %>% 
  set_names(.,.) %>% 
  purrr::map(.x = .,~read.delim(here::here("Project_Datasets","new_peaks",.x)) %>% unlist(use.names = F))

s4 <- lists %>% set_names(.,names(.) %>% str_remove_all("\\.[:print:]*"))
plot(eulerr::euler(s4[-c(3:4)], shape = "ellipse"), quantities = TRUE,main = "Promoter Dream Genes")
plot(eulerr::euler(s4[c(3:4)], shape = "ellipse"), quantities = TRUE,main = "Promoter Genes")

plot(venn(s4), main = "Promoter Genes")
s4 <- lists %>% set_names(.,names(.) %>% str_remove_all("\\.[:print:]*")) %>% 
  map(.x =.,~.x[.x %in% DREAM_targets])
plot(venn(s4), main = "Promoter Genes")
plot(eulerr::euler(s4, shape = "ellipse"), quantities = TRUE, main = "Promoter Genes")

plot(eulerr::euler(s4[c(3:4)], shape = "ellipse"), quantities = TRUE,main = "Promoter Dream Genes")


s4 <- s4[-c(3:4)]
in_cols  = names(s4) 
b_lists <- purrr::map(.x = s4, ~ (s4 %>% unlist() %>% unique()) %in% .x)
Gene_subset <- data.table(Genes = s4 %>% unlist() %>% unique() ) %>% 
  .[,c(in_cols) := b_lists]
# Gene_subset[,DREAM_targets:= if_else(Genes %chin% DREAM_targets,T,F)]
upset(
  Gene_subset, in_cols, width_ratio=0.1,base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=DREAM_targets))+ 
        scale_fill_manual(values=c(
        'TRUE'='pink', 'FALSE'='grey80')
    )
  ),
  min_degree=1,
)
genres_subset = in_cols
arranged = arrange_venn(Gene_subset[1:10,], sets=genres_subset)
m= make_comb_mat(arranged)
set_size(m)

(
  ggplot(arranged)
  + theme_void()
  + coord_fixed()
  + geom_venn_region(Gene_subset, sets=genres_subset, alpha=0.1)
  + geom_point(aes(x=x, y=y, color=region), size=2.5)
    + geom_venn_circle(Gene_subset, sets=genres_subset, size=1.5)
   + geom_venn_label_set(Gene_subset, sets=genres_subset, aes(label=region), outwards_adjust=2.6)
  # + geom_venn_label_region(Gene_subset, sets=genres_subset, aes(label=size), position=position_nudge(y=0.15))
  + scale_color_venn_mix(Gene_subset, sets=genres_subset, guide="none")
  + scale_fill_venn_mix(Gene_subset, sets=genres_subset, guide="none")
)

