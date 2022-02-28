##### BRD4 - LOXL2 higher in Cancer####


if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install('GenomicDataCommons')
library(GenomicDataCommons)
GenomicDataCommons::status()
BRCA_norm = cases() %>% fGenomicDataCommons::ilter(~ project.project_id=='TCGA-BRCA' &
                              samples.sample_type=='Solid Tissue Normal') %>%
    GenomicDataCommons::select(c(default_fields(cases()),'samples.sample_type')) %>%
    response_all() %>% ids()
files() %>% 
    GenomicDataCommons::filter(~ cases.project.project_id=='CPTAC-3') %>% 
    facet("cases.samples.sample_type") %>% aggregations()
map(.x = files() %>% 
          GenomicDataCommons::filter(~ cases.project.project_id=='CPTAC-3') %>% 
        available_fields() %>% str_subset("type"),~ files() %>% 
        GenomicDataCommons::filter(~ cases.project.project_id=='CPTAC-3') %>% 
        facet(.x) %>% aggregations())          
    
files() %>% 
    GenomicDataCommons::filter(~ cases.project.project_id=='CPTAC-3') %>% 
    facet("cases.disease_type") %>% aggregations()                                   
#####
Importing_CommonsData<-function(file_name,GOIs){
    read.table(file_name) %>% 
        mutate(V1 =  str_remove_all(V1,"\\.[:graph:]*$")) %>%
        subset( V1 %in% GOIs)}

GOI <- c(BRD4 = "ENSG00000141867",
         LOXL2 = "ENSG00000134013")
Import_FKPM_caseid <- function(case_ids = NULL,
                               Type = c("primary tumor","Solid Tissue Normal")){
    #case_ids  = BRCA_norm[1:10]
    #Type = "Solid Tissue Normal"
    resp = files() %>% 
        GenomicDataCommons::filter(~ cases.case_id %in% case_ids)%>% 
        GenomicDataCommons::filter( type == 'gene_expression' ) %>%
        GenomicDataCommons::filter( cases.samples.sample_type == Type ) %>%
        GenomicDataCommons::filter( analysis.workflow_type == 'HTSeq - FPKM')  %>%
        manifest()%>%
        ids()
    fnames = lapply(resp,gdcdata)
    purrr::map_dfr(.x = fnames,~Importing_CommonsData(.x,GOI)) %>%
        mutate(Tissue = Type)
}
#set.seed(42)
Solid_Normal_BRCA <-  Import_FKPM_caseid(BRCA_norm,
                                         "Solid Tissue Normal")

Primary_Tim_BRCA <- Import_FKPM_caseid(BRCA_norm,
                                        "primary tumor")
BRCA_TCGA<- rbind(Solid_Normal_BRCA,Primary_Tim_BRCA) %>% 
    inner_join(enframe(GOI, name = "Gene_Name", value = "V1")) %>%
    rename(FPKM = V2) %>%
    mutate(log10_FPKM = log10(FPKM))
BRCA_TCGA %>%    ggplot(aes(x = Tissue, y = log10_FPKM, colour = Tissue))+
    geom_boxplot()+
    geom_point(alpha = 0.1)+
    facet_grid("Gene_Name") + ggtitle("TCGA-BRCA Gene Expression")

P_value_gene <- function(Gene){
    t.test_BRCA_TCGA <- BRCA_TCGA %>% subset(Gene_Name == Gene)
    t.test(x = t.test_BRCA_TCGA %>% subset(Tissue == "primary tumor") %>% pull(log10_FPKM),
           y =  t.test_BRCA_TCGA %>% subset(Tissue == "Solid Tissue Normal") %>% pull(log10_FPKM)) %>%
        pluck("p.value") }
map_dbl(names(GOI),P_value_gene)%>% set_names(names(GOI))
