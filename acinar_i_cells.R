
# check acinar_i cells using the unmodified dataset
# Load required libraries 
required_packages <- c('Matrix', 'dplyr', 'Seurat', 'patchwork', 'SeuratData', 'ggplot2', 'SeuratDisk', 'viridis')
invisible(lapply(required_packages, library, character.only = T))

# obtain directories
# ref_dir with a folder EGAS00001004653, where the adult_pancreas.gz file can loaded
egas_nr <- "EGAS00001004653"
egas_dir <- paste0(ref_dir,egas_nr, "/")

##############################################################

#######################
### EGAS00001004653 ###
#######################

# https://www.sciencedirect.com/science/article/pii/S0016508520353993
# data from
# https://explore.data.humancellatlas.org/projects/b3938158-4e8d-4fdb-9e13-9e94270dde16/project-matrices
# cp: 2/13/25 & ap: 4/6/25 (re-loaded)

##################
# Adult Pancreas #
##################
# read seurat object with annotated cell types
# 43815 genes, 112563 cells
ap <- readRDS(paste0(egas_dir, "adult_pancreas.gz"))
dim(ap)
# [1]  43815 112563

str(ap@meta.data)
# 'data.frame':	112563 obs. of  9 variables:
#  $ nCount_RNA       : num  896 601 811 685 907 ...
#  $ nFeature_RNA     : int  655 461 641 509 671 510 364 368 701 1493 ...
#  $ sample_ID        : chr  "AFES448_FF1" "AFES448_FF1" "AFES448_FF1" "AFES448_FF1" ...
#  $ patient_ID       : chr  "AFES448" "AFES448" "AFES448" "AFES448" ...
#  $ sex              : chr  "male" "male" "male" "male" ...
#  $ age              : chr  "30" "30" "30" "30" ...
#  $ pancreas_location: chr  "head" "head" "head" "head" ...
#  $ procurement_lab  : chr  "Stanford" "Stanford" "Stanford" "Stanford" ...
#  $ Cluster          : Factor w/ 14 levels "Macrophage","Acinar-REG+",..: 5 5 4 4 4 4 4 2 13 13 ...

table(ap@meta.data$patient_ID)
# AFES365 AFES448 AGBR024  TUM_13  TUM_25  TUM_C1 
#   22288   34167   33196    5233    9712    7967 

# add module scores for acinar cell markers
acinar_markers <- c("PRSS1", "PRSS2", "CPA1", "CTRC", "CELA3A", "CELA2A", "PNLIP", "AMY2A", "AMY2B")
acinar_journal_markers <- c("PRSS1","CPA1","CPA2", "RBPJL", "FOXP2")
acinar_markers_valid <- c("PRSS1","CPA1", "RBPJL")
ap <- AddModuleScore(ap, features = list(acinar_markers), name = "Acinar_Marker_Score")
ap <- AddModuleScore(ap, features = list(acinar_journal_markers), name = "Acinar_Journal_Marker_Score")
ap <- AddModuleScore(ap, features = list(acinar_markers_valid), name = "Acinar_Markers_Valid_Score")

############################################################################################################################################################

# create dimplots for all samples combined
p1 <- DimPlot(ap, group.by="Cluster", label=T,label.size=2, repel=T)+ NoLegend()
ap[["Acinar_i"]] <- ap@meta.data$Cluster == "Acinar-i"
p2 <- DimPlot(ap, group.by="Acinar_i", label=T,label.size=2, repel=T, cols= c("gray", "red"))+ NoLegend()
p3 <- FeaturePlot(ap, features = "Acinar_Marker_Score1")&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p4 <- FeaturePlot(ap, features = c("Acinar_Journal_Marker_Score1"))&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p5 <- FeaturePlot(ap, features = c("Acinar_Markers_Valid_Score1"))&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p_patched<-  wrap_plots(p1, p2, p3,p4,p5, ncol=5) + plot_annotation(title= paste0("Cell UMAP & Acinar Markers in all adult pancreas samples")) & theme(plot.title = element_text(size = 14, face = "bold"))
ggsave(paste0(egas_dir, "plots/all_AP_acinar_cells_umap.png"), plot = p_patched, width = 70, height =15, units = "cm")


# create dotplots for all samples combined
ap$acinar_status <- as.character(ap$Cluster)
ap$acinar_status[!grepl("Acinar", ap$Cluster)] <- "Other"

p1 <- DotPlot(ap, features = acinar_markers, group.by = "acinar_status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- DotPlot(ap, features = acinar_journal_markers, group.by = "acinar_status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_patched<-  wrap_plots(p1, p2,ncol=1, guides="collect") + plot_annotation(title= paste0("Acinar Markers in all adult pancreas samples"), tag_levels=list(c("Acinar Markers", "Acinar Journal Markers"))) & theme(plot.title = element_text(size = 14, face = "bold"))
ggsave(paste0(egas_dir, "plots/all_AP_acinar_cells_dotplot.png"), plot = p_patched, width = 25, height =25, units = "cm")

############################################################################################################################################################

# create dimplots for each sample separately
p1 <- DimPlot(ap, group.by="Cluster", split.by="patient_ID", label=T,label.size=2, repel=T)+ NoLegend()
ap[["Acinar_i"]] <- ap@meta.data$Cluster == "Acinar-i"
p2 <- DimPlot(ap, group.by="Acinar_i", split.by="patient_ID",label=T,label.size=2, repel=T, cols= c("gray", "red"))+ NoLegend()
p3 <- FeaturePlot(ap, features = "Acinar_Marker_Score1", split.by="patient_ID")&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p4 <- FeaturePlot(ap, features = c("Acinar_Journal_Marker_Score1"), split.by="patient_ID")&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p5 <- FeaturePlot(ap, features = c("Acinar_Markers_Valid_Score1"), split.by="patient_ID")&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p_patched<-  wrap_plots(p1, p2, p3,p4,p5,nrow=5) + plot_annotation(title= paste0("Cell UMAP & Acinar Markers in each adult pancreas sample")) & theme(plot.title = element_text(size = 14, face = "bold"))
ggsave(paste0(egas_dir, "plots/each_AP_acinar_cells_umap.png"), plot = p_patched, width = 80, height =55, units = "cm")


# create dotplots for each sample separately
p1 <- DotPlot(ap,  split.by="patient_ID",features = acinar_markers, group.by = "acinar_status", cols=viridis(6)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- DotPlot(ap,  split.by="patient_ID", features = acinar_journal_markers, group.by = "acinar_status", cols=viridis(6)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_patched<-  wrap_plots(p1, p2,nrow=2, guides="collect") + plot_annotation(title= paste0("Acinar Markers in each adult pancreas sample"), tag_levels=list(c("Acinar Markers", "Acinar Journal Markers"))) & theme(plot.title = element_text(size = 14, face = "bold"))
ggsave(paste0(egas_dir, "plots/each_AP_acinar_cells_dotplot.png"), plot = p_patched, width = 25, height =55, units = "cm")

############################################################################################################################################################

# calculate statistics (with help of chatGPT)
library(pROC)
library(dplyr)
library(rstatix) # for wilcox_test convenience
# define groups
ap$group3 <- "Other"
ap$group3[grepl("^Acinar", ap$Cluster)] <- "Acinar"
ap$group3[ap$Cluster == "Acinar-i"] <- "Acinar_i"

scores <- c("Acinar_Marker_Score1","Acinar_Journal_Marker_Score1","Acinar_Markers_Valid_Score1")

# 1) Group sizes
table(ap$group3) 
# Acinar Acinar_i    Other 
# 35405    43923    33235 

# 2) Pairwise comparisons (Acinar vs Other, Acinar_i vs Other, Acinar vs Acinar_i)
pairwise_results <- list()
pairs <- list(c("Acinar","Other"), c("Acinar_i","Other"), c("Acinar","Acinar_i"))
for(sc in scores){
  for(pp in pairs){
    a <- ap@meta.data[ap$group3==pp[1], sc]
    b <- ap@meta.data[ap$group3==pp[2], sc]
    wt <- wilcox.test(a, b, alternative="two.sided")
    # effect sizes
    cohend <- (mean(a)-mean(b))/sqrt(((length(a)-1)*sd(a)^2 + (length(b)-1)*sd(b)^2)/(length(a)+length(b)-2))
    cliffs <- function(x,y){ (sum(outer(x,y,">")) - sum(outer(x,y,"<"))) / (length(x)*length(y)) }
    cdelta <- cliffs(a,b)
    r_auc <- roc(c(rep(1,length(a)), rep(0,length(b))), c(a,b), quiet=TRUE)
    pairwise_results[[paste(sc, pp[1], "vs", pp[2])]] <- data.frame(score=sc,
                                                                    group1=pp[1], group2=pp[2],
                                                                    n1=length(a), n2=length(b),
                                                                    mean1=mean(a), sd1=sd(a), mean2=mean(b), sd2=sd(b),
                                                                    cohens_d=cohend, cliffs_delta=cdelta,
                                                                    auroc=as.numeric(auc(r_auc)), p=wt$p.value)}}
pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df
#                                                                        score   group1   group2    n1    n2       mean1       sd1       mean2       sd2      cohens_d cliffs_delta     auroc            p
# Acinar_Marker_Score1 Acinar vs Other                    Acinar_Marker_Score1   Acinar    Other 35405 33235  2.42625546 1.3985815 -0.10789582 0.6786156    2.28319046   0.86631536 0.9331577 0.000000e+00
# Acinar_Marker_Score1 Acinar_i vs Other                  Acinar_Marker_Score1 Acinar_i    Other 43923 33235 -0.05697696 0.6693358 -0.10789582 0.6786156    0.07562034  -0.02236396 0.5111820 9.932991e-08
# Acinar_Marker_Score1 Acinar vs Acinar_i                 Acinar_Marker_Score1   Acinar Acinar_i 35405 43923  2.42625546 1.3985815 -0.05697696 0.6693358    2.34532966   0.86169492 0.9308475 0.000000e+00
# Acinar_Journal_Marker_Score1 Acinar vs Other    Acinar_Journal_Marker_Score1   Acinar    Other 35405 33235  2.81021722 1.3202344 -0.07410696 0.7242081    2.68613274   0.90812852 0.9540643 0.000000e+00
# Acinar_Journal_Marker_Score1 Acinar_i vs Other  Acinar_Journal_Marker_Score1 Acinar_i    Other 43923 33235  0.63534241 0.8815920 -0.07410696 0.7242081    0.86780254   0.61936167 0.8096808 0.000000e+00
# Acinar_Journal_Marker_Score1 Acinar vs Acinar_i Acinar_Journal_Marker_Score1   Acinar Acinar_i 35405 43923  2.81021722 1.3202344  0.63534241 0.8815920    1.97858429   0.76829407 0.8841470 0.000000e+00
# Acinar_Markers_Valid_Score1 Acinar vs Other      Acinar_Markers_Valid_Score1   Acinar    Other 35405 33235  3.72996131 1.7545544 -0.03482292 0.9584076    2.64064004   0.88913301 0.9445665 0.000000e+00
# Acinar_Markers_Valid_Score1 Acinar_i vs Other    Acinar_Markers_Valid_Score1 Acinar_i    Other 43923 33235  0.77917916 1.1925432 -0.03482292 0.9584076    0.74146370   0.51960144 0.7598007 0.000000e+00
# Acinar_Markers_Valid_Score1 Acinar vs Acinar_i   Acinar_Markers_Valid_Score1   Acinar Acinar_i 35405 43923  3.72996131 1.7545544  0.77917916 1.1925432    2.00711021   0.76505213 0.8825261 0.000000e+00


# 3) Per-gene checks: mean, fraction expressing (>0), AUROC (Acinar vs Acinar_i)
genes <- c("PRSS1","PRSS2","CPA1","CPA2","CTRC","CELA3A","CELA2A","PNLIP","AMY2A","AMY2B","RBPJL","FOXP2")
expr <- GetAssayData(ap, slot="data") # log-normalized counts maybe
gene_stats <- lapply(genes, function(g){
  if(!g %in% rownames(expr)) return(NULL)
  v <- as.numeric(expr[g, ])
  a <- v[ap$group3=="Acinar"]
  i <- v[ap$group3=="Acinar_i"]
  other <- v[ap$group3=="Other"]
  df <- data.frame(gene=g,
                   mean_acinar=mean(a), mean_acinar_i=mean(i), mean_other=mean(other),
                   frac_acinar=mean(a>0), frac_acinar_i=mean(i>0), frac_other=mean(other>0))
  # AUROC acinar vs acinar_i for the gene
  roc_obj <- roc(c(rep(1,length(a)), rep(0,length(i))), c(a,i), quiet=TRUE)
  df$auroc_acinar_vs_acinar_i <- as.numeric(auc(roc_obj))
  df})

gene_stats_df <- do.call(rbind, gene_stats)
gene_stats_df
# gene mean_acinar mean_acinar_i mean_other frac_acinar frac_acinar_i frac_other auroc_acinar_vs_acinar_i
# 1   PRSS1    5.317807     0.4168819 0.33805610   0.7866968    0.07276370 0.05581465                0.8714799
# 2   PRSS2    0.000000     0.0000000 0.00000000   0.0000000    0.00000000 0.00000000                0.5000000
# 3    CPA1    4.520975     0.4885747 0.30261645   0.7877418    0.10131366 0.05861291                0.8631042
# 4    CPA2    2.504349     0.2443927 0.14163274   0.6301370    0.06964461 0.04040921                0.7880131
# 5    CTRC    3.365919     0.2476014 0.20316657   0.7142494    0.06395283 0.04913495                0.8368560
# 6  CELA3A    4.274164     0.4315779 0.26059314   0.7843525    0.09161487 0.05421995                0.8607216
# 7  CELA2A    0.411524     0.3935063 0.08765101   0.1938992    0.13371127 0.03788175                0.5225650
# 8   PNLIP    2.445963     0.1535293 0.15098667   0.5089959    0.03902284 0.03529412                0.7396286
# 9   AMY2A    2.899839     0.2672695 0.19502017   0.6083886    0.06226806 0.04525350                0.7772695
# 10  AMY2B    1.564122     0.2783379 0.12991608   0.4832651    0.08997564 0.04570483                0.6990215
# 11  RBPJL    2.507081     2.6294284 0.27963641   0.7063692    0.63101336 0.08454942                0.5490032
# 12  FOXP2    1.109007     1.4097583 0.28430998   0.4674199    0.46647542 0.12155860                0.4442974

# 4) Per-cluster post-hoc tests for a chosen score (Kruskal-Wallis + pairwise)
table(ap$Cluster)
# Macrophage   Acinar-REG+  MUC5B+ Ductal    Acinar-i     Acinar-s      Ductal      Alpha       Beta      Gamma   Delta   Quiescent Stellate  Activated Stellate    Endothelial        Schwann
# 833          3689         225              43923        31716         20610       1607        4246      90      1022    178                 2225                  2119               80 

library(rcompanion)
for(sc in scores){
  kw <- kruskal.test(as.formula(paste(sc, "~ Cluster")), data=ap@meta.data)
  cat("K-W for", sc, ", p=", kw$p.value, "\n")
  pw <- pairwise.wilcox.test(ap@meta.data[[sc]], ap@meta.data$Cluster, p.adjust.method="bonferroni")
  # extract comparisons of interest (Acinar vs Acinar-i etc.)
  print(pw)}

# 5) Per-donor Wilcoxon (get p-values per donor) and p.adjust across donors
per_donor_stats <- do.call(rbind, lapply(unique(ap$patient_ID), function(pid){
  sub <- ap@meta.data[ap$patient_ID==pid, ]
  if(!all(c("Acinar","Acinar_i") %in% unique(sub$group3))) return(NULL)
  a <- sub[sub$group3=="Acinar", "Acinar_Marker_Score1"]
  i <- sub[sub$group3=="Acinar_i", "Acinar_Marker_Score1"]
  p <- wilcox.test(a,i)$p.value
  data.frame(patient_ID=pid, p_raw=p, n_acinar=length(a), n_acinar_i=length(i))}))
per_donor_stats$p_adj <- p.adjust(per_donor_stats$p_raw, method="bonferroni")
per_donor_stats
# patient_ID p_raw n_acinar n_acinar_i p_adj
# 1    AFES448     0     7591      16444     0
# 2    AFES365     0     5886       9860     0
# 3    AGBR024     0    13918       8820     0
# 4     TUM_13     0     1011       2484     0
# 5     TUM_25     0     3654       2810     0
# 6     TUM_C1     0     3345       3505     0

#################################################################################################################################################################

# find markers for each cell type compared to all remaining cell types, report only positive markers
pos_cluster_markers <- FindAllMarkers(ap, only.pos = T)
top10 <- pos_cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
write.csv(top10, paste0(egas_dir, "all_AP_top10_pos_cluster_markers_per_cell_type.csv"))

top10 %>% filter(cluster=="Acinar-i")
# A tibble: 10 × 7
#   p_val avg_log2FC pct.1 pct.2 p_val_adj cluster  gene        
#   <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>    <chr>       
# 1     0       2.09 0.192 0.083         0 Acinar-i ARMC3       
# 2     0       1.97 0.273 0.14          0 Acinar-i ZNF804B     
# 3     0       1.90 0.217 0.128         0 Acinar-i PEX6        
# 4     0       1.90 0.204 0.097         0 Acinar-i CLDN10-AS1  
# 5     0       1.87 0.194 0.108         0 Acinar-i C15orf27    
# 6     0       1.86 0.282 0.145         0 Acinar-i AC010967.2  
# 7     0       1.85 0.253 0.124         0 Acinar-i RP11-624L4.1
# 8     0       1.80 0.295 0.159         0 Acinar-i ANKRD62     
# 9     0       1.78 0.28  0.144         0 Acinar-i RP4-765H13.1
# 10    0       1.76 0.197 0.119         0 Acinar-i AC011306.2  

top10 %>% filter(cluster=="Acinar-s")
# A tibble: 10 × 7
#   p_val avg_log2FC pct.1 pct.2 p_val_adj cluster  gene   
#   <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>    <chr>  
# 1     0       4.05 0.517 0.056         0 Acinar-s PNLIP  
# 2     0       3.86 0.792 0.096         0 Acinar-s PRSS1  
# 3     0       3.82 0.697 0.076         0 Acinar-s SYCN   
# 4     0       3.81 0.781 0.11          0 Acinar-s CELA3B 
# 5     0       3.77 0.721 0.085         0 Acinar-s CTRC   
# 6     0       3.70 0.7   0.08          0 Acinar-s CTRB1  
# 7     0       3.66 0.318 0.035         0 Acinar-s CLPS   
# 8     0       3.66 0.781 0.101         0 Acinar-s PLA2G1B
# 9     0       3.64 0.797 0.112         0 Acinar-s CPA1   
# 10    0       3.53 0.797 0.112         0 Acinar-s PRSS3

top10 %>% filter(cluster=="Acinar-REG+")
# # A tibble: 10 × 7
#   p_val avg_log2FC pct.1 pct.2 p_val_adj cluster     gene  
#   <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>       <chr> 
# 1     0       6.53 1     0.028         0 Acinar-REG+ REG3A 
# 2     0       3.71 0.55  0.089         0 Acinar-REG+ REG3G 
# 3     0       3.62 0.703 0.131         0 Acinar-REG+ REG1B 
# 4     0       2.90 0.785 0.208         0 Acinar-REG+ REG1A 
# 5     0       2.52 0.586 0.15          0 Acinar-REG+ SPINK1
# 6     0       2.09 0.185 0.043         0 Acinar-REG+ RPL13 
# 7     0       2.06 0.243 0.067         0 Acinar-REG+ MT-ND3
# 8     0       2.05 0.17  0.04          0 Acinar-REG+ RPS19 
# 9     0       2.02 0.273 0.07          0 Acinar-REG+ EEF1A1
# 10    0       2.01 0.195 0.049         0 Acinar-REG+ RPL13A
