
# check acinar_i cells using the unmodified dataset

#######################
### EGAS00001004653 ###
#######################

# https://www.sciencedirect.com/science/article/pii/S0016508520353993
# data from
# https://explore.data.humancellatlas.org/projects/b3938158-4e8d-4fdb-9e13-9e94270dde16/project-matrices

##############################################################

# Load required libraries 
required_packages <- c('Matrix', 'dplyr', 'Seurat', 'patchwork', 'SeuratData', 'ggplot2', 'SeuratDisk', 'viridis', 'ggpubr', 'ggrepel', 'metafor', 'tidytext', 'RColorBrewer', 'ggtext')
invisible(lapply(required_packages, library, character.only = T))

# obtain directories
# ref_dir with a folder EGAS00001004653, where the adult_pancreas.gz file can be loaded
egas_nr <- "EGAS00001004653"
egas_dir <- paste0(ref_dir,egas_nr, "/")

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
acinar_markers_valid <- c("PRSS1","CPA1","CPA2", "RBPJL")
ap <- AddModuleScore(ap, features = list(acinar_markers), name = "Acinar_Marker_Score")
ap <- AddModuleScore(ap, features = list(acinar_journal_markers), name = "Acinar_Journal_Marker_Score")
ap <- AddModuleScore(ap, features = list(acinar_markers_valid), name = "Acinar_Markers_Journal_Valid")

ap[["Acinar_i"]] <- ap@meta.data$Cluster == "Acinar-i"
ap$acinar_status <- as.character(ap$Cluster)
ap$acinar_status[!grepl("Acinar", ap$Cluster)] <- "Other"

x_label <- unique(c(acinar_markers, acinar_journal_markers))

############################################################################################################################################################

# create dimplots for all samples combined
pal_col <- colorRampPalette(brewer.pal(3,"Blues"))
col_list <- c(pal_col(length(unique(ap@meta.data$Cluster)[!unique(ap@meta.data$Cluster) %in% "Acinar-i"])), "darkred")
names(col_list) <- c(unique(as.character(ap@meta.data$Cluster))[!unique(ap@meta.data$Cluster) %in% "Acinar-i"], "Acinar-i")
p1 <- DimPlot(ap, group.by="Cluster", label=T,label.size=2, repel=T, cols=col_list)+ NoLegend() + ggtitle("Cell Type: Acinar-i darkred, other in blue colors")
p1[[1]]$layers[[1]]$aes_params$alpha <-  ifelse(ap@meta.data$Acinar_i == "TRUE", 1, .8)
p2 <- DimPlot(ap, group.by="Acinar_i", label=T,label.size=2, repel=T, cols= c("gray", "red"))+ NoLegend()
p3 <- FeaturePlot(ap, features = "Acinar_Marker_Score1")+ ggtitle("MarkerScore") + NoLegend() &
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred")) 
p4 <- FeaturePlot(ap, features = c("Acinar_Journal_Marker_Score1"))+ ggtitle("JournalScore") + NoLegend() &
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred")) 
p5 <- FeaturePlot(ap, features = c("Acinar_Markers_Journal_Valid1")) + ggtitle("ValidScore")&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p_patched<-  wrap_plots(p1, p3, p4, p5, ncol=2, guides="collect") + plot_annotation(title= paste0("UMAP colored by Cell Type and Acinar Marker Expression in all adult pancreas samples")) & theme(plot.title = element_text(size = 14, face = "bold"))
ggsave(paste0(egas_dir, "plots/all_AP_acinar_cells_umap.png"), plot = p_patched, width = 23, height =20, units = "cm")

# create dotplots for all samples combined
p1 <- DotPlot(ap, features = unique(c(acinar_markers, acinar_journal_markers)),cols = c("gray", "blue"), group.by = "acinar_status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Acinar Markers across all adult pancreas samples") + WhiteBackground()
ggsave(paste0(egas_dir, "plots/all_AP_acinar_cells_dotplot.png"), plot = p1, width = 25, height =15, units = "cm")

############################################################################################################################################################

# create dimplots for each sample separately
p1 <- DimPlot(ap, group.by="Cluster", split.by="patient_ID", label=T,label.size=2, repel=T)+ NoLegend()
p2 <- DimPlot(ap, group.by="Acinar_i", split.by="patient_ID",label=T,label.size=2, repel=T, cols= c("gray", "red"))+ NoLegend()
p3 <- FeaturePlot(ap, features = "Acinar_Marker_Score1", split.by="patient_ID")&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p4 <- FeaturePlot(ap, features = c("Acinar_Journal_Marker_Score1"), split.by="patient_ID")&
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred"))
p5 <- FeaturePlot(ap, features = c("Acinar_Markers_Journal_Valid"), split.by="patient_ID")&
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

p1 <- DotPlot(ap, split.by="patient_ID", features = unique(c(acinar_markers, acinar_journal_markers)),cols=viridis(6), group.by = "acinar_status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Acinar Markers across adult pancreas samples")
ggsave(paste0(egas_dir, "plots/each_AP_acinar_cells_dotplot.png"), plot = p1, width = 25, height =35, units = "cm")

############################################################################################################################################################

# calculate statistics (with help of chatGPT)
library(pROC)
library(rstatix) # for wilcox_test convenience
# define groups
ap$group3 <- "Other"
ap$group3[grepl("^Acinar", ap$Cluster)] <- "Acinar"
ap$group3[ap$Cluster == "Acinar-i"] <- "Acinar_i"

# Pairwise comparisons (Acinar vs Other, Acinar_i vs Other, Acinar vs Acinar_i) for
# marker sets
# each gene in marker sets
scores <- c("Acinar_Marker_Score1","Acinar_Journal_Marker_Score1","Acinar_Markers_Journal_Valid1")
genes <- c("PRSS1","PRSS2","CPA1","CPA2","CTRC","CELA3A","CELA2A","PNLIP","AMY2A","AMY2B","RBPJL","FOXP2")

# Group sizes
table(ap$group3) 
# Acinar Acinar_i    Other 
# 35405    43923    33235 

table(ap$Cluster)
# Macrophage   Acinar-REG+  MUC5B+ Ductal    Acinar-i     Acinar-s      Ductal      Alpha       Beta      Gamma   Delta   Quiescent Stellate  Activated Stellate    Endothelial        Schwann
# 833          3689         225              43923        31716         20610       1607        4246      90      1022    178                 2225                  2119               80 

# 1) Cell Level: “If I treat every cell as an independent datapoint, how well can marker X separate these groups?”
pairwise_results <- list()
pairs <- list(c("Acinar","Other"), c("Acinar_i","Other"), c("Acinar","Acinar_i"))
expr <- GetAssayData(ap, slot="data") 
for(g in c(scores, genes)){
  for(pp in pairs){
    if (g %in% scores){
        a <- ap@meta.data[ap$group3==pp[1], g]
        b <- ap@meta.data[ap$group3==pp[2], g]
    }else if (g %in% genes){
      if(!g %in% rownames(expr)) return(NULL)
      a <- ap@meta.data$group3==pp[1]
      b <- ap@meta.data$group3==pp[2]
      a <- expr[g,a]
      b <- expr[g,b]}
    wt <- wilcox.test(a, b, alternative="two.sided")
    # effect sizes
    cohend <- (mean(a)-mean(b))/sqrt(((length(a)-1)*sd(a)^2 + (length(b)-1)*sd(b)^2)/(length(a)+length(b)-2))
    cliffs <- function(x,y){ (sum(outer(x,y,">")) - sum(outer(x,y,"<"))) / (length(x)*length(y)) }
    cdelta <- cliffs(a,b)
    r_auc <- roc(c(rep(1,length(a)), rep(0,length(b))), c(a,b), quiet=TRUE)
    pairwise_results[[paste0(g,"_", pp[1], "_vs_", pp[2])]] <- data.frame(score=g,
                                                                    comparison_groups=paste0(pp[1], "_vs_", pp[2]),
                                                                    type=ifelse(g %in% scores, "Marker_set", "Gene"),
                                                                    group1=pp[1], group2=pp[2],
                                                                    n1=length(a), n2=length(b),
                                                                    mean1=mean(a), sd1=sd(a), mean2=mean(b), sd2=sd(b),
                                                                    frac1=mean(a>0), frac2=mean(b>0),
                                                                    cohens_d=cohend, cliffs_delta=cdelta,
                                                                    ci_low=ci.auc(r_auc)[1], ci_high=ci.auc(r_auc)[3],
                                                                    auroc=as.numeric(auc(r_auc)), p=wt$p.value)}}
pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df$p_adj <- p.adjust(pairwise_df$p, method="bonferroni")
# Reorder by average AUROC so best markers are at top
pairwise_df$score <- reorder(pairwise_df$score, pairwise_df$auroc, FUN=mean)
# add significance label column
pairwise_df$signif <- ifelse(!is.na(pairwise_df$p_adj) & pairwise_df$p_adj < 0.001, "*", "")

write.csv(pairwise_df, paste0(egas_dir, "acinar_auroc_pairwise_across_samples.csv"))
pairwise_df <- read.csv(paste0(egas_dir, "acinar_auroc_pairwise_across_samples.csv"))
pairwise_df[,2:22]
#                             score  comparison_groups       type   group1   group2    n1    n2       mean1       sd1       mean2       sd2      frac1      frac2     cohens_d cliffs_delta    ci_low   ci_high     auroc             p         p_adj signif
#   1           Acinar_Marker_Score1    Acinar_vs_Other Marker_set   Acinar    Other 35405 33235  2.42625546 1.3985815 -0.10789582 0.6786156 0.91286541 0.08692643  2.283190461  0.866315363 0.9310688 0.9352465 0.9331577  0.000000e+00  0.000000e+00      *
#   2           Acinar_Marker_Score1  Acinar_i_vs_Other Marker_set Acinar_i    Other 43923 33235 -0.05697696 0.6693358 -0.10789582 0.6786156 0.19616146 0.08692643  0.075620343 -0.022363955 0.5070947 0.5152693 0.5111820  9.932991e-08  4.171856e-06      *
#   3           Acinar_Marker_Score1 Acinar_vs_Acinar_i Marker_set   Acinar Acinar_i 35405 43923  2.42625546 1.3985815 -0.05697696 0.6693358 0.91286541 0.19616146  2.345329663  0.861694920 0.9289866 0.9327083 0.9308475  0.000000e+00  0.000000e+00      *
#   4   Acinar_Journal_Marker_Score1    Acinar_vs_Other Marker_set   Acinar    Other 35405 33235  2.81021722 1.3202344 -0.07410696 0.7242081 0.94896201 0.16717316  2.686132738  0.908128523 0.9524800 0.9556486 0.9540643  0.000000e+00  0.000000e+00      *
#   5   Acinar_Journal_Marker_Score1  Acinar_i_vs_Other Marker_set Acinar_i    Other 43923 33235  0.63534241 0.8815920 -0.07410696 0.7242081 0.74127450 0.16717316  0.867802542  0.619361668 0.8064478 0.8129139 0.8096808  0.000000e+00  0.000000e+00      *
#   6   Acinar_Journal_Marker_Score1 Acinar_vs_Acinar_i Marker_set   Acinar Acinar_i 35405 43923  2.81021722 1.3202344  0.63534241 0.8815920 0.94896201 0.74127450  1.978584285  0.768294067 0.8815677 0.8867264 0.8841470  0.000000e+00  0.000000e+00      *
#   7  Acinar_Markers_Journal_Valid1    Acinar_vs_Other Marker_set   Acinar    Other 35405 33235  3.72996131 1.7545544 -0.03482292 0.9584076 0.92780681 0.12282233  2.640640040  0.889133010 0.9427647 0.9463683 0.9445665  0.000000e+00  0.000000e+00      *
#   8  Acinar_Markers_Journal_Valid1  Acinar_i_vs_Other Marker_set Acinar_i    Other 43923 33235  0.77917916 1.1925432 -0.03482292 0.9584076 0.64237415 0.12282233  0.741463701  0.519601438 0.7562919 0.7633095 0.7598007  0.000000e+00  0.000000e+00      *
#   9  Acinar_Markers_Journal_Valid1 Acinar_vs_Acinar_i Marker_set   Acinar Acinar_i 35405 43923  3.72996131 1.7545544  0.77917916 1.1925432 0.92780681 0.64237415  2.007110206  0.765052134 0.8799135 0.8851387 0.8825261  0.000000e+00  0.000000e+00      *
#   10                         PRSS1    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  5.31780701 2.9036002  0.33805610 1.4186144 0.78669679 0.05581465  2.158358925  0.747813652 0.8715093 0.8763044 0.8739068  0.000000e+00  0.000000e+00      *
#   11                         PRSS1  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.41688188 1.5262967  0.33805610 1.4186144 0.07276370 0.05581465  0.053229223  0.016266590 0.5063977 0.5098689 0.5081333  1.581546e-19  6.642494e-18      *
#   12                         PRSS1 Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  5.31780701 2.9036002  0.41688188 1.5262967 0.78669679 0.07276370  2.180309065  0.742959772 0.8690961 0.8738637 0.8714799  0.000000e+00  0.000000e+00      *
#   13                         PRSS2    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  0.00000000 0.0000000  0.00000000 0.0000000 0.00000000 0.00000000           NA  0.000000000 0.5000000 0.5000000 0.5000000            NA            NA       
#   14                         PRSS2  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.00000000 0.0000000  0.00000000 0.0000000 0.00000000 0.00000000           NA  0.000000000 0.5000000 0.5000000 0.5000000            NA            NA       
#   15                         PRSS2 Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  0.00000000 0.0000000  0.00000000 0.0000000 0.00000000 0.00000000           NA  0.000000000 0.5000000 0.5000000 0.5000000            NA            NA       
#   16                          CPA1    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  4.52097456 2.4886682  0.30261645 1.2426719 0.78774184 0.05861291  2.124546070  0.744017592 0.8695769 0.8744407 0.8720088  0.000000e+00  0.000000e+00      *
#   17                          CPA1  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.48857465 1.4968064  0.30261645 1.2426719 0.10131366 0.05861291  0.133491232  0.041534291 0.5188666 0.5226677 0.5207671  4.179704e-95  1.755476e-93      *
#   18                          CPA1 Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  4.52097456 2.4886682  0.48857465 1.4968064 0.78774184 0.10131366  2.015013630  0.726208346 0.8606355 0.8655729 0.8631042  0.000000e+00  0.000000e+00      *
#   19                          CPA2    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  2.50434918 2.0135266  0.14163274 0.7106443 0.63013699 0.04040921  1.545958026  0.598080631 0.7963321 0.8017485 0.7990403  0.000000e+00  0.000000e+00      *
#   20                          CPA2  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.24439267 0.9155994  0.14163274 0.7106443 0.06964461 0.04040921  0.123284398  0.029187732 0.5129999 0.5161879 0.5145939  4.780579e-67  2.007843e-65      *
#   21                          CPA2 Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  2.50434918 2.0135266  0.24439267 0.9155994 0.63013699 0.06964461  1.498784246  0.576026275 0.7852430 0.7907833 0.7880131  0.000000e+00  0.000000e+00      *
#   22                          CTRC    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  3.36591913 2.2532971  0.20316657 0.9209985 0.71424940 0.04913495  1.817060316  0.677047360 0.8359319 0.8411155 0.8385237  0.000000e+00  0.000000e+00      *
#   23                          CTRC  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.24760142 0.9746807  0.20316657 0.9209985 0.06395283 0.04913495  0.046678741  0.014258263 0.5054950 0.5087632 0.5071291  3.938499e-17  1.654170e-15      *
#   24                          CTRC Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  3.36591913 2.2532971  0.24760142 0.9746807 0.71424940 0.06395283  1.866192905  0.673711978 0.8342877 0.8394242 0.8368560  0.000000e+00  0.000000e+00      *
#   25                        CELA3A    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  4.27416378 2.3778772  0.26059314 1.1146537 0.78435249 0.05421995  2.139812733  0.746429686 0.8708098 0.8756199 0.8732148  0.000000e+00  0.000000e+00      *
#   26                        CELA3A  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.43157791 1.4041639  0.26059314 1.1146537 0.09161487 0.05421995  0.132807070  0.037073823 0.5167169 0.5203569 0.5185369  7.783257e-83  3.268968e-81      *
#   27                        CELA3A Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  4.27416378 2.3778772  0.43157791 1.4041639 0.78435249 0.09161487  2.020939376  0.721443166 0.8582208 0.8632224 0.8607216  0.000000e+00  0.000000e+00      *
#   28                        CELA2A    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  0.41152398 0.9008435  0.08765101 0.4642879 0.19389917 0.03788175  0.447857380  0.154860538 0.5751231 0.5797374 0.5774303  0.000000e+00  0.000000e+00      *
#   29                        CELA2A  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.39350632 1.0332126  0.08765101 0.4642879 0.13371127 0.03788175  0.365422227  0.098178513 0.5472036 0.5509749 0.5490893  0.000000e+00  0.000000e+00      *
#   30                        CELA2A Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  0.41152398 0.9008435  0.39350632 1.0332126 0.19389917 0.13371127  0.018454005  0.045130019 0.5199608 0.5251692 0.5225650  9.949191e-66  4.178660e-64      *
#   31                         PNLIP    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  2.44596311 2.5158669  0.15098667 0.8153325 0.50899590 0.03529412  1.211795325  0.478772023 0.7366080 0.7421641 0.7393860  0.000000e+00  0.000000e+00      *
#   32                         PNLIP  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.15352932 0.7858379  0.15098667 0.8153325 0.03902284 0.03529412  0.003183586  0.003455673 0.5003830 0.5030727 0.5017278  1.228150e-02  5.158232e-01       
#   33                         PNLIP Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  2.44596311 2.5158669  0.15352932 0.7858379 0.50899590 0.03902284  1.288192183  0.479257274 0.7368798 0.7423774 0.7396286  0.000000e+00  0.000000e+00      *
#   34                         AMY2A    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  2.89983944 2.5034248  0.19502017 0.9374136 0.60838865 0.04525350  1.414195433  0.568974944 0.7817108 0.7872641 0.7844875  0.000000e+00  0.000000e+00      *
#   35                         AMY2A  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.26726947 1.0952528  0.19502017 0.9374136 0.06226806 0.04525350  0.070128979  0.016934623 0.5068771 0.5100576 0.5084673  1.663908e-24  6.988412e-23      *
#   36                         AMY2A Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  2.89983944 2.5034248  0.26726947 1.0952528 0.60838865 0.06226806  1.415016467  0.554538996 0.7744733 0.7800657 0.7772695  0.000000e+00  0.000000e+00      *
#   37                         AMY2B    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  1.56412237 1.7357677  0.12991608 0.6222873 0.48326508 0.04570483  1.086778580  0.443057209 0.7186992 0.7243580 0.7215286  0.000000e+00  0.000000e+00      *
#   38                         AMY2B  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  0.27833785 0.9173504  0.12991608 0.6222873 0.08997564 0.04570483  0.184684801  0.044924239 0.5207168 0.5242074 0.5224621 7.785660e-128 3.269977e-126      *
#   39                         AMY2B Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  1.56412237 1.7357677  0.27833785 0.9173504 0.48326508 0.08997564  0.955549586  0.398042973 0.6960745 0.7019685 0.6990215  0.000000e+00  0.000000e+00      *
#   40                         RBPJL    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  2.50708129 1.7800044  0.27963641 0.9551724 0.70636916 0.08454942  1.545924506  0.631636428 0.8129972 0.8186392 0.8158182  0.000000e+00  0.000000e+00      *
#   41                         RBPJL  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  2.62942843 2.1092449  0.27963641 0.9551724 0.63101336 0.08454942  1.373798756  0.574497132 0.7846436 0.7898535 0.7872486  0.000000e+00  0.000000e+00      *
#   42                         RBPJL Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  2.50708129 1.7800044  2.62942843 2.1092449 0.70636916 0.63101336 -0.062133053 -0.098006357 0.5450821 0.5529242 0.5490032 1.116317e-129 4.688530e-128      *
#   43                         FOXP2    Acinar_vs_Other       Gene   Acinar    Other 35405 33235  1.10900718 1.2770155  0.28430998 0.7959901 0.46741986 0.12155860  0.769721412  0.347542508 0.6705937 0.6769488 0.6737713  0.000000e+00  0.000000e+00      *
#   44                         FOXP2  Acinar_i_vs_Other       Gene Acinar_i    Other 43923 33235  1.40975827 1.5600362  0.28430998 0.7959901 0.46647542 0.12155860  0.873955301  0.377058609 0.6856834 0.6913752 0.6885293  0.000000e+00  0.000000e+00      *
#   45                         FOXP2 Acinar_vs_Acinar_i       Gene   Acinar Acinar_i 35405 43923  1.10900718 1.2770155  1.40975827 1.5600362 0.46741986 0.46647542 -0.208766810 -0.111405276 0.4406324 0.4479623 0.4442974 4.363604e-189 1.832714e-187      *

fp <- ggplot(pairwise_df, aes(y=score, x=auroc, color=comparison_groups, shape=type)) +
  geom_point(size=1) +
  geom_errorbarh(aes(xmin=ci_low, xmax=ci_high), height=0.2) +
  geom_vline(xintercept=0.5, linetype="dashed", color="grey50") +
  geom_text_repel(aes(label=signif), vjust=-0.6, hjust=0.5, size=3, color="black") +
  xlim(0.4, 1.0) +
  scale_shape_manual(values=c("Marker_set"=17, "Gene"=16)) +
  theme(axis.text.y=element_text(size=10, face=c(rep("plain",length(genes)), rep("bold", length(scores)))), legend.text=element_text(size=10)) +
  labs(x="AUROC (95% CI)", y="", color="Comparison", shape="", title="Forest plot for Acinar Genes")
ggsave(paste0(egas_dir, "/plots/forest_auroc_acinar_markers.png"), plot = fp, width = 25, height =20, units = "cm")


# 2) Per-donor Level: “Across independent biological replicates, how consistent is marker X in separating the groups?”
# Acinar Marker Expression Levels across donors, using pooled AUROC (Area under the Receiver Operating Curve) computation, 
# CI by DeLong, pooling by Random Effect Modelling (metafor::rma()) with Restricted Maximum Likelihood (REML), 
# Benjamini-Hochberg adjusted p values for pooled AUROC != 0.5: *p<0.05, **p<0.01, ***p<0.001)
pairwise_results <- list()
pairs <- list(c("Acinar","Other"), c("Acinar_i","Other"), c("Acinar","Acinar_i"))
for (pid in unique(ap@meta.data$patient_ID)){
    for(g in c(scores, genes)){
    for(pp in pairs){
      if (g %in% scores){
        a <- ap@meta.data[ap$patient_ID==pid & ap$group3==pp[1], g]
        b <- ap@meta.data[ap$patient_ID==pid & ap$group3==pp[2], g]
      }else if (g %in% genes){
        if(!g %in% rownames(expr)) next
        a <- ap$patient_ID==pid & ap@meta.data$group3==pp[1]
        b <- ap$patient_ID==pid & ap@meta.data$group3==pp[2]
        a <- expr[g,a]
        b <- expr[g,b]}
      if (length(a) < 2 | length(b) < 2) next  # skip if not enough cells
      wt <- wilcox.test(a, b, alternative="two.sided")
      # effect sizes
      cohend <- (mean(a)-mean(b))/sqrt(((length(a)-1)*sd(a)^2 + (length(b)-1)*sd(b)^2)/(length(a)+length(b)-2))
      cliffs <- function(x,y){ (sum(outer(x,y,">")) - sum(outer(x,y,"<"))) / (length(x)*length(y)) }
      cdelta <- cliffs(a,b)
      r_auc <- roc(c(rep(1,length(a)), rep(0,length(b))), c(a,b), quiet=TRUE)
      pairwise_results[[paste0(pid,"_", g,"_", pp[1], "_vs_", pp[2])]] <- data.frame(patient_ID=pid, score=g,
                                                                            comparison_groups=paste0(pp[1], "_vs_", pp[2]),
                                                                            type=ifelse(g %in% scores, "Marker_set", "Gene"),
                                                                            group1=pp[1], group2=pp[2],
                                                                            n1=length(a), n2=length(b),
                                                                            mean1=mean(a), sd1=sd(a), mean2=mean(b), sd2=sd(b),
                                                                            frac1=mean(a>0), frac2=mean(b>0),
                                                                            cohens_d=cohend, cliffs_delta=cdelta,
                                                                            ci_low=ci.auc(r_auc)[1], ci_high=ci.auc(r_auc)[3],
                                                                            auroc=as.numeric(auc(r_auc)), p=wt$p.value)}}}

per_donor_df <- do.call(rbind, pairwise_results)
per_donor_df$p_adj <- p.adjust(per_donor_df$p, method="bonferroni")

# add significance label column
per_donor_df$signif <- ifelse(!is.na(per_donor_df$p_adj) & per_donor_df$p_adj < 0.001, "*", "")

write.csv(per_donor_df, paste0(egas_dir, "acinar_auroc_pairwise_per_donor.csv"))

#######################################################################

per_donor_df <- read.csv(paste0(egas_dir, "acinar_auroc_pairwise_per_donor.csv"))
scores <- c("Acinar_Marker_Score1","Acinar_Journal_Marker_Score1","Acinar_Markers_Journal_Valid1")
genes <- c("PRSS1","PRSS2","CPA1","CPA2","CTRC","CELA3A","CELA2A","PNLIP","AMY2A","AMY2B","RBPJL","FOXP2")

# create pdf comparing marker genes by donor
res <- list()
pdf(paste0(egas_dir, "/plots/meta_forest_plots.pdf"), width=8, height=8)
for (cg in unique(per_donor_df$comparison_groups)) {
  for (g in c(scores, genes)) {
    # subset per donor values for this gene/marker + comparison
    df <- subset(per_donor_df, score == g & comparison_groups == cg)
    # skip if empty or all NA
    if (nrow(df) == 0 || all(is.na(df$auroc))) next
    # derive SE from CI width
    df$sei <- (df$ci_high - df$ci_low) / 3.92
    # drop rows with missing SE or AUROC
    df <- df[!is.na(df$sei) & !is.na(df$auroc), ]
    
    # only run if there’s at least 2 donors left
    if (nrow(df) < 2) next
    # fit model
    fit <- tryCatch(
      rma(yi = df$auroc, sei = df$sei, method="REML"),
      error = function(e) NULL)
    p_val <- tryCatch(
      rma(yi = df$auroc-0.5, sei = df$sei, method="REML"),
      error = function(e) NULL)
    fit$pval <- p_val$pval
    # only create plots, if res != 0
    if (is.null(fit)) next
    res[[paste0(cg, "_", g)]] <- fit
    # plot
    forest(
      fit,
      slab = df$patient_ID,
      xlab = "AUROC (95% CI)",
      alim = c(0.4, 1),
      main=paste0(" AUROC for ", g, " (", cg, ")"),
      mlab = "Pooled")}}
dev.off()

# create summary of marker gene auroc
meta_summary <- do.call(rbind, lapply(names(res), function(nm) {
  if (is.null(res[[nm]])) return(NULL)
  # Extract the subset of donor-level data used
  if(grepl("Acinar_i",nm)){
    cg <- strsplit(nm, "_")[[1]][1:4] |> paste(collapse="_")
  }else{
    cg <- strsplit(nm, "_")[[1]][1:3] |> paste(collapse="_")}
  g  <- sub(paste0(cg, "_"), "", nm)
  # create df
  data.frame(
    comparison_gene = nm,
    gene = g,
    comparison_group = cg,
    k = res[[nm]]$k, # number of donors
    tau2 = res[[nm]]$tau2, # heterogeneity
    pooled_auc =unname(res[[nm]]$b),
    ci_lb = unname(res[[nm]]$ci.lb),
    ci_ub = unname(res[[nm]]$ci.ub),
    pval = unname(res[[nm]]$pval))}))

# adjust for multiple testing
meta_summary$p_adj <- p.adjust(meta_summary$pval, method="BH")
table(meta_summary$p_adj <0.05)
# FALSE  TRUE
# 5    37
table(meta_summary$p_adj <0.01)
# FALSE  TRUE
# 9    33
table(meta_summary$p_adj <0.001)
# FALSE  TRUE
# 11    31

#######################################################################
# create plots for summary
ggplot(meta_summary, aes(x = pooled_auc, 
                         y = gene, 
                         color = comparison_group)) +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub),
                 height=0.2, position=position_dodge(width=0.5)) +
  geom_vline(xintercept = 0.5, linetype="dashed", color="red") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Pooled AUROC (95% CI)",
    y = "Gene / Marker Set",
    color = "Comparison group",
    title = "Meta-analyzed AUROCs per Donor by gene and cell type")

ggplot(meta_summary, aes(x = pooled_auc, y = gene)) +
  geom_point(size=2) +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height=0.2) +
  geom_vline(xintercept = 0.5, linetype="dashed", color="red") +
  geom_vline(xintercept = 0.4, linewidth=2) +
  geom_vline(xintercept = 1, linewidth=2) +
  theme_minimal(base_size = 12) +
  xlim(0.4,1)+
  labs(x = "Pooled AUROC (95% CI)", y = "Gene / Marker Set", title = "Meta-analyzed AUROCs per Donor by gene and cell type") +
  facet_wrap(~comparison_group)
#######################################################################

# adjust score names
meta_summary$gene <- gsub("Acinar_Marker_Score1","MarkerScore", meta_summary$gene)
meta_summary$gene <- gsub("Acinar_Journal_Marker_Score1","JournalScore", meta_summary$gene)
meta_summary$gene <- gsub("Acinar_Markers_Journal_Valid1","ValidScore", meta_summary$gene)
meta_summary$signif <- ""
meta_summary$signif[meta_summary$p_adj < 0.05] <- "*" 
meta_summary$signif[meta_summary$p_adj < 0.01] <- "**" 
meta_summary$signif[meta_summary$p_adj < 0.001] <- "***" 

canonicals <- c("PRSS1","CPA1","CTRC","CELA3A","CELA2A","PNLIP","AMY2A","AMY2B","RBPJL", "MarkerScore", "ValidScore")
plot_df <- meta_summary %>%
  mutate(
    is_canonical = gene %in% canonicals,
    gene_label = ifelse(is_canonical, paste0("<b>", gene, "</b>"), gene),
    gene_col = ifelse(is_canonical, gene, "Invalid"),
    gene = reorder_within(gene_label, pooled_auc, comparison_group))

valid_genes <- unique(plot_df$gene_col)[!unique(plot_df$gene_col) == "Invalid"]
cols <- RColorBrewer::brewer.pal(length(valid_genes),"Paired")
names(cols) <- valid_genes
cols <- c(cols, "Invalid"="grey")

p1 <- ggplot(plot_df, aes(x = pooled_auc, y = gene, colour = gene_col)) +
  geom_point(size=2.5) +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height=0.2, linewidth=1.5) +
  geom_vline(xintercept = 0.5, linetype="dashed", color="red") +
  xlim(0.45,1) +
  facet_wrap(~comparison_group, scales="free_y") +
  scale_y_reordered() +
  theme(axis.text.y = element_markdown(size=10)) +
  geom_text(aes(label=signif), vjust=-0.6, hjust=0.5, size=3, color="black") +
  scale_color_manual(values=cols, labels=names(cols), guide="none") +
  labs(x="Pooled AUROC (95% CI)", y="Gene / Marker Set", 
       title="Meta-analyzed AUROCs across Donors: Colored by canonical Gene / Acinar Markers (bold)")
ggsave(paste0(egas_dir, "plots/forest_across_donors_colored.png"), plot = p1, width = 35, height =25, units = "cm")

#################################################################################################################################################################

# find markers for each cell type compared to all remaining cell types, report only positive markers
pos_cluster_markers <- FindAllMarkers(ap, only.pos = T)
top10 <- pos_cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
write.csv(top10, paste0(egas_dir, "all_AP_top10_pos_cluster_markers_per_cell_type.csv"))

p2 <- top10 %>% filter(cluster=="Acinar-i")
p2
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

# Select relevant columns for table in plot
p2_display <- p2 %>%
  dplyr::select(gene, avg_log2FC, pct.1, pct.2) %>%
  dplyr::rename(
    "Gene" = gene,
    "avg_log2FC" = avg_log2FC,
    "prop_acinar-i" = pct.1,
    "prop_other" = pct.2)

# Create a ggtexttable
p2_table <- ggtexttable(
  p2_display,
  rows = NULL,
  theme = ttheme("minimal", base_size = 12)) + 
  ggtitle("Top 10 markers for Acinar-i") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# create dimplots for all samples combined, incl. top markers for acinar-i cells
pal_col <- colorRampPalette(brewer.pal(3,"Blues"))
col_list <- c(pal_col(length(unique(ap@meta.data$Cluster)[!unique(ap@meta.data$Cluster) %in% "Acinar-i"])), "darkred")
names(col_list) <- c(unique(as.character(ap@meta.data$Cluster))[!unique(ap@meta.data$Cluster) %in% "Acinar-i"], "Acinar-i")
p1 <- DimPlot(ap, group.by="Cluster", label=T,label.size=3, repel=T, cols=col_list)+ NoLegend() + ggtitle("Cell Type: Acinar-i darkred, other in blue colors") + theme(plot.title = element_text(size = 10, face = "bold"))
p1[[1]]$layers[[1]]$aes_params$alpha <-  ifelse(ap@meta.data$Acinar_i == "TRUE", 1, .8)
p3 <- FeaturePlot(ap, features = "Acinar_Marker_Score1")+ ggtitle("MarkerScore") + theme(plot.title = element_text(size = 11, face = "bold"), legend.text = element_text(face = "bold")) &
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred")) 
p4 <- FeaturePlot(ap, features = c("Acinar_Journal_Marker_Score1"))+ ggtitle("JournalScore") + theme(plot.title = element_text(size = 11, face = "bold"), legend.text = element_text(face = "bold")) &
  scale_color_gradientn(colours = c("gray", "lightblue","pink", "red", "darkred")) 
p_patched<-  wrap_plots(p1, p2_table, p3, p4, ncol = 2) + plot_annotation(title= paste0("UMAP Acinar-i cells, top expression, and marker expression in all adult pancreas samples")) & theme(plot.title = element_text(size = 13, face = "bold", hjust=0.5))
ggsave(paste0(egas_dir, "plots/all_AP_acinar_cells_umap.png"), plot = p_patched, width = 26, height =23, units = "cm")
      

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

sessionInfo()     
# R version 4.4.3 (2025-02-28)
# Platform: x86_64-conda-linux-gnu
# Running under: Ubuntu 24.04.1 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /miniconda_dir/miniconda/envs/pdac/lib/libopenblasp-r0.3.29.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Etc/UTC
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] rstatix_0.7.2                pROC_1.18.5                  ggtext_0.1.2                 tidytext_0.4.3               metafor_4.8-0                numDeriv_2016.8-1.1         
# [7] metadat_1.4-0                ggpubr_0.6.1                 viridis_0.6.5                viridisLite_0.4.2            SeuratDisk_0.0.0.9021        pbmcsca.SeuratData_3.0.0    
# [13] pancreasref.SeuratData_1.0.0 panc8.SeuratData_3.0.2       SeuratData_0.2.2.9001        Matrix_1.7-3                 Seurat_5.0.0                 SeuratObject_5.1.0          
# [19] sp_2.2-0                     ggrepel_0.9.6                ggalluvial_0.12.5            patchwork_1.3.1              ggplot2_3.5.2                RColorBrewer_1.1-3          
# [25] dplyr_1.1.4                 
# 
# loaded via a namespace (and not attached):
#   [1] mathjaxr_1.8-0         rstudioapi_0.17.1      jsonlite_2.0.0         magrittr_2.0.3         spatstat.utils_3.1-4   farver_2.1.2           ragg_1.4.0             vctrs_0.6.5           
# [9] ROCR_1.0-11            spatstat.explore_3.4-3 htmltools_0.5.8.1      broom_1.0.8            janeaustenr_1.0.0      Formula_1.2-5          sctransform_0.4.2      parallelly_1.45.0     
# [17] KernSmooth_2.23-26     tokenizers_0.3.0       htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9             plotly_4.11.0          zoo_1.8-14             igraph_2.1.4          
# [25] mime_0.13              lifecycle_1.0.4        pkgconfig_2.0.3        R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-3     future_1.58.0          shiny_1.11.0          
# [33] digest_0.6.37          tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.5.1          SnowballC_0.7.1        textshaping_1.0.1      labeling_0.4.3         progressr_0.15.1      
# [41] spatstat.sparse_3.1-0  httr_1.4.7             polyclip_1.10-7        abind_1.4-8            compiler_4.4.3         bit64_4.6.0-1          withr_3.0.2            backports_1.5.0       
# [49] carData_3.0-5          fastDummies_1.7.5      ggsignif_0.6.4         MASS_7.3-65            rappdirs_0.3.3         tools_4.4.3            lmtest_0.9-40          httpuv_1.6.16         
# [57] future.apply_1.20.0    goftest_1.2-3          glue_1.8.0             nlme_3.1-168           promises_1.3.3         gridtext_0.1.5         grid_4.4.3             Rtsne_0.17            
# [65] cluster_2.1.8.1        reshape2_1.4.4         generics_0.1.4         hdf5r_1.3.12           gtable_0.3.6           spatstat.data_3.1-6    tidyr_1.3.1            data.table_1.17.6     
# [73] xml2_1.3.8             car_3.1-3              utf8_1.2.6             spatstat.geom_3.4-1    RcppAnnoy_0.0.22       RANN_2.6.2             pillar_1.10.2          stringr_1.5.1         
# [81] spam_2.11-1            RcppHNSW_0.6.0         limma_3.62.2           later_1.4.2            splines_4.4.3          lattice_0.22-7         survival_3.8-3         bit_4.6.0             
# [89] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-2          gridExtra_2.3          scattermore_1.2        statmod_1.5.0          matrixStats_1.5.0     
# [97] stringi_1.8.7          lazyeval_0.2.2         codetools_0.2-20       tibble_3.3.0           cli_3.6.5              uwot_0.2.3             xtable_1.8-4           reticulate_1.42.0     
# [105] systemfonts_1.3.1      dichromat_2.0-0.1      Rcpp_1.0.14            globals_0.18.0         spatstat.random_3.4-1  png_0.1-8              spatstat.univar_3.1-3  parallel_4.4.3        
# [113] presto_1.0.0           dotCall64_1.2          listenv_0.9.1          scales_1.4.0           ggridges_0.5.6         leiden_0.4.3.1         purrr_1.0.4            crayon_1.5.3          
# [121] rlang_1.1.6            cowplot_1.1.3       

