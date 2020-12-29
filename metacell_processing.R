library(metacell) 
library(dplyr) 
library(Seurat)

args = commandArgs(trailingOnly = TRUE)
seurat_obj = args[1]
sin_1 = args[2]
sin_2 = args[3]
pics = args[4]

if(!dir.exists("saved_work")) dir.create("saved_work")
scdb_init("saved_work", force_reinit=T)

if(!dir.exists("figures")) dir.create("figures")
figs_dir = "figures/metacell_figs" 
if(!dir.exists(figs_dir)) dir.create(figs_dir)
kidney = readRDS(seurat_obj)

# ###########################
# #Cell count bar graph
# 
# clusts = levels(kidney@active.ident)
# clusts = clusts[is.na(as.numeric(clusts))]
# cell_counts = c(); 
# for(clust in clusts) {
#   subset = kidney[,kidney@active.ident == clust]
#   cell_counts[clust] = dim(subset)[2]
# }
# 
# cell_counts = sort(cell_counts, decreasing = TRUE)
# clust_cols = ifelse(grepl('+', names(cell_counts), fixed = TRUE), "lightblue", "pink")
# names(clust_cols) = names(cell_counts)
# 
# png(paste0(figs_dir,"/cell_counts.png"), height = 1000, width = 1000)
# par(mar = c(17,10,2,1))
# barplot(cell_counts, col = clust_cols, las = 2,
#         cex.names = 2, cex.axis = 2)
# mtext("Cell count", side = 2, cex = 2.5, line = 6)
# legend("topright", legend = c('singlets', 'doublets'),
#        fill = c('pink','lightblue'), box.lty=0, cex = 2)
# dev.off()
# 
# sin_clusts = names(clust_cols[clust_cols == "pink"])
# sin_clusts_counts = cell_counts[sin_clusts]
# actual_pics = cell_counts[clust_cols == "lightblue"]
# exp_names = combn(sin_clusts, 2, simplify = FALSE) #[1:length(actual_pics)]
# exp_counts = combn(sin_clusts_counts, 2, simplify = FALSE) #[1:length(actual_pics)]
# exp_pics = lapply(exp_counts, sum)
# exp_pics = as.numeric(exp_pics)
# names(exp_pics) = lapply(exp_names, paste, collapse = " + ")
# 
# png(paste0(figs_dir, "/actual_vs_expected_pics.png"), height = 750, width = 1000)
# par(mar = c(13,10,2,1), fig = c(0,0.8,0,1))
# barplot(sort(exp_pics, decreasing = TRUE), names.arg = names(exp_pics), col = 'pink', 
#         las = 2, main = "Expected", cex.names = 1, cex.axis = 2, cex.main = 2)
# mtext("Cell count", side = 2, cex = 2.5, line = 6)
# par(mar = c(13,1,2,1), fig = c(0.8,1,0,1), new = T)
# barplot(actual_pics, names.arg = names(actual_pics), col = 'lightblue',
#         las = 2, main = "Actual", cex.names = 1, cex.axis = 2, cex.main = 2)#, axes = FALSE)
# dev.off()
# 
# std_act_pics = (actual_pics - mean(actual_pics))/sd(actual_pics)
# std_exp_pics = (exp_pics - mean(exp_pics))/sd(exp_pics)
# std_pics = sort(c(std_exp_pics, std_act_pics), decreasing = TRUE)
# 
# png(paste0(figs_dir, "/std_actual_vs_expected_pics.png"), height = 750, width = 1000)
# par(mar = c(13,10,2,1))
# barplot(c(sort(std_exp_pics, decreasing = TRUE), std_act_pics), las = 2,
#         col = c(rep('pink', length(std_exp_pics)), rep('lightblue',length(std_act_pics))))
# mtext("Cell count (standardized)", side = 2, cex = 2.5, line = 6)
# dev.off()

###########################
#Create matrices from seurat objects that can be used by metacell package for singlets, pics and both combined
#Does not include genes that were filtered out when the Seurat object were created

id = "kidney"
id_s = paste0(id, "_singlets")
id_d = paste0(id, "_PIC")

#singlets
singlet_clusters = kidney[,kidney@active.ident == sin_1 | kidney@active.ident == sin_2]
sin_sce = as.SingleCellExperiment(singlet_clusters)
sin_mat = scm_import_sce_to_mat(sin_sce)
scdb_add_mat(id_s, sin_mat)

#PICs
doublet_cluster = kidney[,kidney@active.ident == pics]
pic_sce = as.SingleCellExperiment(doublet_cluster)
pic_mat = scm_import_sce_to_mat(pic_sce)
scdb_add_mat(id_d, pic_mat)

#all cells
scdb_add_mat(id, scm_merge_mats(sin_mat,pic_mat))

##########################
#Remove bad genes and small cells from the 3 matrices

mat = scdb_mat(id)
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T),
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T),
             grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))
mcell_mat_ignore_genes(new_mat_id=id, mat_id=id, bad_genes, reverse=F)
mcell_mat_ignore_small_cells(id, id, 500)
mcell_mat_ignore_genes(new_mat_id=id_s, mat_id=id_s, bad_genes, reverse=F)
mcell_mat_ignore_small_cells(id_s, id_s, 500)
mcell_mat_ignore_genes(new_mat_id=id_d, mat_id=id_d, bad_genes, reverse=F)
mcell_mat_ignore_small_cells(id_d, id_d, 500)

#######################
#Create metacells for singlets and pics separately

mcell_add_gene_stat(gstat_id=id_s, mat_id=id_s, force=T)
mcell_add_gene_stat(gstat_id=id_d, mat_id=id_d, force=T)

mcell_gset_filter_varmean(gset_id=id_s, gstat_id=id_s, T_vm=0.2, force_new=T)
mcell_gset_filter_cov(gset_id = id_s, gstat_id=id_s, T_tot=20, T_top3=3)
mcell_gset_filter_varmean(gset_id=id_d, gstat_id=id_d, T_vm=0.2, force_new=T)
mcell_gset_filter_cov(gset_id = id_d, gstat_id=id_d, T_tot=20, T_top3=3)

scfigs_init(figs_dir)

mcell_plot_gstats(gstat_id=id_s, gset_id=id_s)
mcell_plot_gstats(gstat_id=id_d, gset_id=id_d)

genes_anchors = c('TOP2A', 'MKI67', 'PCNA', 'MCM4', 'CDK1', 'UBE2C')
tab_fn = paste(figs_dir, "lateral_gmods_s.txt", sep="/")
gset_nm = paste0(id_s, "_lateral")
mcell_mat_rpt_cor_anchors(mat_id=id_s, gene_anchors = genes_anchors, cor_thresh = 0.1, gene_anti = c(), tab_fn = tab_fn, sz_cor_thresh = 0.2)
gcor_mat = read.table(tab_fn, header=T)
foc_genes = apply(gcor_mat[, genes_anchors], 1, which.max)
gset = gset_new_gset(sets = foc_genes, desc = "Cell cycle genes")
scdb_add_gset(gset_nm, gset)
mcell_mat_ignore_genes(new_mat_id = gset_nm, mat_id = id_s, ig_genes = names(foc_genes), reverse = T)
mcell_gset_split_by_dsmat(gset_id = gset_nm, mat_id = gset_nm, K = 20, force = T)
gset = scdb_gset(gset_nm)
mcell_plot_gset_cor_mats(gset_id = gset_nm, scmat_id = gset_nm)
lateral_clusts = unique(gset@gene_set[genes_anchors])
mcell_gset_remove_clusts(gset_id = gset_nm, filt_clusts = lateral_clusts, new_id = paste0(id_s, "_lateral_f"), reverse=T)

genes_anchors = c('PCNA')
tab_fn_d = paste(figs_dir, "lateral_gmods_d.txt", sep="/")
gset_nm_d = paste0(id_d, "_lateral")
mcell_mat_rpt_cor_anchors(mat_id=id_d, gene_anchors = genes_anchors, cor_thresh = 0.1, gene_anti = c(), tab_fn = tab_fn_d, sz_cor_thresh = 0.2)
gcor_mat_d = read.table(tab_fn_d, header=T)
foc_genes_d = apply(gcor_mat_d[, genes_anchors], 1, which.max)
gset_d = gset_new_gset(sets = foc_genes_d, desc = "Cell cycle genes")
scdb_add_gset(gset_nm_d, gset_d)
mcell_mat_ignore_genes(new_mat_id = gset_nm_d, mat_id = id_d, ig_genes = names(foc_genes_d), reverse = T)
mcell_gset_split_by_dsmat(gset_id = gset_nm_d, mat_id = gset_nm_d, K = 20)
gset_d = scdb_gset(gset_nm_d)
mcell_plot_gset_cor_mats(gset_id = gset_nm_d, scmat_id = gset_nm_d)
lateral_clusts_d = unique(gset_d@gene_set[genes_anchors])
mcell_gset_remove_clusts(gset_id = gset_nm_d, filt_clusts = lateral_clusts_d, new_id = paste0(id_d, "_lateral_f"), reverse=T)


lateral_gset_id = paste0(gset_nm, "_f")
lateral_gset = scdb_gset(lateral_gset_id)
marker_gset = scdb_gset(id_s)
marker_gset = gset_new_restrict_gset(gset = marker_gset, filt_gset = lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
scdb_add_gset(paste0(id_s, "_filtered"), marker_gset)
message("Using ", length(marker_gset@gene_set), " genes as features") #B68 = 209/194, B/PBMC68 = 318, B64 = 192/187, integrated = 249

lateral_gset_id_d = paste0(gset_nm_d, "_f")
lateral_gset_d = scdb_gset(lateral_gset_id_d)
marker_gset_d = scdb_gset(id_d)
marker_gset_d = gset_new_restrict_gset(gset = marker_gset_d, filt_gset = lateral_gset_d, inverse = T, desc = "cgraph gene markers w/o lateral genes")
scdb_add_gset(paste0(id_d, "_filtered"), marker_gset_d)
message("Using ", length(marker_gset_d@gene_set), " genes as features") #B68 = 140, B/PBMC68 = 153, B64 = 57/66, integrated = 113

mcell_add_cgraph_from_mat_bknn(mat_id=id_s,
                               gset_id = paste0(id_s, "_filtered"),
                               graph_id= id_s,
                               K=100,
                               dsamp=T)
mcell_add_cgraph_from_mat_bknn(mat_id=id_d,
                               gset_id = paste0(id_d, "_filtered"),
                               graph_id= id_d,
                               K=100,
                               dsamp=T)

mcell_coclust_from_graph_resamp(
  coc_id=id_s,
  graph_id=id_s,
  min_mc_size=20,
  p_resamp=0.70, n_resamp=500)
mcell_coclust_from_graph_resamp(
  coc_id=id_d,
  graph_id=id_d,
  min_mc_size=20,
  p_resamp=0.70, n_resamp=500)

mcell_mc_from_coclust_balanced(
  coc_id=id_s,
  mat_id= id_s,
  mc_id= id_s,
  K=30, min_mc_size=30, alpha=2)

mcell_mc_from_coclust_balanced(
  coc_id=id_d,
  mat_id= id_d,
  mc_id= id_d,
  K=20, min_mc_size=30, alpha=1) # changed to k = 20 and alpha = 1 for b64

 mcell_plot_outlier_heatmap(mc_id=id_s, mat_id = id_s, T_lfc=4)
mcell_mc_split_filt(new_mc_id=id_s,
                    mc_id=id_s,
                    mat_id=id_s,
                    T_lfc=4, plot_mats=F)
mcell_plot_outlier_heatmap(mc_id=id_d, mat_id = id_d, T_lfc=4)
mcell_mc_split_filt(new_mc_id=id_d,
                    mc_id=id_d,
                    mat_id=id_d,
                    T_lfc=4, plot_mats=F)

########################
#label metacells

mcell_gset_from_mc_markers(gset_id=paste0(id_s,"_markers"), mc_id=id_s)
mcell_mc_plot_marks(mc_id = id_s, gset_id = paste0(id_s,"_markers"), mat_id = id_s, plot_cells=T)
mc_colorize_default(id_s)
mcell_mc2d_force_knn(mc2d_id = id_s, mc_id = id_s, graph_id = id_s)
mcell_mc2d_plot(mc2d_id = id_s, plot_edges = TRUE)

sin_cl = scdb_mc(id_s)
lfp = log2(sin_cl@mc_fp)

# plt = function(gene1, gene2, lfp, colors)
# {
#   plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
#   text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
# 
# }
# plt("VWF", "ACKR1", lfp = lfp, sin_cl@colors)

# genes1 = c('TRBC2','TCF7','CCR7','CD8A','TOP2A','NME1','CD14',
#            'FCGR3A','GNLY','CLIC3','CST3','ICAM1','VWF','LYVE1','VWF')
# genes2 = c('CD3E','IL7R','S100A4','CD8B','CD8A','IL2RB','LYZ',
#            'MS4A7','NKG7','GZMA','CD79B','VCAM1','CD34','TEK','PECAM1')
# par(mfrow=c(2,2))
# par(mar=c(4,4,1,1))
# for (i in seq_along(genes1)) {
#   plt(gene1 = genes1[i], gene2 = genes2[i], lfp = lfp, colors = sin_cl@colors)
#   abline(h=1, lty=2)
# }

mc_hc = mcell_mc_hclust_confu(mc_id = id_s, graph_id = id_s)
mc_sup = mcell_mc_hierarchy(mc_id = id_s, mc_hc = mc_hc, T_gap = 0.04)
mcell_mc_plot_hierarchy(mc_id = id_s, graph_id = id_s,
                        mc_order = mc_hc$order,
                        sup_mc = mc_sup,
                        width = 1200, height = 2400,
                        min_nmc=2, show_mc_ids = T)

# query_sup_by_mcs = function(sup, mcs)
# {
#   do.call('rbind', lapply(1:length(sup), function(i) { csup = sup[[i]]; data.frame(id=i, n_mcs=length(mcs), n_in=length(intersect(mcs, csup$mcs)), sup_size=length(csup$mcs)) })) %>% filter(n_in > 0) %>% arrange(sup_size)
# }
# query_sup_by_mcs(mc_sup, which(lfp['VWF', ] > 2))
# print(mc_sup[[20]])

mc_colorize_sup_hierarchy(mc_id = id_s, supmc = mc_sup, supmc_key = "supmc.txt") #, gene_key = "colorkey.txt")

mcell_mc_plot_hierarchy(mc_id = id_s,
                        graph_id = id_s,
                        mc_order = mc_hc$order,
                        sup_mc = mc_sup,
                        width=1200, height=2400, min_nmc=2, show_mc_ids= T)
mcell_mc2d_plot(mc2d_id = id_s, plot_edges = TRUE)
mcell_mc_plot_marks(mc_id = id_s, gset_id = paste0(id_s,"_markers"), mat_id = id_s, plot_cells=T)

#Cell subtype markers
mc_tab = read.table("supmc.txt", header = TRUE)
marksdir = "cell_subtype_markers"
if(!dir.exists(marksdir)) dir.create(marksdir)
for(i in 1:nrow(mc_tab)) {
  #sort(apply(lfp[,mc_sup[[i]][["mcs"]]], 1, mean), decreasing = TRUE)
  supid = mc_tab[[i,"supid"]]
  name = gsub("/", "_", mc_tab[[i,"name"]])
  marks = sort(mc_sup[[supid]][["marks"]], decreasing = TRUE)
  marks_gap = sort(mc_sup[[supid]][["marks_gap"]], decreasing = TRUE)
  table = cbind(names(marks), marks, names(marks_gap), marks_gap)
  rownames(table) = c()
  colnames(table) = c( "marks", "avg_lfp", "marks_gap", "avg_lfp")
  write.table(table, file = paste0(marksdir,"/", name ,".tsv"), sep="\t", row.names = FALSE)
}
