#################################################################
#								#
# Run PIC-seq on Biopsy Data		#
#								#
#################################################################
library("metacell")
source("../scripts/metacell_functions.r")
source("../scripts/pic_parser.r")

message("Running PIC-seq on Biopsy Data")
id = "kidney"

############
scdb_init("saved_work/", force_reinit=T)

id_s = paste0(id, "_singlets")
id_d = paste0(id, "_PIC")
sc_mat = scdb_mat(id)
sin_2d = scdb_mc2d(id_s); sin_cl = scdb_mc(id_s); sin_mat = scdb_mat(id_s)
db_cl = scdb_mc(id_d); db_mat = scdb_mat(id_d)

all_cells = union(names(sin_cl@mc), names(db_cl@mc))
umis = as.matrix(sc_mat@mat[,all_cells])
cell_stats = sc_mat@cell_metadata[all_cells,]
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)

supdir = "figures/supp_figures_1-4/"
if(!dir.exists(supdir)) dir.create(supdir)
scfigs_init(supdir)

#sin_stats = sin_mat@cell_metadata[names(sin_cl@mc),]
#sin_comb = with(sin_stats, paste0(sorting.scheme, ".", treatment, ".", timepoint, ".", date))
#names(sin_comb) = rownames(sin_stats)
#db_stats = db_mat@cell_metadata[names(db_cl@mc),]
#db_comb = with(db_stats, paste0(sorting.scheme, ".", treatment, ".", timepoint, ".", date))
#names(db_comb) = rownames(db_stats)
sin_umis = as.matrix(sin_mat@mat[, names(sin_cl@mc)])
sin_n = sweep(sin_umis,2,colSums(sin_umis),"/") * 1000

color_scheme = sin_cl@color_key
color2name = color_scheme$group; names(color2name) = color_scheme$color
name2color = color_scheme$color; names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)

#B68
lin_ord = c("Tcells-CD8", "Tcells-CXCL14",
           "EC-ACKR1", "EC-VWF", "EC-IGFBP5","EC-VCAM1/UBD", "EC-CXCL9")
#B-PBMC68
# lin_ord = c("Tcells-CD8", "Tcells-IL7R", "Tcells-CCL4", "Tcells-GNLY", 
#             "EC-MGP", "EC-CXCL9", "EC-ACKR1", "EC-IGFBP5")
#B64
# lin_ord = c("Tcells-NKG7",
#             "EC-VCAM1/UBD", "EC-HSPA1A", "EC-ACKR1", "EC-ANGPT2", "EC-IGFBP5", "EC-VWF/AQP1",
#             "RBC")
#B64 and B68
#lin_ord = c("Tcells-CXCR4", "Tcells-NKG7", 
#            "EC-VWF/AQP1", "EC-HSPA1A", "EC-IGFBP5", "EC-ANGPT2", "EC-VCAM1/UBD",
#            "RBC", "cont")            

#change range for lin_ord depending on number of t cell and ec subtypes
t_cells = setdiff(names(sin_cl@mc)[ color2name[ sin_cl@colors[ sin_cl@mc]] %in% lin_ord[1:4]], c())
ec_cells = names(sin_cl@mc)[ color2name[ sin_cl@colors[ sin_cl@mc]] %in% lin_ord[5:8]]
cells = union(t_cells, ec_cells)

##################
# define cell cycle modules

anchor_genes = c('TOP2A', 'MKI67', 'PCNA', 'MCM4', 'CDK1', 'UBE2C')

mcell_mat_rpt_cor_anchors(mat_id=id, gene_anchors = anchor_genes, cor_thresh = 0.1, gene_anti = c(),
                          tab_fn = paste0(supdir, "/lateral_gmods.txt"), sz_cor_thresh = 0.2)

gcor_mat = read.delim(paste0(supdir, "/lateral_gmods.txt"), stringsAsFactor=F, row.names=1)
foc_genes = apply(gcor_mat[, gsub("-", ".", anchor_genes)], 1, which.max)
gset = gset_new_gset(sets = foc_genes, desc = "Cell cycle and MHC-II genes")
scdb_add_gset(paste0(id, "_lateral"), gset)
mcell_mat_ignore_genes(new_mat_id = paste0(id, "_lateral"), mat_id = id, ig_genes = names(foc_genes), reverse = T)
mcell_gset_split_by_dsmat(gset_id = paste0(id, "_lateral"), mat_id = paste0(id, "_lateral"), K = 5)
gset = scdb_gset(paste0(id, "_lateral"))
mcell_plot_gset_cor_mats(gset_id = paste0(id, "_lateral"), scmat_id = paste0(id, "_lateral"))
lateral_clusts = unique(gset@gene_set[ anchor_genes])
mcell_gset_remove_clusts(gset_id = paste0(id, "_lateral"), filt_clusts = lateral_clusts, new_id = paste0(id, "_lateral_f"), reverse=T)

############
# Discard irrelevant PIC populations

db_col = "blue"
# clusts = c(paste0("S.", sin_cl@mc), paste0("P.", db_cl@mc))
# names(clusts) = c(names(sin_cl@mc), names(db_cl@mc))
# num_clusts = as.numeric(factor(clusts)); names(num_clusts) = names(clusts)
# cells = names(clusts)

# afp = log2(mc_compute_fp_abs(num_clusts[cells], umis[,cells]))
# clust_title = names(table(clusts))
# colnames(afp) = clust_title
# clust_num = as.numeric(vecsplit(clust_title, "\\.", 2)) #unlist(lapply(sapply(clust_title, strsplit, "\\."), "[[", 2)))
# clust_orig = vecsplit(clust_title, "\\.", 1)

# ga = "Fscn1"; gb = "Trbc2"; a = afp[ga,]; b = afp[gb,];
# gc = "Klrc1"; gd = "Igkc"; c = afp[gc,]; d = afp[gd,];
# ga = "Fcer1g"; gb = "Gzma"; a = afp[ga,]; b = afp[gb,]; 
# db_clusts = grep("P", clust_title, v = T)
# sin_clusts = grep("S", clust_title, v = T)
# good_clusts = names(which(c < 2 & d < 4))
# clust_cols = ifelse(clust_orig == "P", db_col, sin_cl@colors[ clust_num])

levels(color2name) = c(levels(color2name), "PIC")
levels(name2color) = c(levels(name2color), db_col)
color2name[db_col] = "PIC"
name2color["PIC"] = db_col
#good_clusts = names(which(clust_cols == db_col))
good_pics = names(db_cl@mc) #names(clusts)[ clusts %in% good_clusts]

#######################
# feature selection

bad_genes = grep("Gm[0-9].|Mir|-ps|Rpl|Rps|Ig|Jchain", rownames(umis), v=T)
cells = union(t_cells, ec_cells)
#date = "20181211"; tp = "20h" # choose features by analyzing cells from one experiment only
#rel_cells = intersect(cells, rownames(cell_stats)[ cell_stats$timepoint == tp & cell_stats$date == date])

# remove cycling single cells from feature selection
lateral_set = scdb_gset(paste0(id, "_lateral_f"))@gene_set
cc_genes = setdiff(names(lateral_set), "Ldha")
cc_umis = colSums(umis[cc_genes,])
cc_cells = names(which(cc_umis >= 50))

sub_t_cells = setdiff(t_cells, cc_cells)
sub_ec_cells = setdiff(ec_cells, cc_cells)

lr_features = choose_lr_features(id_s, sub_t_cells, sub_ec_cells, bad_genes=bad_genes, cor_n=100, 
                                 must_haves=names(scdb_gset(paste0(id_s, "_filtered"))@gene_set))
mle_features = choose_mle_features(id_s, id_s,  t_cells, ec_cells, bad_genes=c("Ftl1", bad_genes), shared_thresh=1,
                                   existing_list = names(scdb_gset(paste0(id_s, "_filtered"))@gene_set))

##############
# Run PIC-seq

#comb = paste0(cell_stats$timepoint, ".", cell_stats$date); 
#names(comb) = rownames(cell_stats)

numis=1000
ds = .downsamp(umis[,good_pics], numis)
good_pics = intersect(good_pics, colnames(ds))

mle_res = run_pic_seq(id_s, id_s, ds[,good_pics], t_cells, ec_cells,
                      lr_features, mle_features, fname=paste0(supdir, "/FigS2c.png"), bad_genes = bad_genes,
                      comb = NULL, reg = 1e-4, numis = 1000, downsample=F)

###############
# Analyze triplets

triplet_dir = paste0(supdir, "/FigS2f/")
if(!dir.exists(triplet_dir)) dir.create(triplet_dir)

#coc_cells = rownames(cell_stats)[ cell_stats$treatment == "OVA + LPS"]
triplet_res = analyze_triplets(id_s, id_s, id_s, ds[,good_pics], t_cells, ec_cells, mle_res, mle_features, 
                               bad_genes=bad_genes, reg=1e-4, outdir=triplet_dir, downsample=F)

forward_pr = triplet_res$tr_p
reverse_pr = triplet_res$rev_p
tr_res = triplet_res$mle_res

#############
# run mle on simulated PICs

numis = 1000; k = 5000
res = simulate_doublets(id, t_cells, ec_cells, k, comb = NULL, numis = rep(numis, k))

sim_umis = res$sim_umis; sim_info = res$info
sim_cells = names(which(colSums(sim_umis) == numis))
sim_umis = sim_umis[,sim_cells]; sim_info = sim_info[sim_cells,]

sim_alpha = sim_info$alpha.1; names(sim_alpha) = rownames(sim_info)
sim_mle_res = assign_pics_to_singlets(id_s, id_s, sim_umis, t_cells, ec_cells, sim_alpha,
                                      verbose=T, bad_genes = bad_genes, markers = mle_features, reg = 1e-4)

t_confu = table(sin_cl@mc[ as.vector( sim_info$sim.1)], sim_mle_res$a_mc)
t_n = t_confu / rowSums(t_confu)
ec_confu = table(sin_cl@mc[ as.vector( sim_info$sim.2)], sim_mle_res$b_mc)
ec_n = ec_confu / rowSums(ec_confu)

grad = colorRampPalette(c("white", "#FDC51D", "#CA531C", "#951851", "#36277A", "black"))(1000)
t_cls = factor(color2name[ sin_cl@colors[ as.numeric(rownames(t_n))]], levels = lin_ord)
png(paste0(supdir, "/FigS2d.png"), height = 1500, width = 1500)
par(mar = rep(1,4), lwd = 3, fig = c(0.05,1,0.05,1))
image.2(t_n, zlim = c(0,1), col = grad, annotate = "none", hct = t_cls, vct = t_cls); box()
par(fig = c(0.05,1,0,0.05), new = T)
image(matrix(seq_along(t_cls)), axes = F, col = as.character(name2color[ as.vector(sort(t_cls))])); box()
par(fig = c(0,0.05,0.05,1), new	= T)
image(t(seq_along(t_cls)), axes = F, col = as.character(name2color[ as.vector(sort(t_cls))])); box()
dev.off()

ec_cls = factor(color2name[ sin_cl@colors[ as.numeric(rownames(ec_n))]], levels = lin_ord)
png(paste0(supdir, "/FigS2e.png"), height = 1500, width = 1500)
par(mar = rep(1,4), lwd = 3, fig = c(0.05,1,0.05,1))
image.2(ec_n, zlim = c(0,1), col = grad, annotate = "none", hct = ec_cls,vct = ec_cls); box()
par(fig = c(0.05,1,0,0.05), new	= T)
image(matrix(seq_along(ec_cls)), axes = F, col = as.character(name2color[ as.vector(sort(ec_cls))]));box()
par(fig = c(0,0.05,0.05,1), new = T)
image(t(seq_along(ec_cls)), axes = F, col = as.character(name2color[ as.vector(sort(ec_cls))])); box()
dev.off()

png(paste0(supdir, "/FigS2d-e_colorbar.png"), height = 100, width = 1000)
par(mar = c(3,0,0,0))
image(matrix(1:1000), col = grad)
dev.off() 

save.image(".RDataPIC")