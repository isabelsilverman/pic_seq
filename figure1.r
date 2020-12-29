#########################
#			#
# Reproduce Figure 1 for Biopsy Data	#
#			#
#########################

message("Generating Figure 1")
outdir = "figures/figure1"
if(!dir.exists(outdir)) dir.create(outdir)

cells = union(t_cells, ec_cells)

png(paste0(outdir, "/Fig1c.png"), height = 1000, width = 1000)
plot(sin_2d@sc_x[cells], sin_2d@sc_y[cells], pch = 21, bg = sin_cl@colors[ sin_cl@mc[ cells]],
     axes = F, xlab = "", ylab = "", cex = 1.5)
legend("right", legend = sin_cl@color_key$group[1:7], fill = as.character(sin_cl@color_key$color[1:7]),
       box.lty=0, cex = 2)
dev.off()

##################

cells = union(t_cells, ec_cells)
sin_ds = .downsamp(umis[,cells], 500)
good_pics = rownames(mle_res) #[ !(mle_res$forward_triplet | mle_res$reverse_triplet)]

alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
ec_mc = mle_res[good_pics, "b_mc"]; names(ec_mc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc]]; names(parser_t) = good_pics
parser_ec = color2name[ sin_cl@colors[ ec_mc]]; names(parser_ec) = good_pics

t_col = "limegreen"; ec_col =  "firebrick3"

t_clusts = names(table(sin_cl@mc[ t_cells]))
t_ord = t_clusts[order(factor(color2name[ sin_cl@colors[ as.numeric(t_clusts)]], levels = lin_ord))]

ec_clusts = names(table(sin_cl@mc[ ec_cells]))
ec_ord = ec_clusts[order(factor(color2name[ sin_cl@colors[ as.numeric(ec_clusts)]], levels = lin_ord))]

clust_ord = c(t_ord, ec_ord); 
#cells = setdiff(union(t_cells, ec_cells), rownames(cell_stats)[ cell_stats$treatment == "OVA + LPS transwell"])

#tps = factor(cell_stats$timepoint, levels = c("3h", "20h", "48h"))
#g1 = names(tps) = rownames(cell_stats)
#tp_cols = c("lightskyblue1", "lightslateblue", "blue")

umicount = colSums(umis)
ylim = quantile(log2(umicount[union(cells, good_pics)]), c(0,1))
#B68
# t_nms = rev(c('DUSP2', 'CD52', 'NKG7', 'CCL5', 'CXCR4', 'CD3E',
#               'MT1G', 'GPX3', 'FXYD2', "CXCL14"))
# ec_nms = rev(c('VCAM1', 'UBD', 'IGFBP5', 'PLPP1', 'CXCL9', 'LGALS1', 'VWF',
#                'CCL14', 'ACKR1', 'DUSP23'))
#B-PBMC68
t_nms = rev(c('CD8A', 'CD8B', 'NKG7', 'CCL5', 'CCL4', 'GNLY',
              'IL7R', 'LTB', 'CD52', "CXCR4"))
ec_nms = rev(c('ACKR1', 'IFITM3', 'CXCL9', 'CXCL10', 'IGFBP5', 'IFI27', 'MGP',
              'CTGF', 'VWF', 'UBD'))
#64
# t_nms = rev(c('NKG7', 'GNLY', 'CCL4', 'CCL5', 'CD52', 'CST7',
#              'IL32', 'HCST', 'LGALS1', 'S100A4'))
# ec_nms = rev(c('VWF', 'AQP1', 'CLDN5', 'VCAM1', 'UBD', 'ICAM1', 'HSPA1A',
#              'ACKR1', 'ANGPT2', 'IGFBP5'))
#64 and 68
# t_nms = rev(c('NKG7', 'GNLY', 'CCL4', 'CCL5', 'CD52', 'CST7',
#           'IL32', 'S100A4', "JCHAIN", "CXCR4"))
# ec_nms = rev(c('VWF', 'AQP1', 'CLDN5', 'VCAM1', 'UBD', 'ICAM1', 'HSPA1A',
#                'ACKR1', 'ANGPT2', 'IGFBP5'))

png(paste0(outdir, "/Fig1d.png"), height=2000, width=2000)
sin_vec = sin_names[cells]
par(mar = c(0.5,15,0.5,0.5), fig = c(0,1,0,0.45), lwd = 2)
sin_ord = plot_sc_heatmap(id_s, id, t_nms, clusts = sin_vec, good_clusts = lin_ord, cells = cells, annotate=T, normalize=T, lty=1, lwd=2); box()
par(fig = c(0,1,0.45,0.9), new = T)
sin_ord = plot_sc_heatmap(id_s, id, ec_nms, clusts = sin_vec, good_clusts = lin_ord, cells = cells, annotate=T, normalize=T, lty=1, lwd=2); box()
#par(fig=c(0,1,0.85,0.9), new=T)
#image(matrix(as.numeric(tps[sin_ord])), axes = F, col = tp_cols); box()
#par(fig=c(0,1,0.8,0.85), new=T)
#image(matrix(1 * (cell_stats[ sin_ord, "treatment"] == "OVA + LPS")), axes = F, col = c("gray80", "gray20")); box()
par(mar = c(0.5,15,5,0.5), fig=c(0,1,0.9,1), new=T)
plot(log2(umicount[sin_ord]), col = rgb(0.2,0.2,0.2,0.4), pch = 20, xaxs = "i", axes = F, xlab = "", ylab = "", ylim = ylim);
mtext(paste0(length(cells), " single cells"), side = 3, cex = 3, line = 1)
mtext("UMIs(log2)", side = 2, cex = 3, las = 2, line = 3)
axis(2, las = 2); box()
dev.off()

png(paste0(outdir, "/Fig1e.png"), height=200, width=2000)
par(mar = c(0.5,15,0.5,0.5), fig=c(0,1,0,0.5))
image(matrix(seq_along(sin_ord)), col = ifelse(sin_ord %in% t_cells, sin_cl@colors[ sin_cl@mc[ sin_ord]], NA), axes=F)
mtext("T cell identity", side = 2, cex = 3, las = 2, line = 0.5)
par(fig=c(0,1,0.5,1), new=T)
image(matrix(seq_along(sin_ord)), col = ifelse(sin_ord %in% ec_cells, sin_cl@colors[ sin_cl@mc[ sin_ord]], NA), axes=F)
mtext("EC identity", side = 2, cex = 3, las = 2, line = 0.5)
dev.off()

png(paste0(outdir, "/Fig1f.png"), height=2000, width=2000)
db_clusts = interaction(factor(parser_ec, levels = lin_ord),
                        factor(parser_t, levels = lin_ord))
db_vec = as.vector(db_clusts); names(db_vec) = good_pics
par(mar = c(0.5,15,0.5,0.5), fig = c(0,1,0,0.45), lwd = 2)
db_ord = plot_sc_heatmap(id_d, id, t_nms, clusts = db_vec, cells = good_pics, good_clusts = names(table(db_clusts)),
                         annotate=T, normalize=T, lty=1, lwd=2); box()
par(fig = c(0,1,0.45,0.9), new = T)
db_ord = plot_sc_heatmap(id_d, id, ec_nms, clusts = db_vec, cells = good_pics, good_clusts = names(table(db_clusts)),
                         annotate=T, normalize=T, lty=1, lwd=2); box()
#par(fig=c(0,1,0.85,0.9), new=T)
#image(matrix(as.numeric(tps[db_ord])), axes = F, col = tp_cols); box()
par(mar = c(0.5,15,5,0.5), fig=c(0,1,0.9,1), new=T)
plot(log2(umicount[db_ord]), col = rgb(0.2,0.2,0.2,0.4), pch = 20, xaxs = "i", axes = F, xlab = "", ylab = "", ylim = ylim);
mtext(paste0(length(good_pics), " PICs"), side = 3, cex = 3, line = 1)
mtext("UMIs(log2)", side = 2, cex = 3, las = 2, line = 3)
axis(2, las = 2); box()
dev.off()

png(paste0(outdir, "/Fig1g.png"), height=200, width=2000)
par(mar = c(0.5,15,0.5,0.5), fig=c(0,1,0,0.5))
image(matrix(seq_along(db_ord)), col = as.character(name2color[as.character(parser_t[ db_ord])]), axes=F)
mtext("T cell identity", side = 2, cex = 3, las = 2, line = 0.5)
par(fig=c(0,1,0.5,1), new=T)
image(matrix(seq_along(db_ord)), col = as.character(name2color[as.character(parser_ec[ db_ord])]), axes=F)
mtext("EC identity", side = 2, cex = 3, las = 2, line = 0.5)
dev.off()

png(paste0(outdir, "/Fig1h.png"), height=200, width=2000)
par(mar = c(0.5,15,0.5,0.5))
split_count = cbind(alpha[good_pics], 1 - alpha[good_pics])# * umicount[good_pics]
#split_count = split_count / rowSums(split_count)
barplot(t(split_count[db_ord,]), col = c(t_col, ec_col), border = NA, xaxs = "i", space = 0, names.arg = rep("", length(db_ord)), las = 2,
        cex.axis = 2)
mtext(paste("mixing", "factor", sep = "\n"), side = 2, cex = 3, line = 4)
box()
dev.off()
