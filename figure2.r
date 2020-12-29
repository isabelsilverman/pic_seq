########################################
#                       #
# Reproduce Figure 2 for Biopsy Data   #
#                       #
########################################

message("Generating Figure 2")

good_pics = rownames(mle_res)
outdir = "figures/figure2"
if(!dir.exists(outdir)) dir.create(outdir)

supdir = "figures/supp_figures_1-4"

###########
# Compositional changes

#need these for figure 2a,b, c and d
#mixing factor for each pic
alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
#t cell metacell contributing to each pic
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
#ec metacell contributing to each pic
ec_mc = mle_res[good_pics, "b_mc"]; names(ec_mc) = good_pics
#t cell subtype contributing to each pic
parser_t = color2name[ sin_cl@colors[ t_mc]]; names(parser_t) = good_pics
#ec subtype contributing to each pic
parser_ec = color2name[ sin_cl@colors[ ec_mc]]; names(parser_ec) = good_pics

nice_cells = good_pics
# bad_treat = "OVA + LPS transwell"
# bad_cells = rownames(cell_stats)[ cell_stats$treatment == bad_treat]
# t_cells = setdiff(t_cells, bad_cells)
# ec_cells = setdiff(ec_cells, bad_cells)
# nice_cells = setdiff(nice_cells, bad_cells)
# 
# tp_ord = c("3h", "20h", "48h")
# comb = with(cell_stats, paste0(timepoint, "@", sorting.scheme, ".", treatment)); names(comb) = rownames(cell_stats)
 
t_dist = rbind(table(rep("PIC", length(nice_cells)), droplevels(parser_t[nice_cells])), 
 	table(rep("singlets", length(t_cells)), droplevels(sin_names[t_cells])))
t_dist = t_dist[ rowSums(t_dist) > 20,]

# t_dist = t_dist[ rev(order(factor(vecsplit(rownames(t_dist), "@", 1), levels = tp_ord),
# 	factor(vecsplit(rownames(t_dist), "@", 2), 
# 		levels = c("Trbc+.OVA + LPS T only", "Trbc+.OVA + LPS", "doublets.OVA + LPS")))),]
 
ec_dist = rbind(table(rep("PIC", length(nice_cells)), droplevels(parser_ec[nice_cells])), 
         table(rep("singlets", length(ec_cells)), droplevels(sin_names[ec_cells])))
ec_dist = ec_dist[ rowSums(ec_dist) > 20,]

# ec_dist = ec_dist[ rev(order(factor(vecsplit(rownames(ec_dist),	"@", 1), levels = tp_ord),
# 	factor(vecsplit(rownames(ec_dist),"@", 2), 
# 		levels = c("Cd11c+.OVA + LPS DC only", "Cd11c+.OVA + LPS", "doublets.OVA + LPS")))),]
# 
# t_dist = t_dist[, intersect(lin_ord, colnames(t_dist))]
# ec_dist = ec_dist[, intersect(lin_ord, colnames(ec_dist))]

t_n = t_dist / rowSums(t_dist)
ec_n = ec_dist / rowSums(ec_dist)

png(paste0(outdir, "/Fig2a.png"), height=250, width=1000)
par(lwd=6, mar = c(4.5,23.5,3,1))
barplot(t(t_n), horiz = T, col = as.character(name2color[ colnames(t_n)]), space = c(0.3,0.3),
        names.arg = paste0(rownames(t_n)," (",rowSums(t_dist)," cells)"),
        cex.names = 3, las = 2, cex.axis = 2)
mtext("T cell subset distribution", side = 3, cex = 3)
dev.off()

png(paste0(outdir, "/Fig2b.png"), height=250, width=1000)
par(lwd=6, mar = c(4.5,23.5,3,1))
barplot(t(ec_n), horiz = T, col = as.character(name2color[ colnames(ec_n)]), space = c(0.3,0.3), 
        names.arg = paste0(rownames(ec_n)," (",rowSums(ec_dist)," cells)"),
        cex.names = 3, las = 2, cex.axis = 2)
mtext("EC subset distribution", side = 3, cex = 3)
dev.off()

##############
# Compute the expected transcription of PICs

#downsample umi matrix to have 1000 umis for each pic
numis = 1000
ds = .downsamp(umis[,good_pics], numis)
#take only genes that are expressed in singlets
ds_f = ds[setdiff(rownames(sin_cl@e_gc), bad_genes),]
#same set of genes and cells as ds_f but original umis
us = umis[ rownames(ds_f), colnames(ds_f)]
#pics in downsampled umi matrix, should be the same as good_pics
pics = intersect(colnames(ds_f), rownames(mle_res))

#expected expression in pics
exp_us = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[ good_pics, "alpha"], 
	colSums(us), bad_genes = bad_genes)
#expected expression in pics with downsampled umi counts
exp_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[ good_pics, "alpha"], 
	colSums(ds_f), bad_genes = bad_genes)
#expected T cell expression in pics
t_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], rep(1, length(good_pics)),
	colSums(ds_f) * alpha, bad_genes = bad_genes)
#expected EC cell expression in pics
ec_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], rep(0, length(good_pics)),
        colSums(ds_f) * (1 - alpha), bad_genes = bad_genes)
#genes in expected pics
genes = rownames(exp_us)
#observed gene expression summed for all pics
y = rowSums(umis[genes,good_pics])
#expected gene expression summed for all pics
x = rowSums(exp_us[genes,good_pics]);

# generate PICs expected UMIs when modeled as 2T+EC triplets
tr_res = triplet_res$mle_res
triplet_mc = tr_res[pics, c("forward_a1_mc", "forward_a2_mc", "b_mc")]
triplet_alpha = with(tr_res[pics,], cbind(forward_alpha, 1 - forward_alpha) *  alpha); 
exp_tr_us = generate_expected_pics_from_mle(id_s, triplet_mc, triplet_alpha,
        colSums(us), bad_genes = bad_genes)
x2 = rowSums(exp_tr_us[genes, pics])

################
# Figure 2c

reg = 40
#expected gene expression for t cells and ecs individually summed for all pics
sum_t = rowSums(t_n[genes,]); sum_ec = rowSums(ec_n[genes,])
z = log2((10 + sum_ec) / (10 + sum_t))
grad = colorRampPalette(c("limegreen", "gray40", "firebrick3"))(101)
#genes where there is a greater than 2 fold difference between observed and expected expression
fc = 1
disp_genes = names(which(abs(log2((y+reg)/(x+reg))) > fc))
logfc = log2((y+reg)/(x+reg))[disp_genes]
write.table(cbind(disp_genes, logfc), file = paste0(outdir, "/O-E_genes_", fc, "_all.tsv"),sep = '\t', row.names = FALSE)
val = z[disp_genes]; zlim = max(abs(val))
val_n = round((val + zlim) / (2 * zlim) * 100) + 1
#log of expected and observed expression
lx = log2(reg+x); ly = log2(reg+y)
lim = quantile(c(lx,ly), c(0,1))

png(paste0(outdir, "/Fig2c.png"), height=1000, width=1000)
plot(lx, ly, pch = 20, col = "gray", cex = 2,
        xlim = lim, ylim = lim, axes = F, xlab = "", ylab = "")
axis(1); axis(2);
abline(coef = c(fc,1), lty = 2); abline(coef = c(-fc,1), lty = 2)
points(lx[ disp_genes], ly[disp_genes], cex = 3, pch = 21, bg = grad[val_n[disp_genes]])
mtext("Expected total UMIs (log2)", side = 1, cex = 2, line = 2)
mtext("Observed total UMIs (log2)", side = 2, cex = 2, line = 2)
dev.off()

png(paste0(outdir, "/Fig2c_text.png"), height=1000, width=1000)
xy_scatter_genes(x,y, col = c(NA,NA), text = T, reg = reg, fc = fc)
dev.off()

z1 = log2((y+reg)/(x+reg)); z2 = log2((y+reg)/(x2+reg))
png(paste0(supdir, "/FigS4b.png"), height=1000, width=1000)
plot(z1, z2, type = "n", xlab = NULL, ylab = NULL); grid(col = "black"); points(z1, z2, cex = 2, pch = 20, col = "navyblue")
mtext("log2 Observed/Expected (doublets model)", side = 1, cex = 3)
mtext("log2 Observed/Expected (triplets model)", side = 2, cex = 3)
dev.off()
png(paste0(supdir, "/FigS4b_text.png"), height=1000, width=1000)
plot(z1, z2, type = "n"); grid(col = "black"); text(z1, z2, names(z1), col = "navyblue")
dev.off()

#########
# Compute significant changes between PICs observed and expected gene expression

#combination of the two singlet subtypes that make up each pic
pic_comb = interaction(factor(parser_ec, levels = lin_ord),
        factor(parser_t, levels = lin_ord)); names(pic_comb) = names(parser_t)

#combination of singlet subtypes that have more than 30 pics
good_combs = names(which(table(pic_comb) > 30))
#pics with these combinations
analyzed_pics = good_pics[ pic_comb %in% good_combs]
#pics with these combinations labeled by their subtype combination
pic_comb = factor(pic_comb[analyzed_pics], levels = good_combs)

############
# Figure 2d

marker_genes = disp_genes
#marker_genes = c('MGP', 'CTGF', 'ELN', 'BGN', 'TFF3', 'CXCL12')
#marker_genes = c("NAT8", "AQP1", "CRYAB", "GATM", "MT1G", "GPX3", "MT1X", "FXYD2", "STMN1")

#observed mean expression for marker genes in pics for each combination in good_combs
real_m = t(apply(ds_f[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
#expected mean expression for marker genes in pics for each combination in good_combs
exp_m =  t(apply(exp_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
#expected mean expression in t cells for marker genes in pics for each combination in good_combs
t_m =  t(apply(t_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
#expected mean expression in ecs for marker genes in pics for each combination in good_combs
ec_m =  t(apply(ec_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
reg = 0.1
z_m = log2((ec_m + reg) / (t_m + reg));
zlim = max(abs(z_m))

library(Hmisc)

#sum of observed expression for marker genes in pics for each combination in good_combs
m = t(apply(ds_f[marker_genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
#sum of total observed expression in pics for each combination in good_combs
n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
#confidence intervals for each marker gene
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")

circ_frac = 0.1

t_comp = factor(vecsplit(colnames(real_m), "\\.", 2), levels = lin_ord)
ec_comp = factor(vecsplit(colnames(real_m), "\\.", 1), levels = lin_ord)
ord = order(t_comp, ec_comp);

###############
gene = marker_genes[1]
exp_tab = rbind(t_m[gene,ord], ec_m[gene,ord], 0, 0)
real_tab = rbind(0, 0, real_m[gene,ord], 0)
tab = cbind(exp_tab, real_tab)[, rep(seq_along(ord), each = 2) + rep(c(0, length(ord)), length(ord))]
X = barplot(tab)

png(paste0(outdir, "/Fig2d.png"), height=3000, width=2400)
x0_fig = rep(c(0,0.5), each = 3)
x1_fig = rep(c(0.5,1), each = 3)
y0_fig = rep(c(0.1,0.4,0.7), 2)
y1_fig = rep(c(0.4,0.7,1), 2)
for (i in 1:6) {
        gene = marker_genes[i]
        par(fig = c(x0_fig[i], x1_fig[i], y0_fig[i], y1_fig[i]), new=(i>1), lwd=5)
        exp_tab = rbind(t_m[gene,ord], ec_m[gene,ord], 0, 0)
        real_tab = rbind(0, 0, real_m[gene,ord], 0)
        tab = cbind(exp_tab, real_tab)[, rep(seq_along(ord), each = 2) + rep(c(0, length(ord)), length(ord))]
        mtab = round(max(Y[gene,], na.rm=T),2)
        y_ec = -circ_frac * mtab; y_t = y_ec * 2
        ylim = c(0, mtab)
        plot(X, rep(0, length(X)), ylim = ylim, axes=F, xlab = "", ylab = "", type = "n")
        X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), ylim = ylim, las = 2,
                names.arg = rep("", ncol(tab)), space = c(0.5,0), axes = F, add=T)
        mtext(gene, side = 2, cex = 3, line = 0.5)
        obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
        segments(obs_coords, ci.l, y1 = ci.u);
        segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
        axis(2, at = c(0, mtab), las = 2)
        #X = rowMeans(cbind(X[seq(1,length(X),3)], X[seq(2,length(X),3)]))
}

X2 = rowMeans(cbind(X[seq(1,length(X),2)], X[seq(2,length(X),2)]))
par(fig = c(0,0.5,0,0.1), new=T)
plot(X2, rep(0, length(X2)), type="n", axes=F, xlab="", ylab="", ylim = c(-0.3,1.3), xlim = quantile(X,c(0,1)))
segments(X2, 0, y1 = 1)
points(X2, rep(0, length(X2)), pch = 21, bg = as.character(name2color[as.vector(ec_comp[ord])]), cex = 6)
points(X2, rep(1, length(X2)), pch = 21, bg = as.character(name2color[as.vector(t_comp[ord])]), cex = 6)
box()
par(fig = c(0.5,1,0,0.1), new=T)
plot(X2, rep(0, length(X2)), type="n", axes=F, xlab="", ylab="", ylim = c(-0.3,1.3), xlim = quantile(X, c(0,1)))
segments(X2, 0, y1 = 1)
points(X2, rep(0, length(X2)), pch = 21, bg = as.character(name2color[as.vector(ec_comp[ord])]), cex = 6)
points(X2, rep(1, length(X2)), pch = 21, bg = as.character(name2color[as.vector(t_comp[ord])]), cex = 6)
#mtext("Mean expression (down-sampled)", side = 2, cex = 3, line = 2)
box()
dev.off()

#############################
# Alternative to Figure 2d, bar graph for each gene in a separate file

dir.create(paste0(outdir, "/FigS4c/"))
for (gene in marker_genes) {
	png(paste0(outdir, "/FigS4c/", gene, ".png"), height=1000, width=1800)
	par(lwd=5)
	exp_tab = rbind(t_m[gene,ord], ec_m[gene,ord], 0, 0)
	real_tab = rbind(0, 0, real_m[gene,ord], 0)
	tab = cbind(exp_tab, real_tab)[, rep(seq_along(ord), each = 2) + rep(c(0, length(ord)), length(ord))]
	mtab = round(max(Y[gene,], na.rm=T),3)
	y_ec = -circ_frac * mtab; y_t = y_ec * 2
	ylim = c(y_t + y_ec, mtab)
	X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), ylim = ylim, xaxs = "i", las = 2,
		names.arg = rep("", ncol(tab)), space = c(0.5,0), axes = F)
	obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
	segments(obs_coords, ci.l, y1 = ci.u);
	segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
	axis(2, at = c(0, mtab), las = 2)
	X2 = rowMeans(cbind(X[seq(1,length(X),2)], X[seq(2,length(X),2)]))
	segments(X2, y_ec, y1 = y_t)
	points(X2, rep(y_ec, length(X2)), pch = 21, bg = as.character(name2color[ as.vector(ec_comp[ord])]), cex = 10)
	points(X2, rep(y_t, length(X2)), pch = 21, bg = as.character(name2color[ as.vector(t_comp[ord])]), cex = 10)
	mtext(gene, side = 3, cex =3)
	dev.off()
}
