#####################################################
#                       #
# Reproduce Figure 2 with ligand and receptor genes #
#                       #
#####################################################

message("Generating Ligand Receptor Pair Figures")

good_pics = rownames(mle_res)
outdir = "figures/ligand_receptor_figs"
if(!dir.exists(outdir)) dir.create(outdir)

#table = read.csv("../CytokineChemokineLigRec.csv")
#id_p = "CytoChem"
table = read.csv("../PairsLigRecPLUS.csv")
id_p = "LigRec"

pairs = table[, c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")]

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

#################
#Compute the expected transcription of PICs

#downsample umi matrix to have 1000 umis for each pic
numis = 1000
ds = .downsamp(umis[,good_pics], numis)
#take only genes that are expressed in singlets
ds_f = ds[setdiff(rownames(sin_cl@e_gc), bad_genes),]
#same set of genes and cells as ds_f but original umis
us = umis[ rownames(ds_f), colnames(ds_f)]

pairs = pairs[pairs$Ligand.ApprovedSymbol %in% rownames(ds_f) & pairs$Receptor.ApprovedSymbol %in% rownames(ds_f),]
pair_genes = union(pairs$Ligand.ApprovedSymbol, pairs$Receptor.ApprovedSymbol)

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
#observed gene expression for ligands and receptors summed for all pics
y = rowSums(umis[pair_genes,good_pics])
#expected gene expression for ligands and receptors summed for all pics
x = rowSums(exp_us[pair_genes,good_pics]);

##############
# Figure 2c

reg = 40
#expected gene expression for t cells and ecs individually summed for all pics
sum_t = rowSums(t_n[pair_genes,]); sum_ec = rowSums(ec_n[pair_genes,])
z = log2((10 + sum_ec) / (10 + sum_t))
grad = colorRampPalette(c("limegreen", "gray40", "firebrick3"))(101)
fc = 0.5
#genes where there is a greater than 2^fc fold difference between observed and expected expression
disp_genes = names(which(abs(log2((y+reg)/(x+reg))) > fc))
logfc = log2((y+reg)/(x+reg))[disp_genes]
write.table(cbind(disp_genes, logfc), file = paste0(outdir, "/O-E_genes_", fc, "_", id_p, ".tsv"),sep = '\t', row.names = FALSE)
val = z[disp_genes]; zlim = max(abs(val))
val_n = round((val + zlim) / (2 * zlim) * 100) + 1
#log of expected and observed expression
lx = log2(reg+x); ly = log2(reg+y)
lim = quantile(c(lx,ly), c(0,1))

png(paste0(outdir, "/Fig2c_", fc, "_", id_p, ".png"), height=1000, width=1000)
plot(lx, ly, pch = 20, col = "gray", cex = 2,
     xlim = lim, ylim = lim, axes = F, xlab = "", ylab = "")
axis(1); axis(2);
abline(coef = c(fc,1), lty = 2); abline(coef = c(-fc,1), lty = 2)
points(lx[ disp_genes], ly[disp_genes], cex = 3, pch = 21, bg = grad[val_n[disp_genes]])
mtext("Expected total UMIs (log2)", side = 1, cex = 2, line = 2)
mtext("Observed total UMIs (log2)", side = 2, cex = 2, line = 2)
dev.off()

png(paste0(outdir, "/Fig2c_text_", fc, "_", id_p, ".png"), height=1000, width=1000)
xy_scatter_genes(x,y, col = c(NA,NA), text = T, reg = reg, fc = fc)
dev.off()

#################
#Figure 2d

#combination of the two singlet subtypes that make up each pic
pic_comb = interaction(factor(parser_ec, levels = lin_ord),
                       factor(parser_t, levels = lin_ord)); names(pic_comb) = names(parser_t)

#combination of singlet subtypes that have more than 30 pics
good_combs = names(which(table(pic_comb) > 30))
#pics with these combinations
analyzed_pics = good_pics[ pic_comb %in% good_combs]
#pics with these combinations labeled by their subtype combination
pic_comb = factor(pic_comb[analyzed_pics], levels = good_combs)

marker_genes = disp_genes

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

#t cell subtype in each combination of cell subtypes in good_combs 
t_comp = factor(vecsplit(colnames(real_m), "\\.", 2), levels = lin_ord)
#ec subtype in each combination of cell subtypes in good_combs
ec_comp = factor(vecsplit(colnames(real_m), "\\.", 1), levels = lin_ord)
ord = order(t_comp, ec_comp);

dir.create(paste0(outdir, "/Fig2d_", fc, "_", id_p, "/"))
for (gene in marker_genes) {
	png(paste0(outdir, "/Fig2d_", fc, "_", id_p, "/", gene, ".png"), height=1000, width=1800)
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
