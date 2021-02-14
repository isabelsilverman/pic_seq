########################################
#                                      #
# Reproduce Figure 4 for Biopsy Data   #
#                                      #
########################################

message("Generating Figure 4")

#######
# Need this for all parts of figure4

good_pics = rownames(mle_res)
alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
ec_mc = mle_res[good_pics, "b_mc"]; names(ec_mc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc]]; names(parser_t) = good_pics
parser_ec = color2name[ sin_cl@colors[ ec_mc]]; names(parser_ec) = good_pics

outdir = "figures/figure4"
if(!dir.exists(outdir)) dir.create(outdir)

bad_cells = c()
t_cells = setdiff(t_cells, bad_cells)
ec_cells = setdiff(ec_cells, bad_cells)

##########
# 4a and 4b

cells = union(t_cells, ec_cells)
db_cells = good_pics #intersect(good_pics, rownames(cell_stats)[ cell_stats$sorting.scheme == "doublets" & cell_stats$timepoint == "48h"])
sin_cells = cells #intersect(cells, rownames(cell_stats)[ cell_stats$timepoint == "48h" & cell_stats$sorting.scheme %in% c("Trbc+", "Cd11c+")])

# nb_cells =  intersect(sin_cells, rownames(sin_stats)[ sin_stats$treatment == "helminths"])
# pbs_cells = intersect(sin_cells, rownames(sin_stats)[ sin_stats$treatment == "PBS"])
# 
# nb_doublets = intersect(db_cells, rownames(cell_stats)[ cell_stats$treatment == "helminths"])
# pbs_doublets = intersect(db_cells, rownames(cell_stats)[ cell_stats$treatment == "PBS"])

t_clusts = names(table(sin_cl@mc[ t_cells]))
ec_clusts = names(table(sin_cl@mc[ ec_cells]))

obs = c(table(factor(t_mc[db_cells], levels = t_clusts)),
         table(factor(ec_mc[db_cells], levels = ec_clusts)))
exp = c(table(factor(sin_cl@mc[intersect(t_cells, cells)], levels = t_clusts)),
        table(factor(sin_cl@mc[intersect(ec_cells, cells)], levels = ec_clusts)))
obs_n = obs / sum(obs)
exp_n = exp / sum(exp)
reg = 0.01
enr = log2((reg + obs_n) / (reg + exp_n))

# nb_obs = c(table(factor(t_mc[nb_doublets], levels = t_clusts)),
#          table(factor(ec_mc[nb_doublets], levels = ec_clusts)))
# nb_exp = c(table(factor(sin_cl@mc[intersect(t_cells, nb_cells)], levels = t_clusts)),
#         table(factor(sin_cl@mc[intersect(ec_cells, nb_cells)], levels = ec_clusts)))
# nb_obs_n = nb_obs / sum(nb_obs)
# nb_exp_n = nb_exp / sum(nb_exp)
# reg = 0.01
# nb_enr = log2((reg + nb_obs_n) / (reg + nb_exp_n))

clust_ord = names(enr)[ order(factor(color2name[ sin_cl@colors[ as.numeric(names(enr))]], levels = lin_ord), enr)]
t_ord = intersect(clust_ord, names(table(sin_cl@mc[t_cells])))
t_ord = intersect(t_ord, names(which(exp[t_ord] > 5  & obs[t_ord] > 5)))
ec_ord = intersect(clust_ord, names(table(sin_cl@mc[ec_cells])))
ec_ord = intersect(ec_ord, names(which(exp[ec_ord] > 5 & obs[ec_ord] > 5)))
ord = c(t_ord, ec_ord)
t_share = length(t_ord) / length(ord)
IM = enr[ord]
ylim = max(abs(IM)) * c(-1,1)

png(paste0(outdir, "/Fig4a.png"), height = 700, width = 1000)
barplot(enr[t_ord], border = NA, col = sin_cl@colors[ as.numeric(t_ord)], ylim = ylim, names.arg = rep("", length(t_ord)))
dev.off()

png(paste0(outdir, "/Fig4b.png"), height = 700, width = 1000)
barplot(enr[ec_ord], border = NA, col = sin_cl@colors[ as.numeric(ec_ord)], ylim = ylim, names.arg = rep("", length(ec_ord)))
dev.off()

# png(paste0(outdir, "/Fig4c.png"), height = 700, width = 1000)
# barplot(nb_enr[t_ord], border = NA, col = sin_cl@colors[ as.numeric(t_ord)], ylim = ylim, names.arg = rep("", length(t_ord)))
# dev.off()
# 
# png(paste0(outdir, "/Fig4d.png"), height = 700, width = 1000)
# barplot(nb_enr[ec_ord], border = NA, col = sin_cl@colors[ as.numeric(ec_ord)], ylim = ylim, names.arg = rep("", length(ec_ord)))
# dev.off()

#############
# S8a and S8b


# t_dist = rbind(table(rep("PIC", length(good_pics)), color2name[ sin_cl@colors[t_mc[good_pics]]]),
#                table(rep("singlets", length(t_cells)), color2name[ sin_cl@colors[ sin_cl@mc[t_cells]]]))
# t_dist = t_dist[,intersect(lin_ord, colnames(t_dist))]
# t_dist = t_dist[ rowSums(t_dist) > 20,]
# #t_dist = t_dist[ order(factor(vecsplit(rownames(t_dist), "\\.", 2), levels = c("helminths", "PBS")),
# #        factor(vecsplit(rownames(t_dist),"\\.", 1), levels = c("Trbc+", "Cd11c+", "doublets",  "Ag+ Cd11c+", "Ag+ doublets"))),]
# 
# ec_dist = rbind(table(rep("PIC", length(good_pics)), color2name[ sin_cl@colors[ec_mc[good_pics]]]),
#                 table(rep("singlets", length(ec_cells)), color2name[ sin_cl@colors[ sin_cl@mc[ec_cells]]]))
# ec_dist = ec_dist[,intersect(lin_ord, colnames(ec_dist))]
# ec_dist = ec_dist[ rowSums(ec_dist) > 20,]
# #ec_dist = ec_dist[ order(factor(vecsplit(rownames(ec_dist), "\\.", 2), levels = c("helminths", "PBS")),
# #        factor(vecsplit(rownames(ec_dist),"\\.", 1), levels = c("Trbc+", "Cd11c+", "doublets", "Ag+ Cd11c+", "Ag+ doublets"))),]
# 
# t_n = t_dist / rowSums(t_dist); ec_n = ec_dist / rowSums(ec_dist)

# t_melt = melt(t_n)
# t_melt$tp = vecsplit(as.vector(t_melt$Var1), "\\.", 2)
# t_melt$group = vecsplit(as.vector(t_melt$Var1), "\\.", 1)
# t_melt$rep = vecsplit(as.vector(t_melt$Var1), "\\.", 3)

# ec_melt = melt(ec_n)
# ec_melt$tp = vecsplit(as.vector(ec_melt$Var1), "\\.", 2)
# ec_melt$group = vecsplit(as.vector(ec_melt$Var1), "\\.", 1)
# ec_melt$rep = vecsplit(as.vector(ec_melt$Var1), "\\.", 3)

# # t cell
# library(Hmisc)
# t_melt$abs = with(t_melt, t_dist[cbind(as.vector(Var1), as.vector(Var2))])
# t_melt$bar_id = with(t_melt, paste0(Var2, ":", Var1))
# m = with(t_melt, tapply(abs, bar_id, sum))
# tot_df = data.frame(m = m, bar_id = vecsplit(names(m), ":", 2), pop = vecsplit(names(m), ":", 1))
# tot_df$bar_id = as.vector(tot_df$bar_id)
# #tot_df$group = vecsplit(tot_df$bar_id, "#", 1); tot_df$tp = vecsplit(tot_df$bar_id, "#", 2)
# totals = with(tot_df, tapply( m, bar_id, sum))
# tot_df$n = totals[ tot_df$bar_id]
# 
# #treats = c("PBS@48h", "helminths@48h"); groups = c("Trbc+", "doublets")
# #tot_df = tot_df[ vecsplit(as.vector(tot_df$bar_id), "#", 1) %in% groups & vecsplit(as.vector(tot_df$bar_id), "#", 2) %in% treats,]
# 
# conf = with(tot_df, binconf(m,n))
# rownames(conf) = rownames(tot_df)
# #Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,cells2]), db_int[cells2], mean), 3), "*")
# tot_df = cbind(tot_df, conf)
# 
# ylim = c(0,0.63)
# X = t_melt #[t_melt$tp %in% treats & t_melt$group %in% groups,]
# Y = dcast(tot_df[,c("group", "pop", "PointEst")], factor(group) ~ factor(pop, levels = lin_ord), mean)
# Y.m = as.matrix(Y[,-1]); rownames(Y.m) = Y[,1]
# Y = dcast(tot_df[,c("group", "pop", "Lower")], factor(group) ~ factor(pop, levels = lin_ord), mean)
# Y.l = as.matrix(Y[,-1]); rownames(Y.l) = Y[,1]
# Y = dcast(tot_df[,c("group", "pop", "Upper")], factor(group) ~ factor(pop, levels = lin_ord), mean)
# Y.u = as.matrix(Y[,-1]); rownames(Y.u) = Y[,1]
# 
# fill_cols = rep(c("white", "slateblue2"),2)
# X2_1 = X[ X$group %in% rownames(Y.m),]
# 
# Z = barplot(Y.m, beside=T, las = 2, col = fill_cols, space = rep(c(1,0.2,1,0.2), ncol(Y.m) / 2), ylim = ylim) #c(0, max(X2$value) * 1.05))
# dimnames(Z) = dimnames(Y.m)
# X2$title = paste0(X2$Var2, "_", X2$tp)
# X2$coord = Z[cbind(X2$group, X2$title)]
# X2$offset = X2$coord + runif(nrow(X2), -0.4, 0.4)
# X2$side = with(X2, paste0(tp, "#", rep))
# two_sides = names(which(table(X2$side) == ncol(Y.m)))
# X3 = X2[ X2$side %in% two_sides,]
# X3 = X3[ order(X3$group, X3$Var2, X3$side),]
# X3_singlets = X3[ X3$group == "Trbc+",]; X3_doublets = X3[X3$group == "doublets",]
# 
# paired_pvals = rep(NA, nrow(X3_singlets));
# for (i in seq_along(paired_pvals)) {
#         pop = as.vector(X3_singlets[i, "Var2"])
#         sin_cells = names(which(comb == as.vector(X3_singlets[i, "Var1"])))
#         db_cells = names(which(comb == as.vector(X3_doublets[i, "Var1"])))
#         dist = rbind(table(color2name[ sin_cl@colors[ sin_cl@mc[ sin_cells]]] == pop), table(parser_t[db_cells] == pop))
#         paired_pvals[i] = fisher.test(dist)$p.value
# }
# 
# paired_qvals = p.adjust(paired_pvals, "fdr")
# line_col = ifelse(paired_qvals < 0.05, ifelse(X3_doublets$value > X3_singlets$value, "red2", "blue2"), "gray20")
# png(paste0(supdir, "/FigS8a.png"), height=700, width = 1800)
# Z = barplot2(Y.m, beside=T, las = 2, col = "gray80", space = rep(c(1,0.2,1,0.2), ncol(Y.m) / 2), ylim = ylim, axes = F, border = NA,
# 	plot.ci = T, ci.l = Y.l, ci.u = Y.u, ci.lwd=4)
# axis(2); axis(1, at = colMeans(Z), labels = colnames(Z), las = 2)
# segments(X3_singlets$offset, X3_singlets$value, X3_doublets$offset, X3_doublets$value, col = line_col, lwd = ifelse(line_col != "gray20", 3, 1.5))
# with(X2, points(offset, value, pch = ifelse(group == "doublets", 21, 23), cex = 3, bg = name2color[ as.vector(X2$Var2)]))
# dev.off()
# 
# 
# # ec
# ec_melt$abs = with(ec_melt, ec_dist[cbind(as.vector(Var1), as.vector(Var2))])
# ec_melt$bar_id = with(ec_melt, paste0(Var2, ":", group, "#", tp))
# m = with(ec_melt, tapply(abs, bar_id, sum))
# tot_df = data.frame(m = m, bar_id = vecsplit(names(m), ":", 2), pop = vecsplit(names(m), ":", 1))
# tot_df$bar_id = as.vector(tot_df$bar_id)
# tot_df$group = vecsplit(tot_df$bar_id, "#", 1); tot_df$tp = vecsplit(tot_df$bar_id, "#", 2)
# totals = with(tot_df, tapply( m, bar_id, sum))
# tot_df$n = totals[ tot_df$bar_id]
# 
# treats = c("PBS@48h", "helminths@48h"); groups = c("Cd11c+", "doublets") #, "Ag+ Cd11c+", "Ag+ doublets")
# tot_df = tot_df[ vecsplit(as.vector(tot_df$bar_id), "#", 1) %in% groups & vecsplit(as.vector(tot_df$bar_id), "#", 2) %in% treats,]
# 
# conf = with(tot_df, binconf(m,n))
# rownames(conf) = rownames(tot_df)
# tot_df = cbind(tot_df, conf)
# 
# ylim = c(0,0.63)
# X = ec_melt[ec_melt$tp %in% treats & ec_melt$group %in% groups,]
# Y = dcast(tot_df[,c("tp", "group", "pop", "PointEst")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
# Y.m = as.matrix(Y[,-1]); rownames(Y.m) = Y[,1]
# Y = dcast(tot_df[,c("tp", "group", "pop", "Lower")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
# Y.l = as.matrix(Y[,-1]); rownames(Y.l) = Y[,1]
# Y = dcast(tot_df[,c("tp", "group", "pop", "Upper")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
# Y.u = as.matrix(Y[,-1]); rownames(Y.u) = Y[,1]
# 
# X2 = X[ X$group %in% rownames(Y.m) & X$tp %in% treats,]
# Z = barplot(Y.m, beside=T, las = 2, col = fill_cols, space = rep(c(1,0.2,1,0.2), ncol(Y.m) / 2), ylim = ylim) #c(0, max(X2$value) * 1.05))
# dimnames(Z) = dimnames(Y.m)
# X2$coord = Z[cbind(X2$group, paste0(X2$Var2, "_", X2$tp))]
# with(X2, points(coord + runif(nrow(X2), -0.2, 0.2), value, pch = 21, cex = 4, bg = name2color[ as.vector(X2$Var2)]))
# 
# X2$offset = X2$coord + runif(nrow(X2), -0.4, 0.5)
# X2$side = with(X2, paste0(tp, "#", rep))
# two_sides = names(which(table(X2$side) == ncol(Y.m)))
# X3 = X2[ X2$side %in% two_sides,]
# X3 = X3[ order(X3$group, X3$Var2, X3$side),]
# X3_singlets = X3[ X3$group == "Cd11c+",]; X3_doublets = X3[X3$group == "doublets",]
# 
# paired_pvals = rep(NA, nrow(X3_singlets)); 
# for (i in seq_along(paired_pvals)) {
# 	pop = as.vector(X3_singlets[i, "Var2"])
# 	sin_cells = names(which(comb == as.vector(X3_singlets[i, "Var1"])))
# 	db_cells = names(which(comb == as.vector(X3_doublets[i, "Var1"])))
# 	dist = rbind(table(color2name[ sin_cl@colors[ sin_cl@mc[ sin_cells]]] == pop), table(parser_ec[db_cells] == pop))
# 	paired_pvals[i] = fisher.test(dist)$p.value
# }
# 
# paired_qvals = p.adjust(paired_pvals, "fdr")
# line_col = ifelse(paired_qvals < 0.05, ifelse(X3_doublets$value > X3_singlets$value, "red2", "blue2"), "gray20")
# png(paste0(supdir, "/FigS8b.png"), height=700, width = 1800)
# Z = barplot2(Y.m, beside=T, las = 2, col = "gray80", space = rep(c(1,0.2,1,0.2), ncol(Y.m) / 2), ylim = ylim, axes = F, border = NA,
# 	plot.ci = T, ci.l = Y.l, ci.u = Y.u, ci.lwd=4)
# axis(2); axis(1, at = colMeans(Z), labels = colnames(Z), las = 2)
# segments(X3_singlets$offset, X3_singlets$value, X3_doublets$offset, X3_doublets$value, col = line_col, lwd = ifelse(line_col != "gray20", 3, 1.5))
# with(X2, points(offset, value, pch = ifelse(group == "doublets", 21, 23), cex = 3, bg = name2color[ as.vector(X2$Var2)]))
# dev.off()
# 
# ############
# # S8d
# 
# treats = c("helminths@48h"); groups = c("doublets", "Ag+ doublets")
# X = t_melt[t_melt$tp %in% treats & t_melt$group %in% groups,]
# Y = dcast(X[,c(5,2,3)], factor(group, levels = groups) ~ factor(Var2, levels = lin_ord), mean)
# Y2 = as.matrix(Y[,-1]); rownames(Y2) = Y[,1]
# fill_cols = rep(c("gray70", "gray30"),2)
# X2 = X[ X$group %in% rownames(Y2) & X$tp %in% treats,]
# png(paste0(supdir, "/FigS8d.png"), height=700, width=1500)
# par(lwd=6)
# Z = barplot(Y2, beside=T, las = 1, col = fill_cols, ylim = c(0, max(X2$value) * 1.05), space = rep(c(2,0), ncol(Y2)))
# dimnames(Z) = dimnames(Y2)
# X2$coord = Z[cbind(X2$group, as.vector(X2$Var2))]
# with(X2, points(coord + runif(nrow(X2), -0.2, 0.2), value, pch = 21, cex = 4, bg = name2color[ as.vector(X2$Var2)]))
# dev.off()

#################
# 4g and S8c

numis = 1000
ds = .downsamp(umis[,good_pics], numis)
ds_f = ds[setdiff(rownames(sin_cl@e_gc), bad_genes), good_pics]
us = umis[ rownames(ds_f), colnames(ds_f)]

exp_us = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[ good_pics, "alpha"],
        colSums(us), bad_genes = bad_genes)
exp_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[ good_pics, "alpha"],
	colSums(ds_f), bad_genes = bad_genes)
t_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], rep(1, length(good_pics)),
        colSums(ds_f) * alpha, bad_genes = bad_genes)
ec_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], rep(0, length(good_pics)),
	colSums(ds_f) * (1 - alpha), bad_genes = bad_genes)
genes = rownames(exp_us)


sub_pics = good_pics #intersect(good_pics, rownames(cell_stats)[ cell_stats$sorting.scheme == "doublets" & cell_stats$timepoint == "48h"])

#comb = with(cell_stats, paste0(sorting.scheme, ".", treatment, ".", timepoint, ".", replicate)); names(comb) = rownames(cell_stats)
pic_joint = interaction(factor(color2name[sin_cl@colors[t_mc[sub_pics]]], levels = lin_ord),
        factor(color2name[sin_cl@colors[ec_mc[sub_pics]]], levels = lin_ord))
names(pic_joint) = sub_pics

t_ord = intersect(lin_ord, names(table(parser_t))); ec_ord = intersect(lin_ord, names(table(parser_ec)))

#change this to desired T cell subtype
pop = "Tcells-NKG7"; side = 2
cell_thresh = 10
analyzed_pics = intersect(sub_pics, names(parser_t[which(parser_t == pop)]))
pic_comb = factor(pic_joint[analyzed_pics]); 
#good_joint = names(which(rowSums(table(pic_joint[analyzed_pics], as.vector(cell_stats[analyzed_pics, "treatment"])) < 15) == 0))
names(pic_comb) = analyzed_pics
real_m = t(apply(us[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
exp_m =  t(apply(exp_us[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
x2 = rowSums(((real_m - exp_m)^2 )/exp_m)
qvals = p.adjust(1 - pchisq(x2, df = ncol(real_m) - 1), "fdr")
z = log2(real_m / exp_m); 
z[real_m==0] = NA

y = log2(rowSums(us[names(qvals),analyzed_pics]));
genes = names(which(qvals < 1e-6 & apply(us[ names(qvals), analyzed_pics], 1, function(x) sort(x,T)[3]) > 3 & y > 6))

z_reg = log2((real_m + 5) / (exp_m + 5))
IM = z_reg[genes, ]#grep(paste0("^|\\.", pop, "\\.|$"), colnames(z))]
IM = IM[rowSums(!is.na(IM)) > 0,]

# x_min = apply(IM,1,min); x_max = apply(IM,1,max)
# x = ifelse(abs(x_min) > x_max, x_min, x_max)

# library(tglkmeans)
# 
# k = 15
# data = as.data.frame(IM)
# data$id = rownames(data)
# data = data[,c(ncol(data), 1:(ncol(data) - 1))]
# km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
# centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
# x_min = apply(centers,1,min); x_max = apply(centers,1,max)
# centers = centers[ order(abs(x_min) < x_max, max.col(centers)),]
# km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM)
# ord = order(km_clusts)
# cls = cumsum(table(km_clusts)) / length(km_clusts)

# png(paste0(supdir, "/FigS8c.png"), height = nrow(IM) * 12, width=1000)
# par(mar = c(10,15,5,5))
# image.2(IM, balance = T, annotate = "both", hct = km_clusts)
# dev.off()


marker_genes = c("CCL21", "FABP4")
real_m = t(apply(ds_f[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean)); real_m = real_m[,colnames(IM)]
exp_m =  t(apply(exp_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean)); exp_m = exp_m[,colnames(IM)]
t_m =  t(apply(t_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean)); t_m = t_m[,colnames(IM)]
ec_m =  t(apply(ec_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean)); ec_m = ec_m[,colnames(IM)]

library(Hmisc)
m = t(apply(ds_f[marker_genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")
ord = seq_along(n)

png(paste0(outdir, "/FigS4g.png"), height=700, width=2000)
par(mfrow = c(1,2), mar = c(15,5,5,5))
for (gene in marker_genes) {
        exp_tab = rbind(t_m[gene,ord], ec_m[gene,ord], 0, 0)
        real_tab = rbind(0, 0, real_m[gene,ord], 0)
        tab = cbind(exp_tab, real_tab)[, rep(seq_along(ord), each = 2) + rep(c(0,length(ord)), length(ord))]
        mtab = max(Y[gene,], na.rm=T)
        X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), xaxs = "i", las = 2, space = c(1,0), ylim = c(0,mtab), main = gene, cex.main=2)
        obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
        segments(obs_coords, ci.l, y1 = ci.u);
        segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
}
dev.off()

# ########################################
# # 4h and S8e
# 
# ec_pop = "Migec";
# sub_pics = good_pics #intersect(good_pics, rownames(cell_stats)[ cell_stats$timepoint == "48h" & cell_stats$sorting.scheme %in% c("doublets", "Ag+ doublets") & 
# 	#cell_stats$treatment %in% c("PBS","helminths")])
# analyzed_pics = intersect(sub_pics, names(parser_ec)[parser_ec == ec_pop])
# pic_joint = factor(paste0(cell_stats[analyzed_pics, "treatment"], "@", cell_stats[analyzed_pics, "sorting.scheme"]), 
# 	levels = c("PBS@doublets", "helminths@doublets", "helminths@Ag+ doublets"))
# names(pic_joint) = analyzed_pics
# 
# good_joint = names(which(rowSums(table(pic_joint[analyzed_pics], as.vector(cell_stats[analyzed_pics, "treatment"])) < 15) == 0))
# pic_comb = pic_joint
# names(pic_comb) = analyzed_pics
# real_m = t(apply(us[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
# exp_m =  t(apply(exp_us[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
# x2 = rowSums(((real_m - exp_m)^2 )/exp_m)
# qvals = p.adjust(1 - pchisq(x2, df = ncol(real_m) - 1), "fdr")
# z = log2(real_m / exp_m); 
# z[real_m==0] = NA
# 
# y = log2(rowSums(us[names(qvals),analyzed_pics]));
# genes = union(names(which(qvals < 1e-5 & apply(us[ names(qvals), analyzed_pics], 1, function(x) sort(x,T)[3]) > 3 & y > 5)), c())
# 
# z_reg = log2((real_m + 5) / (exp_m + 5))
# IM = z_reg[genes, ]#grep(paste0("^|\\.", pop, "\\.|$"), colnames(z))]
# IM = IM[rowSums(!is.na(IM)) > 0,]
# 
# x_min = apply(IM,1,min); x_max = apply(IM,1,max)
# x = ifelse(abs(x_min) > x_max, x_min, x_max)
# 
# library(tglkmeans)
# 
# k = 12
# data = as.data.frame(IM)
# data$id = rownames(data)
# data = data[,c(ncol(data), 1:(ncol(data) - 1))]
# km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
# #gc = km$cluster$clust; names(gc) = km$cluster$id
# #cls = cumsum(table(gc)) / length(gc)
# centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
# x_min = apply(centers,1,min); x_max = apply(centers,1,max)
# centers = centers[ order(abs(x_min) < x_max, max.col(centers)),]
# km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM)
# ord = order(km_clusts)
# cls = cumsum(table(km_clusts)) / length(km_clusts)
# 
# png(paste0(supdir, "/FigS8e.png"), height = nrow(IM) * 12, width=1000)
# par(mar = c(5,15,3,3))
# image.2(IM, balance = T, annotate = "both", hct = km_clusts)
# dev.off()


# genes = c("Ccl22", "Ccl17", "Cd40", "Ebi3", "Dll4")
# real_m = t(apply(ds_f[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
# exp_m =  t(apply(exp_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
# t_m =  t(apply(t_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
# ec_m =  t(apply(ec_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
# 
# library(Hmisc)
# m = t(apply(ds_f[genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
# n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
# Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")
# ord = seq_along(n)
# 
# png(paste0(outdir, "/Fig4h.png"), height=700, width=2500)
# par(lwd=2, mfrow = c(1,length(genes)))
# for (gene in genes) {
#         exp_tab = rbind(t_m[gene,], ec_m[gene,], 0, 0)
#         real_tab = rbind(0, 0, real_m[gene,], 0)
#         tab = cbind(exp_tab, real_tab)[, rep(seq_len(ncol(real_m)), each = 2) + rep(c(0,ncol(real_m)), ncol(real_m))]
#         mtab = max(Y[gene,], na.rm=T)
#         X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), xaxs = "i", las = 2, space = c(1,0), ylim = c(0,mtab), main = gene, cex.main=2)
#         obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
#         segments(obs_coords, ci.l, y1 = ci.u);
# 	segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
# }
# dev.off()

#########
# # S8f

# genes = c("Serpinb9", "Serpinb9b", "Fabp4", "Fabp5", "Ccl9", "Nrp2",
# 	"Bcl2a1d", "Bcl2l1", "Tnfsf9", "Cd86")
# 
# real_m = t(apply(ds_f[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
# exp_m =  t(apply(exp_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
# t_m =  t(apply(t_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
# ec_m =  t(apply(ec_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
# 
# library(Hmisc)
# m = t(apply(ds_f[genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
# n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
# Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")
# ord = seq_along(n)
# 
# dir.create(paste0(supdir, "/FigS8f/"))
# for (gene in genes) {
# 	png(paste0(supdir, "/FigS8f/", gene, ".png"), height=700, width=1200)
# 	par(lwd=2)
#         exp_tab = rbind(t_m[gene,], ec_m[gene,], 0, 0)
#         real_tab = rbind(0, 0, real_m[gene,], 0)
# 	tab = cbind(exp_tab, real_tab)[, rep(seq_len(ncol(real_m)), each = 2) + rep(c(0,ncol(real_m)), ncol(real_m))]
#         mtab = max(Y[gene,], na.rm=T)
#         X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), xaxs = "i", las = 2, space = c(1,0), ylim = c(0,mtab))
#         obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
# 	segments(obs_coords, ci.l, y1 = ci.u);
#         segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
# 	dev.off()
# }
