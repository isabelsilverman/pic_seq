library("Biostrings")

good_pics = rownames(mle_res)
samples = c("B1", "B2")
name = paste(samples, collapse = "_")
comb_clones = c()
for(samp in samples) {
  i = 1
  file = paste0("vdj_clones/", samp, "/", i)
  while(file.exists(file)) {
    clones_table = readDNAStringSet(file)
    clone_list = vecsplit(clones_table@ranges@NAMES, "_", 1)
    clones = clone_list[-c(1,2)]
    names(clones) = rep(paste0(clone_list[1], "_", samp), length(clones))
    comb_clones = c(comb_clones, clones)
    i = i + 1
    file = paste0("vdj_clones/", samp, "/", i)
  }
}

sin_clones = intersect(names(sin_cl@mc), comb_clones)
sin_clones = comb_clones[comb_clones %in% sin_clones]
num_sin_clones = length(sin_clones)
num_sin = c(length(names(sin_cl@mc)) - num_sin_clones, num_sin_clones)

pic_clones = intersect(good_pics, comb_clones)
pic_clones = comb_clones[comb_clones %in% pic_clones]
num_pic_clones = length(pic_clones)
num_pic = c(length(good_pics) - num_pic_clones, num_pic_clones)

figdir = "figures/vdj_figs/"
if(!dir.exists(figdir)) dir.create(figdir)

labels = c("nonclones","clones")
png(paste0(figdir, name, "_clones_vs_nonclones.png"), height = 300, width = 1000)
par(mar = c(1,7,2,5), fig = c(0,0.5,0,1))
pie(num_sin, labels = labels, col = c("pink", "deeppink"),
    main = "Singlets", cex = 2, cex.main = 2)
par(fig = c(0.5,1,0,1), new = T)
pie(num_pic, labels = labels, col = c("lightblue", "darkblue"),
    main = "PICs", cex = 2, cex.main = 2)
dev.off()

png(paste0(figdir, name, "_clones_vs_nonclones_bar.png"), height = 500, width = 500)
barplot(c(num_sin_clones, num_pic_clones), names.arg = c("singlets", "PICs"), col = c("deeppink", "darkblue"),
        main = "Number of clones in singlets vs. PICs", cex.names = 2)
dev.off()

sin_clones_mc = color2name[sin_cl@colors[sin_cl@mc[sin_clones]]]; names(sin_clones_mc) = sin_clones
sin_dist = table(droplevels(sin_clones_mc))

t_mc = mle_res[pic_clones, "a_mc"]; names(t_mc) = pic_clones
parser_t = color2name[ sin_cl@colors[ t_mc]]; names(parser_t) = pic_clones
ec_mc = mle_res[pic_clones, "b_mc"]; names(ec_mc) = pic_clones
parser_ec = color2name[ sin_cl@colors[ ec_mc]]; names(parser_ec) = pic_clones

pic_clones_joint = interaction(factor(parser_ec, levels = lin_ord),
                        factor(parser_t, levels = lin_ord)); names(pic_clones_joint) = names(parser_t)
pic_dist = table(droplevels(pic_clones_joint))

png(paste0(figdir, name, "_clones_celltypes.png"), height = 300, width = 1000)
par(mar = c(1,10,2,10), fig = c(0,0.5,0,1))
pie(sin_dist, col = as.character(name2color[names(sin_dist)]),
    radius = 1, cex = 1.5, main = "Singlets", cex.main = 2)
par( fig = c(0.5,1,0,1), new = T)
col = colorRampPalette(c("lightblue", "darkblue"))(length(pic_dist))
pie(pic_dist, radius = 1, cex = 1, col = col,
    main = "PICs", cex.main = 2)
dev.off()
