
# Calcul du score pseudobulk ( cf AddModuleScore)


score_pseudobulk = function(counts, genes_list){

    # expression moyenne des gènes sur la population totale des cellules
    count_mean = Matrix::rowMeans(counts, na.rm = TRUE)
    count_mean = count_mean[!is.na(count_mean)]
    count_mean = count_mean[order(count_mean)]

    # cut into n=24 intervals
    count_cut = ggplot2::cut_number(count_mean + rnorm(length(count_mean))/1e+30,
                                    n = 24, labels = FALSE, right = FALSE)
    names(count_cut) = names(count_mean)
    
    # get the interval number for each gene in the signature
    bin_num = c()
    for (g in genes_list){
        n = as.vector(count_cut[g])
        
        bin_num = append(bin_num, n)
    }

    # get the name of 100 random genes from each of the selected interval
    ctr_genes = list()
    for (i in bin_num){
        g = sample(count_cut[count_cut == i])
        names_g = names(g)
        ctr_genes = append(ctr_genes, names_g)
    }
    ctr_genes = unique(ctr_genes)

    ####################################
    # score genes contrôles par cellule (cluster)
    score_ctr = Matrix::colMeans(counts[rownames(counts) %in% ctr_genes,], na.rm = TRUE)

    # score genes de la liste 
    score_genes = Matrix::colMeans(counts[rownames(counts) %in% genes_list,], na.rm = TRUE)

    # score signature
    score_sig = score_genes - score_ctr

    return(score_sig)
}

# attribution du score aux cellules du cluster
ExpandScore = function(score, so, resolution, name_score){ 
        num_res = as.numeric(levels(so@meta.data[,resolution])[so@meta.data[,resolution]])
    for (c in sort(unique(so@meta.data[,resolution]))){
        cell_score_cluster = as.vector(score[paste("cluster_", c, sep="")])
        num_res[num_res == c] = cell_score_cluster
    }
    so@meta.data[,name_score] = num_res
    return(so)
}



# plot:

CoordClusterLabel = function(so, resolution){
    df_coord = so@reductions$umap@cell.embeddings
    label_x = c()
    label_y = c()
    for( c in sort(unique(so@meta.data[, resolution]))){
        coord = df_coord[so@meta.data[, resolution] == c,]
        label_x = append(label_x, median(coord[, "UMAP_1"]))
        label_y = append(label_y, median(coord[, "UMAP_2"]))
    }
    names(label_x) = sort(unique(so@meta.data[, resolution]))
    names(label_y) = sort(unique(so@meta.data[, resolution]))
    return(list(label_x, label_y))
}



PlotPseudobulk = function(so, resolution, name_signature, title =""){
    require(ggplot2)
    require(RColorBrewer)
    require(ggpubr)

    coord_label_cluster = CoordClusterLabel(so, resolution)
    df_coord = as.data.frame(so@reductions$umap@cell.embeddings)
    score = so@meta.data[, name_signature]

    # création fonction ramp_palette 
    palette_purd = brewer.pal(9, "PuRd")
    ramp_palette = colorRampPalette(palette_purd)
    # liste des couleurs
    color_score = ramp_palette(length(unique(score)))

    # add legend 
    fp = FeaturePlot(so, features = name_signature,
                label = TRUE, label.size = 7) + 
                scale_color_gradientn(colours = brewer.pal(n = 9, name = "PuRd" ))


    legend = ggpubr::get_legend(fp)
    
    plt = ggplot(df_coord, aes(x= UMAP_1, 
                                y = UMAP_2, 
                                fill= as.factor(score) , 
                                group = as.factor(so@meta.data[, resolution]))) + 
        stat_bin_hex(bins = 50) + scale_fill_manual(values = color_score) + 
        NoLegend() + 
        annotate("text", x=coord_label_cluster[[1]], 
                        y=coord_label_cluster[[2]], 
                        label=names(coord_label_cluster[[1]]), size=8) +
                        
        ggtitle(title)
        
    return(cowplot::plot_grid(plt, legend, ncol = 2, rel_widths = c(1, 0.15)))
}

# plot from list of so
PlotPseudobulkFromList = function(list_so, resolution, sig_name, list_title){

                                # création list_grob
                                list_grob = list()
                                for (i in 1:length(list_so)){
                                    so = list_so[[i]]
                                    p = PlotPseudobulk(so, resolution, sig_name, title = list_title[i]);
                                    p = ggplotify::as.grob(p)
                                    list_grob = append(list_grob, list(p))
                                }

                                # multiple plots
                                cowplot::plot_grid(plotlist = list_grob, 
                                                    ncol = 2, scale = 0.9)

                                }
