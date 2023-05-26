
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


# HEATMAP
# data frame des scores par cluster
MakeDFHeatmap = function(list_set, resolution, name_sig){

    try(if(length(names(list_set)) != length(list_set))
            stop("List's sets must have a name.", call. = FALSE))

    set_name = names(list_set)
    
    list_score = lapply(list_set, function(x, y){
        score_sc = x@meta.data[, c(resolution, name_sig)]
        score_sc = score_sc[order(score_sc[, resolution]), ]
        score_cluster = data.frame(unique(score_sc[, name_sig]))
        return(score_cluster)
        })
    df = data.frame(row.names = paste("cluster_", sort(unique(list_set[[1]]@meta.data[, resolution])), sep = ""))
    for ( c in names(list_set)){
        df = cbind(df, list_score[c])
    }
    colnames(df) = names(list_score)
    return(df)
}

# log2 ratio des comparaisons

LogRatioDE = function(df_score){

    combinaisons = combn(names(df_score), 2) # matrice 2xcomb(5,2)
    names_score = names(df_score)
    names_ratio = c()
    for (i in 1: ncol(combinaisons)){
        comb = combinaisons[,i]
        score_comb = log2(df_score[, comb[2]] / df_score[, comb[1]])
        name_comb = paste(comb[2], "/", comb[1], sep = " ")
        names_ratio = append(names_ratio, name_comb)
        df_score = cbind(df_score, score_comb)  
    }
    names(df_score) = c(names_score, names_ratio)
    return(df_score)
}

HeatmapLogRatio = function(df_score){
    df_ratio = df_score[, grep("/", colnames(df_score))]

    min_max = c(min(df_ratio), max(df_ratio))
    min_max_neg = min_max < 0
    if(sum(min_max_neg) == 1){
        breaks = c(-max(abs(min_max)), -max(abs(min_max))/2, 0, max(abs(min_max))/2, max(abs(min_max)))
        print("1er")
        }else{
        breaks = seq(min_max[1], min_max[2], length.out = 5)
        print("2eme")
        }

    purd = brewer.pal(9, "PuRd")
    bugn = brewer.pal(9, "BuGn")

    pal_pos = purd[c(1, 5, 9)]
    pal_neg = bugn[c(9,5)]
    pal = c(pal_neg, pal_pos)
    col_fun = circlize::colorRamp2(breaks, pal)

    hm_logratio_score = Heatmap(as.matrix(df_ratio), 
        heatmap_legend_param = list(title = "Log2(ratio)"),
        column_names_gp = gpar(fontsize = 8, fontface = "bold"),
        column_names_rot = 45, 
        col = col_fun,
        column_title = "Log2 score ratio",
        column_title_gp = gpar(fontface = "bold"))
    
    return(hm_logratio_score)
}




