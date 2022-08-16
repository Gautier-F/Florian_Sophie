library(rjson)
library(Seurat)
library(stringr)


# téléchargement des signatures en json sur le site du Broad

read_signature = function(js_file){
    require(rjson)
    temp = fromJSON(file = js_file)   
    sig = temp[[1]]$geneSymbols
    name = names(temp)
    return(list(name, sig ) )
}

# construction de la liste de signatures
files = list.files("Signature_broad/")
list_sig = list()
names_list = c()
for( i in seq_along( files ) ){
    temp = read_signature( paste( "Signature_broad/", files[i] , sep = "") )
    sig = list(c(temp[[2]]))
    names_list[i] = temp[[1]][1]
    list_sig = append( list_sig, sig )
}
names(list_sig) = names_list

# SC de souris => rectifier les noms des gènes
for (i in seq_along(list_sig)){
    list_sig[[i]] = str_to_title(list_sig[[i]])
}

 
for ( rds in list.files("RDS/")){
    so = readRDS(paste("RDS/", rds, sep=""))
    so = AddModuleScore(object = so, features = list_sig, name = names_list)
    saveRDS(so, file = paste("RDS/", rds, sep=""))
}

dir.create("Signatures_plot")

for (rds in list.files("RDS/")){
    name_dir = gsub(".rds", "", rds)
    if (!dir.exists(paste("Signatures_plot/",name_dir, sep = ""))){
        dir.create(paste("Signatures_plot/",name_dir, sep = ""))
    }
    so = readRDS(paste("RDS/", rds, sep=""))
    for ( sig in colnames(so@meta.data)[30:49]){
        fp = FeaturePlot(so, features = sig, min.cutoff = "q1")
        png(paste("Signatures_plot/",name_dir, "/", sig, ".png", sep = ""))
        print(fp)
        dev.off()
    }
}

