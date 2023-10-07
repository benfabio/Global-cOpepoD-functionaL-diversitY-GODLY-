
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 06/09/23: R script analyse the covariance between the copepod FD indices and the various variables depicting the functioning of the BCP (productivity related variables and export related variables) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load the final data.frame mae in Script#5.1; find some cool statistical analyses to perform (clustergram, PCA, MDS...)
# - perform cool statistical analyses

### Latest update: 12/09/23 (for preparing OSM24 abstract)

### ------------------------------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("ncdf4")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("gplots")
library("viridis")
library("parallel")
library("ggpubr")
library("ggthemes")
library("marmap")
library("FactoMineR")
library("ggcorrplot")

world <- map_data("world") # coastlines for maps
world2 <- map_data("world2") # coastlines for maps

setwd("/net/kryo/work/fabioben/GODLY/data") # working dir on kryo
fd <- get(load("table_mean_ann_FD_indices_baseline+BCP_06.09.23.RData"))

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1°) Option: clustergram - select variables of interest, log-transform, normalize, center/scale to variance, cluster cells (lat or lon+lat) based on Ward + variables of interest

rownames(fd) <- fd$cell_id

# Retain variables of interest: don't take: 
# NPP (too many NaNs, keep NPP v2 --> sensitivity analyses
# PICO/NANO/MICRO --> keep more complex and biologically accurate PFTs (Diato/Dino/Hapto etc.)
# FRic (Faith instead)
# Jne (ßratio instead) 
# Biovolume/Slope/Export (too many NaN) --> sensitivity analyses

names2keep <- colnames(fd)[c(4:8,10:12,14,16,21:25,29,30,32:35,39)] ; names2keep
data4pca <- fd[,c("x","y",names2keep)]

# Adjust colnames
colnames(data4pca)[12] <- "RaoQ"
colnames(data4pca)[20] <- "NPP"
colnames(data4pca)[24] <- "Slope"
# colnames(data4pca)

# Normalize SR
data4pca$SR <- data4pca$SR/max(data4pca$SR)

# Log transform
data4pca[,c(13:22)] <- log10(data4pca[,c(13:22)])
# Remove NAs
data4pca <- na.omit(data4pca) 
# summary(data4pca); dim(data4pca)


### Option A) Cluster on full global grid (all 35299 1°x1° cells)
# Center and scale values
scaled.data4pca <- data4pca
scaled.data4pca[,c(3:length(scaled.data4pca))] <- base::scale(scaled.data4pca[,c(3:length(scaled.data4pca))], center = T, scale = T)
# summary(scaled.data4pca)

# Perform clustering based on Euclidean distance matrices
# Cluster variables - try different linkages
vars_clust_avg <- hclust(dist(t(scaled.data4pca[,c(3:length(scaled.data4pca))])), method = "average")
#vars_clust_cpt <- hclust(dist(t(scaled.data4pca[,c(3:length(scaled.data4pca))])), method = "complete")
vars_clust_ward <- hclust(dist(t(scaled.data4pca[,c(3:length(scaled.data4pca))])), method = "ward.D2")
# plot(vars_clust_avg)
# plot(vars_clust_cpt)
# plot(vars_clust_ward)

# Cluster space based on z-scores of the variables
space_clust <- hclust(dist(scaled.data4pca[,c(3:length(scaled.data4pca))]), method = "ward.D2") # takes a while because of the dimensionality of the data
# plot(space_clust)
klusters <- cutree(space_clust, k = 5)
scaled.data4pca$k <- as.numeric(klusters) # str(scaled.data4pca$k)

# Plot heatmap
# ?heatmap
# Define more appealing color divergent palette
pale <- rev(brewer.pal(11,"RdBu")) #; pale
pale[c(1:3)] <- "#2166ac"
pale[c(9:11)] <- "#b2182b"
heatmap(x = as.matrix(scaled.data4pca[,c(3:24)]),
        Rowv = as.dendrogram(space_clust),
        Colv = as.dendrogram(vars_clust_avg),
        scale = "none",
        col = pale, # red = low; blue = high
        cexCol = 0.8)

# And don't forget to map associated clusters
map.k <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(k)), data = scaled.data4pca) +
    scale_fill_manual(name = "Clusters", values = rev(c("#1a9850","#a6d96a","#80cdc1","#fee08b","#fdae61"))) +
    geom_contour(colour = "black", binwidth = 1, size = .4, aes(x = x, y = y, z = k), data = scaled.data4pca) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# Fancier version of the map: style like in Benedetti et al. (2022)
# Rotate longitudes
scaled.data4pca$x2 <- scaled.data4pca$x
scaled.data4pca[scaled.data4pca$x < 0,"x2"] <- (scaled.data4pca[scaled.data4pca$x < 0,"x"])+360
# summary(scaled.data4pca$x2)

map.k.fancy <- ggplot() + geom_tile(aes(x = x2, y = y, fill = factor(k)), data = scaled.data4pca) +
    geom_contour(colour = "grey25", binwidth = 1, size = .4, aes(x = x2, y = y, z = k), data = scaled.data4pca) +
    scale_fill_manual(name = "Clusters", values = rev(c("#1a9850","#a6d96a","#80cdc1","#fee08b","#fdae61"))) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                             panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    coord_map("mollweide", orientation = c(90,-180,0))

# And plot latitude FD indices values across clusters to identify them back in the heatmap
data4box <- scaled.data4pca[,c(2,3:12,25)]
colnames(data4box)[1] <- "Latitude"
data4box$Latitude <- abs(data4box$Latitude)
data4box <- melt(data4box, id.vars = "k")
# head(data4box)
boxes <- ggplot(data = data4box) + geom_boxplot(aes(x = factor(k), y = value, fill = factor(k))) +
    scale_fill_manual(name = "Clusters", values = rev(c("#4575b4","#74add1","#e0f3f8","#f46d43","#fee08b"))) +
    geom_hline(yintercept = 0, linetype = "dashed") + xlab("Cluster") + ylab("") + theme_bw() +
    facet_wrap(.~ factor(variable), ncol = 4, scales = "free_y")

# Save plots
setwd("/net/kryo/work/fabioben/GODLY/plots")
#ggsave(plot = map.k, filename = "map_clusters_mean_ann_FD+BCP_ward_08.09.23.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = map.k.fancy, filename = "map_clusters_mean_ann_FD+BCP_ward_08.09.23_v3.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = map.k.fancy, filename = "map_clusters_mean_ann_FD+BCP_ward_08.09.23_v3.jpg", dpi = 300, width = 7, height = 4)
ggsave(plot = boxes, filename = "boxplots_clusters_mean_ann_FD_ward_08.09.23.pdf", dpi = 300, width = 8, height = 7)

jpeg("heatmap_mean_ann_FD+BCP_avg+ward_08.09.23_v2.jpeg", height = 6, width = 6, units = 'in', res = 600)
heatmap(x = as.matrix(scaled.data4pca[,c(3:24)]),
        Rowv = as.dendrogram(space_clust),
        Colv = as.dendrogram(vars_clust_avg),
        scale = "none",
        col = pale, # 
        cexCol = 0.8)
dev.off()



### Option B) Cluster latitudinal averages


### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 2°) Perform PCA with BCP varibales as quanti sup
colnames(data4pca)
pca <- PCA(X = data4pca[,c(3:length(data4pca))], scale.unit = T, graph = F, ncp = 5, quanti.sup = c(11:22))
# summary(pca)
#                        Dim.1   Dim.2   Dim.3   
# Variance               6.899   1.572   0.791 
# % of var.             68.992  15.724   7.913 
# Cumulative % of var.  68.992  84.715  92.629
plot.PCA(pca, axes = c(1,2), choix = "var") 
plot.PCA(pca, axes = c(2,3), choix = "var") 


### Perform PCA only on tropical grid cells (0-30) and extra tropicals ones
tropical.pca <- PCA(X = data4pca[which(abs(data4pca$y) <= 30),c(3:length(data4pca))], scale.unit = T, graph = F, ncp = 5, quanti.sup = c(11:22))
# summary(tropical.pca)
#                         Dim.1 Dim.2   Dim.3   
#Variance               6.194   2.041   0.933  
#% of var.              61.938  20.410  9.332 
#Cumulative % of var.   61.938  82.349  91.68
plot.PCA(tropical.pca, axes = c(1,2), choix = "var") 
# plot.PCA(tropical.pca, axes = c(2,3), choix = "var") 


extra.pca <- PCA(X = data4pca[which(abs(data4pca$y) > 30),c(3:length(data4pca))], scale.unit = T, graph = F, ncp = 5, quanti.sup = c(11:22))
# summary(extra.pca)
#                       Dim.1   Dim.2   Dim.3 
# Variance              6.172   1.654   1.090  
# % of var.             61.722  16.536  10.89
# Cumulative % of var.  61.722  78.258  89.15

plot.PCA(extra.pca, axes = c(1,2), choix = "var") 
plot.PCA(extra.pca, axes = c(2,3), choix = "var") 


### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------