
## Script to create clusters for time data based on differentially expressed genes.
## Ashfaq Ali (ashfaq.ali@nbis.se)
## Activate R environment
renv::init()
renv::snapshot()
library(limma)
library(tidyverse)
library(ggthemes)


## Load data

load("res.rda")
load("normalised.rda")

### Extract the significant proteins

results <- list(resf1, resf2, resf3, resf4, resf5, resf6)
prot_list <- lapply(results, function(x) rownames(subset(x, (adj.P.Val< 0.05) & abs(logFC) > 0) ))
prot_list <- unique(unlist(prot_list) )


# read the data

prot_int <- data_prot_filt %>% 
            as.data.frame(row.names = rownames(data_prot_filt)) %>% 
            rownames_to_column(var = "protein")


sample_info <- data.frame(samples = colnames(data_prot_filt), time = rep(1:4,each=4))


# Summarise intesities by timepoint  
prot_int_mean <- prot_int %>% 
    # convert to long format
    pivot_longer(cols = sample_info$samples, names_to = "samples", values_to = "intensity")  %>% # transform data to long format
    # join with sample info table
    full_join(sample_info, by = ("samples"),copy =TRUE) %>% # Join to import sample annotations
    # for each protein create a group
    group_by(protein) %>% # 
    # scale the intensities column
    mutate(cts_scaled = (intensity - mean(intensity))/sd(intensity)) %>% 
    # for each protein and timpoint
    group_by(protein, time) %>%
    # calculate the mean (scaled) intensity
    summarise(mean_cts_scaled = mean(cts_scaled), 
              nrep = dplyr::n()) %>% 
    ungroup() 

hclust_matrix <- data_prot_filt[prot_list,] %>% t() %>% 
    # apply scaling to each column of the matrix (proteins)
    scale() %>% 
    # transpose back so proteins are as rows again
    t()
prot_dist <- dist(hclust_matrix) ## cal

prot_hclust <- hclust(prot_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
png("./Cluster_dand.png", res = 300, width =10, height = 10, units = "cm")
plot(prot_hclust, labels = FALSE)
abline(h = 4.6, col = "dred", lwd = 2)
dev.off()

prot_cluster <- cutree(prot_hclust, h = 4, 12 ) %>% 
    # turn the named vector into a tibble
    enframe() %>% 
    # rename some of the columns
    rename(protein = name, cluster = value)

#head(prot_cluster)

## Save summary data to file
prot_int_cluster <- prot_int_mean %>% 
    inner_join(prot_cluster, by = "protein")

writexl::write_xlsx(prot_int_cluster,"./Clusters_Time.xlsx")

### Plot the clusters 

plot <- prot_int_cluster %>% 
    ggplot(aes(time, mean_cts_scaled, color=factor(cluster))) +
    geom_line(aes(group = protein))+
    geom_smooth(stat = "smooth", size=2)+
    theme_pander() +
    labs(
        colour = "Cluster ID"
    )  +
    theme(
        #legend.position = c(.9, .7),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    ) +
    scale_y_continuous("Mean Cluster Intensity (Scaled)") +
    scale_x_continuous("Time point 1 to 4") +
    
    ggtitle("Trajectories of protein clusters")

ggsave(plot, file=paste("./",
                        "Clusters", ".png", sep=''), 
       scale=2,
       units = "cm", height = 15, width = 20)

# save plots as .pdf
ggsave(plot, file=paste("./",
                        "Clusters", ".pdf", sep=''), 
       scale=2,
       units = "cm", height = 15, width = 20)

## Create a separate plot for individual clusters
plot1 <- prot_int_cluster %>% 
    ggplot(aes(time, mean_cts_scaled, color=factor(cluster))) +
    geom_line(aes(group = protein))+
    geom_smooth(stat = "smooth", size=2)+
    theme_pander() +
    labs(
        colour = "Cluster ID"
    )  +
    theme(
        #legend.position = c(.9, .7),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    ) +
    scale_y_continuous("Mean Cluster Intensity (Scaled)") +
    scale_x_continuous("Time point 1 to 4") +
    
    ggtitle("Trajectories of protein clusters")+    facet_grid(cols  = vars(cluster))

ggsave(plot1, file=paste("./",
                        "Clusters_invidual", ".png", sep=''), 
       scale=2,
       units = "cm", height = 15, width = 40)

# save plots as .pdf
ggsave(plot1, file=paste("./",
                        "Clusters_individual", ".pdf", sep=''), 
       scale=2,
       units = "cm", height = 15, width = 40)
