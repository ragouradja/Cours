


```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

```{r}
load("~/cours/Cours/genomics/TP_Genomics/Data_HCA/expression_matrix_31T_3N_samples.RData")
```


```{r}
dim(exp)
```
54675 probes (sondes = gene)


```{r}

head(exp, 5)
```


```{r}
most_variant = apply(exp, 1, sd)
```

```{r}
head(most_variant,5)
```


```{r}
all_var = sort(most_variant)
```

```{r}
most_variant_500 = tail(all_var, 500)
```


```{r}
names_variant = names(most_variant_500)
```

```{r}
mat = exp[names_variant,]
mat
```

```{r}
mean_mat = apply(mat, 1, mean)
```



```{r}
mat_center = mat - mean_mat
```


```{r}
dist_mat = dist(t(mat_center))
```


```{r}
cluster = hclust(dist_mat, method = "ward.D")
```


```{r}
plot(cluster)
```

groupe de tumeur
consensus : 
```{r}
library(ConsensusClusterPlus)
```

```{r}
pdf("plot.pdf")
cons_mat_plot = ConsensusClusterPlus(d = mat_center, clusterAlg = "hc", distance = "euclidean", innerLinkage = 'ward.D', finalLinkage = "ward.D", maxK = 10, pFeature = 0.8, pItem = 0.8)
dev.off()
```
5 groupes selon le graph (endroit où ça casse)



```{r}
col_group = cons_mat_plot[[6]]$consensusClass
name_col_group = names(col_group)
```


```{r}
load("~/cours/Cours/genomics/TP_Genomics/Data_HCA/Clinical_annotations.RData")
```


```{r}
annot
```


```{r}
expGroupBis = rep(NA, 65)
expGroupBis
```

```{r}
new_annot = cbind(annot , expGroupBis)
```


```{r}
new_annot[name_col_group,14] = col_group
col_group
```

```{r}
temp_dtf = data.frame(annot, row.names = annot$Sample.ID)
temp_dtf = cbind(temp_dtf, expGroupBis)
temp_dtf
```


```{r}
temp_dtf[name_col_group,14] = col_group
```


```{r}
final_annot = cbind(annot,temp_dtf$expGroupBis)
```

```{r}
rownames(temp_dtf) = rownames(annot)
```


```{r}
temp_dtf
```

```{r}
table(temp_dtf[,"Molecular.group"])
```

```{r}
table(temp_dtf[,"expGroupBis"])
```
```{r}
chisq.test(temp_dtf$Molecular.group, temp_dtf$expGroupBis)
```
Les deux distributions sont différentes car la pvalue est inférieure à 5% (1.278e-15)


```{r}
library(ComplexHeatmap)
```




```{r}
Heatmap(mat_center,cluster_columns = cluster)
```


```{r}
pdf("heatmap_dendo.pdf")
Heatmap(cluster$merge)
dev.off()
```


```{r}
consensus_tree= cons_mat_plot[[6]]$consensusTree
```

attention, utiliser les cluster comme option et mat center
```{r}
pdf("heatmap_consensus.pdf")
Heatmap(mat_center, cluster_columns = consensus_tree)
dev.off()
```
```{r}
ifelse(temp_dtf[1], 
```


```{r}
dtf_annotation = data.frame(cbind(temp_dtf$Molecular.group, temp_dtf$Pathological.diagnosis))
colnames(dtf_annotation) = c("Molecular.group","Pathological.diagnosis")
dtf_annotation

dtf_annotation <- annot[match(colnames(mat),annot$Sample.ID), c("Molecular.group","Pathological.diagnosis")]
dtf_annotation
```


```{r}
top_annot = HeatmapAnnotation(df = dtf_annotation)
```


```{r}
pdf("heatmap_annot.pdf")
Heatmap(mat_center, top_annotation = top_annot, cluster_columns = cluster)
dev.off()
```
# PCA
prendre les 4 premieres
```{r}
comp = prcomp(t(mat_center), center = F, scale = F)
#mycol = as.numeric(annot[match(colnames(mat),annot$Sample.ID), "expGroupBis")])
summary(comp)
plot(comp$x[,1],comp$x[,2])
barplot(summary(comp)$importance["Proportion of variance",], names.arg = colnames(comp$x))
```

