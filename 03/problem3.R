# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("limma")
# BiocManager::install("leukemiasEset")


library("leukemiasEset")
library("shiny")
library("e1071")
library("limma")
library("gplots")

# a) load DataSet 
df <- data("leukemiasEset")

# b) extract expressionMatrix
e_matrix = exprs(leukemiasEset)
# c) Phenodata
p_data <- pData(leukemiasEset)
print(summary(p_data))

# d) + f)
heatmap.2(e_matrix)


# e)

design <- model.matrix(~ 0 + factor(p_data$LeukemiaType))
colnames(design) <- levels(p_data$LeukemiaType)

fit <- lmFit(leukemiasEset, design)

# None vs any leukemia
contrast.matrix <- makeContrasts(NoL - (ALL + AML + CLL + CML)/4, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")
results <- decideTests(fit2)
vennDiagram(results)

# Acute vs chronic 
contrast.matrix2 <- makeContrasts((ALL + AML)/2 - (CLL + CML)/2, levels = design)

fit3 <- contrasts.fit(fit, contrast.matrix2)
fit3 <- eBayes(fit3)

topTable(fit3, adjust="BH")
results <- decideTests(fit3)
vennDiagram(results)
