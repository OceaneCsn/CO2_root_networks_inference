library(stringr)
library(DIANE)
library(tictoc)
library(ggplot2)
library(ggpubr)
library(patchwork)


############ This script perfoms GRN inference of the root CO2 response
############ In arabidopsis under the N*CO2 combinatorial design
############ It was done on a server allow parallel computing, and will
############ be much longer on a personal computer. It is thus recommended
############ Not to re-run these steps, but rather directly load the 
############ file in "Data/network_N_CO2.RData", saved at the end of this
############ script.


data <- read.csv("Data/raw_expression.csv", h = T, row.names = "Gene")
annotate_condition <- function(condition) {
  co2 <- substr(condition, 1, 1)
  nitrate <- substr(condition, 2, 2)
  iron <- substr(condition, 3, 3)
  if (co2 == 'c')
    res <- "ACO2"
  else
    res <- "ECO2"
  if (nitrate == 'n')
    res <- paste(res, "Low-Nitrate", sep = '.')
  else
    res <- paste(res, "High-Nitrate", sep = '.')
  if (iron == 'f')
    res <- paste(res, "Iron-starvation", sep = '.')
  else
    res <- paste(res, "Iron-supply", sep = '.')
  return(res)
}

colnames(data) <-
  paste0(
    sapply(str_split_fixed(colnames(data), '_', 2)[, 1],
           annotate_condition),
    '_',
    str_split_fixed(colnames(data), '_', 2)[, 2]
  )

########## normalization with tmm and filtering < 10 avg counts

tcc_object <-
  DIANE::normalize(data, str_split_fixed(colnames(data), '_', 2)[, 1], iteration = FALSE)
threshold = 10 * ncol(data)
tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
normalized_counts <- TCC::getNormalizedData(tcc_object)
normalized_counts <-
  normalized_counts[, order(colnames(normalized_counts))]


###################################### DEA

fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions =
                                   str_split_fixed(colnames(data), '_', 2)[, 1])

topTags <-
  DIANE::estimateDEGs(
    fit,
    reference = "ACO2.Low-Nitrate.Iron-supply",
    perturbation = "ECO2.Low-Nitrate.Iron-supply",
    p.value = 0.05,
    lfc = 0
  )
genes <- topTags$table$genes


#removing iron starvation
normalized_counts <-
  normalized_counts[, !str_detect(colnames(normalized_counts), 'starvation')]


###################### Gene grouping

data("regulators_per_organism")

regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])

DIANE::draw_PCA(data = normalized_counts)

r <- group_regressors(genes = genes, normalized.count = 
                        normalized_counts[genes,], regressors = regressors, 
                      corr_thr = 0.95)

visNetwork::visNetwork(r$correlated_regressors_graph$nodes,
                       r$correlated_regressors_graph$edges)

data <- r$counts
grouped_genes <- r$grouped_genes
grouped_regressors <- r$grouped_regressors


#################### Inference

tic()
set.seed(999)
mat <- DIANE::network_inference(data, str_split_fixed(colnames(data), '_', 2)[,1], 
                                regressors = grouped_regressors, 
                                targets = grouped_genes, nTrees = 4000, nCores = 6)
toc()

nGenes = length(grouped_genes)
nRegulators = length(grouped_regressors)
density = 0.01
get_nEdges(density, nGenes, nRegulators)

#testing


DIANE::estimate_test_edges_time(mat, normalized_counts = data, density = density,
                                nGenes = nGenes,
                                nRegulators = nRegulators, nShuffle = 1000,
                                nTrees = 4000, verbose = FALSE, nCores = 32)

tic()
set.seed(999)
res <- DIANE::test_edges(mat, normalized_counts = data, density = density,
                         nGenes = nGenes,
                         nRegulators = nRegulators, nShuffle = 1000,
                         nTrees = 4000, verbose = FALSE, nCores = 32)
toc()
######################## results

res$fdr_nEdges_curve

res$pvalues_distributions + xlim(0,0.1)

net <- DIANE::network_from_tests(res$links, fdr = 0.05) 

net_data <- network_data(net, regulators_per_organism[["Arabidopsis thaliana"]], 
                         gene_annotations$`Arabidopsis thaliana`)

draw_discarded_edges(res$links, net_data)
draw_network(net_data$nodes, net_data$edges)
DIANE::draw_network_degrees(net_data$nodes, net)
ranking <- net_data$nodes[order(-net_data$nodes$degree),]

network_N_CO2 <- list(edges_testing_results = res, 
                      network = net, 
                      network_data = net_data,
                      ranking = ranking)
save(network_N_CO2, file = "Results/network_N_CO2.RData")




##### writing tables of the targets of three genes of interest (DIANE is required)

# myb15 <- describe_node(network_N_CO2$network, "AT3G23250")$targets
# edf3 <- describe_node(network_N_CO2$network, "AT3G25730")$targets
# wox11 <- describe_node(network_N_CO2$network, "AT3G03660")$targets
# 
# 
# write.table(na.omit(get_gene_information(myb15, organism = "Arabidopsis thaliana")), 
#             file = "Results/CiblesMYB15.tsv", quote = F, sep = '\t')
# 
# write.table(na.omit(get_gene_information(wox11, organism = "Arabidopsis thaliana")), 
#             file = "Results/CiblesWOX11.tsv", quote = F, sep = '\t')
# 
# write.table(na.omit(get_gene_information(edf3, organism = "Arabidopsis thaliana")), 
#             file = "Results/CiblesEDF3.tsv", quote = F, sep = '\t')
