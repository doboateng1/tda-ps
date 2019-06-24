
##############################################################################################################################################

rm(list = ls())


setwd("C:/Users/doforib/Desktop/Attack trials")

source("Resil_Motiffunc.R")
source("Resilience_func.R")



###################
#=== LIBRARIES ===#
###################

library(igraph)
library(TDA)
library(rio)


###################
#=== CONSTANTS ===#
###################

frac = seq(0,0.8, by = 0.1)


data4 = import("branches.tsv", format = "csv") #=== Branch info only ===# 

new_data = as.matrix(data4[,c(2,3)])
gg = graph_from_edgelist(new_data, directed = T)


igraph.options(vertex.size=7, edge.arrow.size=0.5, edge.color = "black",vertex.label=NA) 
plot(gg, main = "nesta_case118_ieee(study-01)")



#=== Edge weights ===#

data5 = import("result-1.tsv", format = "csv")

weights = as.numeric(data5[1, c(361:546)]) #mas-flow in each line

E(gg)$weight = weights



##############################################################################################################################################

################
#=== MOTIFS ===#
################
 
#=== 1. NODE ===#

# Degree #
NodeM_Deg = Resilience_3Motifs(attack = "node", directed = "TRUE", graph_type = "weighted", attack_type = "degree", gg, frac)
rownames(NodeM_Deg) = strrep(" ", 1 : nrow(NodeM_Deg))
NodeM_DegMatrix = as.matrix(NodeM_Deg)

# Betweenness Centrality #
NodeM_Betw = Resilience_3Motifs(attack = "node", directed = "TRUE", graph_type = "weighted", attack_type = "V_betweeness", gg, frac)
rownames(NodeM_Betw) = strrep(" ", 1 : nrow(NodeM_Betw))
NodeM_BetwMatrix = as.matrix(NodeM_Betw)

# Strength #
NodeM_Str = Resilience_3Motifs(attack = "node", directed = "TRUE", graph_type = "weighted", attack_type = "strength", gg, frac)
rownames(NodeM_Str) = strrep(" ", 1 : nrow(NodeM_Str))
NodeM_StrMatrix = as.matrix(NodeM_Str)


#=== 2. EDGE ===#

# Edge Betweenness Centrality #
EdgeM_Betw = Resilience_3Motifs(attack = "edge", directed = "TRUE", graph_type = "weighted", attack_type = "E_betweeness", gg, frac)
rownames(EdgeM_Betw) = strrep(" ", 1 : nrow(EdgeM_Betw))
EdgeM_BetwMatrix = as.matrix(EdgeM_Betw)

# Plain weights #
EdgeM_Wei = Resilience_3Motifs(attack = "edge", directed = "TRUE", graph_type = "weighted", attack_type = "weight_hier", gg, frac)
rownames(EdgeM_Wei) = strrep(" ", 1 : nrow(EdgeM_Wei))
EdgeM_WeiMatrix = as.matrix(EdgeM_Wei)



##############################################################################################################################################


#############
#=== TDA ===#
#############


# Degree #
NodeTDA_Deg = Resilience_Attacks2(attack = "node", graph_type = "weighted", analysis = "TDA", type = "degree", gg, frac)
NodeTDA_DegMatrix = as.matrix(NodeTDA_Deg)
colnames(NodeTDA_DegMatrix) = c("fr", "Betti0", "Betti1","Wass01")

NodeTDA_Betw = Resilience_Attacks2(attack = "node", graph_type = "weighted", analysis = "TDA", type = "V_betweenness", gg, frac)
NodeTDA_BetwMatrix = as.matrix(NodeTDA_Betw)
colnames(NodeTDA_BetwMatrix) = c("fr", "Betti0", "Betti1","Wass01")

NodeTDA_Str = Resilience_Attacks2(attack = "node", graph_type = "weighted", analysis = "TDA", type = "strength", gg, frac)
NodeTDA_StrMatrix = as.matrix(NodeTDA_Str)
colnames(NodeTDA_StrMatrix) = c("fr", "Betti0", "Betti1","Wass01")




# E_Betw #
EdgeTDA_Betw= Resilience_Attacks2(attack = "edge", graph_type = "weighted", analysis = "TDA", type = "E_betweenness", gg, frac)
EdgeTDA_BetwMatrix = as.matrix(EdgeTDA_Betw)
colnames(EdgeTDA_BetwMatrix) = c("fr", "Betti0", "Betti1","Wass01")

EdgeTDA_Wei= Resilience_Attacks2(attack = "edge", graph_type = "weighted", analysis = "TDA", type = "weight_hier", gg, frac)
EdgeTDA_WeiMatrix = as.matrix(EdgeTDA_Wei)
colnames(EdgeTDA_WeiMatrix) = c("fr", "Betti0", "Betti1","Wass01")


##############################################################################################################################################



rm(list = ls())


setwd("C:/Users/doforib/Desktop/Attack trials")

source("Resil_Motiffunc.R")
source("Resilience_func.R")


###################
#=== LIBRARIES ===#
###################

library(igraph)
library(TDA)


###################
#=== CONSTANTS ===#
###################

frac = seq(0,0.3, by = 0.1)

nd = 100 # node
ed = 200# edge


edge_weights = runif(ed)

##############################################################################################################################################

#####################
#=== GRAPH TYPES ===#
#####################

# 1. Erdos-Renyi #

ER_G           = sample_gnm(nd, ed, directed = F, loops = FALSE)
E(ER_G)$weight = edge_weights


##############################################################################################################################################

#######################
#=== 3-NODE MOTIFS ===#
#######################

#=== 1. NODE ===#

# Degree #
NodeM_Deg = Resilience_3Motifs(attack = "node", directed = "FALSE", graph_type = "weighted", attack_type = "degree", ER_G, frac)
rownames(NodeM_Deg) = strrep(" ", 1 : nrow(NodeM_Deg))
NodeM_DegMatrix = as.matrix(NodeM_Deg)

# Betweenness Centrality #
NodeM_Betw = Resilience_3Motifs(attack = "node", directed = "FALSE", graph_type = "weighted", attack_type = "V_betweeness", ER_G, frac)
rownames(NodeM_Betw) = strrep(" ", 1 : nrow(NodeM_Betw))
NodeM_BetwMatrix = as.matrix(NodeM_Betw)

# Strength #
NodeM_Str = Resilience_3Motifs(attack = "node", directed = "FALSE", graph_type = "weighted", attack_type = "strength", ER_G, frac)
rownames(NodeM_Str) = strrep(" ", 1 : nrow(NodeM_Str))
NodeM_StrMatrix = as.matrix(NodeM_Str)


#=== 2. EDGE ===#

# Edge Betweenness Centrality #
EdgeM_Betw = Resilience_3Motifs(attack = "edge", directed = "FALSE", graph_type = "weighted", attack_type = "E_betweeness", ER_G, frac)
rownames(EdgeM_Betw) = strrep(" ", 1 : nrow(EdgeM_Betw))
EdgeM_BetwMatrix = as.matrix(EdgeM_Betw)

# Plain weights #
EdgeM_Wei = Resilience_3Motifs(attack = "edge", directed = "FALSE", graph_type = "weighted", attack_type = "weight_hier", ER_G, frac)
rownames(EdgeM_Wei) = strrep(" ", 1 : nrow(EdgeM_Wei))
EdgeM_WeiMatrix = as.matrix(EdgeM_Wei)



##############################################################################################################################################
#============================================================================================================================================#
##############################################################################################################################################


#######################
#=== 4-NODE MOTIFS ===#
#######################

#=== 1. NODE ===#

# Degree #
NodeM_Deg = Resilience_4Motifs(attack = "node", directed = "FALSE", graph_type = "weighted", attack_type = "degree", ER_G, frac)
rownames(NodeM_Deg) = strrep(" ", 1 : nrow(NodeM_Deg))
NodeM_DegMatrix = as.matrix(NodeM_Deg)

# Betweenness Centrality #
NodeM_Betw = Resilience_4Motifs(attack = "node", directed = "FALSE", graph_type = "weighted", attack_type = "V_betweeness", ER_G, frac)
rownames(NodeM_Betw) = strrep(" ", 1 : nrow(NodeM_Betw))
NodeM_BetwMatrix = as.matrix(NodeM_Betw)

# Strength #
NodeM_Str = Resilience_4Motifs(attack = "node", directed = "FALSE", graph_type = "weighted", attack_type = "strength", ER_G, frac)
rownames(NodeM_Str) = strrep(" ", 1 : nrow(NodeM_Str))
NodeM_StrMatrix = as.matrix(NodeM_Str)


#=== 2. EDGE ===#

# Edge Betweenness Centrality #
EdgeM_Betw = Resilience_4Motifs(attack = "edge", directed = "FALSE", graph_type = "weighted", attack_type = "E_betweeness", ER_G, frac)
rownames(EdgeM_Betw) = strrep(" ", 1 : nrow(EdgeM_Betw))
EdgeM_BetwMatrix = as.matrix(EdgeM_Betw)

# Plain weights #
EdgeM_Wei = Resilience_4Motifs(attack = "edge", directed = "FALSE", graph_type = "weighted", attack_type = "weight_hier", ER_G, frac)
rownames(EdgeM_Wei) = strrep(" ", 1 : nrow(EdgeM_Wei))
EdgeM_WeiMatrix = as.matrix(EdgeM_Wei)


##############################################################################################################################################


#############
#=== TDA ===#
#############


# Degree #
NodeTDA_Deg = Resilience_Attacks2(attack = "node", graph_type = "weighted", analysis = "TDA", type = "degree", ER_G, frac)
NodeTDA_DegMatrix = as.matrix(NodeTDA_Deg)
colnames(NodeTDA_DegMatrix) = c("fr", "Betti0", "Betti1","Wass01")

NodeTDA_Betw = Resilience_Attacks2(attack = "node", graph_type = "weighted", analysis = "TDA", type = "V_betweenness", ER_G, frac)
NodeTDA_BetwMatrix = as.matrix(NodeTDA_Betw)
colnames(NodeTDA_BetwMatrix) = c("fr", "Betti0", "Betti1","Wass01")

NodeTDA_Str = Resilience_Attacks2(attack = "node", graph_type = "weighted", analysis = "TDA", type = "strength", ER_G, frac)
NodeTDA_StrMatrix = as.matrix(NodeTDA_Str)
colnames(NodeTDA_StrMatrix) = c("fr", "Betti0", "Betti1","Wass01")




# E_Betw #
EdgeTDA_Betw= Resilience_Attacks2(attack = "edge", graph_type = "weighted", analysis = "TDA", type = "E_betweenness", ER_G, frac)
EdgeTDA_BetwMatrix = as.matrix(EdgeTDA_Betw)
colnames(EdgeTDA_BetwMatrix) = c("fr", "Betti0", "Betti1","Wass01")

EdgeTDA_Wei= Resilience_Attacks2(attack = "edge", graph_type = "weighted", analysis = "TDA", type = "weight_hier", ER_G, frac)
EdgeTDA_WeiMatrix = as.matrix(EdgeTDA_Wei)
colnames(EdgeTDA_WeiMatrix) = c("fr", "Betti0", "Betti1","Wass01")



##############################################################################################################################################











