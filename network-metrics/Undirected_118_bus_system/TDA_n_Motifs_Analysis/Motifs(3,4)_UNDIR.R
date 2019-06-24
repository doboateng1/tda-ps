
rm(list = ls())

setwd("C:/Users/doforib/Desktop/Attack trials/Undirected_graph_analysis")

library(igraph)


# # pleae read 118_graph.gml and result.csv file before running program
#write_graph(case_118_network,"118_graph.gml", format = "gml")

new_case118_result1 = read.csv("result-1.csv")[,-2]

############################################################################################################################################

# UNDIRECTED GRAPH - SEQUENTIAL ATTACKS #
#: Constant max flow

data4 = read.csv("edge_label.csv")
data4 = as.matrix(data4)


weights1 = new_case118_result1[1,360:545]#=== max per-unit flow in line/transformer i ===# 186
weights2 = new_case118_result1[1,546:599]#=== max per-unit generation at generator i ===# 54
weights3 = new_case118_result1[1,600:698]#=== max per-unit consumption at load i ===# (some loads don't exist in the count) 99



weights1C = new_case118_result1[-1,360:545]#=== distinct per-unit flow in line/transformer i ===# 186
weights2C = new_case118_result1[-1,546:599]#=== distinct per-unit generation at generator i ===# 54
weights3C = new_case118_result1[-1,600:698]#=== distinct per-unit consumption at load i ===# (some loads don't exist in the count) 99


#=== Undirected graph ===#
mmh = graph_from_edgelist(data4[,1:2], directed = F)
#E(mmh)$label = data4[,3]
#E(mmh)$weight = weights1
igraph.options(vertex.size=7, edge.arrow.size=0.2, edge.color = "black",vertex.label=NA) 
plot(mmh, main = "graph")



data5 = new_case118_result1[,2:119]#== bus exists or not ===#
data6 = new_case118_result1[,120:305]#== flow line exists or not ===#
data7 = new_case118_result1[,306:359]#== whether load is in service or not ===#


##############################################################################################################################

#=== buses disconnected ===#
dd = matrix(NA, 119,118)
dob = data5

for(i in 1:119){
  
  mm = which(dob[(i+1),] == "FALSE")
  
  if(length(mm)>1){
    dd[i,1:length(mm)] = mm
  }else{ dd[i,length(mm)] = mm}
  
}


#=== flows disconnected by branch removal ===#
ff = matrix(NA, 119,186)
fob = data6

for(i in 1:119){
  
  mm = which(fob[(i+1),] == "FALSE")
  
  if(length(mm)>1){
    ff[i,1:length(mm)] = mm
  }else{ ff[i,length(mm)] = mm}
  
}



#=== loads disconnected by branch removal ===#
gg = matrix(NA, 119,54)
gob = data7

for(i in 1:119){
  
  mm = which(gob[(i+1),] == "FALSE")
  
  if(length(mm)>1){
    gg[i,1:length(mm)] = mm
  }else{ gg[i,length(mm)] = mm}
  
}


###############################################################################################################
##===========================================================================================================##
###############################################################################################################


n2 = c()
V_1 = V_2 = V_3 = V_4 = V_5 = V_6 = Tot_M = Tot_V = c()
C_V1 = C_V2 = C_V3 = C_V4 = C_V5 = C_V6 = c()


#=== 3-node motifs ===#
dob = data5

network_org = mmh
network = network_org


for(i in 1:81){
  
  m2 = motifs(network, 3)
  m2[is.na(m2)]  =  0
  
  n02 = count_motifs(network, 3)
  n2[i] =  n02
  
  V_1[i] = m2[3]  
  V_2[i] = m2[4] 
  
  
  Tot_V[i] = sum(m2)
  
  
  network = delete.vertices(network_org, c(which(dob[(i+1),] == "FALSE")) )
  
  
}


n22 = n2[1]

C_V1 = V_1/n22
C_V2 = V_2/n22


resultsN= data.frame(c(0:80),Tot_V,V_1,V_2,C_V1,C_V2)
write.csv(resultsN, "Motifs_UNDIR_3N_80.csv")

###############################################################################################################
##===========================================================================================================##
###############################################################################################################

n2 = c()
V_1 = V_2 = V_3 = V_4 = V_5 = V_6 = Tot_M = Tot_V = c()
C_V1 = C_V2 = C_V3 = C_V4 = C_V5 = C_V6 = c()



#=== 4-node motifs ===#
dob = data5

network_org = mmh
network = network_org


for(i in 1:81){
  
  m2 = motifs(network, 4)
  m2[is.na(m2)]  =  0
  
  n02 = count_motifs(network, 4)
  n2[i] =  n02
  
  V_1[i] = m2[5]  
  V_2[i] = m2[7] 
  V_3[i] = m2[8] 
  V_4[i] = m2[9]   
  V_5[i] = m2[10] 
  V_6[i] = m2[11]  
  
  Tot_V[i] = sum(m2)
  
  network = delete.vertices(network_org, c(which(dob[(i+1),] == "FALSE")) )
  
  
}


n22 = n2[1]

C_V1 = V_1/n22 
C_V2 = V_2/n22 
C_V3 = V_3/n22  
C_V4 = V_4/n22 
C_V5 = V_5/n22
C_V6 = V_6/n22


resultsN = data.frame(c(0:80),Tot_V,V_1,V_2,V_3,V_4,V_5,V_6,C_V1,C_V2,C_V3,C_V4,C_V5,C_V6)
write.csv(resultsN, "Motifs_UNDIR_4N_80.csv")





##############################################################################################################################

# DIRECTED GRAPH - SEQUENTIAL ATTACKS #
#: Constant max flow

data4 = read.csv("edge_label.csv")
data4 = as.matrix(data4)


weights1 = new_case118_result1[1,360:545]#=== max per-unit flow in line/transformer i ===# 186
weights2 = new_case118_result1[1,546:599]#=== max per-unit generation at generator i ===# 54
weights3 = new_case118_result1[1,600:698]#=== max per-unit consumption at load i ===# (some loads don't exist in the count) 99



weights1C = new_case118_result1[-1,360:545]#=== distinct per-unit flow in line/transformer i ===# 186
weights2C = new_case118_result1[-1,546:599]#=== distinct per-unit generation at generator i ===# 54
weights3C = new_case118_result1[-1,600:698]#=== distinct per-unit consumption at load i ===# (some loads don't exist in the count) 99


#=== Directed graph ===#
mmh = graph_from_edgelist(data4[,1:2], directed = T)
E(mmh)$label = data4[,3]
E(mmh)$weight = weights1
igraph.options(vertex.size=7, edge.arrow.size=0.2, edge.color = "black",vertex.label=NA) 
plot(mmh, main = "graph")



data5 = new_case118_result1[,2:119]#== bus exists or not ===#
data6 = new_case118_result1[,120:305]#== flow line exists or not ===#
data7 = new_case118_result1[,306:359]#== whether load is in service or not ===#

##############################################################################################################################



n2 = c()
V_1 = V_2 = V_3 = V_4 = V_5 = V_6 = Tot_M = Tot_V = c()
C_V1 = C_V2 = C_V3 = C_V4 = C_V5 = C_V6 = c()


#=== 3-node motifs ===#
dob = data5

network_org = mmh
network = network_org


for(i in 1:4){
  
  m2 = motifs(network, 3)
  m2[is.na(m2)]  =  0
  
  n02 = count_motifs(network, 3)
  n2[i] =  n02
  
  V_1[i] = m2[3]  
  V_2[i] = m2[4] 
  
  
  Tot_V[i] = sum(m2)
  
  
  network = delete.vertices(network_org, c(which(dob[(i+1),] == "FALSE")) )
  
  
}


n22 = n2[1]

C_V1 = V_1/n22
C_V2 = V_2/n22


resultsN= data.frame( Tot_V,V_1,V_2,C_V1,C_V2)



#=== 4-node motifs ===#
dob = data5

network_org = mmh
network = network_org


for(i in 1:4){
  
  m2 = motifs(network, 4)
  m2[is.na(m2)]  =  0
  
  n02 = count_motifs(network, 4)
  n2[i] =  n02
  
  V_1[i] = m2[5]  
  V_2[i] = m2[7] 
  V_3[i] = m2[8] 
  V_4[i] = m2[9]   
  V_5[i] = m2[10] 
  V_6[i] = m2[11]  
  
  Tot_V[i] = sum(m2)
  
  
  network = delete.vertices(network_org, c(which(dob[(i+1),] == "FALSE")) )
  
  
}


n22 = n2[1]

C_V1 = V_1/n22 
C_V2 = V_2/n22 
C_V3 = V_3/n22  
C_V4 = V_4/n22 
C_V5 = V_5/n22
C_V6 = V_6/n22


resultsN = data.frame(Tot_V,V_1,V_2,V_3,V_4,V_5,V_6,C_V1,C_V2,C_V3,C_V4,C_V5,C_V6)




##############################################################################################################################

