
rm(list = ls())

setwd("C:/Users/doforib/Desktop/Attack trials/Undirected_graph_analysis")


###################
#=== LIBRARIES ===#
###################

library(igraph)
library(TDA)
library(rio)


new_case118_result1 = read.csv("result-1.csv")[,-2]


data4 = import("branches.tsv", format = "csv") #=== Branch info only ===# 

new_data = as.matrix(data4[,c(2,3)])
gg = graph_from_edgelist(new_data, directed = F)


igraph.options(vertex.size=7, edge.arrow.size=0.5, edge.color = "black",vertex.label=NA) 
plot(gg, main = "nesta_case118_ieee(study-01)")





weights1 = as.numeric(new_case118_result1[1,360:545])#=== max per-unit flow in line/transformer i ===# 186
weights2 = new_case118_result1[1,546:599]#=== max per-unit generation at generator i ===# 54
weights3 = new_case118_result1[1,600:698]#=== max per-unit consumption at load i ===# (some loads don't exist in the count) 99



weights1C = new_case118_result1[-1,360:545]#=== distinct per-unit flow in line/transformer i ===# 186
weights2C = new_case118_result1[-1,546:599]#=== distinct per-unit generation at generator i ===# 54
weights3C = new_case118_result1[-1,600:698]#=== distinct per-unit consumption at load i ===# (some loads don't exist in the count) 99




#=== Undirected graph ===#
mmh = gg
E(mmh)$weight = as.numeric(abs(1/weights1))
igraph.options(vertex.size=7, edge.arrow.size=0.2, edge.color = "black",vertex.label=NA) 
plot(mmh, main = "graph")



data5 = new_case118_result1[,2:119]#== bus exists or not ===#
data6 = new_case118_result1[,120:305]#== flow line exists or not ===#
data7 = new_case118_result1[,306:359]#== whether load is in service or not ===#






##############################################################################################################################

#=== Constant weights/distances ===#

dob = data5

network_org = mmh
network = network_org


#=== TDA necessities ===#

cap = 3 # dimension cap
delta = 0.10
filt_len = 100

betti_0 = betti_1 = wasser_dist_01 = c()


for(i in 1:10){
  
  
  A1 =  as_adjacency_matrix(network, attr="weight")
  A2 = as.matrix(A1)
  
  A2[A2==0] = 999
  diag(A2)=0
  
  
  d = length(V(network))
  print(i)
  
  # writing data into file M.txt
  cat(d,file='M.txt',append=F,sep = '\n')
  cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') 
  cat(A2,file='M.txt',append = T) 
  
  system('perseusWin.exe distmat M.txt Moutput')
  
  # read Betti numbers from file Moutput_betti.txt
  
  betti_data = as.matrix(read.table('Moutput_betti.txt'))
  betti_index = setdiff(0:filt_len,betti_data[,1])
  
  for (k in betti_index) 
    if (k < length(betti_data[ ,1])) 
    {
      betti_data = rbind(betti_data[1:k, ], betti_data[k,], betti_data[(k+1):length(betti_data[,1]), ])
      betti_index = betti_index + 1
    } else
      betti_data = rbind(betti_data[1:k,], betti_data[k,])
  
  betti_0 = rbind(betti_0, betti_data[,2])
  betti_1 = rbind(betti_1, betti_data[,3])
  
  
  # read birth and death times for each dimension
  
  # dim = 0
  persist_data = as.matrix(read.table('Moutput_0.txt'))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P = cbind(rep(0, nrow(persist_data)), persist_data)
  
  # dim = 1
  if (file.info('Moutput_1.txt')$size>0)
  { 
    persist_data = as.matrix(read.table('Moutput_1.txt', blank.lines.skip = T))
    persist_data[persist_data[,2] == -1, 2] = filt_len + 1
    persist_data = persist_data/(filt_len + 1)
    P = rbind(P, cbind(rep(1, nrow(persist_data)), persist_data))
    
  }
  
  if (i == 1) P_org = P  
  
  wasser_dist_01 = c(wasser_dist_01, wasserstein(P_org, P, dimension = c(0,1)))
  
  network = delete.vertices(network_org, c(which(dob[(i+1),] == "FALSE")) )
  #E(network)$weight = weights1C[(1+i),] #== necessary for TDA and not motif analysis ===#
  
}

norm_const0 = norm(as.matrix(betti_0[1,]), type = '2')
betti_0_DistS = as.matrix(dist(betti_0))/norm_const0

norm_const1 = norm(as.matrix(betti_1[1,]),type = '2')
betti_1_DistS = as.matrix(dist(betti_1))/norm_const1

resultsN = data.frame(betti_0_DistS[1,], betti_1_DistS[1,], wasser_dist_01)




##############################################################################################################################

#: Changing flows/ edge weights

dob = data5

network_org = mmh
E(network_org) = abs(1/weights1C[1,])

network = network_org


#=== TDA necessities ===#

cap = 3 # dimension cap
delta = 0.10
filt_len = 100

betti_0 = betti_1 = wasser_dist_01 = dgf = c()
nodes_rem = matrix(NA, 10, 118)


for(i in 1:10){
  
  
  A1 =  as_adjacency_matrix(network, attr="weight")
  A2 = as.matrix(A1)
  
  A2[A2==0] = 999
  diag(A2)=0
  
  
  d = length(V(network))
  print(i)
  
  # writing data into file M.txt
  cat(d,file='M.txt',append=F,sep = '\n')
  cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') 
  cat(A2,file='M.txt',append = T) 
  
  system('perseusWin.exe distmat M.txt Moutput')
  
  # read Betti numbers from file Moutput_betti.txt
  
  betti_data = as.matrix(read.table('Moutput_betti.txt'))
  betti_index = setdiff(0:filt_len,betti_data[,1])
  
  for (k in betti_index) 
    if (k < length(betti_data[ ,1])) 
    {
      betti_data = rbind(betti_data[1:k, ], betti_data[k,], betti_data[(k+1):length(betti_data[,1]), ])
      betti_index = betti_index + 1
    } else
      betti_data = rbind(betti_data[1:k,], betti_data[k,])
  
  betti_0 = rbind(betti_0, betti_data[,2])
  betti_1 = rbind(betti_1, betti_data[,3])
  
  
  # read birth and death times for each dimension
  
  # dim = 0
  persist_data = as.matrix(read.table('Moutput_0.txt'))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P = cbind(rep(0, nrow(persist_data)), persist_data)
  
  # dim = 1
  if (file.info('Moutput_1.txt')$size>0)
  { 
    persist_data = as.matrix(read.table('Moutput_1.txt', blank.lines.skip = T))
    persist_data[persist_data[,2] == -1, 2] = filt_len + 1
    persist_data = persist_data/(filt_len + 1)
    P = rbind(P, cbind(rep(1, nrow(persist_data)), persist_data))
    
  }
  #which(weights1C!=0)
  
  if (i == 1) P_org = P  
  
  wasser_dist_01 = c(wasser_dist_01, wasserstein(P_org, P, dimension = c(0,1)))
  
  network = network_org
  E(network)$weight = 1/(abs(as.numeric(weights1C[(i+1),])) + 0.5)
  network = delete.vertices(network, c(which(dob[(i+1),] == "FALSE")) )
  
  dgf = which(dob[(i+1),] == "FALSE")
  
  if(length(dgf) >1){
      nodes_rem[i,1:length(dgf)] = dgf
  } else if(length(dgf) == "0"){
    nodes_rem[i,1] = 0
  }else{nodes_rem[i,1] = dgf}
  
  
    
}

norm_const0 = norm(as.matrix(betti_0[1,]), type = '2')
betti_0_DistS = as.matrix(dist(betti_0))/norm_const0

norm_const1 = norm(as.matrix(betti_1[1,]),type = '2')
betti_1_DistS = as.matrix(dist(betti_1))/norm_const1

resultsN = data.frame(betti_0_DistS[1,], betti_1_DistS[1,], wasser_dist_01)

nodes_rem

##############################################################################################################################






