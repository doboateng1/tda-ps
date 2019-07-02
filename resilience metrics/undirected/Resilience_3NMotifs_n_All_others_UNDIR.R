

rm(list = ls())

start.time = Sys.time()


############################################################################################################################

setwd("C:/Users/Owner/OneDrive/Desktop/NREL")



library(igraph)
library(TDA)


# # pleae read 118_graph.gml and result.csv file before running program
#write_graph(case_118_network,"118_graph.gml", format = "gml")
new_case118_result1 = read.csv("result-1.csv")[,-2]

input_graph = read.graph("118_graph.gml",format=c("gml")) # 118_graph.gml is built based on graph.dot/graph.svg #


#=== LIBRARIES ===#



#=== FUNCTION ===#

output_graph_under_attack_f = function(case_118_network = input_graph, initial_status_row = 2, sequence_row){
  
  # 1, import graph.dot graph #
  case_118_network_name_is_label = case_118_network # corresponding to input_graph
  V(case_118_network_name_is_label)$name = V(case_118_network_name_is_label)$label
  edgelist_with_label = get.edgelist(case_118_network_name_is_label)
  edgelist_with_label_branch = cbind(edgelist_with_label, edge_attr(case_118_network_name_is_label)$label)
  colnames(edgelist_with_label_branch) = c("From_Bus","To_Bus","Branch")
  edgelist_with_label_branch = cbind(edgelist_with_label_branch, paste("F_",edgelist_with_label_branch[,3],sep = ""))
  colnames(edgelist_with_label_branch) = c("From_Bus","To_Bus","Branch","F_Branch")
  
  # 2, import initial status from result-1.tsv when sequence = 0 #
  new_case118_result1_seq0 = new_case118_result1[initial_status_row,] # where new_case118_result1 is result-1.tsv; corresponding to initial_status
  new_case118_result1_seq0_Flow_direction_info = new_case118_result1_seq0[which(colnames(new_case118_result1)=="F_1"):which(colnames(new_case118_result1)=="F_186")] # len is 186
  initial_direction_sign = sign(new_case118_result1_seq0_Flow_direction_info)
  
  # 3, assign the F_sign when sequence = 0 to edgelist_with_label_branch generated from graph.dot #
  sign_vector = vector(length = dim(edgelist_with_label_branch)[1])
  for (i in 1:dim(edgelist_with_label_branch)[1]) {
    sign_vector[i] = initial_direction_sign[1,
                                            paste("F_",which(names(initial_direction_sign) == edgelist_with_label_branch[i,4]),sep = "")]
  }
  
  edgelist_with_label_branch_sign = cbind(edgelist_with_label_branch, sign_vector)
  colnames(edgelist_with_label_branch_sign) = c("From_Bus","To_Bus","Branch","F_Branch","F_Sign")
  
  # 4, groundtruth initial edgelist info with branch label, bus label, and flow sign; 
  #that is change the direction with negative compared with the edgelist from graph.dot #
  edgelist_with_label_branch_sign_df = as.data.frame(edgelist_with_label_branch_sign)
  
  
  initial_from_to_bus = data.frame(From_Bus_initial = rep(0,dim(edgelist_with_label_branch)[1]),
                                   To_Bus_Initial = rep(0,dim(edgelist_with_label_branch)[1]),Branch = rep(0,dim(edgelist_with_label_branch)[1]),
                                   Branch_num = rep(0,dim(edgelist_with_label_branch)[1]),F_Sig_zero = rep(0,dim(edgelist_with_label_branch)[1]))
  zero_flow_from_to_bus = data.frame(From_Bus_initial = rep(0,4),
                                     To_Bus_Initial = rep(0,4),Branch = rep(0,4),
                                     Branch_num = rep(0,4))
  a = 1
  for (j in 1:dim(edgelist_with_label_branch)[1]) {
    if(edgelist_with_label_branch_sign_df[j,5] == 1){
      initial_from_to_bus[j,1] = edgelist_with_label_branch_sign[j,1]
      initial_from_to_bus[j,2] = edgelist_with_label_branch_sign[j,2]
      initial_from_to_bus[j,3] = edgelist_with_label_branch_sign[j,4]
      initial_from_to_bus[j,4] = edgelist_with_label_branch_sign[j,3]
      initial_from_to_bus[j,5] = edgelist_with_label_branch_sign[j,5]
    }else if(edgelist_with_label_branch_sign_df[j,5] == -1){
      initial_from_to_bus[j,1] = edgelist_with_label_branch_sign[j,2]
      initial_from_to_bus[j,2] = edgelist_with_label_branch_sign[j,1]
      initial_from_to_bus[j,3] = edgelist_with_label_branch_sign[j,4]
      initial_from_to_bus[j,4] = edgelist_with_label_branch_sign[j,3]
      initial_from_to_bus[j,5] = edgelist_with_label_branch_sign[j,5]
    }else{
      initial_from_to_bus[j,1] = edgelist_with_label_branch_sign[j,1]
      initial_from_to_bus[j,2] = edgelist_with_label_branch_sign[j,2]
      initial_from_to_bus[j,3] = edgelist_with_label_branch_sign[j,4]
      initial_from_to_bus[j,4] = edgelist_with_label_branch_sign[j,3]
      initial_from_to_bus[j,5] = edgelist_with_label_branch_sign[j,5]
      zero_flow_from_to_bus[a,1] = edgelist_with_label_branch_sign[j,1]
      zero_flow_from_to_bus[a,2] = edgelist_with_label_branch_sign[j,2]
      zero_flow_from_to_bus[a,3] = edgelist_with_label_branch_sign[j,4]
      zero_flow_from_to_bus[a,4] = edgelist_with_label_branch_sign[j,3]
      a = a+1
    }
  }
  
  
  # 5, generate seq_zero_case118_ieee_network with correct edge direction when sequence = 0 and then remove the edges with F_=0 when sequence = 0
  seq_zero_case118_ieee_network = graph_from_edgelist(as.matrix(initial_from_to_bus[,c(1:2)]))
  # IGRAPH dc75189 DN-- 118 186 --#
  seq_zero_case118_ieee_network = seq_zero_case118_ieee_network %>% set_edge_attr("name", value = initial_from_to_bus[,3])
  # delete the edges with F_ - 0 #
  seq_zero_case118_ieee_network = delete_edges(seq_zero_case118_ieee_network, zero_flow_from_to_bus[,3])
  # IGRAPH dc75189 DN-- 118 182 --#
  
  # 6, assign name, label, color, feature in case_118_network to seq_zero_case118_ieee_network
  four_attrs_case_118_network = data.frame(name = vertex_attr(case_118_network)$name,
                                           label = vertex_attr(case_118_network)$label,
                                           color = vertex_attr(case_118_network)$color,
                                           feature = vertex_attr(case_118_network)$feature)
  
  
  four_attrs_seq_zero_case118_ieee_network = data.frame(name = rep(0,length(vertex_attr(seq_zero_case118_ieee_network)$name)),
                                                        label = vertex_attr(seq_zero_case118_ieee_network)$name,
                                                        color = rep(0,length(vertex_attr(seq_zero_case118_ieee_network)$name)),
                                                        feature = rep(0,length(vertex_attr(seq_zero_case118_ieee_network)$name)))
  
  for (kk in 1:dim(four_attrs_seq_zero_case118_ieee_network)[1]) {
    tmp_label = which(vertex_attr(case_118_network)$label == vertex_attr(seq_zero_case118_ieee_network)$name[kk])
    four_attrs_seq_zero_case118_ieee_network[kk,1] = vertex_attr(case_118_network)$name[tmp_label]
    four_attrs_seq_zero_case118_ieee_network[kk,3] = vertex_attr(case_118_network)$color[tmp_label]
    four_attrs_seq_zero_case118_ieee_network[kk,4] = vertex_attr(case_118_network)$feature[tmp_label]
  }
  
  four_attrs_seq_zero_case118_ieee_network$name = four_attrs_seq_zero_case118_ieee_network$name #paste("bus",four_attrs_seq_zero_case118_ieee_network$name,sep = "")
  four_attrs_seq_zero_case118_ieee_network_dataframe = four_attrs_seq_zero_case118_ieee_network
  four_attrs_seq_zero_case118_ieee_network_mat = as.matrix(four_attrs_seq_zero_case118_ieee_network_dataframe)
  all.equal(four_attrs_seq_zero_case118_ieee_network_mat[,2],vertex_attr(seq_zero_case118_ieee_network)$name) #TRUE
  
  V(seq_zero_case118_ieee_network)$feature = as.numeric(four_attrs_seq_zero_case118_ieee_network_mat[,4])
  V(seq_zero_case118_ieee_network)$name = four_attrs_seq_zero_case118_ieee_network_mat[,1]
  V(seq_zero_case118_ieee_network)$label = four_attrs_seq_zero_case118_ieee_network_mat[,2]
  V(seq_zero_case118_ieee_network)$color = four_attrs_seq_zero_case118_ieee_network_mat[,3]
  # IGRAPH dc75189 DN-- 118 182 -- #
  
  initial_from_to_bus$F_Sig_zero = as.numeric(initial_from_to_bus$F_Sig_zero)
  initial_from_to_bus$Branch_f = paste("f_",initial_from_to_bus[,4],sep = "") 
  
  # 7, example on sequence = 17; below 19 corresponding to variable - sequence_row
  F_sign_seq_under_attack = sign(new_case118_result1[sequence_row,c(which(colnames(new_case118_result1) == "F_1") : which(colnames(new_case118_result1) == "F_186"))])
  
  # 8, obtain the nodes' labels with b_ is false when sequence = 19
  result_seq_under_attack_nodes_info = new_case118_result1[sequence_row,c(which(colnames(new_case118_result1) == "b_1") : which(colnames(new_case118_result1) == "b_118"))]
  removed_nodes_label = names(result_seq_under_attack_nodes_info)[which(result_seq_under_attack_nodes_info*1 == 0)]
  
  # rename edges' names for seq_zero_case118_ieee_network use f_ not F_#
  length(edge_attr(seq_zero_case118_ieee_network)$name)
  edge_label_num = as.numeric(gsub("F_", "", edge_attr(seq_zero_case118_ieee_network)$name)) # extract the number in edge_attr(seq_zero_case118_iee_network)$name; since i want to rebuild the f_ for edges' names
  E(seq_zero_case118_ieee_network)$name = paste("f_",edge_label_num, sep="")
  # rename complete #
  
  # 9, generage comb_seq0_seq_under_attack_stage1 which combine initial_from_to_bus and F_sign with sequence = under attack
  F_sign_seq_under_attack_mapping = vector(length = dim(initial_from_to_bus)[1])
  for (ii in 1:dim(initial_from_to_bus)[1]) {
    tmp = which(initial_from_to_bus[,3] == names(F_sign_seq_under_attack)[ii])
    F_sign_seq_under_attack_mapping[tmp] = as.numeric(F_sign_seq_under_attack[ii])
  }
  
  comb_seq0_seq_under_attack_stage1 = cbind(initial_from_to_bus, F_sign_seq_under_attack_mapping)
  
  # Firstly, find the edges with F_ = 0 when sequence = 17(under attack) - stage 1#
  edges_with_zero_flow_label = comb_seq0_seq_under_attack_stage1[which(comb_seq0_seq_under_attack_stage1$F_sign_seq_under_attack_mapping == 0),"Branch_f"] 
  
  # Secondly, remove the rows with F_ = 0 when sequence = 17 - generate comb_seq0_seq_under_attack_stage2 #
  comb_seq0_seq_under_attack_stage2 = comb_seq0_seq_under_attack_stage1[-which(comb_seq0_seq_under_attack_stage1$F_sign_seq_under_attack_mapping == 0),]
  rownames(comb_seq0_seq_under_attack_stage2) = c(1:dim(comb_seq0_seq_under_attack_stage2)[1])
  
  # Thirdly, calcualte the difference between seq17 and initial which try to find the which edges change the direction - generage comb_seq0_seq_under_attack_stage3 #
  comb_seq0_seq_under_attack_stage3 = comb_seq0_seq_under_attack_stage2
  comb_seq0_seq_under_attack_stage3$diff = comb_seq0_seq_under_attack_stage3$F_sign_seq_under_attack_mapping - comb_seq0_seq_under_attack_stage3$F_Sig_zero
  
  # Fourthly, find out the edges with direction changed - generate comb_seq0_seq_under_attack_stage4 #
  #-------------fixed already-----------------------------------#
  # Debug for comb_seq0_seq_under_attack_stage4 - 06/21 08:33am #
  # if there is not edge direction changed i.e., all diff =0 #
  # add condition - whether all difference equal to 0 i.e., not edge direction changed#
  if(!all(comb_seq0_seq_under_attack_stage3$diff == 0)){
    
    comb_seq0_seq_under_attack_stage4 = comb_seq0_seq_under_attack_stage3[-which(comb_seq0_seq_under_attack_stage3$diff==0),]
    rownames(comb_seq0_seq_under_attack_stage4) = c(1:dim(comb_seq0_seq_under_attack_stage4)[1])
    
    # Fifthly, build the updated (new when seq = 17) from_to_bus dataframe after sequence = 17 attack - generate comb_seq0_seq_under_attack_stage5 #
    comb_seq0_seq_under_attack_stage5 = data.frame(From_Bus = rep(0,dim(comb_seq0_seq_under_attack_stage4)[1]), 
                                                   To_Bus = rep(0,dim(comb_seq0_seq_under_attack_stage4)[1]), Branch_F = rep(0,dim(comb_seq0_seq_under_attack_stage4)[1]), Branch_f = rep(0,dim(comb_seq0_seq_under_attack_stage4)[1]))
    
    comb_seq0_seq_under_attack_stage5$From_Bus = comb_seq0_seq_under_attack_stage4$To_Bus_Initial
    comb_seq0_seq_under_attack_stage5$To_Bus = comb_seq0_seq_under_attack_stage4$From_Bus_initial
    comb_seq0_seq_under_attack_stage5$Branch_F = comb_seq0_seq_under_attack_stage4$Branch
    comb_seq0_seq_under_attack_stage5$Branch_f = comb_seq0_seq_under_attack_stage4$Branch_f
    
    # change the label and name for seq_zero_case118_ieee_network #
    tmp_label_store = V(seq_zero_case118_ieee_network)$label
    V(seq_zero_case118_ieee_network)$label = V(seq_zero_case118_ieee_network)$name
    V(seq_zero_case118_ieee_network)$name = tmp_label_store
    # change over #
    
    # 10, create the edgelist with "edge name" for sequence = 0 from seq_zero_case118_ieee_network #
    initial_graph_edgelist_with_Branch_f = cbind(get.edgelist(seq_zero_case118_ieee_network), edge_attr(seq_zero_case118_ieee_network)$name) # - shape is  182*3
    
    # 11, find the final short version of initial from_to bus info and Branch_f WITHOUT F_Sig_zero = 0 scenario #
    label_without_f_sig_zero = which(comb_seq0_seq_under_attack_stage4$F_Sig_zero != 0 )
    df_with_edge_changed_with_F_signotzero = comb_seq0_seq_under_attack_stage4[label_without_f_sig_zero, c("From_Bus_initial", "To_Bus_Initial","Branch_f")]
    rownames(df_with_edge_changed_with_F_signotzero) = c(1:dim(df_with_edge_changed_with_F_signotzero)[1]) #
    
    # 12, remove direction changed edges from initial_graph_edgelist_with_Branch_f - generate edgelist_change_stage1 #
    edgelist_change_stage1_label = (!duplicated(rbind(initial_graph_edgelist_with_Branch_f, as.matrix(df_with_edge_changed_with_F_signotzero)), fromLast = T))[1:dim(initial_graph_edgelist_with_Branch_f)[1]]
    edgelist_change_stage1 = initial_graph_edgelist_with_Branch_f[edgelist_change_stage1_label,]
    
    # 13, combine delted verion above i.e., stage 1, with comb_seq0_seq_under_attack_stage5 # 
    final_bus_fromto_Branch_f = rbind(edgelist_change_stage1, as.matrix(comb_seq0_seq_under_attack_stage5[,c("From_Bus", "To_Bus","Branch_f")]))
    
    # 14, remove the edges with F_  = 0 when sequence = 17 #
    final_verion_edgelist_with_Branch_f = final_bus_fromto_Branch_f[(!final_bus_fromto_Branch_f[,3] %in% edges_with_zero_flow_label),]
    
    # 15, create igraph graphs from data frames #
    final_verion_graph = graph_from_data_frame(as.data.frame(final_verion_edgelist_with_Branch_f[,c(1:3)]), directed=TRUE, vertices=V(seq_zero_case118_ieee_network)$name)
  }else if(all(comb_seq0_seq_under_attack_stage3$diff == 0)){
    seq_zero_case118_ieee_network = graph_from_edgelist(as.matrix(initial_from_to_bus[,c(1:2)]))
    seq_zero_case118_ieee_network = seq_zero_case118_ieee_network %>% set_edge_attr("name", value = initial_from_to_bus[,3])
    
    #------#
    # repeat 6, assign name, label, color, feature in case_118_network to seq_zero_case118_ieee_network (another version)
    four_attrs_case_118_network = data.frame(name = vertex_attr(case_118_network)$name,
                                             label = vertex_attr(case_118_network)$label,
                                             color = vertex_attr(case_118_network)$color,
                                             feature = vertex_attr(case_118_network)$feature)
    
    
    four_attrs_seq_zero_case118_ieee_network = data.frame(name = rep(0,length(vertex_attr(seq_zero_case118_ieee_network)$name)),
                                                          label = vertex_attr(seq_zero_case118_ieee_network)$name,
                                                          color = rep(0,length(vertex_attr(seq_zero_case118_ieee_network)$name)),
                                                          feature = rep(0,length(vertex_attr(seq_zero_case118_ieee_network)$name)))
    
    for (kk in 1:dim(four_attrs_seq_zero_case118_ieee_network)[1]) {
      tmp_label = which(vertex_attr(case_118_network)$label == vertex_attr(seq_zero_case118_ieee_network)$name[kk])
      four_attrs_seq_zero_case118_ieee_network[kk,1] = vertex_attr(case_118_network)$name[tmp_label]
      four_attrs_seq_zero_case118_ieee_network[kk,3] = vertex_attr(case_118_network)$color[tmp_label]
      four_attrs_seq_zero_case118_ieee_network[kk,4] = vertex_attr(case_118_network)$feature[tmp_label]
    }
    
    four_attrs_seq_zero_case118_ieee_network$name = four_attrs_seq_zero_case118_ieee_network$name #paste("bus",four_attrs_seq_zero_case118_ieee_network$name,sep = "")
    four_attrs_seq_zero_case118_ieee_network_dataframe = four_attrs_seq_zero_case118_ieee_network
    four_attrs_seq_zero_case118_ieee_network_mat = as.matrix(four_attrs_seq_zero_case118_ieee_network_dataframe)
    all.equal(four_attrs_seq_zero_case118_ieee_network_mat[,2],vertex_attr(seq_zero_case118_ieee_network)$name) #TRUE
    
    V(seq_zero_case118_ieee_network)$feature = as.numeric(four_attrs_seq_zero_case118_ieee_network_mat[,4])
    V(seq_zero_case118_ieee_network)$name = four_attrs_seq_zero_case118_ieee_network_mat[,1]
    V(seq_zero_case118_ieee_network)$label = four_attrs_seq_zero_case118_ieee_network_mat[,2]
    V(seq_zero_case118_ieee_network)$color = four_attrs_seq_zero_case118_ieee_network_mat[,3]
    # IGRAPH dc75189 DN-- 118 186 -- #
    #------#
    
    # rename edges' names for seq_zero_case118_ieee_network use f_ not F_#
    length(edge_attr(seq_zero_case118_ieee_network)$name)
    edge_label_num = as.numeric(gsub("F_", "", edge_attr(seq_zero_case118_ieee_network)$name)) # extract the number in edge_attr(seq_zero_case118_iee_network)$name; since i want to rebuild the f_ for edges' names
    E(seq_zero_case118_ieee_network)$name = paste("f_",edge_label_num, sep="")
    # rename complete #
    
    # change the name and label for seq_zero_case118_ieee_network #
    tmp_label_store = V(seq_zero_case118_ieee_network)$label
    V(seq_zero_case118_ieee_network)$label = V(seq_zero_case118_ieee_network)$name
    V(seq_zero_case118_ieee_network)$name = tmp_label_store
    # change complete #
    
    initial_graph_edgelist_with_Branch_f_case2 = cbind(get.edgelist(seq_zero_case118_ieee_network), edge_attr(seq_zero_case118_ieee_network)$name)
    final_verion_edgelist_with_Branch_f = initial_graph_edgelist_with_Branch_f_case2[(!initial_graph_edgelist_with_Branch_f_case2[,3] %in% edges_with_zero_flow_label),]
    
    # fixed (already) for colnames can not work for dataframe with only one row - 06/23 #
    if(!is.null(dim(final_verion_edgelist_with_Branch_f))){
      colnames(final_verion_edgelist_with_Branch_f) = c("From_Bus", "To_Bus","Branch_f")}
    else{
      final_verion_edgelist_with_Branch_f = data.frame(t(unlist(final_verion_edgelist_with_Branch_f)))
      colnames(final_verion_edgelist_with_Branch_f) = c("From_Bus", "To_Bus","Branch_f")
    }
    final_verion_graph = graph_from_data_frame(as.data.frame(final_verion_edgelist_with_Branch_f[,c(1:3)]), directed=TRUE, vertices=V(seq_zero_case118_ieee_network)$name)
  }
  # 16, delete the nodes with b_ = false #
  final_version_graph_after_delete_false_nodes = delete_vertices(final_verion_graph, removed_nodes_label)
  
  # 17, assign various attributes from seq_zero_case118_ieee_network to final_version_graph_after_delete_false_nodes
  vertex_attr(final_version_graph_after_delete_false_nodes)$busname = rep(NA, length(vertex_attr(final_version_graph_after_delete_false_nodes)$name))
  vertex_attr(final_version_graph_after_delete_false_nodes)$feature = rep(NA, length(vertex_attr(final_version_graph_after_delete_false_nodes)$name))
  
  for (uu in 1:length(vertex_attr(final_version_graph_after_delete_false_nodes)$name)) {
    tmp_label_ii = vertex_attr(seq_zero_case118_ieee_network)$name %in% vertex_attr(final_version_graph_after_delete_false_nodes)$name[uu]
    vertex_attr(final_version_graph_after_delete_false_nodes)$busname[uu] = vertex_attr(seq_zero_case118_ieee_network)$label[tmp_label_ii]
    vertex_attr(final_version_graph_after_delete_false_nodes)$feature[uu] = vertex_attr(seq_zero_case118_ieee_network)$feature[tmp_label_ii]
  }
  
  # assign Branch_F to final_version_graph_after_delete_false_nodes #
  edge_num_with_Branch_f = as.numeric(gsub("f_", "", edge_attr(final_version_graph_after_delete_false_nodes)$Branch_f))
  edge_attr(final_version_graph_after_delete_false_nodes)$Branch_F = paste("F_",edge_num_with_Branch_f, sep="")
  
  # assign F_ to edge weight #
  edge_attr(final_version_graph_after_delete_false_nodes)$weight = rep(NA, length(edge_attr(final_version_graph_after_delete_false_nodes)$Branch_f))
  for (mm in 1:length(edge_attr(final_version_graph_after_delete_false_nodes)$Branch_f)) {
    sequnce_under_attack_Flow_info = new_case118_result1[sequence_row,c(which(colnames(new_case118_result1) == "F_1") : which(colnames(new_case118_result1) == "F_186"))]
    tmp_weight_label = names(sequnce_under_attack_Flow_info) %in% edge_attr(final_version_graph_after_delete_false_nodes)$Branch_F[mm]
    edge_attr(final_version_graph_after_delete_false_nodes)$weight[mm] = abs(as.numeric(sequnce_under_attack_Flow_info[tmp_weight_label]))
  }
  
  return(final_version_graph_after_delete_false_nodes) # the output graph with parallel edges
}



############################################################################################################################

#=== MOTIFS ===#

t = 2
fnm = 81

n1 = n2 = numeric(fnm)
m2 = matrix(NA, nrow = fnm, ncol = 16)
T_1 = T_2 = V_1 = V_2 = V_3 = V_4 = V_5 = V_6 = Tot_T = Tot_V = c()
C_T1 = C_T2 = C_V1 = C_V2 = C_V3 = C_V4 = C_V5 = C_V6 = c()

for(i in 1:fnm){
  
  
  graph_structure_under_sequence_row_attack = simplify(output_graph_under_attack_f(case_118_network = input_graph, initial_status_row = 2, sequence_row = t))
  graph_structure_under_sequence_row_attack_undirected = as.undirected(graph_structure_under_sequence_row_attack)
  gg_dir = graph_structure_under_sequence_row_attack_undirected
  
  
  ########################
  #=== Motif analysis ===#
  ########################
  
  
  
  #------ Motif size = 3 -------#
  m1 = motifs(gg_dir, 3)
  m1[is.na(m1)]  =  0
  
  n01 = count_motifs(gg_dir, 3)
  n1[i] =  n01
  
  T_1[i] = m1[3]  
  T_2[i] = m1[4] 
 
  
  Tot_T[i] = sum(m1[3],m1[4])
  
  
  
  #------ Motif size = 4 -------#
  m2 = motifs(gg_dir, 4)
  m2[is.na(m2)]  =  0
  
  n02 = count_motifs(gg_dir, 4)
  n2[i] =  n02
  
  V_1[i] = m2[5]  
  V_2[i] = m2[7] 
  V_3[i] = m2[8] 
  V_4[i] = m2[9]   
  V_5[i] = m2[10] 
  V_6[i] = m2[11]  
  
  Tot_V[i] = sum(m2[5],m2[7],m2[8],m2[9],m2[10],m2[11])
  
  t = t + 1
  
}
# 
# n12 = n1[1]
# 
# n22 = n2[1]

# Motif concentration #

C_T1 = T_1/n1; C_T1[is.na(C_T1)] = 0
C_T2 = T_2/n1; C_T2[is.na(C_T2)] = 0

C_V1 = V_1/n2; C_V1[is.na(C_V1)] = 0
C_V2 = V_2/n2; C_V2[is.na(C_V2)] = 0
C_V3 = V_3/n2; C_V3[is.na(C_V3)] = 0 
C_V4 = V_4/n2; C_V4[is.na(C_V4)] = 0
C_V5 = V_5/n2; C_V5[is.na(C_V5)] = 0
C_V6 = V_6/n2; C_V6[is.na(C_V6)] = 0


resultsN = data.frame(c(0: (fnm - 1)), Tot_T, Tot_V, T_1, T_2, C_T1, C_T2, V_1, V_2, V_3, V_4, V_5, V_6, C_V1, C_V2, C_V3, C_V4, C_V5, C_V6 )
write.csv(resultsN, "Motifs_UNDIR_118-Bus.csv")


#=== TDA ===#



t = 2
cap = 3 # dimension cap
delta = 0.10
filt_len = 10


wasser_dist_01 = c()

for(i in 1:fnm){
  
  
  graph_structure_under_sequence_row_attack = simplify(output_graph_under_attack_f(case_118_network = input_graph, initial_status_row = 2, sequence_row = t))
  graph_structure_under_sequence_row_attack_undirected = as.undirected(graph_structure_under_sequence_row_attack)
  gg_dir = graph_structure_under_sequence_row_attack_undirected
  E(gg_dir)$weight  =  (E(gg_dir)$weight)/max(E(gg_dir)$weight)
  
  
  ######################
  #=== TDA analysis ===#
  ######################
  
  
  
  A1 =  get.adjacency(gg_dir, attr="weight")
  A2 = as.matrix(A1)
  
  A2[A2==0] = 999
  diag(A2)=0
  
  
  d = length(V(gg_dir))
  print(i)
  
  # writing data into file M.txt
  cat(d,file='M.txt',append=F,sep = '\n')
  cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') 
  cat(A2,file='M.txt',append = T) 
  
  system('perseusWin.exe distmat M.txt Moutput')
  
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
  
  t = t + 1
  
}

resultsN = data.frame(c(0: (fnm - 1)), wasser_dist_01)
write.csv(resultsN, "TDA_UNDIR_118.csv")

#=== Giant component ===#


t = 2

GC = c()


for(i in 1:fnm){
  
  
  graph_structure_under_sequence_row_attack = simplify(output_graph_under_attack_f(case_118_network = input_graph, initial_status_row = 2, sequence_row = t))
  graph_structure_under_sequence_row_attack_undirected = as.undirected(graph_structure_under_sequence_row_attack)
  gg_dir = graph_structure_under_sequence_row_attack_undirected
  
  
  
  
  #########################
  #=== Giant component ===#
  #########################
  
  
  #------ Small world properties  -------#
  compts = components(gg_dir)    #largest Connected components/Giant component 
  GC[i] = max(compts$csize) 
  
  
  t = t +1
  
  
}


resultsN = data.frame(c(0: (fnm - 1)), GC)
write.csv(resultsN, "GCC_UNDIR_118.csv")


#=== Connectivity Loss ===#

t = 2

fin = c()

for(i in 1:fnm){
  
  
  graph_structure_under_sequence_row_attack = simplify(output_graph_under_attack_f(case_118_network = input_graph, initial_status_row = 2, sequence_row = t))
  graph_structure_under_sequence_row_attack_undirected = as.undirected(graph_structure_under_sequence_row_attack)
  gg_dir = graph_structure_under_sequence_row_attack_undirected
  
  
  
  
  ###########################
  #=== Connectivity Loss ===#
  ###########################
  
  
  dist2 <- distances(gg_dir)
  dist2[dist2 == Inf] <- 0
  dist2[dist2 > 0] <- 1
  tot2 <- sum(dist2)
  if(i == 1){tot = tot2}
  fin[i] <- tot - tot2
  
  
  t = t +1
  
}


fin = fin/tot


resultsN = data.frame(c(0: (fnm - 1)), fin)
write.csv(resultsN, "LC_UNDIR_118.csv")

#=== Small-worldness property ===#

t = 2

CLUS = DIAM = APL = c()


for(i in 1:fnm){
  
  
  graph_structure_under_sequence_row_attack = simplify(output_graph_under_attack_f(case_118_network = input_graph, initial_status_row = 2, sequence_row = t))
  graph_structure_under_sequence_row_attack_undirected = as.undirected(graph_structure_under_sequence_row_attack)
  gg_dir = graph_structure_under_sequence_row_attack_undirected
  
  
  
  
  ##################################
  #=== Small-worldness property ===#
  ##################################
  
  CLUS[i] = transitivity(gg_dir) # Clustering Coefficient
  DIAM[i] = diameter(gg_dir) # Diameter
  APL[i] = average.path.length(gg_dir) # Average path length
  
  
  
  t = t +1
  
}


resultsN = data.frame(c(0: (fnm - 1)), APL, CLUS, DIAM)
write.csv(resultsN, "SWProp_UNDIR_118.csv")

############################################################################################################################

#=== Load served (BLACKOUT SIZE/ Load shed = 1 - Load served) ===#

Load_matrix = new_case118_result1[which(colnames(new_case118_result1)=="L_1"):which(colnames(new_case118_result1)=="L_118")] 
LoadShs_mat = apply(Load_matrix, 1, function(x){sum(x)})
LS_value    = LoadShs_mat/LoadShs_mat[1]
write.csv(LS_value, "LS_Value.csv")




end.time = Sys.time()
time.taken = end.time - start.time

