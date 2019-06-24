###############################################################################################################################################################################

########################################################################
#===== Ensure that Perseus is within same folder as the main code =====#
########################################################################


Resilience_3Motifs = function(attack = c("node", "edge"),  directed = c("TRUE", "FALSE"), graph_type = c("weighted", "unweighted") ,  attack_type = c("degree", "V_betweeness", "strength" , "weight_hier", "E_betweeness"), net_val, frac_val){
  
  frac_len = length(frac_val)
  network = net_val
  
  network_org = network
  
  deg = degree(network_org)
  d_ord = length(deg)
  deg_org = order(deg, decreasing = T)
  
  bet = betweenness(network_org)
  b_ord = length(bet)
  bet_org = order(bet, decreasing = T)
  
  e_bet = edge_betweenness(network_org)
  eb_ord = length(e_bet)
  ebet_org = order(e_bet, decreasing = T)
  
  if(graph_type == "weighted"){
    
    str_vals = strength(network_org)
    s_ord = length(str_vals)
    str_ord = order(str_vals, decreasing = T)
    
    ed_val = E(network_org)$weight
    ed_ord = length(ed_val)
    edge_org = order(ed_val, decreasing = T)
  }else{
    str_vals = ed_val = c()
  }
  
  
  n1 = n2 = numeric(frac_len)
  V_1 = V_2 = V_3 = V_4 = V_5 = V_6 = Tot_M = Tot_V = numeric(frac_len)
  C_V1 = C_V2 = C_V3 = C_V4 = C_V5 = C_V6 = numeric(frac_len)
 
  
  ########################
  ##### NODE ATTACKS #####
  ########################

  if(attack == "node"){
    
   
    if(directed == "TRUE"){
      
      tt = matrix(NA, nrow = (frac_len), ncol = 16)
     
      #------ ORIGINAL GRAPH MOTIFS -------#
      m_OR = triad.census(network_org)
      Tot_OR_M =  sum(m_OR)
      
      
      #------ AFTER ATTACK CALCULATIONS -------#
      
      for (i in 1:(frac_len)){ 
        
        
        if (attack_type == "degree") #===== Degree - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            nodes_to_delete = V(network_org)[deg_org[1:round(d_ord*frac_val[i])]]
            network = delete_vertices(network_org, nodes_to_delete)
            
          }
          
        }
        
        
        else if (attack_type == "V_betweeness") #===== Betweeness centrality - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            nodes_to_delete = V(network_org)[bet_org[1:round(b_ord*frac_val[i])]]
            network = delete_vertices(network_org, nodes_to_delete)    
          }
          
        }
        
        else if (attack_type == "strength") #===== Strength - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            nodes_to_delete = V(network_org)[str_ord[1:round(s_ord*frac_val[i])]]
            network = delete_vertices(network_org, nodes_to_delete)    
          }
          
        }
        
        
        #------ Motif size = 3 -------#
        m2 = triad.census(network)
        Tot_M[i] =  sum(m2)
        
        tt[i,] = m2
        
        
      }
      
      
      tt_full = rbind(m_OR, tt)
      C_Motifs = apply(tt_full, 2, function(x){x/Tot_OR_M})
      
      resultsN = data.frame(c(0, frac_val), c(Tot_OR_M,Tot_M), tt_full, C_Motifs)
      colnames(resultsN) = c("frac", "Tot_T", paste0("T_", 1:16), paste0("C_T", 1:16))
     
    }
    
    else if(directed == "FALSE"){
      
      
      
      for (i in 1:(frac_len + 1)){ 
        
        
      m2 = motifs(network, 3)
      m2[is.na(m2)]  =  0
      
      n02 = count_motifs(network, 3)
      n2[i] =  n02
      
      V_1[i] = m2[3]  
      V_2[i] = m2[4] 
     
      
      Tot_V[i] = sum(m2)
      
      
      
      if (attack_type == "degree") #===== Degree - Based Attacks =====# 
      {
        if (i <= frac_len)
        {
          nodes_to_delete = V(network_org)[deg_org[1:round(d_ord*frac_val[i])]]
          network = delete_vertices(network_org, nodes_to_delete)
          
        }
        
      }
      
      
      else if (attack_type == "V_betweeness") #===== Betweeness centrality - Based Attacks =====# 
      {
        if (i <= frac_len)
        {
          nodes_to_delete = V(network_org)[bet_org[1:round(b_ord*frac_val[i])]]
          network = delete_vertices(network_org, nodes_to_delete)    
        }
        
      }
      
      else if (attack_type == "strength") #===== Strength - Based Attacks =====# 
      {
        if (i <= frac_len)
        {
          nodes_to_delete = V(network_org)[str_ord[1:round(s_ord*frac_val[i])]]
          network = delete_vertices(network_org, nodes_to_delete)    
        }
        
      }
      
        
      }
      
    
    n22 = n2[1]
    
    C_V1 = V_1/n22
    C_V2 = V_2/n22
    
    
      resultsN= data.frame( c(0,frac_val),Tot_V,V_1,V_2,C_V1,C_V2)
      colnames(resultsN) = c("frac", "Tot_T", paste0("T_", 1:2), paste0("C_T", 1:2))
      
    }
    
    

  }
  
  
  
  ########################
  ##### EDGE ATTACKS #####
  ########################
  
  if(attack == "edge"){
    
    
    if(directed == "TRUE"){
      
      tt = matrix(NA, nrow = (frac_len), ncol = 16)
      
      #------ ORIGINAL GRAPH MOTIFS -------#
      m_OR = triad.census(network_org)
      Tot_OR_M =  sum(m_OR)
      
      
      #------ AFTER ATTACK CALCULATIONS -------#
      
      for (i in 1:(frac_len)){ 
        
        
        if (attack_type == "weight_hier") #===== Plain weight - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            edges_to_delete = E(network_org)[edge_org[1:round(ed_ord*frac_val[i])]]
            network = delete_edges(network_org, edges_to_delete)
            
          }
          
        }
        
        else if (attack_type == "E_betweeness") #===== Betweeness centrality - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            edges_to_delete = E(network_org)[ebet_org[1:round(eb_ord*frac_val[i])]]
            network = delete_edges(network_org, edges_to_delete)    
          }
          
        }
        
        
        
        #------ Motif size = 3 -------#
        m2 = triad.census(network)
        Tot_M[i] =  sum(m2)
        
        tt[i,] = m2
        
        
      }
      
      
      tt_full = tt
      C_Motifs = apply(tt_full, 2, function(x){x/Tot_OR_M})
      
      resultsN = data.frame(frac_val, Tot_M, tt_full, C_Motifs)
      colnames(resultsN) = c("frac", "Tot_T", paste0("T_", 1:16), paste0("C_T", 1:16))
      
    }
    
    else if(directed == "FALSE"){
      
      
      
      for (i in 1:(frac_len)){ 
        
        
        m2 = motifs(network, 3)
        m2[is.na(m2)]  =  0
        
        n02 = count_motifs(network, 3)
        n2[i] =  n02
        
        V_1[i] = m2[3]  
        V_2[i] = m2[4] 
        
        
        Tot_V[i] = sum(m2)
        
        
        if (attack_type == "weight_hier") #===== Plain Weight - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            edges_to_delete = E(network_org)[edge_org[1:round(ed_ord*frac_val[i])]]
            network = delete_edges(network_org, edges_to_delete)
            
          }
          
        }
        
        else if (attack_type == "E_betweeness") #===== Betweeness centrality - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            edges_to_delete = E(network_org)[ebet_org[1:round(eb_ord*frac_val[i])]]
            network = delete_edges(network_org, edges_to_delete)    
          }
          
        }
        
 
        
      }
      
      
      n22 = n2[1]
      
      C_V1 = V_1/n22
      C_V2 = V_2/n22
      
      
      resultsN= data.frame( frac_val,Tot_V,V_1,V_2,C_V1,C_V2)
      colnames(resultsN) = c("frac", "Tot_T", paste0("T_", 1:2), paste0("C_T", 1:2))
      
    }
    
    
    
  }
  
  return(resultsN)
  
}


Resilience_4Motifs = function(attack = c("node", "edge"),  directed = c("TRUE", "FALSE"), graph_type = c("weighted", "unweighted") ,  attack_type = c("degree", "V_betweeness", "strength" , "weight_hier", "E_betweeness"), net_val, frac_val){
  
  frac_len = length(frac_val)
  network = net_val
  
  network_org = network
  
  deg = degree(network_org)
  d_ord = length(deg)
  deg_org = order(deg, decreasing = T)
  
  bet = betweenness(network_org)
  b_ord = length(bet)
  bet_org = order(bet, decreasing = T)
  
  e_bet = edge_betweenness(network_org)
  eb_ord = length(e_bet)
  ebet_org = order(e_bet, decreasing = T)
  
  if(graph_type == "weighted"){
    
    str_vals = strength(network_org)
    s_ord = length(str_vals)
    str_ord = order(str_vals, decreasing = T)
    
    ed_val = E(network_org)$weight
    ed_ord = length(ed_val)
    edge_org = order(ed_val, decreasing = T)
  }else{
    str_vals = ed_val = c()
  }
  
  
  n1 = n2 = numeric(frac_len)
  V_1 = V_2 = V_3 = V_4 = V_5 = V_6 = Tot_M = Tot_V = numeric(frac_len)
  C_V1 = C_V2 = C_V3 = C_V4 = C_V5 = C_V6 = numeric(frac_len)
  
  
  ########################
  ##### NODE ATTACKS #####
  ########################
  
  if(attack == "node"){
    
    
    if(directed == "TRUE"){
      print("Code search in progress. Thank you for your patience!")
    }
    
    else if(directed == "FALSE"){
      
      for (i in 1:(frac_len)){ 
        
        
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
        
        
        
        if (attack_type == "degree") #===== Degree - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            nodes_to_delete = V(network_org)[deg_org[1:round(d_ord*frac_val[i])]]
            network = delete_vertices(network_org, nodes_to_delete)
            
          }
          
        }
        
        
        else if (attack_type == "V_betweeness") #===== Betweeness centrality - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            nodes_to_delete = V(network_org)[bet_org[1:round(b_ord*frac_val[i])]]
            network = delete_vertices(network_org, nodes_to_delete)    
          }
          
        }
        
        else if (attack_type == "strength") #===== Strength - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            nodes_to_delete = V(network_org)[str_ord[1:round(s_ord*frac_val[i])]]
            network = delete_vertices(network_org, nodes_to_delete)    
          }
          
        }
        
        
      }
      
      
      n22 = n2[1]
      
      C_V1 = V_1/n22 
      C_V2 = V_2/n22 
      C_V3 = V_3/n22  
      C_V4 = V_4/n22 
      C_V5 = V_5/n22
      C_V6 = V_6/n22
      
      
      resultsN = data.frame(frac_val,Tot_V,V_1,V_2,V_3,V_4,V_5,V_6,C_V1,C_V2,C_V3,C_V4,C_V5,C_V6)
      colnames(resultsN) = c("frac", "Tot_M", paste0("M_", 1:6), paste0("C_M", 1:6))
      
    }
    
    
    
  }
  
  
  
  ########################
  ##### EDGE ATTACKS #####
  ########################
  
  if(attack == "edge"){
    
    
    if(directed == "TRUE"){
      print("Code search in progress. Thank you for your patience!")
    }
    
    else if(directed == "FALSE"){
      
      for (i in 1:(frac_len)){ 
        
        
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
        
        
        
        if (attack_type == "weight_hier") #===== Plain Weight - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            edges_to_delete = E(network_org)[edge_org[1:round(ed_ord*frac_val[i])]]
            network = delete_edges(network_org, edges_to_delete)
            
          }
          
        }
        
        else if (attack_type == "E_betweeness") #===== Betweeness centrality - Based Attacks =====# 
        {
          if (i <= frac_len)
          {
            edges_to_delete = E(network_org)[ebet_org[1:round(eb_ord*frac_val[i])]]
            network = delete_edges(network_org, edges_to_delete)    
          }
          
        }
        
        
        
      }
      
      n22 = n2[1]
      
      C_V1 = V_1/n22 
      C_V2 = V_2/n22 
      C_V3 = V_3/n22  
      C_V4 = V_4/n22 
      C_V5 = V_5/n22
      C_V6 = V_6/n22
      
      
      resultsN = data.frame(frac_val,Tot_V,V_1,V_2,V_3,V_4,V_5,V_6,C_V1,C_V2,C_V3,C_V4,C_V5,C_V6)
      colnames(resultsN) = c("frac", "Tot_M", paste0("M_", 1:6), paste0("C_M", 1:6))
      
    }
    
    
    
  }
  
  return(resultsN)
  
}


###############################################################################################################################################################################
