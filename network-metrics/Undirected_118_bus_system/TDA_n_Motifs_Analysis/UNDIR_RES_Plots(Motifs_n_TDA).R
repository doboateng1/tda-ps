
rm(list = ls())

setwd("C:/Users/doforib/Desktop/Attack trials/Undirected_graph_analysis")


##############
#=== DATA ===#
##############

res1 = read.csv("TDA_Max_flow_weight_function_80.csv")
resX = res1[,2]
res1 = (as.data.frame(apply(res1, 2, function(x){1 - (x/max(x))})))[,-1]
res1$V_Remv = resX

res2 = read.csv("TDA_Changing_flow_weight_function_80.csv")
res2 = (as.data.frame(apply(res2, 2, function(x){1 - (x/max(x))})))[,-1]
res2$V_Remv = resX

res3 = read.csv("Motifs_UNDIR_3N_80Reliab.csv")
res3 = subset(res3, select = c(V_Remv,DRt))


res4 = read.csv("Motifs_UNDIR_4N_80Reliab.csv")
res4 = subset(res4, select = c(V_Remv,DRt))


res5 = read.csv("Motifs_UNDIR_3N_80.csv")
res6 = read.csv("Motifs_UNDIR_4N_80.csv")


########################################################################################################################

setwd("C:/Users/doforib/Desktop/Attack trials/Undirected_graph_analysis/plots")



########################################################################################################################################################################################
#######################################################################################################################################################

#######################################################################################################################################################
# CONSTANT WEIGHT #
#######################################################################################################################################################

#=== All + TDA together ===#

par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(resX, res3$DRt, type = "l", xlab = "Nodes removed",lty=6, lwd=4, main = "Constant weight (Max_flow) TDA and Motif analysis of nesta_case_118_study-02",
     ylab =  "Robustness metric",col = "purple", ylim = c(0,1), xlim = c(0,80))                                 #Wass01#
#lines(resX, data41,         type = "l",          lty = 1,lwd=2, col = "red")
lines(resX, res4$DRt,      type = "l",          lty = 2, lwd=2, col = "cyan")
lines(resX, res1$Betti0,   type = "l",          lty = 1, lwd=3, col = "green")
lines(resX, res1$Betti1,   type = "l", pch = 2, lty = 2, lwd=2, col = "orange")
lines(resX, res1$Wass01,   type = "o", pch = 21,lty = 3, lwd=2, col = "red")

grid();
#par(xpd = T)
legend("bottomleft", bty = "n", inset = 0.005,
       c(expression(paste(Delta, MR[3])), expression(paste(Delta, MR[4])) , expression(paste(1,"-",Delta ,beta[0])), expression(paste(1,"-",Delta ,beta[1])),expression(paste(1," - ",Delta ,W[2](D[0],D[P])))),
       cex = 0.8,lty=c(6,2,1,2,3),pch = c(NA, NA,NA,2,21),lwd=c(4,2,3,2,2), col=c("purple", "cyan", "green", "orange","red")) 

par(mar = c(5, 4, 4, 2) + 0.1)



#######################################################################################################################################################
# DYNAMIC WEIGHT #
#######################################################################################################################################################

#=== All + TDA together ===#

par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(resX, res3$DRt, type = "l", xlab = "Nodes removed",lty=6, lwd=4, main = "Changing weight (Flow) TDA and Motif analysis of nesta_case_118_study-02",
     ylab =  "Robustness metric",col = "purple", ylim = c(0,1), xlim = c(0,80))                                 #Wass01#
#lines(resX, data41,         type = "l",          lty = 1,lwd=2, col = "red")
lines(resX, res4$DRt,      type = "l",          lty = 2, lwd=2, col = "cyan")
lines(resX, res2$Betti0,   type = "l",          lty = 1, lwd=3, col = "green")
lines(resX, res2$Betti1,   type = "l", pch = 2, lty = 2, lwd=2, col = "orange")
lines(resX, res2$Wass01,   type = "o", pch = 21,lty = 3, lwd=2, col = "red")

grid();
#par(xpd = T)
legend("bottomleft", bty = "n", inset = 0.005,
       c(expression(paste(Delta, MR[3])), expression(paste(Delta, MR[4])) , expression(paste(1,"-",Delta ,beta[0])), expression(paste(1,"-",Delta ,beta[1])),expression(paste(1," - ",Delta ,W[2](D[0],D[P])))),
       cex = 0.8,lty=c(6,2,1,2,3),pch = c(NA, NA,NA,2,21),lwd=c(4,2,3,2,2), col=c("purple", "cyan", "green", "orange","red")) 

par(mar = c(5, 4, 4, 2) + 0.1)





#######################################################################################################################################################
# ALL 3-N Motif CONCENTRATION #
#######################################################################################################################################################

#=== All + TDA together ===#

par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(resX, res5$C_V1, type = "l", xlab = "Nodes removed",lty=6, lwd=4, main = "3-Node motif concentration",
     ylab =  "Concentration",col = "purple", ylim = c(0,1), xlim = c(0,80))                                 #Wass01#
#lines(resX, data41,         type = "l",          lty = 1,lwd=2, col = "red")
lines(resX, res5$C_V2,      type = "l",          lty = 2, lwd=2, col = "red")

grid();
#par(xpd = T)
legend("topright", bty = "n", inset = 0.005,
       c(expression(paste(T[1])), expression(paste(T[2]))),
       cex = 0.8,lty=c(6,2),pch = c(NA, NA),lwd=c(4,2), col=c("purple", "red")) 

par(mar = c(5, 4, 4, 2) + 0.1)



#######################################################################################################################################################
# ALL 4-N Motif CONCENTRATION #
#######################################################################################################################################################

#=== All + TDA together ===#

par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(resX, res6$C_V1, type = "l", xlab = "Nodes removed",lty=6, lwd=4, main = "4-Node motif concentration",
     ylab =  "Concentration",col = "purple", ylim = c(0,0.7), xlim = c(0,80))                                 #Wass01#
lines(resX, res6$C_V2,   type = "l",          lty = 1,lwd=2, col = "red")
lines(resX, res6$C_V3,   type = "l",          lty = 2, lwd=2, col = "cyan")
lines(resX, res6$C_V4,   type = "l",          lty = 1, lwd=3, col = "green")
lines(resX, res6$C_V5,   type = "l", pch = 2, lty = 2, lwd=2, col = "orange")
lines(resX, res6$C_V6,   type = "o", pch = 21,lty = 3, lwd=2, col = "red")

grid();
#par(xpd = T)
legend("topright", bty = "n", inset = 0.005,
       c(expression(paste(M[1])), expression(paste(M[2])) , expression(paste(M[3])),expression(paste(M[4])),expression(paste(M[5])),expression(paste(M[6]))),
       cex = 0.8,lty=c(6,2,1,2,3),pch = c(NA, NA,NA,2,21),lwd=c(4,2,3,2,2), col=c("purple", "cyan", "green", "orange","red")) 

par(mar = c(5, 4, 4, 2) + 0.1)


########################################################################################################################




