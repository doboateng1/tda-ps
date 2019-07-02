
rm(list = ls())

#setwd("C:/Users/doforib/Desktop/Attack trials/Undirected_graph_analysis")

setwd("C:/Users/Owner/OneDrive/Desktop/NREL/results")


#####################################################################################################################

##############
#=== DATA ===#
##############

res1 = read.csv("TDA_DIR_118.csv")
resX = res1[,2]
res1 = (as.data.frame(apply(res1, 2, function(x){1 - (x/max(x))})))[,-1]
colnames(res1) = c("V_Remv", "Wass01")
res1$V_Remv = resX


res2 = read.csv("TDA_UNDIR_118.csv")
res2 = (as.data.frame(apply(res2, 2, function(x){1 - (x/max(x))})))[,-1]
colnames(res2) = c("V_Remv", "Wass01")
res2$V_Remv = resX



res3 = read.csv("Motifs_DIR_3N_80Reliab.csv")
res3 = subset(res3, select = -X)

res4 = read.csv("Motifs_UNDIR_3N_80Reliab.csv")
res4 = subset(res4, select = -X)


res5 = read.csv("Motifs_UNDIR_4N_80Reliab.csv")
res5 = subset(res5, select = -X)




res5 = read.csv("Motifs_UNDIR_3N_80.csv")
res6 = read.csv("Motifs_UNDIR_4N_80.csv")


########################################################################################################################


setwd("C:/Users/Owner/OneDrive/Desktop/NREL/results/plots")



########################################################################################################################################################################################
#######################################################################################################################################################

#######################################################################################################################################################
# DIRECTED #
#######################################################################################################################################################


par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(resX, res1$Wass01, type = "l", xlab = "Number of node(s) removed",lty=6, lwd=4, main = "Directed: TDA(Wass) vs. 3N-joint motif reliability/survival probability",
     ylab =  "Robustness metric",col = "purple", ylim = c(0,1), xlim = c(0,80))                                 #Wass01#
lines(resX, res3$DRt,      type = "l",          lty = 2, lwd=2, col = "cyan")

grid();
#par(xpd = T)
legend("topright", bty = "n", inset = 0.005,
       c(expression(paste(1," - ",Delta ,W[2](D[0],D[P]))), expression(paste(Delta, MR))),
       cex = 1.4,lty=c(6,2),pch = c(NA, NA),lwd=c(4,2), col=c("purple","green")) 

par(mar = c(5, 4, 4, 2) + 0.1)




par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(resX, res1$Wass01, type = "l", xlab = "Number of node(s) removed",lty=6, lwd=4, main = "Directed: TDA(Wass) vs. 3N-motifs survival probability",
     ylab =  "Robustness metric",col = "brown", ylim = c(0,1), xlim = c(0,80))                                 #Wass01#
lines(resX, res3$S_V1,                           lty = 4, lwd=2, col = "red")
lines(resX, res3$S_V2,      type = "l",          lty = 6, lwd=2, col = "blue")
lines(resX, res3$S_V4,      type = "s",          lty = 1, lwd=2, col = "cyan")
lines(resX, res3$S_V5,      type = "s",          lty = 3, lwd=2, col = "green")
lines(resX, res3$S_V6,      type = "l",          lty = 1, lwd=2, col = "black")
lines(resX, res3$S_V9,      type = "o",          lty = 3, lwd=2, col = "orange")

grid();
#par(xpd = T)
legend("topright", bty = "n", inset = 0.005,
       c(expression(paste(1," - ",Delta ,W[2](D[0],D[P]))), expression(paste(Delta, S[T[1]])),  expression(paste(Delta, S[T[2]])), expression(paste(Delta, S[T[4]])), expression(paste(Delta, S[T[5]])), expression(paste(Delta, S[T[6]])), expression(paste(Delta, S[T[9]]))),
       cex = 1.4,lty=c(6,4,6,1,3,1,4), lwd=c(4,2,2,2,2,2,3), col=c("brown", "red", "blue", "cyan", "green", "black", "orange")) 

par(mar = c(5, 4, 4, 2) + 0.1)




#######################################################################################################################################################
# UNDIRECTED #
#######################################################################################################################################################


par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(resX, res2$Wass01, type = "l", xlab = "Number of node(s) removed",lty=6, lwd=4, main = "Undirected: TDA(Wass) vs. 3N-joint motif reliability/survival probability",
     ylab =  "Robustness metric",col = "purple", ylim = c(0,1), xlim = c(0,80))                                 #Wass01#
lines(resX, res4$DRt,      type = "l",          lty = 2, lwd=2, col = "red")

grid();
#par(xpd = T)
legend("topright", bty = "n", inset = 0.005,
       c(expression(paste(1," - ",Delta ,W[2](D[0],D[P]))), expression(paste(Delta, MR))),
       cex = 1.4,lty=c(6,2),lwd=c(4,2), col=c("purple","red")) 

par(mar = c(5, 4, 4, 2) + 0.1)




par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(resX, res2$Wass01, type = "l", xlab = "Number of node(s) removed",lty=6, lwd=4, main = "Undirected: TDA(Wass) vs. 3N-motifs survival probability",
     ylab =  "Robustness metric",col = "brown", ylim = c(0,1), xlim = c(0,80))                                 #Wass01#
lines(resX, res4$S_T1,      type = "s",          lty = 4, lwd=2, col = "red")
lines(resX, res4$S_T2,      type = "s",          lty = 6, lwd=2, col = "blue")

grid();
#par(xpd = T)
legend("topright", bty = "n", inset = 0.005,
       c(expression(paste(1," - ",Delta ,W[2](D[0],D[P]))), expression(paste(Delta, S[T[1]])),  expression(paste(Delta, S[T[2]])), expression(paste(Delta, S[T[4]])), expression(paste(Delta, S[T[5]])), expression(paste(Delta, S[T[6]])), expression(paste(Delta, S[T[9]]))),
       cex = 1.4,lty=c(6,4,6,1,3,1,4), lwd=c(4,2,2,2,2,2,3), col=c("brown", "red", "blue", "cyan", "green", "black", "orange")) 

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
