
rm(list = ls())

#setwd("C:/Users/doforib/Desktop/Attack trials/Undirected_graph_analysis")

setwd("C:/Users/Owner/OneDrive/Desktop/NREL/results")

#####################################################################################################################

#################
#=== LIBRARY ===#
#################

library(survival) 


##############
#=== DATA ===#
##############

dataM = read.csv("Motifs_UNDIR_118-Bus.csv")

#=== 3N ===#
dm01 = dataM[which(colnames(dataM)=="T_1"):which(colnames(dataM)=="T_2")]
dm01 = cbind(dataM[,2], dm01)

#####################################################################################################################

#######################################################
#=============== Survival curves of T1 ===============#
#######################################################


D1_V1          =  abs(diff(dm01$T_1, lag = 1)) # Death #
Time1_V1       =  dm01[-1,1]     # delete first one #
df1            =  data.frame(D1_V1,Time1_V1)
dm1            =  df1[rep(rownames(df1), df1$D1_V1), ]
colnames(dm1)  =  c("Death", "Time") # Making column name same #
dm1$event1     =  1



#========= Non-parametric Survival Model ==========#
NPsV1          =  survfit(Surv(Time, event1) ~ 1,  type="kaplan-meier", conf.type="log", data=dm1) 

#========= Parametric Survival Model ==========#
PsV1           =  survreg(Surv(Time, event1) ~ 1,  data=dm1,dist="exponential")
lamda1         = exp(-PsV1$coefficients[[1]])
t              = c(Time1_V1)
S_V1           = exp(-lamda1*t) 


#========= Plot nonparametric vs parametreic =========#
# plot(NPsV1, conf.int=FALSE, lwd=2, lty=1, xlab="Time", ylab="Survival Probability", main="V1", xlim = c(0,80), ylim = c(0,1))
# lines(S_V1, col="red",lwd=2,lty=1)
# legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    


# Kaplan-Meier estimate (KME) #
# MLE - Max likelihood estimate (Exponential model) #

#######################################################
#=============== Survival curves of T2 ===============#
#######################################################


D1_V2          = abs(diff(dm01$T_2, lag = 1)) # Death
Time1_V2       = dm01[-1,1]           
df2            = data.frame(D1_V2,Time1_V2)
dm2            = df2[rep(rownames(df2), df2$D1_V2), ]
colnames(dm2)  = c("Death", "Time") #Making column name same
dm2$event1     = 1


#========= Non-parametric Survival Model ==========#
NPsV2          = survfit(Surv(Time, event1) ~ 1,  type="kaplan-meier", conf.type="log", data=dm2)

#========= Parametric Survival Model ==========#
PsV2           = survreg(Surv(Time, event1) ~ 1,  data=dm2, dist="exponential")
lamda2         = exp(-PsV2$coefficients[[1]])# 0.03107279
S_V2           = exp(-lamda2*t) 


#========= Plot nonparametric vs parametreic =========#
# plot(NPsV2,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability",main="V2", xlim = c(0,80), ylim = c(0,1))
# lines(S_V2,col="red",lwd=2,lty=1)
# legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))




####################################################################
###################### Dependent Comoponent ########################
####################################################################

t = c(0, Time1_V1)  

DRt = numeric(length(t))

for (i in (1:length(t))) { #i=1
  
  PAc =  ((1-exp(-lamda1*t[i]))*(1-exp(-lamda2*t[i])))#*(1-exp(-lamda3*t[i]))* (1-exp(-lamda4*t[i]))*(1-exp(-lamda5*t[i]))*(1-exp(-lamda6*t[i])))
  
  Mn = min((1-exp(-lamda1*t[i])),(1-exp(-lamda2*t[i])))#,(1-exp(-lamda3*t[i])), (1-exp(-lamda4*t[i])),(1-exp(-lamda5*t[i])),(1-exp(-lamda6*t[i])))
  
  DRt[i] =  1-sqrt(Mn*PAc)
  
}


S_T1 = exp(-lamda1*t) 
S_T2 = exp(-lamda2*t)  


D_d1 = data.frame(t,S_T1,S_T2,DRt)
colnames(D_d1) = c("V_Remv", "S_T1", "S_T2","DRt")

write.csv(D_d1, "Motifs_UNDIR_3N_80Reliab.csv")



##############
#=== DATA ===#
##############

##############
#=== DATA ===#
##############

#=== TDA ===#
res1 = read.csv("TDA_DIR_118.csv")
resX = res1[,2]
res1 = (as.data.frame(apply(res1, 2, function(x){1 - (x/max(x))})))[,-1]
colnames(res1) = c("V_Remv", "Wass01")
res1$V_Remv = resX


res2 = read.csv("TDA_UNDIR_118.csv")
res2 = (as.data.frame(apply(res2, 2, function(x){1 - (x/max(x))})))[,-1]
colnames(res2) = c("V_Remv", "Wass01")
res2$V_Remv = resX

#=== LOAD SERVED ===#
LS_vals = read.csv("LS_Value.csv")[2:82,]
LS_vals[,1] = 0:80
colnames(LS_vals) = c("Time", "LS")

#=== GCC ===#
GCC_vals = read.csv("GCC_UNDIR_118.csv")[,-1]
GCC_vals = apply(GCC_vals, 2, function(x){x/max(x)})
GCC_vals[,1] = resX
colnames(GCC_vals) = c("Time", "GCC")

#=== CONNECTIVITY LOSS ===#
LC_vals = read.csv("LC_UNDIR_118.csv")[,-1]
LC_vals = apply(LC_vals, 2, function(x){1-(x/max(x))})
LC_vals[,1] = resX
colnames(LC_vals) = c("Time", "LC")


#=== FLOW SERVED ===#
FL_vals = read.csv("FL_value.csv")[,-1]
FL_vals = FL_vals[-1]



#=== SMALL WORLDNESS PROPERTY ===#
res3 = read.csv("SWProp_UNDIR_118.csv")

APL_vals = res3[,c(2,3)]
APL_vals = apply(APL_vals, 2, function(x){x/max(x)})
APL_vals[,1] = resX
colnames(APL_vals) = c("Time", "APL")


CLUS_vals = res3[,c(2,4)]
CLUS_vals = apply(CLUS_vals, 2, function(x){(x/max(x))})
CLUS_vals[,1] = resX
colnames(CLUS_vals) = c("Time", "CLUS")

DIAM_vals = res3[,c(2,5)]
DIAM_vals = apply(DIAM_vals, 2, function(x){(x/max(x))})
DIAM_vals[,1] = resX
colnames(DIAM_vals) = c("Time", "DIAM")


###################################################################################################

par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(NPsV1,conf.int=FALSE, xlab="Number of node(s) removed",ylab="Load served/ Robustness metric", xlim = c(0,80), ylim = c(0,1),
     main="Undirected: Robustness metrics (3-N motifs,TDA only) & Load served",  lty = 1, lwd = 2, col = "purple")
lines(NPsV2,conf.int=FALSE,   type = "o", pch = 7, lty = 1, lwd = 2, col = "red")
lines(c(0:80),LS_vals$LS,     type = "o", pch = 21,lty = 1, lwd = 1, col = "orange")
lines(res2$V_Remv, res2$Wass01, type = "l",        lty = 6, lwd = 3, col = "brown")
grid();
#par(xpd = T)
legend("bottomleft", bty = "n", inset = 0.005,
       c(expression(paste(T[1])), expression(paste(T[2])),expression(paste(Load," ",served)),expression(paste(1," - ",Delta ,W[2](D[0],D[P])))),
       cex = 1.2,lty=c(1,1,1,6),pch = c(NA, 7,21,NA),lwd=c(2,2,1,3), col=c("purple","red","orange", "brown")) 

par(mar = c(5, 4, 4, 2) + 0.1)





par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(res2$V_Remv, res2$Wass01, xlab="Number of node(s) removed",ylab="Load served/Robustness metric", xlim = c(0,80), ylim = c(0,1),
     main="Undirected: Robustness metrics (3-N motifs) & Load served", type = "l",lty = 1, lwd = 2, col = "purple")
lines(t,LS_vals$LS,     type = "o", pch = 5, lty = 1, lwd = 1, col = "orange")
lines(t,DRt,            type = "l",          lty = 6, lwd = 4, col = "brown")
lines(t,GCC_vals[,2],   type = "l",          lty = 2, lwd = 1, col = "red")
lines(t,LC_vals[,2],    type = "o", pch = 1, lty = 1, lwd = 2, col = "green")
lines(t,APL_vals[,2],   type = "l",          lty = 3, lwd = 1, col = "blue")

grid();
#par(xpd = T)
legend("bottomleft", bty = "n", inset = 0.005,
       c(expression(paste(1," - ",Delta ,W[2](D[0],D[P]))),expression(paste(Load," ",served)),expression(paste(MR)), expression(paste(1, "-", Delta, S)), expression(paste(1,"-",Delta, LC)), expression(paste(Delta, APL))),
       cex = 0.9,lty=c(1,1,6,2,1,3),pch = c(NA,5,NA,NA,1,NA),lwd=c(2,1,4,1,2,1), col=c("purple","orange", "brown","red", "green", "blue")) 

par(mar = c(5, 4, 4, 2) + 0.1)













#####################################################################################################################

#=== 4N ===#
dm01 = dataM[which(colnames(dataM)=="V_1"):which(colnames(dataM)=="V_6")]
dm01 = cbind(dataM[,2], dm01)

#######################################################
#=============== Survival curves of V1 ===============#
#######################################################


D1_V1          =  abs(diff(dm01$V_1, lag = 1)) # Death #
Time1_V1       = dm01[-1,1]     # delete first one #
df1            = data.frame(D1_V1,Time1_V1)
dm1            =  df1[rep(rownames(df1), df1$D1_V1), ]
colnames(dm1)  =  c("Death", "Time") # Making column name same #
dm1$event1     = 1



#========= Non-parametric Survival Model ==========#
NPsV1          =  survfit(Surv(Time, event1) ~ 1,  type="kaplan-meier", conf.type="log", data=dm1) 

#========= Parametric Survival Model ==========#
PsV1           =  survreg(Surv(Time, event1) ~ 1,  data=dm1,dist="exponential")
lamda1         = exp(-PsV1$coefficients[[1]])
t              = c(Time1_V1)
S_V1           = exp(-lamda1*t) 


#========= Plot nonparametric vs parametreic =========#
# plot(NPsV1, conf.int=FALSE, lwd=2, lty=1, xlab="Time", ylab="Survival Probability", main="V1", xlim = c(0,80), ylim = c(0,1))
# lines(S_V1, col="red",lwd=2,lty=1)
# legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    


# Kaplan-Meier estimate (KME) #
# MLE - Max likelihood estimate (Exponential model) #

#######################################################
#=============== Survival curves of V2 ===============#
#######################################################


D1_V2          = abs(diff(dm01$V_2, lag = 1)) # Death
Time1_V2       = dm01[-1,1]           
df2            = data.frame(D1_V2,Time1_V2)
dm2            = df2[rep(rownames(df2), df2$D1_V2), ]
colnames(dm2)  = c("Death", "Time") #Making column name same
dm2$event1     = 1


#========= Non-parametric Survival Model ==========#
NPsV2          = survfit(Surv(Time, event1) ~ 1,  type="kaplan-meier", conf.type="log", data=dm2)

#========= Parametric Survival Model ==========#
PsV2           = survreg(Surv(Time, event1) ~ 1,  data=dm2, dist="exponential")
lamda2         = exp(-PsV2$coefficients[[1]])# 0.03107279
S_V2           = exp(-lamda2*t) 


#========= Plot nonparametric vs parametreic =========#
# plot(NPsV2,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability",main="V2", xlim = c(0,80), ylim = c(0,1))
# lines(S_V2,col="red",lwd=2,lty=1)
# legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    
# 

#######################################################
#=============== Survival curves of V3 ===============#
#######################################################

D1_V3          = abs(diff(dm01$V_3, lag = 1)) # Death
Time1_V3       = dm01[-1,1]    
df3            = data.frame(D1_V3,Time1_V3)
dm3            = df3[rep(rownames(df3), df3$D1_V3), ]
colnames(dm3)  = c("Death", "Time") #Making column name same
dm3$event1     = 1


#========= Non-parametric Survival Model ==========#
NPsV3          =  survfit(Surv(Time, event1) ~ 1,  type="kaplan-meier", conf.type="log", data=dm3)

#========= Parametric Survival Model ==========#
PsV3           =  survreg(Surv(Time, event1) ~ 1,  data=dm3, dist="exponential")
lamda3         = exp(-PsV3$coefficients[[1]]) #lemda=exp(-Intercept)
S_V3           = exp(-lamda3*t) 


#========= Plot nonparametric vs parametreic =========#
# plot(NPsV3,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V3",xlim = c(0,80), ylim = c(0,1))
# lines(S_V3,col="red",lwd=2,lty=1)
# legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    


#######################################################
#=============== Survival curves of V4 ===============#
#######################################################

D1_V4          = abs(diff(dm01$V_4, lag = 1)) # Death
Time1_V4       = dm01[-1,1]   
df4            = data.frame(D1_V4,Time1_V4)
dm4            = df4[rep(rownames(df4), df4$D1_V4), ]
colnames(dm4)  = c("Death", "Time") # Making column name same
dm4$event1     = 1

#========= Non-parametric Survival Model ==========#
NPsV4          =  survfit(Surv(Time, event1) ~ 1,  type="kaplan-meier", conf.type="log", data=dm4)

#========= Parametric Survival Model ==========#
PsV4           =  survreg(Surv(Time, event1) ~ 1,  data=dm4, dist="exponential")
lamda4         = exp(-PsV4$coefficients[[1]]) #lemda=exp(-Intercept)
S_V4           = exp(-lamda4*t) 


#========= Plot nonparametric vs parametreic =========#
# plot(NPsV4,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V4",xlim = c(0,80), ylim = c(0,1))
# lines(S_V4,col="red",lwd=2,lty=1)
# legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    



#######################################################
#=============== Survival curves of V5 ===============#
#######################################################

D1_V5          = abs(diff(dm01$V_5, lag = 1)) # Death
Time1_V5       = dm01[-1,1]   
df5            = data.frame(D1_V5,Time1_V5)
dm5            = df5[rep(rownames(df5), df5$D1_V5), ]
colnames(dm5)  = c("Death", "Time") # Making column name same
dm5$event1     = 1

#========= Non-parametric Survival Model ==========#
NPsV5          =  survfit(Surv(Time, event1) ~ 1,  type="kaplan-meier", conf.type="log", data=dm5)

#========= Parametric Survival Model ==========#
PsV5           =  survreg(Surv(Time, event1) ~ 1,  data=dm5, dist="exponential")
lamda5         = exp(-PsV5$coefficients[[1]]) #lemda=exp(-Intercept)
S_V5           = exp(-lamda5*t) 


#========= Plot nonparametric vs parametreic =========#
# plot(NPsV5,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V5",xlim = c(0,80), ylim = c(0,1))
# lines(S_V5,col="red",lwd=2,lty=1)
# legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    



#######################################################
#=============== Survival curves of V6 ===============#
#######################################################

D1_V6          = abs(diff(dm01$V_6, lag = 1)) # Death
Time1_V6       = dm01[-1,1]   
df6            = data.frame(D1_V6,Time1_V6)
dm6            = df6[rep(rownames(df6), df6$D1_V6), ]
colnames(dm6)  = c("Death", "Time") # Making column name same
dm6$event1     = 1

#========= Non-parametric Survival Model ==========#
NPsV6          =  survfit(Surv(Time, event1) ~ 1,  type="kaplan-meier", conf.type="log", data=dm6)

#========= Parametric Survival Model ==========#
PsV6           =  survreg(Surv(Time, event1) ~ 1,  data=dm6, dist="exponential")
lamda6         = exp(-PsV6$coefficients[[1]]) #lemda=exp(-Intercept)
S_V6           = exp(-lamda6*t) 


#========= Plot nonparametric vs parametreic =========#
# plot(NPsV6,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V6",xlim = c(0,80), ylim = c(0,1))
# lines(S_V6,col="red",lwd=2,lty=1)
# legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    


####################################################################
###################### Dependent Comoponent ########################
####################################################################

t = c(0, Time1_V1)  

DRt = numeric(length(t))

for (i in (1:length(t))) { #i=1
  
  PAc =  ((1-exp(-lamda1*t[i]))*(1-exp(-lamda2*t[i]))*(1-exp(-lamda3*t[i]))* (1-exp(-lamda4*t[i]))*(1-exp(-lamda5*t[i]))*(1-exp(-lamda6*t[i])))
  
  Mn = min((1-exp(-lamda1*t[i])),(1-exp(-lamda2*t[i])),(1-exp(-lamda3*t[i])), (1-exp(-lamda4*t[i])),(1-exp(-lamda5*t[i])),(1-exp(-lamda6*t[i])))
  
  DRt[i] =  1-sqrt(Mn*PAc)
  
}


S_V1 = exp(-lamda1*t) 
S_V2 = exp(-lamda2*t)  
S_V3 = exp(-lamda3*t) 
S_V4 = exp(-lamda4*t) 
S_V5 = exp(-lamda5*t) 
S_V6 = exp(-lamda6*t) 



D_d1 = data.frame(t,S_V1,S_V2,S_V3,S_V4,S_V5,S_V6,DRt)
colnames(D_d1) = c("V_Remv", "S_V1", "S_V2", "S_V3", "S_V4", "S_V5", "S_V6", "DRt")

write.csv(D_d1, "Motifs_UNDIR_4N_80Reliab.csv")



##################################################################################################################

##############
#=== DATA ===#
##############

#=== TDA ===#
res1 = read.csv("TDA_DIR_118.csv")
resX = res1[,2]
res1 = (as.data.frame(apply(res1, 2, function(x){1 - (x/max(x))})))[,-1]
colnames(res1) = c("V_Remv", "Wass01")
res1$V_Remv = resX


res2 = read.csv("TDA_UNDIR_118.csv")
res2 = (as.data.frame(apply(res2, 2, function(x){1 - (x/max(x))})))[,-1]
colnames(res2) = c("V_Remv", "Wass01")
res2$V_Remv = resX

#=== LOAD SERVED ===#
LS_vals = read.csv("LS_Value.csv")[2:82,]
LS_vals[,1] = 0:80
colnames(LS_vals) = c("Time", "LS")

#=== GCC ===#
GCC_vals = read.csv("GCC_UNDIR_118.csv")[,-1]
GCC_vals = apply(GCC_vals, 2, function(x){x/max(x)})
GCC_vals[,1] = resX
colnames(GCC_vals) = c("Time", "GCC")

#=== CONNECTIVITY LOSS ===#
LC_vals = read.csv("LC_UNDIR_118.csv")[,-1]
LC_vals = apply(LC_vals, 2, function(x){1-(x/max(x))})
LC_vals[,1] = resX
colnames(LC_vals) = c("Time", "LC")


#=== FLOW SERVED ===#
FL_vals = read.csv("FL_value.csv")[,-1]
FL_vals = FL_vals[-1,]
FL_vals[,1] = resX
colnames(FL_vals) = c("Time", "FL")



#=== SMALL WORLDNESS PROPERTY ===#
res3 = read.csv("SWProp_UNDIR_118.csv")

APL_vals = res3[,c(2,3)]
APL_vals = apply(APL_vals, 2, function(x){x/max(x)})
APL_vals[,1] = resX
colnames(APL_vals) = c("Time", "APL")


CLUS_vals = res3[,c(2,4)]
CLUS_vals = apply(CLUS_vals, 2, function(x){(x/max(x))})
CLUS_vals[,1] = resX
colnames(CLUS_vals) = c("Time", "CLUS")

DIAM_vals = res3[,c(2,5)]
DIAM_vals = apply(DIAM_vals, 2, function(x){(x/max(x))})
DIAM_vals[,1] = resX
colnames(DIAM_vals) = c("Time", "DIAM")



###################################################################################################

par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(NPsV1,conf.int=FALSE, xlab="Number of node(s) removed",ylab="Load served/Robustness metric", xlim = c(0,80), ylim = c(0,1),
     main="Undirected: Robustness metrics (4-N motifs, TDA only) & Load served",  lty = 1, lwd = 2, col = "purple")
lines(NPsV2,conf.int=FALSE,   type = "o", pch = 7, lty = 1, lwd = 2, col = "red")
lines(NPsV3,conf.int=FALSE,   type = "l",          lty = 4, lwd = 2, col = "cyan")
lines(NPsV4,conf.int=FALSE,   type = "l",          lty = 1, lwd = 3, col = "green")
lines(NPsV5,conf.int=FALSE,   type = "l",          lty = 2, lwd = 2, col = "blue")
lines(NPsV6,conf.int=FALSE,   type = "s",          lty = 3, lwd = 2, col = "black")
lines(c(0:80),LS_vals$LS,     type = "o", pch = 21,lty = 1, lwd = 1, col = "orange")
lines(res2$V_Remv, res2$Wass01, type = "l",        lty = 6, lwd = 3, col = "brown")
grid();
#par(xpd = T)
legend("bottomleft", bty = "n", inset = 0.005,
       c(expression(paste(M[1])), expression(paste(M[2])) , expression(paste(M[3])),expression(paste(M[4])),expression(paste(M[5])),expression(paste(M[6])),expression(paste(Load," ",served)),expression(paste(1," - ",Delta ,W[2](D[0],D[P])))),
       cex = 1.2,lty=c(1,1,4,1,2,3,1,6),pch = c(NA, 7,NA,NA,NA,NA,21,NA),lwd=c(2,2,2,3,2,2,1,3), col=c("purple","red", "cyan", "green", "blue","black", "orange", "brown")) 

par(mar = c(5, 4, 4, 2) + 0.1)



par(mar = c(5, 4, 4, 2) + 0.1)

par(xpd = F)

plot(res2$V_Remv, res2$Wass01, xlab="Number of node(s) removed",ylab="Load served/Robustness metric", xlim = c(0,80), ylim = c(0,1),
     main="Undirected: Robustness metrics (4-N motifs) & Load served", type = "l",lty = 1, lwd = 2, col = "purple")
lines(t,LS_vals$LS,     type = "o", pch = 5, lty = 1, lwd = 1, col = "orange")
lines(t,DRt,            type = "l",          lty = 6, lwd = 3, col = "brown")
lines(t,GCC_vals[,2],   type = "l",          lty = 2, lwd = 1, col = "red")
lines(t,LC_vals[,2],    type = "o", pch = 1, lty = 1, lwd = 2, col = "green")
lines(t,APL_vals[,2],   type = "l",          lty = 3, lwd = 3, col = "blue")

grid();
#par(xpd = T)
legend("bottomleft", bty = "n", inset = 0.005,
       c(expression(paste(1," - ",Delta ,W[2](D[0],D[P]))),expression(paste(Load," ",served)),expression(paste(MR)), expression(paste(1, "-", Delta, S)), expression(paste(1,"-",Delta, LC)), expression(paste(Delta, APL))),
       cex = 0.9,lty=c(1,1,6,2,1,3),pch = c(NA,5,NA,NA,1,NA),lwd=c(2,1,3,1,2,3), col=c("purple","orange", "brown","red", "green", "blue")) 

par(mar = c(5, 4, 4, 2) + 0.1)


plot(t, CLUS_vals[,2], type = "o", col = "blue", main = "Normalized clustering coefficient", xlab = "Number of node(s) removed", ylab = "Normalized change")
plot(t, DIAM_vals[,2], type = "o", col = "red", main = "Normalized diameter", xlab = "Number of node(s) removed", ylab = "Normalized change")


##################################################################################################

dataCOR = cbind(LS_vals[,2], res2$Wass01, DRt, LC_vals[,2], GCC_vals[,2], APL_vals[,2])
colnames(dataCOR) = c("LS", "Wass", "MR", "LC","GCC","APL")


#######################
#=== SCATTER PLOTS ===#
#######################

pairs(dataCOR)


#####################
#=== CORRELATION ===#
#####################

library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(PerformanceAnalytics)


#==============================================================================#

my_data <- dataCOR
chart.Correlation(my_data, histogram=TRUE, pch=19)


#==============================================================================#


M = cor(dataCOR)
corrplot(M, type = "upper", order = "hclust",
         col = brewer.pal(n = 8, name = "RdYlBu"))


source("http://www.sthda.com/upload/rquery_cormat.r")

require("corrplot")
rquery.cormat(dataCOR)






#==============================================================================#
