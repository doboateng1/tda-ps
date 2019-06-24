
rm(list = ls())

setwd("C:/Users/doforib/Desktop/Attack trials/Undirected_graph_analysis")


#####################################################################################################################

##############
#=== DATA ===#
##############

dm01 = read.csv("Motifs_UNDIR_4N_80.csv")

#################
#=== LIBRARY ===#
#################

library(survival) 

#####################################################################################################################

#######################################################
#=============== Survival curves of V1 ===============#
#######################################################


D1_V1          =  abs(diff(dm01$V_1, lag = 1)) # Death #
Time1_V1       = dm01[-1,2]     # delete first one #
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
plot(NPsV1, conf.int=FALSE, lwd=2, lty=1, xlab="Time", ylab="Survival Probability", main="V1", xlim = c(0,80), ylim = c(0,1))
lines(S_V1, col="red",lwd=2,lty=1)
legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    


# Kaplan-Meier estimate (KME) #
# MLE - Max likelihood estimate (Exponential model) #

#######################################################
#=============== Survival curves of V2 ===============#
#######################################################


D1_V2          = abs(diff(dm01$V_2, lag = 1)) # Death
Time1_V2       = dm01[-1,2]           
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
plot(NPsV2,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability",main="V2", xlim = c(0,80), ylim = c(0,1))
lines(S_V2,col="red",lwd=2,lty=1)
legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    


#######################################################
#=============== Survival curves of V3 ===============#
#######################################################

D1_V3          = abs(diff(dm01$V_3, lag = 1)) # Death
Time1_V3       = dm01[-1,2]    
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
plot(NPsV3,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V3",xlim = c(0,80), ylim = c(0,1))
lines(S_V3,col="red",lwd=2,lty=1)
legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    


#######################################################
#=============== Survival curves of V4 ===============#
#######################################################

D1_V4          = abs(diff(dm01$V_4, lag = 1)) # Death
Time1_V4       = dm01[-1,2]   
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
plot(NPsV4,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V4",xlim = c(0,80), ylim = c(0,1))
lines(S_V4,col="red",lwd=2,lty=1)
legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    



#######################################################
#=============== Survival curves of V5 ===============#
#######################################################

D1_V5          = abs(diff(dm01$V_5, lag = 1)) # Death
Time1_V5       = dm01[-1,2]   
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
plot(NPsV5,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V5",xlim = c(0,80), ylim = c(0,1))
lines(S_V5,col="red",lwd=2,lty=1)
legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    



#######################################################
#=============== Survival curves of V6 ===============#
#######################################################

D1_V6          = abs(diff(dm01$V_6, lag = 1)) # Death
Time1_V6       = dm01[-1,2]   
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
plot(NPsV6,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V6",xlim = c(0,80), ylim = c(0,1))
lines(S_V6,col="red",lwd=2,lty=1)
legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    


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


