

rm(list = ls())

setwd("C:/Users/doforib/Desktop/Attack trials/Yuzhou_code")


#####################################################################################################################

##############
#=== DATA ===#
##############

dm01 = read.csv("Motifs_UNDIR_3N_80.csv")

#################
#=== LIBRARY ===#
#################

library(survival) 

#####################################################################################################################

#=============== Survival curves of V1 ===============#

D1_V1          =  abs(diff(dm01$V_1, lag = 1)) # Death #
Time1_V1       =  dm01[-1,2]    # delete first one #
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
plot(NPsV2,conf.int=FALSE,lwd=2,lty=1,xlab="Time",ylab="Survival Probability", main="V2", xlim = c(0,80), ylim = c(0,1))
lines(S_V2,col="red",lwd=2,lty=1)
legend('topright',c( "KME","MLE"), lty=c(1,1), lwd=c(2,2), col=c("black","red"))    



#####################################################################################################################


####################################################################
###################### Dependent Comoponent ########################
####################################################################

t = c(0, Time1_V1)      #Time1_V1 = dm01$Time[-1]
DRt = numeric(length(t))

for (i in (1:length(t))) { 
  
  PAc =  ((1-exp(-lamda1*t[i]))*(1-exp(-lamda2*t[i])))#*(1-exp(-lamda3*t[i]))* (1-exp(-lamda4*t[i])))*(1-exp(-lamda5*t[i])))
  
  Mn = min((1-exp(-lamda1*t[i])),(1-exp(-lamda2*t[i])))#,(1-exp(-lamda3*t[i])), (1-exp(-lamda4*t[i])))#,(1-exp(-lamda5*t[i])))
  
  DRt[i] =  1-sqrt(Mn*PAc)
  
}

S_V1 = exp(-lamda1*t) 
S_V2 = exp(-lamda2*t)  
# S_V3 = exp(-lamda3*t) 
# S_V4 = exp(-lamda4*t) 
#S_V5 = exp(-lamda5*t) 


D_d1 = data.frame(t,S_V1,S_V2,DRt)#S_V3,S_V4,DRt)
colnames(D_d1) = c("V_Remv", "S_V1", "S_V2","DRt")

write.csv(D_d1, "Motifs_UNDIR_3N_80Reliab.csv")


