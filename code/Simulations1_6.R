rm(list=ls())
library(lme4);library(ggplot2);library(cowplot); library(gridExtra); library(grid);library(patchwork)

iter1 = 1000
iter2 = 1


###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 45 genotypes; 4 reps each (180)
#effect:: keeps same effects for the simulation

df = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,900)))
  b = as.factor(sort(rep(0:719,5)))  
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:900),]
  mm = mm[,-c(1,5:183)]
  v = sample(sqrt(c(seq(1,100,by=.5), seq(101, 1500,by=5))),3,replace=T)
  
  vL = v[1]
  vG = v[2]
  vE = 5
  

  u.L = matrix((c(rnorm(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm(540,0,vG)),nrow=540)
  
  L = as.factor(sort(rep(0:2,900)))
  G = as.factor(sort(rep(0:539,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:543)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm(2700,0,vE),nrow=2700)

    y2 = rowSums(y) + e
    
    #expectations
    def$vQ[i] = var(rowSums(y))
    def$h.exp[i] = (var(y[,1]) + var(y[,2])) / (var(y[,1]) + var(y[,2]) + var(e))
    def$hm.exp[i] = (var(y[,1])) / (var(y[,1]) + var(y[,2]) + var(e))
    def$vL[i] = var(y[,1])
    def$vG[i] = var(y[,2])
    def$vE[i] = var(e)
    def$vP[i] = var(y2)
    
    data = data.frame(Y=y2, L=L, G=G)
    
    #random effect model
    model <- lmer(Y ~(1|L/G),data,verbose=F,
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model <- update(model,start = getME(model,c("theta","fixef")))
    model <- update(model,start = getME(model,c("theta","fixef")))
    
    model.G <- lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    
    def$c.vm[i] = def$vm[i] * ((2 * 900) / (539 * 5))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] =  def$c.vq[i] - def$vg[i]
    def$diff[i] =  def$vq[i] - def$vg[i]
    
  }
  
  df = rbind(df,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

save(list="df",file="~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_5R.Rdata")

load("~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_5R.Rdata")
#df = df[which((df$hm.c - df$hm.exp) < .25 & (df$hm.c - df$hm.exp) > -.25),]

(b1 = ggplot(data=as.data.frame(df)) +
  geom_point(aes(x = df[,4], y = df[,14]),alpha=0.5,color="darkred") +
  geom_point(aes(x = df[,4], y = df[,19]),alpha=0.5,color="darkblue") +
  geom_abline(slope = 1,intercept = 0) +
  ggtitle(bquote("One Locus")) + 
  xlab(bquote(H[M]^2)) +ylab(bquote(hat(H)[M]^2 ~ " and "~ hat(H)[M~"*"]^2))+
  theme_bw() +   theme(plot.title = element_text(size = 10),
                       axis.title = element_text(size = 10),
                       axis.text = element_text(size=10))+
  xlim(c(0,1)) + ylim(c(0,2)))

###################################################################################
###################################################################################
###################################################################################
# 2 locus: 3 genotypes each ; 60 reps each (360)
# 45 genotypes; 4 reps each (180)
#effect:: keeps same effects for the simulation

df5 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,900)))
  b = as.factor(c(rep(0,900),rep(c(rep(1,300),rep(2,300),rep(3,300)),3)))
  c= as.factor(sort(rep(0:719,5)))
  mm = model.matrix(~(a + b) + c)
  
  mm = mm[-c(1:900),]
  mm = mm[,-c(1,8:186)]
  
  v = sample(sqrt(c(seq(1,100,by=.5), seq(101, 1500,by=5))),4,replace=T)
  
  vL1 = v[1]
  vL2 = v[2]
  vG = v[3]
  vE = 5
  
  u.L1 = matrix((c(rnorm(3,0,vL1))),nrow=3)
  u.L2 = matrix((c(rnorm(3,0,vL2))),nrow=3)
  u.G = matrix(c(rnorm(540,0,vG)),nrow=540)
  
  L1 = as.factor(sort(rep(0:2,900)))
  L2 = as.factor(rep(c(rep(0,300),rep(1,300),rep(2,300)),3))
  G = as.factor(sort(rep(0:539,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L1
  y2 = mm[,c(4:6)] %*% u.L2
  y3 = mm[,c(7:546)] %*% u.G
  y = cbind(y1,y2,y3)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm(2700,0,vE),nrow=2700)
    
    y2 = rowSums(y) + e
    
    #expectations
    def$vQ[i] = var(rowSums(y))
    def$h.exp[i] = (var(y[,1]) + var(y[,2]) +var(y[,3])) / (var(y[,1]) + var(y[,2]) +var(y[,3]) + var(e))
    def$hm.exp[i] = (var(y[,1]) + var(y[,2])) / (var(y[,1]) + var(y[,2]) +var(y[,3]) + var(e))
    def$vL1[i] = var(y[,1])
    def$vL2[i] = var(y[,2])
    def$vG[i] = var(y[,3])
    def$vE[i] = var(e)
    def$vP[i] = var(y2)
    
    data = data.frame(Y=y2, L1=L1, L2=L2, G=G)
    data$L12 = interaction(data$L1, data$L2) 
    data$L12G = interaction(data$L1, data$L2,data$G,drop=T) 

    #random effect model
    model = lmer(Y ~(1|L1) + (1|L2) + (1|L12) + (1|L12G),data,verbose=F,
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model <- update(model,start = getME(model,c("theta","fixef")))
    model <- update(model,start = getME(model,c("theta","fixef")))
    
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    def$v12g[i] = VarCorr(model)[[1]][1]
    def$v12[i] = VarCorr(model)[[2]][1]
    def$v2[i] = VarCorr(model)[[3]][1]
    def$v1[i] = VarCorr(model)[[4]][1]
    def$ve[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$c.v1[i] = def$v1[i] * ((2 * 900) / (539 * 5))
    def$c.v2[i] = def$v2[i] * ((2 * 900) / (539 * 5))
    def$c.v12[i] = def$v12[i] * (((2 + 2 + 4) * 300) / (539 * 5))
    
    def$vq[i] = def$v1[i] + def$v2[i] + def$v12[i] + def$v12g[i]
    def$c.vq[i] = def$c.v1[i] + def$c.v2[i] + def$c.v12[i] + def$v12g[i]
    
    def$diff.c[i] =  def$c.vq[i] - def$vg[i]
    def$diff[i] =  def$vq[i] - def$vg[i]
    
    def$h.c[i] = (def$c.v1[i] + def$c.v2[i] + def$c.v12[i] + def$v12g[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.v1[i] + def$c.v2[i] + def$c.v12[i]) / (def$vg[i] + def$ve[i])

    def$h.o[i] = (def$v1[i] + def$v2[i] + def$v12[i] + def$v12g[i]) / (def$vg[i] + def$ve[i])
    def$hm.o[i] = (def$v1[i] + def$v2[i] + def$v12[i]) / (def$vg[i] + def$ve[i])

  }
  
  df5 = rbind(df5,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

save(list="df5",file="~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_2L_5R.Rdata")

load("~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_2L_5R.Rdata")
#df5 = df5[which((df5$hm.c - df5$hm.exp) < .25 & (df5$hm.c - df5$hm.exp) > -.25),]

hist(df5$diff.c)

(b5 =ggplot(data=as.data.frame(df5)) + 
  geom_point(aes(x = df5[,4], y = df5[,26]),alpha=0.5,color="darkred") +
  geom_point(aes(x = df5[,4], y = df5[,24]),alpha=0.5,color="darkblue") +
  geom_abline(slope = 1,intercept = 0) +
  ggtitle(bquote("Two Loci") ) + 
  xlab(bquote(H[M]^2)) +ylab("")+
  theme_bw() +   theme(plot.title = element_text(size = 10),
                       axis.title = element_text(size = 10),
                       axis.text = element_text(size=10),
                       axis.text.y = element_blank(),
                       axis.title.y = element_blank())+
  xlim(c(0,1)) + ylim(c(0,2)))

###################################################################################
###################################################################################
###################################################################################
# 3 locus: 3 genotypes each ; 60 reps each (270)
# 54 genotypes; 5 reps each (180)
#effect:: keeps same effects for the simulation

df8 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,900)))
  b = as.factor(c(rep(0,900),rep(c(rep(1,300),rep(2,300),rep(3,300)),3)))
  c = as.factor(c(rep(0,900),rep(c(rep(1,100),rep(2,100),rep(3,100)),9)))
  d = as.factor(sort(rep(0:719,5)))
  mm = model.matrix(~(a+b + c) + d)
  
  
  mm = mm[-c(1:900),]
  mm = mm[,-c(1,11:189)]
  
  v = sample(sqrt(c(seq(1,100,by=.5), seq(101, 1500,by=5))),5,replace=T)
  
  vL1 = v[1]
  vL2 = v[2]
  vL3 = v[3]
  vG = v[4]
  vE = 5
  
  u.L1 = matrix((c(rnorm(3,0,vL1))),nrow=3)
  u.L2 = matrix((c(rnorm(3,0,vL2))),nrow=3)
  u.L3 = matrix((c(rnorm(3,0,vL3))),nrow=3)
  u.G = matrix(c(rnorm(540,0,vG)),nrow=540)
  
  L1 = as.factor(sort(rep(0:2,900)))
  L2 = as.factor(rep(c(rep(0,300),rep(1,300),rep(2,300)),3))
  L3 = as.factor(rep(c(rep(0,100),rep(1,100),rep(2,100)),9))
  G = as.factor(sort(rep(0:539,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L1
  y2 = mm[,c(4:6)] %*% u.L2
  y3 = mm[,c(7:9)] %*% u.L3
  y4 = mm[,c(10:549)] %*% u.G
  y = cbind(y1,y2,y3,y4)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm(2700,0,vE),nrow=2700)
    
    y2 = rowSums(y) + e
    
    #expectations
    def$vQ[i] = var(rowSums(y))
    def$h.exp[i] = (var(y[,1]) + var(y[,2]) + var(y[,3]) + var(y[,4])) / (var(y[,1]) + var(y[,2]) + var(y[,3]) + var(y[,4]) + var(e))
    def$hm.exp[i] = (var(y[,1]) + var(y[,2]) + var(y[,3])) / (var(y[,1]) + var(y[,2]) + var(y[,3]) + var(y[,4]) + var(e))
    def$vL1[i] = var(y[,1])
    def$vL2[i] = var(y[,2])
    def$vL3[i] = var(y[,3])
    def$vG[i] =  var(y[,4])
    def$vE[i] = var(e)
    def$vP[i] = var(y2)
    
    data = data.frame(Y=y2, L1=L1, L2=L2, L3=L3, G=G)
    data$L12 = interaction(data$L1, data$L2) 
    data$L13 = interaction(data$L1, data$L3) 
    data$L23 = interaction(data$L2, data$L3) 
    data$L123 = interaction(data$L1, data$L2, data$L3) 
    data$L123G = interaction(data$L1, data$L2, data$L3, data$G,drop=T) 
    
    #random effect model
    model = lmer(Y ~(1|L1) + (1|L2) + (1|L3) + (1|L12) + (1|L13) + (1|L23) + (1|L123) + (1|L123G),data,verbose=F,
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model <- update(model,start = getME(model,c("theta","fixef")))
    model <- update(model,start = getME(model,c("theta","fixef")))
    
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    
    def$v123g[i] = VarCorr(model)[[1]][1]
    def$v123[i] = VarCorr(model)[[2]][1]
    def$v23[i] = VarCorr(model)[[3]][1]
    def$v13[i] = VarCorr(model)[[4]][1]
    def$v12[i] = VarCorr(model)[[5]][1]
    def$v3[i] = VarCorr(model)[[6]][1]
    def$v2[i] = VarCorr(model)[[7]][1]
    def$v1[i] = VarCorr(model)[[8]][1]
    def$ve[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$c.v1[i] = def$v1[i] * ((2 * 900) / (539 * 5))
    def$c.v2[i] = def$v2[i] * ((2 * 900) / (539 * 5))
    def$c.v3[i] = def$v3[i] * ((2 * 900) / (539 * 5))
    def$c.v12[i] = def$v12[i] * (((2 + 2 + 4) * 300) / (539 * 5))
    def$c.v13[i] = def$v13[i] * (((2 + 2 + 4) * 300) / (539 * 5))
    def$c.v23[i] = def$v23[i] * (((2 + 2 + 4) * 300) / (539 * 5))
    def$c.v123[i] = def$v123[i] * (((2 + 2 + 2 + 4 + 4 + 4 + 8) * 100) / (539 * 5))
    
    def$vq[i] = def$v1[i] + def$v2[i] + def$v3[i] + 
      def$v12[i] + def$v13[i] + def$v23[i] + 
      def$v123[i] + def$v123g[i]
    
    def$c.vq[i] = def$c.v1[i] + def$c.v2[i] + def$c.v3[i] + 
      def$c.v12[i] + def$c.v13[i] + def$c.v23[i] +
      def$c.v123[i] + def$v123g[i]
    
    def$diff.c[i] =  def$c.vq[i] - def$vg[i]
    def$diff[i] =  def$vq[i] - def$vg[i]
    
    def$h.c[i] = (def$c.v1[i] + def$c.v2[i] + def$c.v3[i] + 
                    def$c.v12[i] + def$c.v13[i] + def$c.v23[i] +
                    def$c.v123[i] + def$v123g[i]) / 
      (def$vg[i] + def$ve[i])
    
    def$hm.c[i] = (def$c.v1[i] + def$c.v2[i] + def$c.v3[i] + 
                     def$c.v12[i] + def$c.v13[i] + def$c.v23[i] +
                     def$c.v123[i]) /  
      (def$vg[i] + def$ve[i])
    
    
    def$h.o[i] = (def$v1[i] + def$v2[i] + def$v3[i] + 
                    def$v12[i] + def$v13[i] + def$v23[i] + 
                    def$v123[i] + def$v123g[i]) / 
      (def$vg[i] + def$ve[i])
    
    def$hm.o[i] = (def$v1[i] + def$v2[i] + def$v3[i] + 
                     def$v12[i] + def$v13[i] + def$v23[i] +
                     def$v123[i]) / 
      (def$vg[i] + def$ve[i])
  }
  
  df8 = rbind(df8,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

save(list="df8",file="~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_3L_5R.Rdata")

load("~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_3L_5R.Rdata")
#df8 = df8[which((df8$hm.c - df8$hm.exp) < .2& (df8$hm.c - df8$hm.exp) > -.2),]
colnames(df8)[33]

(b8 =ggplot(data=as.data.frame(df8)) + 
  geom_point(aes(x = df8[,4], y = df8[,35]),alpha=0.5,color="darkred") +
  geom_point(aes(x = df8[,4], y = df8[,33]),alpha=0.5,color="darkblue") +
  geom_abline(slope = 1,intercept = 0) +
  ggtitle(bquote("Three Loci") )  + 
  xlab(bquote(H[M]^2)) +ylab(bquote(hat(H)[M]^2 ~ " and "~ hat(H)[M~"*"]^2))+
  theme_bw() +   theme(plot.title = element_text(size = 10),
                       axis.title = element_text(size = 10),
                       axis.text = element_text(size=10),
                       axis.text.y = element_blank(),
                       axis.title.y = element_blank())+
  xlim(c(0,1)) + ylim(c(0,2)))


###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 45 genotypes; 4 reps each (180)
#effect:: keeps same effects for the simulation

df11 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(c(rep(0,900),rep(1,675), rep(2,1350), rep(3,675))))
  b= as.factor(sort(rep(0:719,5)))
  mm = model.matrix(~a + b)
  mm = mm[-c(1:900),]
  mm = mm[,-c(1,5:183)]
  v = sample(sqrt(c(seq(1,100,by=.5), seq(101, 1500,by=5))),3,replace=T)
  
  vL = v[1]
  vG = v[2]
  vE = 5
  
  u.L = matrix((c(rnorm(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm(540,0,vG)),nrow=540)
  
  L = as.factor(sort(c(rep(0,675), rep(1,1350), rep(2,675))))
  G = as.factor(sort(rep(0:539,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:543)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm(2700,0,vE),nrow=2700)
    
    y2 = rowSums(y) + e
    
    #expectations
    def$vQ[i] = var(rowSums(y))
    def$h.exp[i] = (var(y[,1]) + var(y[,2])) / (var(y[,1]) + var(y[,2]) + vE**2)
    def$hm.exp[i] = (var(y[,1])) / (var(y[,1]) + var(y[,2]) + var(e))
    def$vL[i] = var(y[,1])
    def$vG[i] = var(y[,2])
    def$vE[i] = var(e)
    def$vP[i] = var(y2)
    
    data = data.frame(Y=y2, L=L, G=G)
    
    #random effect model
    model = lmer(Y ~(1|L/G),data,verbose=F,
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model <- update(model,start = getME(model,c("theta","fixef")))
    model <- update(model,start = getME(model,c("theta","fixef")))
    
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    ssgm = sum((c(135,270,135) - mean(c(135,270,135)))^2)
    
    def$c.vm[i] = def$vm[i] * ((1 / 539) * (((540 * 2)/3) - (ssgm/540)))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] =  def$c.vq[i] - def$vg[i]
    def$diff[i] =  def$vq[i] - def$vg[i]
    
  }
  
  df11 = rbind(df11,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

save(list="df11",file="~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_5R_121_ASV.Rdata")

load("~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_5R_121_ASV.Rdata")
#df11 = df11[which(df11[,20] < 25 & df11[,20] > -25),]



(b13 = ggplot(data=as.data.frame(df11)) +
  geom_point(aes(x = df11[,4], y = df11[,14]),alpha=0.5,color="darkred") +
  geom_point(aes(x = df11[,4], y = df11[,19]),alpha=0.5,color="darkblue") +
  geom_abline(slope = 1,intercept = 0) +
  ggtitle(bquote("One Locus - Unbalanced")) + 
  xlab(bquote(H[M]^2)) +ylab(bquote(hat(H)[M]^2 ~ " and "~ hat(H)[M~"*"]^2))+
  theme_bw() +   theme(plot.title = element_text(size = 10),
                       axis.title = element_text(size = 10),
                       axis.text = element_text(size=10))+
  xlim(c(0,1)) + ylim(c(0,2)))


###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 45 genotypes; 4 reps each (180)
#effect:: keeps same effects for the simulation

df2 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,900)))
  b = as.factor(sort(rep(0:719,5)))  
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:900),]
  mm = mm[,-c(1,5:183)]
  v = sample(sqrt(c(seq(1,100,by=.5), seq(101, 1500,by=5))),3,replace=T)
  
  vL = v[1]
  vG = v[2]
  vE = 5
  
  u.L = matrix((c(rnorm(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm(540,0,vG)),nrow=540)
  
  L = as.factor(sort(rep(0:2,900)))
  G = as.factor(sort(rep(0:539,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:543)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm(2700,0,vE),nrow=2700)
    
    y2 = rowSums(y) + e
    
    x = sample(c(1:1800), 1800 * 0.9, replace = F)
    
    #expectations
    def$vQ[i] = var(rowSums(y[x,]))
    def$h.exp[i] = (var(y[x,1]) + var(y[x,2])) / (var(y[x,1]) + var(y[x,2]) + var(e[x]))
    def$hm.exp[i] = (var(y[x,1])) / (var(y[x,1]) + var(y[x,2]) + var(e[x]))
    def$vL[i] = var(y[x,1])
    def$vG[i] = var(y[x,2])
    def$vE[i] = var(e[x])
    def$vP[i] = var(y2[x])
    
    data = data.frame(Y=y2, L=L, G=G)
    
    
    data = data[x,]
    data$G = droplevels(data$G)   
    data$L = droplevels(data$L)   
    
    #random effect model
    model = lmer(Y ~(1|L/G),data,verbose=F,
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model <- update(model,start = getME(model,c("theta","fixef")))
    model <- update(model,start = getME(model,c("theta","fixef")))
    
    
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    length(unique(data$G[which(data$L == 0)]))
    length(unique(data$G[which(data$L == 1)]))
    length(unique(data$G[which(data$L == 2)]))
    
    
    ssgm = sum((c(length(unique(data$G[which(data$L == 0)])),length(unique(data$G[which(data$L == 1)])),length(unique(data$G[which(data$L == 2)]))) - mean(c(length(unique(data$G[which(data$L == 0)])),length(unique(data$G[which(data$L == 1)])),length(unique(data$G[which(data$L == 2)])))))^2)
    ng = length(unique(data$G))
    
    def$c.vm[i] = def$vm[i] *  ((1 / (ng-1)) * (((ng * 2)/3) - (ssgm/ng)))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] =  def$c.vq[i] - def$vg[i]
    def$diff[i] =  def$vq[i] - def$vg[i]
    
  }
  
  df2 = rbind(df2,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

save(list="df2",file="~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_5R_10miss_ASV.Rdata")

load("~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_5R_10miss_ASV.Rdata")
#df2 = df2[which(df2[,20] < 25 & df2[,20] > -25),]


(b23 = ggplot(data=as.data.frame(df2)) +
    geom_point(aes(x = df2[,4], y = df2[,14]),alpha=0.5,color="darkred") +
    geom_point(aes(x = df2[,4], y = df2[,19]),alpha=0.5,color="darkblue") +
    geom_abline(slope = 1,intercept = 0) +
    ggtitle(bquote("One Locus - 10% Missing")) + 
    xlab(bquote(H[M]^2)) + ylab("")+
    theme_bw() +   theme(plot.title = element_text(size = 10),
                         axis.title = element_text(size = 10),
                         axis.text = element_text(size=10))+
    xlim(c(0,1)) + ylim(c(0,2)))



###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 45 genotypes; 4 reps each (180)
#effect:: keeps same effects for the simulation

df3 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,900)))
  b = as.factor(sort(rep(0:719,5)))  
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:900),]
  mm = mm[,-c(1,5:183)]
  v = sample(sqrt(c(seq(1,100,by=1), seq(101, 1500,by=5))),3,replace=T)
  
  vL = v[1]
  vG = v[2]
  vE = 5
  
  u.L = matrix((c(rnorm(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm(540,0,vG)),nrow=540)
  
  L = as.factor(sort(rep(0:2,900)))
  G = as.factor(sort(rep(0:539,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:543)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm(2700,0,vE),nrow=2700)
    
    y2 = rowSums(y) + e
    
    x = sample(c(1:1800), 1800 * 0.667, replace = F)
    
    #expectations
    def$vQ[i] = var(rowSums(y[x,]))
    def$h.exp[i] = (var(y[x,1]) + var(y[x,2])) / (var(y[x,1]) + var(y[x,2]) + var(e[x]))
    def$hm.exp[i] = (var(y[x,1])) / (var(y[x,1]) + var(y[x,2]) + var(e[x]))
    def$vL[i] = var(y[x,1])
    def$vG[i] = var(y[x,2])
    def$vE[i] = var(e[x])
    def$vP[i] = var(y2[x])
    
    data = data.frame(Y=y2, L=L, G=G)
    
    data = data[x,]
    #random effect model
    model = lmer(Y ~(1|L/G),data,verbose=F,
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model <- update(model,start = getME(model,c("theta","fixef")))
    model <- update(model,start = getME(model,c("theta","fixef")))
    
    
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol)))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    ssgm = sum((c(length(unique(data$G[which(data$L == 0)])),length(unique(data$G[which(data$L == 1)])),length(unique(data$G[which(data$L == 2)]))) - mean(c(length(unique(data$G[which(data$L == 0)])),length(unique(data$G[which(data$L == 1)])),length(unique(data$G[which(data$L == 2)])))))^2)
    ng = length(unique(data$G))
    
    def$c.vm[i] = def$vm[i] *  ((1 / (ng-1)) * (((ng * 2)/3) - (ssgm/ng)))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] =  def$c.vq[i] - def$vg[i]
    def$diff[i] =  def$vq[i] - def$vg[i]
    
  }
  
  df3 = rbind(df3,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

save(list="df3",file="~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_5R_33miss_ASV.Rdata")

load("~/Box/Public/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_5R_33miss_ASV.Rdata")
#df3 = df3[which(df3[,20] < 25 & df3[,20] > -25),]

(b33 = ggplot(data=as.data.frame(df3)) +
    geom_point(aes(x = df3[,4], y = df3[,14]),alpha=0.5,color="darkred") +
    geom_point(aes(x = df3[,4], y = df3[,19]),alpha=0.5,color="darkblue") +
    geom_abline(slope = 1,intercept = 0) +
    ggtitle(bquote("One Locus - 33% Missing")) + 
    xlab(bquote(H[M]^2)) + ylab("")+
    theme_bw() +   theme(plot.title = element_text(size = 10),
                         axis.title = element_text(size = 10),
                         axis.text = element_text(size=10))+
    xlim(c(0,1)) + ylim(c(0,2)))

jpeg(filename = "Summary.jpg",width = 24, height=16,units="cm",pointsize = 10,quality=100,res=300)
(b1 + b5 + b8) / (b13 + b23 + b33) + plot_annotation(tag_levels = "A",tag_suffix = ".")
dev.off()
