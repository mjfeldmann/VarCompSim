rm(list=ls())
library(lme4);library(ggplot2);library(cowplot); library(gridExtra); library(grid); library(tidybayes)
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

iter1 = 1000
iter2 = 1

###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 90 genotypes; 2 reps each (180)
#effect:: keeps same effects for the simulation
df2 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,600)))
  b= as.factor(sort(rep(0:1199,2)))
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:600),]
  mm = mm[,-c(1,5:303)]
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(900,0,vG)),nrow=900)
  
  L = as.factor(sort(rep(0:2,600)))
  G = as.factor(sort(rep(0:899,2)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:903)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm2(1800,0,vE),nrow=1800)
    
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
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model <- update(model,start = getME(model,c("theta","fixef")))
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    def$c.vm[i] = def$vm[i] * ((2 * 300) / (899 ))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] = ( def$c.vq[i] - def$vg[i]) / def$vg[i]
    def$diff[i] = ( def$vq[i] - def$vg[i] )/ def$vg[i] 
    
  }
  
  df2 = rbind(df2,(def))
  setTxtProgressBar(pb,j)
}
close(pb)
colnames(df2)


#save(list="df2",file="Sim_1L_2R.Rdata")

#load("~/Box/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_2R.Rdata")


###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 90 genotypes; 2 reps each (180)
#effect:: keeps same effects for the simulation
df5 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,1500)))
  b= as.factor(sort(rep(0:1199,5)))
  mm = model.matrix(~a + b)
  
  dim(mm)
  
  mm = mm[-c(1:1500),]
  mm = mm[,-c(1)]
  mm = mm[,-c(4:302)]
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(900,0,vG)),nrow=900)
  
  L = as.factor(sort(rep(0:2,1500)))
  G = as.factor(sort(rep(0:899,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:903)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm2(4500,0,vE),nrow=4500)
    
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
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model <- update(model,start = getME(model,c("theta","fixef")))
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    def$c.vm[i] = def$vm[i] * ((2 * 300) / (899 ))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] = ( def$c.vq[i] - def$vg[i]) / def$vg[i]
    def$diff[i] = ( def$vq[i] - def$vg[i] )/ def$vg[i] 
    
  }
  
  df5 = rbind(df5,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

#save(list="df5",file="Sim_1L_2R.Rdata")

#load("~/Box/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_2R.Rdata")


###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 18 genotypes; 10 reps each (180)
#effect:: keeps same effects for the simulation

df10 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,3000)))
  b= as.factor(sort(rep(0:1199,10)))
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:3000),]
  mm = mm[,-c(1)]
  mm = mm[,-c(4:302)]
  
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(900,0,vG)),nrow=900)
  
  L = as.factor(sort(rep(0:2,3000)))
  G = as.factor(sort(rep(0:899,10)))

  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:903)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm2(9000,0,vE),nrow=9000)
    
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
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model <- update(model,start = getME(model,c("theta","fixef")))
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    def$c.vm[i] = def$vm[i] * ((2 * 300) / (899 ))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] =  (def$c.vq[i] - def$vg[i]) / def$vg[i] 
    def$diff[i] =  (def$vq[i] - def$vg[i]) / def$vg[i]
    
  }
  
  df10 = rbind(df10,(def))
  setTxtProgressBar(pb,j)
}
close(pb)


#save(list="df10",file="Sim_1L_10R.Rdata")

#load("~/Box/Desktop/Reading/Strawberry_Research/2_Genetic_Variance_Estimate/R_Data/Sim_1L_10R.Rdata")

###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 18 genotypes; 10 reps each (180)
#effect:: keeps same effects for the simulation

df20 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,6000)))
  b= as.factor(sort(rep(0:1199,20)))
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:6000),]
  mm = mm[,-c(1)]
  mm = mm[,-c(4:302)]
  
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(900,0,vG)),nrow=900)
  
  L = as.factor(sort(rep(0:2,6000)))
  G = as.factor(sort(rep(0:899,20)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:903)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm2(18000,0,vE),nrow=18000)
    
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
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model <- update(model,start = getME(model,c("theta","fixef")))
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    def$c.vm[i] = def$vm[i] * ((2 * 300) / (899 ))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] =  (def$c.vq[i] - def$vg[i]) / def$vg[i] 
    def$diff[i] =  (def$vq[i] - def$vg[i]) / def$vg[i]
    
  }
  
  df20 = rbind(df20,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 18 genotypes; 10 reps each (180)
#effect:: keeps same effects for the simulation
50/20
6000*2.5
df50 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,15000)))
  b= as.factor(sort(rep(0:1199,50)))
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:15000),]
  mm = mm[,-c(1)]
  mm = mm[,-c(4:302)]
  
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(900,0,vG)),nrow=900)
  
  L = as.factor(sort(rep(0:2,15000)))
  G = as.factor(sort(rep(0:899,50)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:903)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm2(45000,0,vE),nrow=45000)
    
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
                 control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model <- update(model,start = getME(model,c("theta","fixef")))
    model.G = lmer(Y ~ (1|G),data,verbose=F,
                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol=formals(isSingular)$tol))); model.G <- update(model.G,start = getME(model.G,c("theta","fixef")))
    
    def$vgm[i] = VarCorr(model)[[1]][1]
    def$vm[i] = VarCorr(model)[[2]][1]
    def$ve.l[i] = attr(VarCorr(model), "sc")**2
    def$vg[i] = VarCorr(model.G)[[1]][1]
    
    def$h.obs[i] = (def$vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    
    def$hm.obs[i] = (def$vm[i]) / (def$vg[i] + def$ve[i])
    
    def$c.vm[i] = def$vm[i] * ((2 * 300) / (899 ))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] =  (def$c.vq[i] - def$vg[i]) / def$vg[i] 
    def$diff[i] =  (def$vq[i] - def$vg[i]) / def$vg[i]
    
  }
  
  df50 = rbind(df50,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

df2$rg = factor("2")
df5$rg = factor("5")
df10$rg = factor("10")
df20$rg = factor("20")
df50$rg = factor("50")

df.rg = rbind(df2, df5, df10, df20, df50)

save(list="df.rg",file="bias_rg.RData")


r1 = ggplot(data=df.rg) + 
  geom_boxplot(aes(x=rg, y = diff),fill="gray",outlier.shape = NA,width=0.75) + geom_jitter(aes(x=rg, y = diff),alpha=0.35,color="darkred",width = 0.2) +
  #geom_boxplot(aes(x=rg, y = diff.c),fill="gray",outlier.shape = NA,width=0.75) + geom_jitter(aes(x=rg, y = diff.c),alpha=0.35,color="darkblue",width = 0.2) +
  xlab(bquote(r[G])) + ylab(bquote(Bias[G]^AMV)) + ylim(c(-.015,.6)) +
  geom_hline(yintercept = median(df.rg$diff),lty=2,color="darkred") +
  #geom_hline(yintercept = median(df.rg$diff.c),lty=2,color="darkblue") +
  theme_bw() +  theme(plot.title = element_text(size = 10),
                      axis.title = element_text(size = 10),
                      axis.text = element_text(size=10))

r2 = ggplot(data=df.rg) + 
  #geom_boxplot(aes(x=rg, y = diff),fill="gray",outlier.shape = NA,width=0.75) + geom_jitter(aes(x=rg, y = diff),alpha=0.35,color="darkred",width = 0.2) +
  geom_boxplot(aes(x=rg, y = diff.c),fill="gray",outlier.shape = NA,width=0.75) + geom_jitter(aes(x=rg, y = diff.c),alpha=0.35,color="darkblue",width = 0.2) +
  xlab(bquote(r[G])) + ylab(bquote(Bias[G]^ASV)) + ylim(c(-.005,.005)) +
  #geom_hline(yintercept = median(df.rg$diff),lty=2,color="darkred") +
  geom_hline(yintercept = median(df.rg$diff.c),lty=2,color="darkblue") +
  theme_bw() +  theme(plot.title = element_text(size = 10),
                      axis.title = element_text(size = 10),
                      axis.text = element_text(size=10))

gA <- ggplotGrob(r1)
gB <- ggplotGrob(r2)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
jpeg("BiasG_rg.jpg",width=4, height=8,units="in", pointsize = 10,quality=100,res=300)
grid.arrange(gA, gB, ncol=1)
dev.off()

