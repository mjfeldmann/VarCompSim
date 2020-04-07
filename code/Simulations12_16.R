rm(list=ls())
library(lme4);library(ggplot2);library(cowplot); library(gridExtra); library(grid)
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

iter1 = 1000
iter2 = 1

###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 90 genotypes; 2 reps each (180)
#effect:: keeps same effects for the simulation
df900 = c()
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
  
  df900 = rbind(df900,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 90 genotypes; 2 reps each (180)
#effect:: keeps same effects for the simulation
df450 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,750)))
  b= as.factor(sort(rep(0:599,5)))
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:750),]
  mm = mm[,-c(1)]
  mm = mm[,-c(4:152)]
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(450,0,vG)),nrow=450)
  
  L = as.factor(sort(rep(0:2,750)))
  G = as.factor(sort(rep(0:449,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:453)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm2(2250,0,vE),nrow=2250)
    
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
    
    def$c.vm[i] = def$vm[i] * ((2 * 150) / (449 ))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] = ( def$c.vq[i] - def$vg[i]) / def$vg[i]
    def$diff[i] = ( def$vq[i] - def$vg[i] )/ def$vg[i] 
    
  }
  
  df450 = rbind(df450,(def))
  setTxtProgressBar(pb,j)
}
close(pb)


###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 90 genotypes; 2 reps each (180)
#effect:: keeps same effects for the simulation
df1800 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,3000)))
  b= as.factor(sort(rep(0:2399,5)))
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:3000),]
  mm = mm[,-c(1)]
  mm = mm[,-c(4:602)]
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(1800,0,vG)),nrow=1800)
  
  L = as.factor(sort(rep(0:2,3000)))
  G = as.factor(sort(rep(0:1799,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:1803)] %*% u.G
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
    
    def$c.vm[i] = def$vm[i] * ((2 * 600) / (1799 ))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] = ( def$c.vq[i] - def$vg[i]) / def$vg[i]
    def$diff[i] = ( def$vq[i] - def$vg[i] )/ def$vg[i] 
    
  }
  
  df1800 = rbind(df1800,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 90 genotypes; 2 reps each (180)
#effect:: keeps same effects for the simulation
df3600 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,6000)))
  b= as.factor(sort(rep(0:4799,5)))
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:6000),]
  mm = mm[,-c(1)]
  mm = mm[,-c(4:1202)]
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(3600,0,vG)),nrow=3600)
  
  L = as.factor(sort(rep(0:2,6000)))
  G = as.factor(sort(rep(0:3599,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:3603)] %*% u.G
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
    
    def$c.vm[i] = def$vm[i] * ((2 * 1200) / (3599 ))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] = ( def$c.vq[i] - def$vg[i]) / def$vg[i]
    def$diff[i] = ( def$vq[i] - def$vg[i] )/ def$vg[i] 
    
  }
  
  df3600 = rbind(df3600,(def))
  setTxtProgressBar(pb,j)
}
close(pb)

###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 90 genotypes; 2 reps each (180)
#effect:: keeps same effects for the simulation

df7200 = c()
pb <- txtProgressBar(min=0,max=iter1,style=3)
for(j in 1:iter1){
  #design
  a = as.factor(sort(rep(0:3,12000)))
  b= as.factor(sort(rep(0:9599,5)))
  mm = model.matrix(~a + b)
  
  mm = mm[-c(1:12000),]
  mm = mm[,-c(1)]
  mm = mm[,-c(4:2402)]
  #v = sample(sqrt(seq(1,100,by=1)),3,replace=T)
  
  vL = sqrt(75)
  vG = 5
  vE = 5
  
  u.L = matrix((c(rnorm2(3,0,vL))),nrow=3)
  u.G = matrix(c(rnorm2(7200,0,vG)),nrow=7200)
  
  L = as.factor(sort(rep(0:2,12000)))
  G = as.factor(sort(rep(0:7199,5)))
  
  #Generate data
  y1 = mm[,c(1:3)] %*% u.L
  y2 = mm[,c(4:7203)] %*% u.G
  y = cbind(y1,y2)
  
  def = data.frame(iter=c(1:iter2))
  
  for(i in 1:iter2){
    #add Error
    e = matrix(rnorm2(36000,0,vE),nrow=36000)
    
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
    
    def$c.vm[i] = def$vm[i] * ((2 * 2400) / (7199))
    def$vq[i] = def$vm[i] + def$vgm[i]
    def$c.vq[i] = def$c.vm[i] + def$vgm[i]
    
    def$h.c[i] = (def$c.vm[i] + def$vgm[i]) / (def$vg[i] + def$ve[i])
    def$hm.c[i] = (def$c.vm[i]) / (def$vg[i] + def$ve[i])
    
    def$diff.c[i] = ( def$c.vq[i] - def$vg[i]) / def$vg[i]
    def$diff[i] = ( def$vq[i] - def$vg[i] )/ def$vg[i] 
    
  }
  
  df7200 = rbind(df7200,(def))
  setTxtProgressBar(pb,j)
}
close(pb)


df450$ng = factor("450")
df900$ng = factor("900")
df1800$ng = factor("1800")
df3600$ng = factor("3600")
df7200$ng = factor("7200")

df.ng = rbind(df450, df900, df1800, df3600, df7200)

save(list="df.ng",file="bias_ng.RData")

load("~/bias_ng.RData")
load("~/bias_rg.RData")

library(ggplot2); library(patchwork); library(gridExtra)

###########################
###########################
###########################
(n1 = ggplot(data=df.ng) + 
   geom_jitter(aes(x=ng, y = (vm/vL)-1),alpha=0.35,color="darkred",width = 0.2) +
   #geom_boxplot(aes(x=ng, y = diff.c),fill="gray",outlier.shape = NA,width=0.75) + geom_jitter(aes(x=ng, y = diff.c),alpha=0.35,color="darkblue",width = 0.2) +
   xlab(bquote(italic(n)[G])) + ylab(bquote(frac(hat(theta)[M]^AMV - s[M]^2, s[M]^2))) + ylim(c(-1.5,2.5)) +
   theme_bw() +  theme(plot.title = element_text(size = 10),
                       axis.title = element_text(size = 10),
                       axis.text = element_text(size=10)))

(n2 = ggplot(data=df.ng) + 
    geom_jitter(aes(x=ng, y = (c.vm/vL)-1),alpha=0.35,color="darkblue",width = 0.2) +
    xlab(bquote(italic(n)[G])) + ylab(bquote(frac(hat(theta)[M]^ASV - s[M]^2, s[M]^2))) + ylim(c(-1.5,2.5)) +
    theme_bw() +  theme(plot.title = element_text(size = 10),
                        axis.title = element_text(size = 10),
                        axis.text = element_text(size=10)))


jpeg("BiasG_ng.jpg",width=8, height=4,units="in", pointsize = 10,quality=100,res=300)
n1 + n2 + plot_annotation(tag_levels = "A",tag_suffix = ".")
dev.off()

(r1 = ggplot(data=df.rg) + 
    geom_jitter(aes(x=rg, y = (vm/vL)-1),alpha=0.35,color="darkred",width = 0.2) +
    #geom_boxplot(aes(x=ng, y = diff.c),fill="gray",outlier.shape = NA,width=0.75) + geom_jitter(aes(x=ng, y = diff.c),alpha=0.35,color="darkblue",width = 0.2) +
    xlab(bquote(italic(r)[G])) + ylab(bquote(frac(hat(theta)[M]^AMV - s[M]^2, s[M]^2))) + 
    ylim(c(-1.5,2.5)) +
    #geom_hline(yintercept = median(df.ng$diff.c),lty=2,color="darkblue") +
    theme_bw() +  theme(plot.title = element_text(size = 10),
                        axis.title = element_text(size = 10),
                        axis.text = element_text(size=10)))

(r2 = ggplot(data=df.rg) + 
    geom_jitter(aes(x=rg, y = (c.vm/vL)-1),alpha=0.35,color="darkblue",width = 0.2) +
    xlab(bquote(italic(r)[G])) + ylab(bquote(frac(hat(theta)[M]^ASV - s[M]^2, s[M]^2))) + 
    ylim(c(-1.5,2.5)) +
    theme_bw() +  theme(plot.title = element_text(size = 10),
                        axis.title = element_text(size = 10),
                        axis.text = element_text(size=10)))

jpeg("BiasG_rg.jpg",width=8, height=4,units="in", pointsize = 10,quality=100,res=300)
r1 + r2 + plot_annotation(tag_levels = "A",tag_suffix = ".")
dev.off()

median(df.ng$c.vm[which(df.ng$ng=="450")] /df.ng$vL[which(df.ng$ng=="450")])
median(df.ng$c.vm[which(df.ng$ng=="900")] /df.ng$vL[which(df.ng$ng=="900")])
median(df.ng$c.vm[which(df.ng$ng=="1800")] /df.ng$vL[which(df.ng$ng=="1800")])
median(df.ng$c.vm[which(df.ng$ng=="3600")] /df.ng$vL[which(df.ng$ng=="3600")])
median(df.ng$c.vm[which(df.ng$ng=="7200")] /df.ng$vL[which(df.ng$ng=="7200")])

median(df.ng$vm[which(df.ng$ng=="450")] /df.ng$vL[which(df.ng$ng=="450")])
median(df.ng$vm[which(df.ng$ng=="900")] /df.ng$vL[which(df.ng$ng=="900")])
median(df.ng$vm[which(df.ng$ng=="1800")] /df.ng$vL[which(df.ng$ng=="1800")])
median(df.ng$vm[which(df.ng$ng=="3600")] /df.ng$vL[which(df.ng$ng=="3600")])
median(df.ng$vm[which(df.ng$ng=="7200")] /df.ng$vL[which(df.ng$ng=="7200")])

median(df.ng$vm/df.ng$vL)
median(df.rg$vm/df.rg$vL)
