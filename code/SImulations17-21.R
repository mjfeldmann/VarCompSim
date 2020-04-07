rm(list=ls())
library(lme4);library(ggplot2);library(cowplot); library(gridExtra); library(grid);library(patchwork)
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

iter1 = 1000
iter2 = 1

###################################################################################
###################################################################################
###################################################################################
# 1 locus: 3 genotypes; 60 reps each (180)
# 45 genotypes; 4 reps each (180)
#effect:: keeps same effects for the simulation

df = c()
for(k in 1:5){
  l = c(2,5,sqrt(75),sqrt(225), sqrt(1450))
  for(j in 1:iter1){
    #design
    a = as.factor(sort(rep(0:3,750)))
    b= as.factor(sort(rep(0:599,5)))
    mm = model.matrix(~a + b)
    
    mm = mm[-c(1:750),]
    mm = mm[,-c(1)]
    mm = mm[,-c(4:152)]
    
    vL = l[k]
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
      print(paste0(k," - ", j))
    }
  
  df = rbind(df,(def))
}}

df$hm.exp2 = as.factor(round(df$hm.exp,2))

save(list="df",file="bias_hm.RData")
load("~/bias_hm.RData")
median((df$vm - df$vL)/df$vL)
(a1 = ggplot(data=as.data.frame(df)) +
  geom_jitter(aes(x = df$hm.exp2, y = (df$vm - df$vL)/df$vL),alpha=0.5,width=0.2,col="darkred") +
  xlab(bquote(italic(H)[M]^2)) + ylab(bquote(frac(hat(theta)[M]^AMV - s[M]^2, s[M]^2)))+
  theme_bw() +   theme(plot.title = element_text(size = 10),
                       axis.title = element_text(size = 10),
                       axis.text = element_text(size=10)) +
  ylim(c(-1.5,2.5)))
(b1 = ggplot(data=as.data.frame(df)) +
  geom_jitter(aes(x = df$hm.exp2, y = (df$c.vm - df$vL)/df$vL),alpha=0.5,width=0.2,col="darkblue") +
  xlab(bquote(italic(H)[M]^2)) +ylab(bquote(frac(hat(theta)[M]^ASV - s[M]^2, s[M]^2)))+
  theme_bw() +   theme(plot.title = element_text(size = 10),
                       axis.title = element_text(size = 10),
                       axis.text = element_text(size=10)) +
  ylim(c(-1.5,2.5)))
median((df$vm - df$vL)/df$vL)
median((df$c.vm - df$vL)/df$vL)

jpeg("BiasG_hm.jpg",width=8, height=4,units="in", pointsize = 10,quality=100,res=300)
a1 + b1 + plot_annotation(tag_levels = "A",tag_suffix = ".")
dev.off()
