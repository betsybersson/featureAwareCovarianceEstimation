set.seed(123)
###########################
source("helpers.R")
source("functions.R")
source("Schiavon_SISGaussian.R")
library(parallel)
library(dplyr)
library(readxl)
library(RColorBrewer)
###########################

###########################
## parameters to create fake data
## chemical class info
chemnames = c("EHDPP","TCEP","TCPP","TDCPP","TPP","2IPPDPP", 
              "24DIPPDPP","4tBPDPP","B4TBPP","DEP","DiBP","BBP",   
              "DBP","DEHP","DINP","TOTM","BPA","EPB",
              "MPB","PPB","TCS")
chemclass = c(rep("OPE",9), rep("phthalate",7),rep("phenol",5))
method.names = c("WB","Dust")

g = length(unique(chemclass)) # number of education groups
p1 = 21 # number of chemicals
p2 = length(method.names) # number of sources

n = 73 # sample size
p = p1*p2# feature dimension

# block covariance matrix for 2 groups
Sigma = matrix(.2,ncol = p, nrow = p)
Sigma[1:p1,1:p1] = eye(p1)*(1-.8) + .8
Sigma[(p1+1):p,(p1+1):p] = eye(p1)*(1-.5) + .5


Y = rmatnorm(matrix(0,ncol = p,nrow=n),V =Sigma)
Y = scale(Y)

###########################

###########################
## meta covariates
meta.covs = read_excel("VPandHPVforBetsy.xlsx") |> 
  select(3,4) 
colnames(meta.covs) = c("VP","prodVol") 
meta.covs = meta.covs |> 
  mutate(prodVol = as.numeric(as.factor(prodVol))-1,
         VP.pascal = VP * 133.322) |> select(-VP) |> 
  as.matrix() 

X = cbind(getMatDesignMat(p1,p2),
          rbind(getCatDesignMat(chemclass),
                getCatDesignMat(chemclass)),  # once for each source
          rbind(meta.covs,meta.covs))  # once for each source
which.cov.group = c(rep(1,p1),rep(2,p2),
                    rep(3,length(unique(chemclass))),
                    4,5)
###########################

###########################
# run GS
tic = Sys.time()
out.face.cusp  = FACE_cusp_GS(Y,X,
                            k = p,
                            S = 20000,
                            burnin = 10000,
                            my.seed = 123,
                            alpha = floor(p/3),
                            a.theta = 1/2, b.theta = 1/2,
                            which.cov.group = which.cov.group)
toc.face.cusp = Sys.time() - tic
face.cusp = qr.solve(matrix(colMeans(out.face.cusp$cov.inv),ncol = p)) 
###########################

###########################
## for plotting

# cor mat for all k
cor.hat.stein.collect = cov2cor(face.cusp)

# get color limits so color scale is uniform across plots
univ.col.lims = c(min(0,min(cor.hat.stein.collect)),
                  1)
col.brewer.pal = c(brewer.pal(9,"GnBu")[2],"#FFFFFF", brewer.pal(9,"YlOrRd")) #revs


ope.last.ind = max(which(chemclass=="OPE"))
phthal.last.ind = max(which(chemclass=="phthalate"))

plot.tesie.lines = function(){
  # line break representing the two sources 
  abline(h=p/2+.5,col="black", lwd = 3) 
  abline(v=p/2+.5,col="black", lwd = 3)
  # line break representing OPEs
  abline(v = ope.last.ind + .5)
  abline(v = p/2 + ope.last.ind + .5)
  abline(h = p - (ope.last.ind) + .5)
  abline(h = p - (p/2 + ope.last.ind) + .5)
  # line break representing phthalate
  abline(v = phthal.last.ind + .5)
  abline(v = p/2 + phthal.last.ind + .5)
  abline(h = p - (phthal.last.ind) + .5)
  abline(h = p - (p/2 + phthal.last.ind) + .5)
  # label
  # y-axis labels
  axis(side=2,at=c(p - ope.last.ind/2,(p - ope.last.ind)/2),
       labels=rep("OPE",2),tick=FALSE,
       las=3, cex.axis=1) 
  axis(side=2,at=c(p - ope.last.ind/2 - (phthal.last.ind - ope.last.ind),
                   p - ope.last.ind/2 - (phthal.last.ind - ope.last.ind) - p/2),
       labels=rep("Phthalate",2),tick=FALSE,
       las=3, cex.axis=1) 
  axis(side=2,at=c(p - phthal.last.ind - sum(chemclass=="phenol")/2 + .5,
                   sum(chemclass=="phenol")/2 + .5),
       labels=rep("Phenol",2),tick=FALSE,
       las=3, cex.axis=1) 
  axis(2, at = c(p/2/2,p-p/2/2),
       labels=rev(method.names),
       tick=FALSE,las=3,
       cex.axis=1.8,line=1.3) 
  # x-axis labels (more straightforward)
  axis(side=1,at=c(ope.last.ind/2,(p + ope.last.ind)/2),
       labels=rep("OPE",2),tick=FALSE,
       las=1, cex.axis=1) 
  axis(side=1,at=c(ope.last.ind + sum(chemclass=="phthalate")/2 + .5,
                   (p + ope.last.ind + phthal.last.ind)/2 + .5),
       labels=rep("Phthalate",2),tick=FALSE,
       las=1, cex.axis=1) 
  axis(side=1,at=c(p/2 - sum(chemclass=="phenol")/2,
                   p - sum(chemclass=="phenol")/2),
       labels=rep("Phenol",2),tick=FALSE,
       las=1, cex.axis=1) 
  axis(1, at = c(p/2/2,p-p/2/2),
       labels=method.names,
       tick=FALSE,las=1,
       cex.axis=1.8,line=1.3) 
}

# plot
image(x = 1:p, y = 1:p, z = cor.hat.stein.collect[,p:1],
      col = col.brewer.pal,
      xaxt="n",yaxt="n",xlab="",ylab="",
      cex = 3, cex.main = 1.4,
      zlim = univ.col.lims,
      main = paste0("Posterior Correlation (FACE)"))
plot.tesie.lines()
cor.legend.inds = c(1,3,6,9,11)
legend(
  x = "bottomleft", inset  = 0.052, # position
  legend = round(seq(univ.col.lims[1], 
                     univ.col.lims[2], 
                     length.out = length(col.brewer.pal)), 2)[cor.legend.inds], # only include 5 valus
  fill   = col.brewer.pal[cor.legend.inds],  # representative colors
  cex    = 1, 
  bg = "white")


####################################
library(MASS)

alpha = 0.05
####################################

####################################
###### UQ MLE

### analyze tesie sample covariance
sample.cov = cov(Y)

## get some pvalues
MLE.pvals = matrix(0,ncol = p, nrow = p)
for ( j in 1:p){ 
  for ( k in 1:p){
    MLE.pvals[j,k] = cor.test(Y[,j],Y[,k], 
                              alternative = "two.sided")$p.value
  }
}
MLE.pvals = MLE.pvals[upper.tri(MLE.pvals)]
test.incl0.BY = p.adjust(MLE.pvals, method = 'BY', 
                         n = length(MLE.pvals))
test.incl0.BY = ifelse(test.incl0.BY > alpha,1,0)


Sig.Inclu0 = matrix(-1,nrow=p,ncol=p)
Sig.Inclu0[upper.tri(Sig.Inclu0)] = test.incl0.BY
Sig.Inclu0[lower.tri(Sig.Inclu0)] = t(Sig.Inclu0)[lower.tri(Sig.Inclu0)] 
Sig.Inclu0[upper.tri(Sig.Inclu0)] = -1

image(x = 1:p, y = 1:p, z = Sig.Inclu0[,p:1],
      xaxt="n",yaxt="n",xlab="",ylab="",
      col = c("Gray","White","Blue"),
      cex = 3, cex.main = 1.4,
      main = "Sample Correlation: Element-wise test of 0 corr.")
plot.tesie.lines()
legend(3,12,legend = c("Reject","Don't Reject"),
       fill=c("White","Blue"),bty="y")
####################################

# ####################################
# ###### UQ SIS
# 
# which.uppertri = matrix(1:(p^2),ncol = p)[upper.tri(eye(p))]
# uppertri.mcmc = t(apply(schiavon_out$covSamples,3,function(j)c(j)[which.uppertri]))
# 
# SIS.CIs = apply(uppertri.mcmc,2,function(j)quantile(j,c(alpha/2,1-alpha/2)))
# test.incl0 = apply(SIS.CIs,2,function(j)(0>=j[1])&(0<=j[2])) |> as.numeric()
# test.incl0 |> sum()
# which(test.incl0==0)
# 
# 
# Sig.Inclu0 = matrix(-1,nrow=p,ncol=p)
# Sig.Inclu0[upper.tri(Sig.Inclu0)] = test.incl0
# Sig.Inclu0[lower.tri(Sig.Inclu0)] = t(Sig.Inclu0)[lower.tri(Sig.Inclu0)] 
# Sig.Inclu0[upper.tri(Sig.Inclu0)] = -1
# 
# image(x = 1:p, y = 1:p, z = Sig.Inclu0[,p:1],
#       xaxt="n",yaxt="n",xlab="",ylab="",
#       col = c("Gray","White","Blue"),
#       cex = 3, cex.main = 1.4,
#       main = "SIS: Element-wise test of 0 corr.")
# plot.tesie.lines()
# legend(3,12,legend = c("Reject","Don't Reject"),
#        fill=c("White","Blue"),bty="y")
# 
# ####################################


####################################
###### UQ FACE

which.uppertri = matrix(1:(p^2),ncol = p)[upper.tri(eye(p))]
uppertri.mcmc = out.face.cusp$cov[,which.uppertri]

FACE.CIs = apply(uppertri.mcmc,2,function(j)quantile(j,c(alpha/2,1-alpha/2)))
test.incl0 = apply(FACE.CIs,2,function(j)(0>=j[1])&(0<=j[2])) |> as.numeric()

Sig.Inclu0 = matrix(-1,nrow=p,ncol=p)
Sig.Inclu0[upper.tri(Sig.Inclu0)] = test.incl0
Sig.Inclu0[lower.tri(Sig.Inclu0)] = t(Sig.Inclu0)[lower.tri(Sig.Inclu0)] 
Sig.Inclu0[upper.tri(Sig.Inclu0)] = -1

image(x = 1:p, y = 1:p, z = Sig.Inclu0[,p:1],
      xaxt="n",yaxt="n",xlab="",ylab="",
      col = c("Gray","White","Blue"),
      cex = 3, cex.main = 1.4,
      main = "FACE: Element-wise test of 0 corr.")
plot.tesie.lines()
legend(3,12,legend = c("Reject","Don't Reject"),
       fill=c("White","Blue"),bty="y")

####################################

