set.seed(123)
###########################
source("helpers.R")
source("functions.R")
source("Schiavon_SISGaussian.R")
library(parallel)
library(dplyr)
library(readxl)
# Plot
library(reshape2)
library(ggplot2)
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
Y.orig = Y

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
## MCMC parameters
mcmc.s = 20
mcmc.burnin = 10

###########################


###########################
# run prediction loop
n.lod.set =  c(seq(from=5,to=15,by=5), 30, 45, 60)
method.names = c("MR","MR.I","S.IMP","SIS")
grand.mse = matrix(NA,nrow = length(n.lod.set), 
                   ncol = length(method.names)) 
colnames(grand.mse) = method.names
rownames(grand.mse) = paste0("n",n.lod.set)
for (n.ind in 1:length(n.lod.set)){
  
  ## reset Y
  Y = Y.orig
  
  ## get lod values (n.lod from each column)
  n.lod = n.lod.set[n.ind] ## how many of the smallest values do we want to remove
  ind.lod = apply(Y,2,function(j) order(j)[1:n.lod])
  lod = apply(Y,2,function(j) mean(j[order(j)[n.lod:(n.lod+1)]])) # avg between last missing value and next smallest value
  
  ## change ind.lod values to the lod
  Y.na = Y
  for ( j in 1:p ){
    Y[ind.lod[,j],j] = lod[j]/sqrt(2)
    Y.na[ind.lod[,j],j] = NA
  }
  # get binary lod matrix
  which.lod = matrix(0,nrow = n,ncol=p)
  for (j in 1:p){
    which.lod[ind.lod[,j],j] = 1
  }
  row.ind.lod = sort(unique(c(ind.lod)))
  ## extract true test data
  test.out = (Y.orig*which.lod)[row.ind.lod,] |> c()
  
  just.imputedvals = (Y*which.lod)[row.ind.lod,] |> c()
  grand.mse[n.ind,"S.IMP"] = loss_sq(test.out,just.imputedvals)
  
  
  ## face.i mse
  face.ylod  = FACE_cusp_GS(Y,X=NA,
                          k = p,
                          S = mcmc.s,
                          burnin = mcmc.burnin,
                          my.seed = 123,
                          alpha = floor(p/3),
                          a.theta = 1/2, b.theta = 1/2,
                          which.lod = which.lod,
                          row.ind.lod = row.ind.lod,
                          lod = lod)$y.lod
  
  # mse
  grand.mse[n.ind,"MR.I"] = loss_sq(test.out,colMeans(face.ylod))
  # print update to screen
  print(paste0("Finished running FACE.I for n.lod = ", n.lod))  
  print(paste0("MSE = ", grand.mse[n.ind,"MR.I"]))  
  print(paste0("---------------------------------"))
  
  
  ## face.chemclass mse
  face.ylod  = FACE_cusp_GS(Y,X,
                           k = p,
                           S = mcmc.s,
                           burnin = mcmc.burnin,
                           my.seed = 123,
                           alpha = floor(p/3),
                           a.theta = 1/2, b.theta = 1/2,
                           which.lod = which.lod,
                           row.ind.lod = row.ind.lod,
                           lod = lod,
                           which.cov.group = which.cov.group)$y.lod

  # mse
  grand.mse[n.ind,"MR"] = loss_sq(test.out,colMeans(face.ylod))
  # print update to screen
  print(paste0("Finished running FACE for n.lod = ", n.lod))  
  print(paste0("MSE = ", grand.mse[n.ind,"MR"]))  
  print(paste0("---------------------------------"))
  
  
  ## sis mse
  schiavon_out = Mcmc_SIS(Y,X,
                          as = 1, bs = 1, alpha = round(p/3),
                          a_theta = 2, b_theta = 2,
                          nrun = mcmc.s, burn = mcmc.burnin,
                          thin = 1,
                          start_adapt = mcmc.s, kmax = p + 1, # don't adapt
                          my_seed =  123,
                          output = c("ylodMean"),
                          which.lod = which.lod,
                          row.ind.lod = row.ind.lod,
                          lod = lod)$ylodMean
  
  # mse
  grand.mse[n.ind,"SIS"] = loss_sq(test.out,schiavon_out)
  # print update to screen
  print(paste0("Finished running SIS for n.lod = ", n.lod))
  print(paste0("MSE = ", grand.mse[n.ind,"SIS"]))
  print(paste0("---------------------------------"))
  
}

grand.mse = rbind(rep(0,ncol(grand.mse)),grand.mse)
rownames(grand.mse)[1] = "n0"

mse.melt = melt(grand.mse) |> 
  mutate(Var1 = as.character(Var1)) |> 
  mutate(n.miss = as.numeric(substr(Var1, 2, nchar(Var1)))) |> 
  mutate(value = (sqrt(value))) |> 
  mutate(pct.detect = (n-n.miss)/n*100)


mse.melt |> 
  ggplot(aes(x = pct.detect, y = value, color = Var2)) +
  geom_point() + geom_line() +
  labs(x = "% Detected (per variable)",#expression(n[test]~" (per column)"),
       y = "RMSE",
       color = "Method",
       title = "Imputation RMSE Comparison") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_blank(), # Ensure white background
    axis.line = element_line(color = "black"), # Add axis lines
    strip.text = element_text(face = "bold"), # Bold facet labels
    panel.spacing = unit(.05, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(color = "black", size = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5,vjust=1.8),
    text = element_text(family = "Times",size = 18)) +
  scale_color_manual(values = c("MR" = "red", "S.IMP" = "blue",
                                "SIS" = "purple", "MR.I" = "orange"))





