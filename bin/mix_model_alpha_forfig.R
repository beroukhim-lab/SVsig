#######write optimizing alpha function in R#########
#takes in model1, model2, mfull binned events matrix and returns optimal alphas over three parameters#

library(nloptr)
library(cbinom)
library(data.table)
library(Matrix)
library(R.matlab)

#adj.base <- "/xchip/beroukhimlab/Kiran/adjancencies"
#SVsig.base <- "/xchip/beroukhimlab/Kiran/git/2dmodel/SVsig"
SVsig.base <- "/Users/shu/SVsig_labcopy"

#c.base <- "/xchip/beroukhimlab/Kiran/complex"

#read matlab with model1, model2, and mfull
#complex weighted events
mat <- readMat(file.path(SVsig.base, "debug_alpha.mat"))

#simple unweighted events
#mat <- readMat(file.path(SVsig.base, "debug_alpha_simple.mat"))


#turn annot.tiles into a logical matrix and then into a series of indices
annot.tiles <- lapply(1:4, function(i) { 
  m <- as(mat$annot.tiles[, ,i], "lMatrix")
  which(m == TRUE)
})

#double break join 
(model1 <- as(mat$model1, "CsparseMatrix"))
#break invasion
(model2 <- as(mat$model2, "CsparseMatrix"))


#read in mfull and pa
#mfull <- fread(file.path(adj.base, "20190426mfull.txt"))
#pa <- fread(file.path(adj.base, "20190426pa.txt"))

#turn into matrices
#mfull <- simplify2array(mfull)
#pa <- simplify2array(pa)

#turn into sparse matrices
#(mfs <- as(mfull, "CsparseMatrix"))
#(pas <- as(pa, "CsparseMatrix"))


#find non zero elements in mfs s 




#alpha <- 0.5

#####
##constrained optimization wrt BIC
##optimize parameters for short, long, and interchromosomal events all at the same time

#for simple events replicate Ofer's poisson approximation
eval_f2 <- function(alpha) {
  
  mix_model <- Matrix(0, nrow(mat$bins), nrow(mat$bins), sparse=TRUE)
  
  for (i in 1:4)  { 
    
    
    mix_model[annot.tiles[[i]]] <-  alpha[i]*model1[annot.tiles[[i]]] + (1- alpha[i])*model2[annot.tiles[[i]]]      
    
  }
  
  #normalize because multiplying by alphas means that probabilities no longer sum to 1
  mix_model <- mix_model/sum(mix_model)       
  
  idx_nnz <- intersect(which(mix_model > 0), which(mat$mfull > 0))
  
  nume = sum(mat$mfull)
  
  ll <- sum(log(dbinom(mat$mfull[idx_nnz], nume, mix_model[idx_nnz])))
  
  BIC = -2*ll + log(nume)*length(alpha)
  
  print(BIC)
  
}

eval_f1 <- function(alpha) { 
  
  mix_model <- Matrix(0, nrow(mat$bins), nrow(mat$bins), sparse=TRUE)
  
  for (i in 1:4)  { 
    
    
    mix_model[annot.tiles[[i]]] <-  alpha[i]*model1[annot.tiles[[i]]] + (1- alpha[i])*model2[annot.tiles[[i]]]      
    
  }
  
  #normalize because multiplying by alphas means that probabilities no longer sum to 1
  mix_model <- mix_model/sum(mix_model)       
  
  #idxz <- intersect(which(mix_model == 0), which(mat$mfull == 0))
  idxnnz <- which(mat$mfull > 0)
  #try this including all the zero tiles and see what happens to the alphas
  
  nume = sum(mat$mfull)
  
  #is this understanding the inputs that come in outside of the function?
  #print(nume)
  #yes
  
  #loglikelihood of continous binomial
  #ll <- sum(dcbinom(mat$mfull[idx_nnz], nume, mix_model[idx_nnz], log = TRUE))
  
  #calculate the value of density function for each value of the matrix that is not BOTH zero in mix_model and mfull
  lf <- dcbinom(mat$mfull[idxnnz], nume, mix_model[idxnnz])
  #lf <- dcbinom(as.vector(mat$mfull), nume, as.vector(mix_model))
  
  #which(lf ==  0)
  print(paste0("...removing ", length(which(lf == 0)), " zero values from likelihood function"))
  #removing 0 values (possibly areas of positive selection?!)
  lf <- lf[lf!=0]
  #calculate log likelihood
  ll <- sum(log(lf))
  
  
  
  #BIC
  #length alpha is the number of parameters  is just 3 here 
  
  BIC = -2*ll + log(nume)*length(alpha)
  #print(length(alpha))
  print(BIC)
  
  return(BIC)
  
}


res1 <- nloptr( x0= c(runif(1), runif(1), runif(1), runif(1)),
                eval_f=eval_f1,
                lb = c(0, 0, 0, 0),
                ub = c(1, 1, 1, 1),
                opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-6, "ftol_rel" = 1e-6)
)
print(res1)



#using optimal alphas create the final mix model matrix

opt_alpha <- res1$solution
#minimum of 20  runs 
opt_alpha <- stability_complex[[19]]$alphas

mix_modelf <- Matrix(0, nrow(mat$bins), nrow(mat$bins), sparse=TRUE)

for (i in 1:4)  { 
  
  mix_modelf[annot.tiles[[i]]] <-  opt_alpha[i]*model1[annot.tiles[[i]]] + (1- opt_alpha[i])*model2[annot.tiles[[i]]]      
  
}

#normalize because multiplying by alphas means that probabilities no longer sum to 1
mix_modelf <- mix_modelf/sum(mix_modelf)       



#save optimal alphas and final mix model matrix as matlab files for import into matlab
writeMat(file.path(SVsig.base, "cbinom_alpha.mat"), mix_model = mix_modelf, opt_alpha = opt_alpha )



#test stability of alpha optimization for complex events 
stability_complex <- lapply(1:20, function(i) { 
  res <- nloptr( x0= c(runif(1), runif(1), runif(1), runif(1)),
                 eval_f=eval_f1,
                 lb = c(0, 0, 0, 0),
                 ub = c(1, 1, 1, 1),
                 opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-6, "ftol_rel" = 1e-6))
  
  return(list(alphas  = res$solution,
              BIC = res$objective)) })

#saveRDS(stability_complex, file = file.path(c.base, "stability_complex.rds"))

#########create figure########

sc <- sapply(stability_complex, function(x) x$alphas)
bc <- sapply(stability_complex, function(x) x$BIC)

pdf(file = file.path(c.base, "20200210stability_complex_BIC.pdf"))
hist(bc, main = "Minimum BIC across 20 Iterations of Constrained Optimization", xlab = "Minimum BIC", col = "blue")
dev.off()

pdf(file = file.path(c.base, "20200210stability_complex_alpha.pdf"))
boxplot(t(sc), names = c("Short", "Long", "Interchromosomal"))
dev.off()



##############for checking ofer's simple events#######################
#test stability of alpha function



stability_ofer <- lapply(1:20, function(i) { 
  res <- nloptr( x0= c(runif(1), runif(1), runif(1)),
                 eval_f=eval_f2,
                 lb = c(0, 0, 0),
                 ub = c(1, 1, 1),
                 opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-6, "ftol_rel" = 1e-6))
  
  return(list(alphas  = res$solution,
              BIC = res$objective)) })

saveRDS(stability_ofer, file = file.path(adj.base, "stability_ofer.rds"))

s <- sapply(stability_ofer, function(x) x$alphas)
b <- sapply(stability_ofer, function(x) x$BIC)

pdf(file = file.path(adj.base, "figs/stability_ofer_BIC.pdf"))
hist(b, main = "Minimum BIC across 20 Iterations of Constrained Optimization", xlab = "Minimum BIC", col = "blue")
dev.off()

pdf(file = file.path(adj.base, "figs/stability_ofer_alpha.pdf"))
boxplot(t(s), names = c("Short", "Long", "Interchromosomal"))
dev.off()

####################3D plots to compare data to each model###############
library(plotly)

#read in bins
bins <- mat$bins
#create bin_id
bins <- cbind(bins, 1:nrow(bins))


#order bins
bins <- bins[order(bins[,4]),]

#create matrix in terms of number of breakpoints




#for short events, long events, and interchromosomal events create surface plot




#mean and variance for short, long and interchromosomal
#short  <- sapply(stability_ofer, function(x) x$alphas)[1,]
#long <- sapply(stability_ofer, f
#inter <- sapply(stability_ofer, function(x) x$alphas[3,])

#mean(short)
#mean(long)
#mean(inter)

#sd(short)
#sd(long)
#sd(inter)

#create overlapping histogram stability plot
pdf(file = file.path(adj.base, "/figs/alpha_stability_boxplot.pdf"))
#hist(short, col =  rgb(0,0,1,0.5))
#hist(long, col = rgb(1, 0.5, 0, 0.5))
#hist(inter, col = rgb(1, 0, 0), add  = TRUE)
boxplot(short, long, inter, col = c("blue", "orange", "red"), main = "Stability of Alpha Parameter", xlab = "Rearrangment Type", names = c("short", "long", "interchromosomal"))


dev.off()


#########CODE GRAVEYARD###############3

#do this for each type of events (short, long, interchromsomal)
#stable <- lapply(1:10, function(i) {

#model1 <- numeric()
#model2 <- numeric()

#for (i in 1:3) {

#idx_annot <- which(annot.tiles[[i]] == TRUE)


#model1 <- model1 + mat$model1[idx_annot] 
#model2 <- model2 + mat$model2[idx_annot]

})
#do this for all three alphas at the same time
#eval_f0 <- function(alpha) {

#create mix_model matrix which is model1*alpha + model2*(1 - alpha)
#background model according to alpha choosen
#mix_model <-  model1 * alpha + model2*(1 - alpha) 

#must both have 
#idx_nnz <- intersect(idx_annot, which(mat$mfull > 0)) 

nume = sum(mat$mfull[idx_nnz])
#is this understanding the inputs that come in outside of the function?
print(nume)

#loglikelihood of continous binomial
ll <- sum(dcbinom(mat$mfull[idx_nnz], nume, mix_model[idx_nnz], log = TRUE))


#BIC
#num param is just 3 
return(-2*ll + log(nume)*mat$num.param[1,1])

}

#how stable is this solution?
stable <- sapply(1:10, function(i) { 
  
  res1 <- nloptr( x0= runif(1),
                  eval_f=eval_f0,
                  lb = 0,
                  ub = 1,
                  opts = list("algorithm"="NLOPT_LN_COBYLA")
  )
  return(res1$solutions)})

###test to make sure its a global minimum





##create plots of continous binomial vs binomial
N = 1e6
p = 0.0001
#simulate discrete binomial
binom <- rbinom(N, 10e3, p)


#dbinom
pdf(file = file.path(adj.base, "figs/binomialcurvs.pdf"))
f <- function(x) dcbinom(x, N, p[3])
h <- Vectorize(f)
curve(h, to = 30, from = 0, n = 31, xlim = NULL, ylab = NULL)
dev.off()





##############compare model1, model2, and mfull from different models##########33

#upper triangle
#conditional fragility function
mat1 <- readMat(file.path(SVsig.base, "debug_alpha_simple.mat"))

#symmetric both sides
#regular fragility function
mat2 <- readMat(file.path(SVsig.base, "88hitsalpha.mat"))

#upper triangle of the matrix
mat2$mfull00[lower.tri(mat2$mfull00)] <- 0 
diag(mat2$mfull00) <- diag(mat2$mfull00)/2


identical(mat1$mfull, mat2$mfull00)
#okay so they are both identical


#do the same for the probability matrices
mat2$p[lower.tri(mat2$p)] <- 0 
diag(mat2$p) <- diag(mat2$p)/2

mat2$p.mult[lower.tri(mat2$p)] <- 0 
diag(mat2$p.mult) <- diag(mat2$p.mult)/2
identical(mat1$model2, mat2$p.mult)
#these are identical


#something is happening in calculating model1 that is wrong/inconsistent!!! 
#aha!!! its no. annot!!! 





