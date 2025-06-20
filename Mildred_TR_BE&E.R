######TEST RUN FOR MILDRED########

library("lme4")
library("glmm")
library("swdpwr")
library("MASS")
library("glmmTMB")
library("meta")
library("tidyverse")


#SW schematic
SWdesmat <- function(Ts) {
  Xsw <- matrix(data=0, ncol = (Ts), nrow = (Ts-1))
  for(i in 1:(Ts-1)) {
    Xsw[i,(i+1):Ts] <- 1
  }
  return(Xsw)
}
#Batched SW schematic
batchSWscheme <- function(olap, S, Ts,  Kseq){
  #Creates a batched SW design schematic with
  #S batches of Ts-period stepped wedge designs, 
  #with overlap between batches of olap periods
  #and K clusters in each sequence
  
  g<- Ts- olap #gap between each batch and the next
  
  SWbatch <- matrix(data=NA, nrow=S*(Ts-1), ncol=Ts + (S-1)*(Ts-olap))  
  
  for(i in 1:S){
    SWbatch[((i-1)*(Ts-1) + 1):(i*(Ts-1)) , 
            ((i-1)*g +1):((i-1)*g + Ts)] <- SWdesmat(Ts)
    
  }
  
  return( SWbatch[sort(rep(1:nrow(SWbatch), Kseq)), ])
  
  
}

#################################################################################################
LongitudinalCRT_onedataset_conts <- 
  function(DesMatrix, Kseq, m, Batch, Teffs, Rxeff, ICC, CAC, corrstruct = 0, treateffvar = 0){
    #A function to generate one dataset from a longitudinal CRT with design matrix
    #given by DesMatrix
    #Kseq: number of clusters assigned to each sequence of DesMatrix
    #m: number of observations in each cluster in each period
    #Batch: vector indicating which row of DesMatrix belongs to which batch
    #Teffs: the vector of time effects
    #Rxeff: the treatment effect
    #ICC: the intracluster correlation
    #CAC: the cluster autocorrelation: exchangeable denoted by CAC = 1. 
    #corrstruct = type of correlation structure
    # corrstruct = 0 indicates either block exchangeable
    # corrstruct = 1 indicates discrete-time decay
    # If CAC = 1, then whatever the value of corrstruct, the model will be exchangeable
    #treateffvar is the variance of the treatment effect (allowed to vary across batches)
    # treateffvar = 0 indicates no treatment effect variation across batches
    # treateffvar > 0 indicates there is treatment effect variation across batches
    
    
    #Assume total variance of 1
    sigma_eps2 = 1-ICC # error variance
    sigmaA2 = ICC # variance of cluster random effect
    
    
    
    #Teffs is a matrix of time effects, of dimension 
    # Ttot columns = total number of periods in the entire design
    # Nseq rows = total number of sequences in the design
    #Some of the elements of Teffs may be NA, when no clusters are observed in that
    #sequence in that period.
    Ttot = ncol(Teffs)
    Nseq = nrow(Teffs)
    #This should be of the same dimension as DesMatrix.
    #Generate a vector for the design and a vector for the time effects
    fulldesmat <- DesMatrix[sort(rep(1:nrow(DesMatrix), Kseq)), ]
    Xvec <- rep(as.vector(t(fulldesmat)), each = m)
    
    fulltimemat <- Teffs[sort(rep(1:nrow(Teffs), Kseq)), ]
    Tvec <- rep(as.vector(t(fulltimemat)), each = m)
    
    #Generate the data
    #Error terms:
    epsi = rnorm(Nseq*Kseq*Ttot*m,mean=0,sd=sqrt(sigma_eps2)) # epsilon (error term)
    
    #Random effects for clusters
    
    if(corrstruct == 1) {
      CPrandomeffs = mvrnorm(n=Nseq*Kseq,mu = rep(0,Ttot),
                             Sigma = sigmaA2 * (CAC ^ as.matrix(dist(1:Ttot))) ) 
    }
    else {
      #CPrandomeffs = rnorm(Nseq*Kseq*Ttot,mean=0, sd=sqrt(sigmaG2))
      CPrandomeffs = mvrnorm(n=Nseq*Kseq,mu = rep(0,Ttot), 
                             Sigma = sigmaA2 *(CAC*matrix(data = 1, nrow=Ttot, ncol=Ttot) 
                                               +(1-CAC)*diag(1, nrow=Ttot)) )
    }
    
    CPrandeffsi <- rep(t(CPrandomeffs), each=m)
    
    clusterVi <- rep(seq(1:(Nseq*Kseq)), each=Ttot*m)
    clusterf = factor(clusterVi)
    
    timeVi <- rep(seq(1:Ttot), each=m)
    timeVi <- rep(timeVi, times=(Kseq*Nseq))
    timef = factor(timeVi)
    
    clusbytime <- rep(1:(Ttot*Kseq*Nseq), each = m)
    clusbytimef = factor(clusbytime)
    
    #indicator for batch
    batch <- rep(Batch, each = m*Ttot*Kseq)
    batchf = factor(batch)
    
    #varying the treatment effect across batches:
    nbatch = as.integer(unique(Batch)[length(unique(Batch))])
    rand_TE = rnorm(nbatch, Rxeff, sqrt(treateffvar)) #generating one value for each batch.
    final_Rxeff = rep(rand_TE, each = Ttot*m*Kseq*Nseq/nbatch) #Replicating each value so that each batch receive its respective treatment effect.
    
    #Put everything together to get the outcomes:
    YY1i = final_Rxeff*Xvec + Tvec + epsi + CPrandeffsi
    Ydf1 = data.frame(YY1i, clusterf, timef, clusbytimef, Xvec, batchf)
    Ydf1 <- Ydf1[!is.na(Xvec),]
    
    return(Ydf1)  
    
    
  }

##################################################################################
#Second step
BatchedPowerSim_conts <- function(DesMatrix, Kseq, m, Batch, Teffs, Rxeff, ICC, CAC, corrstruct, treateffvar){
  #A function to generate and analyse a single dataset
  #DesMatrix = the design matrix
  #Kseq = number of clusters assigned to each sequence
  #m = number of observations in each cluster in each period
  #Teffs = thee matrix of time effects
  #Rxeff = the treatment effect
  #ICC = intracluster correlation
  #CAC = cluster autocorrelation
  #corrstruct = type of correlation structure
  # corrstruct = 0 indicates either block exchangeable
  # corrstruct = 1 indicates discrete-time decay
  # If CAC = 1, then whatever the value of corrstruct, the model will be exchangeable
  #treateffvar is the variance of the treatment effect (allowed to vary across batches)
  # treateffvar = 0 indicates no treatment effect variation across batches
  # treateffvar > 0 indicates there is treatment effect variation across batches
  
  #Simulate the data:
  temp <- LongitudinalCRT_onedataset_conts(DesMatrix, Kseq, m, Batch, Teffs, Rxeff, ICC, CAC, corrstruct, treateffvar)
  
  nbatch = as.integer(unique(Batch)[length(unique(Batch))]) #Getting number of batch.
  
  #Type of linear mixed models:
  #A: shared period effects across batches & no treatment effect heterogeneity across batches
  #B: separate period effects across batches & no treatment effect heterogeneity across batches
  #C: shared period effects across batches & treatment effect heterogeneity across batches
  #D: separate period effects across batches & treatment effect heterogeneity across batches
  
  #tryCatch status explanation:
  #0 indicates no warning and error.
  #1 indicates warning related to non-convergence issue.
  #2 indicates other type of warning that are less harmful.
  #3 indicates error.
  
  if (CAC == 1){ #If satisfied, then fit exchangeable model
    #Model A:
    mod = NULL
    mod_A = tryCatch({
      mod = lmer(YY1i~Xvec +timef + (1|clusterf), temp)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_A_res = mod_A[[1]] #Getting the result from the model
    mod_A_status = mod_A[[2]] #Getting the tryCatch status
    
    if (is.null(mod_A_res)) {
      final_mod_A = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_A_status)
    }
    
    else {
      final_mod_A = c(TE = fixef(mod_A_res)[2], SE_TE =sqrt(vcov(mod_A_res)[2,2]), AIC = AIC(mod_A_res), BIC = BIC(mod_A_res), numb_success = 1, mod_A_status)
    }
    
    
    #Model B:
    mod = NULL
    mod_B = tryCatch({
      mod = lmer(YY1i~Xvec +  timef*batchf + (1|clusterf), temp)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_B_res = mod_B[[1]] #Getting the result from the model
    mod_B_status = mod_B[[2]] #Getting the tryCatch status
    
    if (is.null(mod_B_res)) {
      final_mod_B = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_B_status)
    }
    
    else {
      final_mod_B = c(TE = fixef(mod_B_res)[2], SE_TE = sqrt(vcov(mod_B_res)[2,2]), AIC = AIC(mod_B_res), BIC = BIC(mod_B_res), numb_success = 1, mod_B_status)
    }
    
    
    #Model C:
    mod = NULL
    mod_C = tryCatch({
      mod = lmer(YY1i~Xvec +timef + (1|clusterf) + (-1+Xvec| batchf), temp)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_C_res = mod_C[[1]] #Getting the result from the model
    mod_C_status = mod_C[[2]] #Getting the tryCatch status
    
    if (is.null(mod_C_res)) {
      final_mod_C = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_C_status)
    }
    
    else {
      final_mod_C = c(TE = fixef(mod_C_res)[2], SE_TE = sqrt(vcov(mod_C_res)[2,2]), AIC = AIC(mod_C_res), BIC = BIC(mod_C_res), numb_success = 1, mod_C_status)
    }
    
    #Model D:
    mod = NULL
    mod_D = tryCatch({
      mod = lmer(YY1i~Xvec + timef*batchf + (1|clusterf) + (-1+Xvec| batchf), temp)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_D_res = mod_D[[1]] #Getting the result from the model
    mod_D_status = mod_D[[2]] #Getting the tryCatch status
    
    if (is.null(mod_D_res)) {
      final_mod_D = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA,  mod_D_status)
    }
    
    else {
      final_mod_D = c(TE = fixef(mod_D_res)[2], SE_TE = sqrt(vcov(mod_D_res)[2,2]), AIC = AIC(mod_D_res), BIC = BIC(mod_D_res), numb_success = 1, mod_D_status)
    }
    
    #Creating a matrix to contain the meta-analysis results. There are 3 columns in this matrix to contain: treatment effect, standard error of treatment effect, and status of the trycatch. 
    #Remember trycatch status: 
    #0 implies no warning and error
    #1 implies non-convergence
    #2 implies warning other than non-convergence
    #3 implies error
    
    batchresults <- matrix(data=NA, nrow = nbatch, ncol = 3)
    
    #Loop to collect the treatment effect across batches
    mod = NULL
    for (i in 1:nbatch) {
      mod = NULL
      fit_batch <- tryCatch({
        mod = lmer(YY1i~Xvec +timef + (1|clusterf), 
                   temp[temp$batchf == i,])
        list(mod, 0)
      },
      warning = function(war){
        print(paste("MY_WARNING1: ", war))
        if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
          return(list(NULL, 1))
        }
        else {
          return(list(mod, 2))
        }    
      },
      error = function(err){
        print(paste("MY_ERROR1: ", err))
        return(list(mod, 3))    
      })
      
      fit_batch_res = fit_batch[[1]]
      fit_batch_status = fit_batch[[2]]
      
      if (is.null(fit_batch_res)) {
        batchresults[i,1] <- NA
        batchresults[i,2] <- NA
        batchresults[i,3] <- fit_batch_status
      }
      else {
        batchresults[i,1] <- lme4::fixef(fit_batch_res)[2]
        batchresults[i,2] <- vcov(fit_batch_res)[2,2]
        batchresults[i,3] <- fit_batch_status
      }
    }
    
    # if (sum(is.na(batchresults[,1])) > 0) { #If the batchresults contains any NA then make all 
    #   #treatment effect and SE of treatment effect to be NA so the meta-analysis will just produce NA.
    #   batchresults[, 1] = NA
    #   batchresults[, 2] = NA
    # }
    
    #Pooling the treatment effects using meta-analysis
    
    m.gen.fixed = metagen(TE = batchresults[,1],
                          seTE = sqrt(batchresults[,2]),
                          sm = "MD",
                          random = FALSE)
    
    m.gen.random1 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "classic" )
    
    m.gen.random2 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK" )
    
    m.gen.random3 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "KR" )
    
    m.gen.random4 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "classic" )
    
    m.gen.random5 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "HK" )
    
    m.gen.random6 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "KR" )
    
    m.gen.random7 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "se")
    
    m.gen.random8 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "IQWiG6")
    
    
    m.gen.random9 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "se")
    
    m.gen.random10 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "IQWiG6")
    
    #Number of time models inside the meta-analysis procedure managed to converge. Suppose batch = 5, if 4 out of 5 models ran successfully
    #then number_success = 4 instead of 5.
    numb_success_meta = sum(!(is.na(batchresults[,1])))
    
    final_output = c(mod_A = final_mod_A, 
                     mod_B = final_mod_B,
                     mod_C = final_mod_C,
                     mod_D = final_mod_D,
                     me.fixed.TE = m.gen.fixed$TE.common, me.fixed.SE_TE = m.gen.fixed$seTE.common, me.fixed.lower =  m.gen.fixed$lower.common, me.fixed.upper = m.gen.fixed$upper.common, me.fixed.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.classic.TE = m.gen.random1$TE.random, me.random.reml.classic.SE_TE =  m.gen.random1$seTE.random, me.random.reml.classic.lower = m.gen.random1$lower.random, me.random.reml.classic.upper = m.gen.random1$upper.random, me.random.reml.classic.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.TE = m.gen.random2$TE.random, me.random.reml.HK.SE_TE = m.gen.random2$seTE.random, me.random.reml.HK.lower =  m.gen.random2$lower.random, me.random.reml.HK.upper = m.gen.random2$upper.random, me.random.reml.HK.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.KR.TE = m.gen.random3$TE.random, me.random.reml.KR.SE_TE = m.gen.random3$seTE.random, me.random.reml.KR.lower =  m.gen.random3$lower.random, me.random.reml.KR.upper = m.gen.random3$upper.random, me.random.reml.KR.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.classic.TE = m.gen.random4$TE.random, me.random.SJ.classic.SE_TE = m.gen.random4$seTE.random, me.random.SJ.classic.lower = m.gen.random4$lower.random, me.random.SJ.classic.upper = m.gen.random4$upper.random, me.random.SJ.classic.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.TE = m.gen.random5$TE.random, me.random.SJ.HK.SE_TE = m.gen.random5$seTE.random, me.random.SJ.HK.lower = m.gen.random5$lower.random, me.random.SJ.HK.upper = m.gen.random5$upper.random, me.random.SJ.HK.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.KR.TE = m.gen.random6$TE.random, me.random.SJ.KR.SE_TE = m.gen.random6$seTE.random, me.random.SJ.KR.lower = m.gen.random6$lower.random, me.random.SJ.KR.upper = m.gen.random6$upper.random, me.random.SJ.KR.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.SEadhoc.TE = m.gen.random7$TE.random, me.random.reml.HK.SEadhoc.SE_TE = m.gen.random7$seTE.random, me.random.reml.HK.SEadhoc.lower =  m.gen.random7$lower.random, me.random.reml.HK.SEadhoc.upper = m.gen.random7$upper.random, me.random.reml.HK.SEadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.IQadhoc.TE = m.gen.random8$TE.random, me.random.reml.HK.IQadhoc.SE_TE = m.gen.random8$seTE.random, me.random.reml.HK.IQadhoc.lower =  m.gen.random8$lower.random, me.random.reml.HK.IQadhoc.upper = m.gen.random8$upper.random, me.random.reml.HK.IQadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.SEadhoc.TE = m.gen.random9$TE.random, me.random.SJ.HK.SEadhoc.SE_TE = m.gen.random9$seTE.random, me.random.SJ.HK.SEadhoc.lower = m.gen.random9$lower.random, me.random.SJ.HK.SEadhoc.upper = m.gen.random9$upper.random, me.random.SJ.HK.SEadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.IQadhoc.TE = m.gen.random10$TE.random, me.random.SJ.HK.IQadhoc.SE_TE = m.gen.random10$seTE.random, me.random.SJ.HK.IQadhoc.lower = m.gen.random10$lower.random, me.random.SJ.HK.IQadhoc.upper = m.gen.random10$upper.random, me.random.SJ.HK.IQadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ",")
    )
                     
    
    return(final_output)
    
  }
  
  if (corrstruct != 1 & CAC != 1){ #If satisfied, then fit block exchangeable model
    
    #Model A: 
    mod = NULL
    mod_A = tryCatch({
      mod = lmer(YY1i~Xvec +timef + (1|clusterf) +(1|clusbytimef), temp)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_A_res = mod_A[[1]] #Getting the result from the model
    mod_A_status = mod_A[[2]] #Getting the tryCatch status
    
    if (is.null(mod_A_res)) {
      final_mod_A = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_A_status)
    }
    
    else {
      final_mod_A = c(TE = fixef(mod_A_res)[2], SE_TE =sqrt(vcov(mod_A_res)[2,2]), AIC = AIC(mod_A_res), BIC = BIC(mod_A_res), numb_success = 1, mod_A_status)
    }
    
    
    #Model B:
    mod = NULL
    mod_B = tryCatch({
      mod = lmer(YY1i~Xvec +timef*batchf + (1|clusterf) +(1|clusbytimef), temp)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_B_res = mod_B[[1]] #Getting the result from the model
    mod_B_status = mod_B[[2]] #Getting the tryCatch status
    
    if (is.null(mod_B_res)) {
      final_mod_B = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_B_status)
    }
    
    else {
      final_mod_B = c(TE = fixef(mod_B_res)[2], SE_TE = sqrt(vcov(mod_B_res)[2,2]), AIC = AIC(mod_B_res), BIC = BIC(mod_B_res), numb_success = 1, mod_B_status)
    }
    
    
    #Model C:
    mod = NULL
    mod_C = tryCatch({
      mod = lmer(YY1i~Xvec +timef + (1|clusterf) +(1|clusbytimef) + (-1+Xvec| batchf), temp)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_C_res = mod_C[[1]] #Getting the result from the model
    mod_C_status = mod_C[[2]] #Getting the tryCatch status
    
    if (is.null(mod_C_res)) {
      final_mod_C = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_C_status)
    }
    
    else {
      final_mod_C = c(TE = fixef(mod_C_res)[2], SE_TE = sqrt(vcov(mod_C_res)[2,2]), AIC = AIC(mod_C_res), BIC = BIC(mod_C_res), numb_success = 1, mod_C_status)
    }
    
    #Model D:
    mod = NULL
    mod_D = tryCatch({
      mod = lmer(YY1i~Xvec +timef*batchf + (1|clusterf) +(1|clusbytimef) + (-1+Xvec| batchf), temp)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_D_res = mod_D[[1]] #Getting the result from the model
    mod_D_status = mod_D[[2]] #Getting the tryCatch status
    
    if (is.null(mod_D_res)) {
      final_mod_D = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_D_status)
    }
    
    else {
      final_mod_D = c(TE = fixef(mod_D_res)[2], SE_TE = sqrt(vcov(mod_D_res)[2,2]), AIC = AIC(mod_D_res), BIC = BIC(mod_D_res), numb_success = 1, mod_D_status)
    }
    
    #Creating a matrix to contain the meta-analysis results. There are 3 columns in this matrix to contain: treatment effect, standard error of treatment effect, and status of the trycatch. 
    #Remember trycatch status: 
    #0 implies no error or warning
    #1 implies non-convergence
    #2 implies warning other than non-convergence
    #3 implies error
    batchresults <- matrix(data=NA, nrow = nbatch, ncol = 3)
    
    #Loop to collect the treatment effect across batches
    mod = NULL
    for (i in 1:nbatch) {
      mod = NULL
      fit_batch <- tryCatch({
        mod = lmer(YY1i~Xvec +timef + (1|clusterf) +(1|clusbytimef), 
                   temp[temp$batchf == i,])
        list(mod, 0)
      },
      warning = function(war){
        print(paste("MY_WARNING1: ", war))
        if(grepl("Model failed to converge", war$message, ignore.case = TRUE)) {
          return(list(NULL, 1))
        }
        else {
          return(list(mod, 2))
        }    
      },
      error = function(err){
        print(paste("MY_ERROR1: ", err))
        return(list(mod, 3))    
      })
      
      fit_batch_res = fit_batch[[1]]
      fit_batch_status = fit_batch[[2]]
      
      if (is.null(fit_batch_res)) {
        batchresults[i,1] <- NA
        batchresults[i,2] <- NA
        batchresults[i,3] <- fit_batch_status
      }
      else {
        batchresults[i,1] <- lme4::fixef(fit_batch_res)[2]
        batchresults[i,2] <- vcov(fit_batch_res)[2,2]
        batchresults[i,3] <- fit_batch_status
      }
    }
    
    #if (sum(is.na(batchresults[,1])) > 0) { #If the batchresults contains any NA then make all 
    #  #treatment effect and SE of treatment effect to be NA so the meta-analysis will just produce NA.
    #  batchresults[, 1] = NA
    #  batchresults[, 2] = NA
    #}
    
    #Pooling the treatment effects using meta-analysis
    
    m.gen.fixed = metagen(TE = batchresults[,1],
                          seTE = sqrt(batchresults[,2]),
                          sm = "MD",
                          random = FALSE)
    
    m.gen.random1 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "classic" )
    
    m.gen.random2 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK" )
    
    m.gen.random3 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "KR" )
    
    m.gen.random4 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "classic" )
    
    m.gen.random5 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "HK" )
    
    m.gen.random6 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "KR" )
    
    m.gen.random7 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "se")
    
    m.gen.random8 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "IQWiG6")
    
    
    m.gen.random9 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "se")
    
    m.gen.random10 <- metagen(TE = batchresults[,1],
                              seTE = sqrt(batchresults[,2]),
                              sm = "MD",
                              common = FALSE,
                              method.tau = "SJ",
                              method.random.ci = "HK",
                              adhoc.hakn.ci = "IQWiG6")
    
    #Number of time models inside the meta-analysis procedure managed to converge. Suppose batch = 5, if 4 out of 5 models ran successfully
    #then number_success = 4 instead of 5.
    numb_success_meta = sum(!(is.na(batchresults[,1])))
    
    final_output = c(mod_A = final_mod_A, 
                     mod_B = final_mod_B,
                     mod_C = final_mod_C,
                     mod_D = final_mod_D,
                     me.fixed.TE = m.gen.fixed$TE.common, me.fixed.SE_TE = m.gen.fixed$seTE.common, me.fixed.lower =  m.gen.fixed$lower.common, me.fixed.upper = m.gen.fixed$upper.common, me.fixed.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.classic.TE = m.gen.random1$TE.random, me.random.reml.classic.SE_TE =  m.gen.random1$seTE.random, me.random.reml.classic.lower = m.gen.random1$lower.random, me.random.reml.classic.upper = m.gen.random1$upper.random, me.random.reml.classic.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.TE = m.gen.random2$TE.random, me.random.reml.HK.SE_TE = m.gen.random2$seTE.random, me.random.reml.HK.lower =  m.gen.random2$lower.random, me.random.reml.HK.upper = m.gen.random2$upper.random, me.random.reml.HK.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.KR.TE = m.gen.random3$TE.random, me.random.reml.KR.SE_TE = m.gen.random3$seTE.random, me.random.reml.KR.lower =  m.gen.random3$lower.random, me.random.reml.KR.upper = m.gen.random3$upper.random, me.random.reml.KR.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.classic.TE = m.gen.random4$TE.random, me.random.SJ.classic.SE_TE = m.gen.random4$seTE.random, me.random.SJ.classic.lower = m.gen.random4$lower.random, me.random.SJ.classic.upper = m.gen.random4$upper.random, me.random.SJ.classic.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.TE = m.gen.random5$TE.random, me.random.SJ.HK.SE_TE = m.gen.random5$seTE.random, me.random.SJ.HK.lower = m.gen.random5$lower.random, me.random.SJ.HK.upper = m.gen.random5$upper.random, me.random.SJ.HK.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.KR.TE = m.gen.random6$TE.random, me.random.SJ.KR.SE_TE = m.gen.random6$seTE.random, me.random.SJ.KR.lower = m.gen.random6$lower.random, me.random.SJ.KR.upper = m.gen.random6$upper.random, me.random.SJ.KR.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.SEadhoc.TE = m.gen.random7$TE.random, me.random.reml.HK.SEadhoc.SE_TE = m.gen.random7$seTE.random, me.random.reml.HK.SEadhoc.lower =  m.gen.random7$lower.random, me.random.reml.HK.SEadhoc.upper = m.gen.random7$upper.random, me.random.reml.HK.SEadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.IQadhoc.TE = m.gen.random8$TE.random, me.random.reml.HK.IQadhoc.SE_TE = m.gen.random8$seTE.random, me.random.reml.HK.IQadhoc.lower =  m.gen.random8$lower.random, me.random.reml.HK.IQadhoc.upper = m.gen.random8$upper.random, me.random.reml.HK.IQadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.SEadhoc.TE = m.gen.random9$TE.random, me.random.SJ.HK.SEadhoc.SE_TE = m.gen.random9$seTE.random, me.random.SJ.HK.SEadhoc.lower = m.gen.random9$lower.random, me.random.SJ.HK.SEadhoc.upper = m.gen.random9$upper.random, me.random.SJ.HK.SEadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.IQadhoc.TE = m.gen.random10$TE.random, me.random.SJ.HK.IQadhoc.SE_TE = m.gen.random10$seTE.random, me.random.SJ.HK.IQadhoc.lower = m.gen.random10$lower.random, me.random.SJ.HK.IQadhoc.upper = m.gen.random10$upper.random, me.random.SJ.HK.IQadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ",")
    )
    
    
    return(final_output)
  }
  
  if (corrstruct == 1 & CAC != 1){ #If satisfied, then fit discrete time decay model
    
    #Model A: 
    mod = NULL
    mod_A = tryCatch({
      mod = glmmTMB(YY1i ~ Xvec + timef + ar1(timef + 0 | clusterf), REML=T, 
                    data = temp, family = gaussian)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("convergence", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_A_res = mod_A[[1]] #Getting the result from the model
    mod_A_status = mod_A[[2]] #Getting the tryCatch status
    
    if (is.null(mod_A_res)) {
      final_mod_A = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_A_status)
    }
    
    else {
      final_mod_A = c(TE = fixef(mod_A_res)$cond[2], SE_TE =sqrt(vcov(mod_A_res)$cond[2,2]), AIC = AIC(mod_A_res), BIC = BIC(mod_A_res), numb_success = 1, mod_A_status)
    }
    
    
    #Model B:
    mod = NULL
    mod_B = tryCatch({
      mod = glmmTMB(YY1i ~ Xvec + timef*batchf + ar1(timef + 0 | clusterf), REML=T, 
                    data = temp, family = gaussian)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("convergence", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_B_res = mod_B[[1]] #Getting the result from the model
    mod_B_status = mod_B[[2]] #Getting the tryCatch status
    
    if (is.null(mod_B_res)) {
      final_mod_B = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_B_status)
    }
    
    else {
      final_mod_B = c(TE = fixef(mod_B_res)$cond[2], SE_TE = sqrt(vcov(mod_B_res)$cond[2,2]), AIC = AIC(mod_B_res), BIC = BIC(mod_B_res), numb_success = 1, mod_B_status)
    }
    
    
    #Model C:
    mod = NULL
    mod_C = tryCatch({
      mod = glmmTMB(YY1i ~ Xvec + timef + ar1(timef + 0 | clusterf)
                    + (Xvec + 0| batchf), REML=T, 
                    data = temp, family = gaussian)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("convergence", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_C_res = mod_C[[1]] #Getting the result from the model
    mod_C_status = mod_C[[2]] #Getting the tryCatch status
    
    if (is.null(mod_C_res)) {
      final_mod_C = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_C_status)
    }
    
    else {
      final_mod_C = c(TE = fixef(mod_C_res)$cond[2], SE_TE = sqrt(vcov(mod_C_res)$cond[2,2]), AIC = AIC(mod_C_res), BIC = BIC(mod_C_res), numb_success = 1, mod_C_status)
    }
    
    #Model D:
    mod = NULL
    mod_D = tryCatch({
      mod = glmmTMB(YY1i ~ Xvec + timef*batchf + ar1(timef + 0 | clusterf)
                    + (Xvec + 0| batchf), REML=T, 
                    data = temp, family = gaussian)
      list(mod, 0)
    },
    warning = function(war){
      print(paste("MY_WARNING1: ", war))
      if(grepl("convergence", war$message, ignore.case = TRUE)) {
        return(list(NULL, 1))
      }
      else {
        return(list(mod, 2))
      }    
    },
    error = function(err){
      print(paste("MY_ERROR1: ", err))
      return(list(mod, 3))    
    })
    
    mod_D_res = mod_D[[1]] #Getting the result from the model
    mod_D_status = mod_D[[2]] #Getting the tryCatch status
    
    if (is.null(mod_D_res)) {
      final_mod_D = c(TE = NA, SE_TE = NA, AIC = NA, BIC = NA, numb_success = NA, mod_D_status)
    }
    
    else {
      final_mod_D = c(TE = fixef(mod_D_res)$cond[2], SE_TE = sqrt(vcov(mod_D_res)$cond[2,2]), AIC = AIC(mod_D_res), BIC = BIC(mod_D_res), numb_success = 1, mod_D_status)
    }
    
    #Creating a matrix to contain the meta-analysis results. There are 3 columns in this matrix to contain: treatment effect, standard error of treatment effect, and status of the trycatch. 
    #Remember trycatch status: 
    #0 implies no error or warning
    #1 implies non-convergence
    #2 implies warning other than non-convergence
    #3 implies error
    batchresults <- matrix(data=NA, nrow = nbatch, ncol = 3)
    
    #Loop to collect the treatment effect across batches
    mod = NULL
    for (i in 1:nbatch) {
      mod = NULL
      fit_batch <- tryCatch({
        mod = glmmTMB(YY1i ~ Xvec + timef + ar1(timef + 0 | clusterf), REML=T, 
                      data = temp[temp$batchf == i,], family = gaussian)
        list(mod, 0)
      },
      warning = function(war){
        print(paste("MY_WARNING1: ", war))
        if(grepl("convergence", war$message, ignore.case = TRUE)) {
          return(list(NULL, 1))
        }
        else {
          return(list(mod, 2))
        }    
      },
      error = function(err){
        print(paste("MY_ERROR1: ", err))
        return(list(mod, 3))    
      })
      
      fit_batch_res = fit_batch[[1]]
      fit_batch_status = fit_batch[[2]]
      
      if (is.null(fit_batch_res)) {
        batchresults[i,1] <- NA
        batchresults[i,2] <- NA
        batchresults[i,3] <- fit_batch_status
      }
      else {
        batchresults[i,1] <- fixef(fit_batch_res)$cond[2]
        batchresults[i,2] <- vcov(fit_batch_res)$cond[2,2]
        batchresults[i,3] <- fit_batch_status
      }
    }
    
    #if (sum(is.na(batchresults[,1])) > 0) { #If the batchresults contains any NA then make all 
    #  #treatment effect and SE of treatment effect to be NA so the meta-analysis will just produce NA.
    #  batchresults[, 1] = NA
    #  batchresults[, 2] = NA
    #}
    
    #Pooling the treatment effects using meta-analysis
    
    m.gen.fixed = metagen(TE = batchresults[,1],
                          seTE = sqrt(batchresults[,2]),
                          sm = "MD",
                          random = FALSE)
    
    m.gen.random1 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "classic" )
    
    m.gen.random2 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK" )
    
    m.gen.random3 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "KR" )
    
    m.gen.random4 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "classic" )
    
    m.gen.random5 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "HK" )
    
    m.gen.random6 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "KR" )
    
    m.gen.random7 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "se")
    
    m.gen.random8 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "REML",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "IQWiG6")
    
    
    m.gen.random9 <- metagen(TE = batchresults[,1],
                             seTE = sqrt(batchresults[,2]),
                             sm = "MD",
                             common = FALSE,
                             method.tau = "SJ",
                             method.random.ci = "HK",
                             adhoc.hakn.ci = "se")
    
    m.gen.random10 <- metagen(TE = batchresults[,1],
                              seTE = sqrt(batchresults[,2]),
                              sm = "MD",
                              common = FALSE,
                              method.tau = "SJ",
                              method.random.ci = "HK",
                              adhoc.hakn.ci = "IQWiG6")
    
    #Number of time models inside the meta-analysis procedure managed to converge. Suppose batch = 5, if 4 out of 5 models ran successfully
    #then number_success = 4 instead of 5.
    numb_success_meta = sum(!(is.na(batchresults[,1])))
    
    final_output = c(mod_A = final_mod_A, 
                     mod_B = final_mod_B,
                     mod_C = final_mod_C,
                     mod_D = final_mod_D,
                     me.fixed.TE = m.gen.fixed$TE.common, me.fixed.SE_TE = m.gen.fixed$seTE.common, me.fixed.lower =  m.gen.fixed$lower.common, me.fixed.upper = m.gen.fixed$upper.common, me.fixed.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.classic.TE = m.gen.random1$TE.random, me.random.reml.classic.SE_TE =  m.gen.random1$seTE.random, me.random.reml.classic.lower = m.gen.random1$lower.random, me.random.reml.classic.upper = m.gen.random1$upper.random, me.random.reml.classic.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.TE = m.gen.random2$TE.random, me.random.reml.HK.SE_TE = m.gen.random2$seTE.random, me.random.reml.HK.lower =  m.gen.random2$lower.random, me.random.reml.HK.upper = m.gen.random2$upper.random, me.random.reml.HK.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.KR.TE = m.gen.random3$TE.random, me.random.reml.KR.SE_TE = m.gen.random3$seTE.random, me.random.reml.KR.lower =  m.gen.random3$lower.random, me.random.reml.KR.upper = m.gen.random3$upper.random, me.random.reml.KR.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.classic.TE = m.gen.random4$TE.random, me.random.SJ.classic.SE_TE = m.gen.random4$seTE.random, me.random.SJ.classic.lower = m.gen.random4$lower.random, me.random.SJ.classic.upper = m.gen.random4$upper.random, me.random.SJ.classic.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.TE = m.gen.random5$TE.random, me.random.SJ.HK.SE_TE = m.gen.random5$seTE.random, me.random.SJ.HK.lower = m.gen.random5$lower.random, me.random.SJ.HK.upper = m.gen.random5$upper.random, me.random.SJ.HK.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.KR.TE = m.gen.random6$TE.random, me.random.SJ.KR.SE_TE = m.gen.random6$seTE.random, me.random.SJ.KR.lower = m.gen.random6$lower.random, me.random.SJ.KR.upper = m.gen.random6$upper.random, me.random.SJ.KR.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.SEadhoc.TE = m.gen.random7$TE.random, me.random.reml.HK.SEadhoc.SE_TE = m.gen.random7$seTE.random, me.random.reml.HK.SEadhoc.lower =  m.gen.random7$lower.random, me.random.reml.HK.SEadhoc.upper = m.gen.random7$upper.random, me.random.reml.HK.SEadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.reml.HK.IQadhoc.TE = m.gen.random8$TE.random, me.random.reml.HK.IQadhoc.SE_TE = m.gen.random8$seTE.random, me.random.reml.HK.IQadhoc.lower =  m.gen.random8$lower.random, me.random.reml.HK.IQadhoc.upper = m.gen.random8$upper.random, me.random.reml.HK.IQadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.SEadhoc.TE = m.gen.random9$TE.random, me.random.SJ.HK.SEadhoc.SE_TE = m.gen.random9$seTE.random, me.random.SJ.HK.SEadhoc.lower = m.gen.random9$lower.random, me.random.SJ.HK.SEadhoc.upper = m.gen.random9$upper.random, me.random.SJ.HK.SEadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ","),
                     me.random.SJ.HK.IQadhoc.TE = m.gen.random10$TE.random, me.random.SJ.HK.IQadhoc.SE_TE = m.gen.random10$seTE.random, me.random.SJ.HK.IQadhoc.lower = m.gen.random10$lower.random, me.random.SJ.HK.IQadhoc.upper = m.gen.random10$upper.random, me.random.SJ.HK.IQadhoc.numb_success = numb_success_meta, paste(batchresults[,3], collapse = ",")
    )
    
    return(final_output)
  } 
}



##################################################################################
#Functions for continuous outcomes
BatchedPowerSim_conts_wrap <- 
  function(nrep, Ts, nbatch, olap, Kseq, m, TimeEffsInd, Rxeff, ICC, CAC, corrstruct, treateffvar){
    #Wrapper function forthe BatchedPowerSim function. 
    #nrep = number of replications in the simulation study
    #Ts = number of periods in each batch of the design
    #olap = number of periods of overlap between successive batches
    #Kseq = number of clusters assigned to each sequence
    #m = number of observations in each cluster in each period
    #TimeEffsInd = 0 if time effects are shared across batches; 
    #            = 1 if separate time effects for each batch
    #Rxeff = the treatment effect
    #ICC = intracluster correlation
    #CAC = cluster autocorrelation
    #corrstruct = type of correlation structure
    # corrstruct = 0 indicates either block exchangeable
    # corrstruct = 1 indicates discrete-time decay
    # If CAC = 1, then whatever the value of corrstruct, the model will be exchangeable
    #treateffvar is the variance of the treatment effect (allowed to vary across batches)
    # treateffvar = 0 indicates no treatment effect variation across batches
    # treateffvar > 0 indicates there is treatment effect variation across batches
    
    #Generate the design matrix
    DesMatrix <- batchSWscheme(olap, nbatch, Ts, 1)
    Batch <- rep(seq(1:nbatch), each = (Ts-1))
    
    
    if(TimeEffsInd == 0){
      #Non-null time effects, shared across batches
      TimeEffs <- 0.1*matrix(data=seq(1:ncol(DesMatrix)), nrow=nrow(DesMatrix), 
                             ncol=ncol(DesMatrix), byrow=TRUE)
      TimeEffs <- TimeEffs*(!is.na(DesMatrix))
    }
    else if(TimeEffsInd == 1){
      #Non-null time effects, NOT shared across batches
      TimeEffs <- matrix(data=rnorm(n=ncol(DesMatrix)*nbatch), nrow=nbatch, ncol=ncol(DesMatrix) )
      TimeEffs <- TimeEffs[sort(rep(1:nrow(TimeEffs), (Ts-1))), ]
      TimeEffs <- TimeEffs*(!is.na(DesMatrix))   
    }
    
    output <- replicate(nrep, BatchedPowerSim_conts(DesMatrix, Kseq, m, Batch, TimeEffs, Rxeff, ICC, CAC, corrstruct, treateffvar))
    
    #Getting the tryCatch() status
    output_status = output[seq(from = 6, to = nrow(output), by = 6),] 
    
    #Getting the output without tryCatch() status so the matrix can be converted into numeric.
    output = output[-seq(from = 6, to = nrow(output), by = 6), ] 
    
    #Now convert it into numeric.
    class(output) = "numeric"
    
    #4 models with correct within-cluster correlation structure have been fit to the simulated datasets:
    #A: shared period effects across batches & no treatment effect heterogeneity across batches
    #B: separate period effects across batches & no treatment effect heterogeneity across batches
    #C: shared period effects across batches & treatment effect heterogeneity across batches
    #D: separate period effects across batches & treatment effect heterogeneity across batches
    
    #For each, want the following:
    #- bias
    #-Monte Carlo SE of bias estimate
    #- empirical standard error
    #- MSE
    #- Average Model SE
    #- Coverage
    #- Number of replication with converged model and no errors (nrep_conv)
    #- Number of models inside the procedures across all replications that managed to run without an issue (non-harmful warning, no error, or managed to converge)
    #- Number of warnings
    #- Number of errors
    
    #NOTE: nrep_conv and numb_success are the same for non meta-analysis models, but will be different for the meta-analysis models.
    #For meta-anaylsis models, nrep_conv returns the number of replication where the final result of meta-analysis is available (non-missing).
    #However, numb_conv returns the number of times the models inside the meta-analysis procedure converged/has no errors.
    #Suppose, there are 5 batches and 3 replications. Then there are 15 (5 models one for each batch repeated 3 times) models that need to be run for the meta-analysis.
    #If 12 models out of 15 managed to run without an issue, then numb_conv = 12. But if only 1 replication out of 3 replications has full data (no error or convergence issue)
    #then nrep_conv = 1. 
    
    #NOTE: Due to the same issue as before, numb_warnings and numb_errors will be higher for meta-analysis models.
    
    myresults <- data.frame(model = rep("A", 15), Bias = 1, Bias_MC_SE = 1, Emp_SE = 1, MSE = 1, Avg_ModSE = 1, Coverage = 1, nrep_conv = 1, numb_success = 1, numb_warnings = 1, numb_errors = 1)
    
    myresults = myresults %>%
      mutate(model = c("mod_A", "mod_B", "mod_C", "mod_D", "MA.fixed", "MA.random.reml.classic", "MA.random.reml.HK", "MA.random.reml.KR", "MA.random.SJ.classic", "MA.random.SJ.HK", "MA.random.SJ.KR", "MA.random.reml.HK.SEadhoc", "MA.random.reml.HK.IQadhoc",  "MA.random.SJ.HK.SEadhoc", "MA.random.SJ.HK.IQadhoc"))
    
    #Increment of 4 because each model contains: treatment effect estimate (TEE), SE of TEE, AIC, BIC, and number of success.
    seq = seq(from = 1, to = nrow(output), 5) 
    
    for (i in 1:15) {
      nrep_conv = nrep - sum(is.na(output[seq[i], ]))
      
      #Bias
      myresults[i, 2] = mean(output[seq[i], ], na.rm=TRUE)-Rxeff
      
      #Bias monte carlo SE
      myresults[i, 3] = sqrt(sum((output[seq[i],]-mean(output[seq[i],], na.rm = T))^2, na.rm = T)/(nrep_conv*(nrep_conv-1)))
      
      #Empirical SE
      myresults[i, 4] =  sqrt(sum((output[seq[i],]-mean(output[seq[i],], na.rm = T))^2, na.rm = T)/(nrep_conv-1))
      
      #MSE
      myresults[i, 5] = (sum((output[seq[i],]-Rxeff)^2, na.rm = T)/(nrep_conv))
      
      #Average model SE
      myresults[i, 6] = sqrt(sum((output[seq[i]+1,])^2, na.rm = T)/(nrep_conv))
      
      #Coverage
      
      if (myresults[i, "model"] %in% c("mod_A", "mod_B", "mod_C", "mod_D")) { #This is for the linear mixed model approach.
        myresults[i, 7] = sum((output[seq[i],] - (-qnorm((1-0.95)/2)*output[seq[i]+1,])) <= Rxeff & Rxeff <= (output[seq[i],] + (-qnorm((1-0.95)/2)*output[seq[i]+1,])), na.rm = T) / nrep_conv
      }
      
      if (!(myresults[i, "model"] %in% c("mod_A", "mod_B", "mod_C", "mod_D"))) { #This is for the meta-analysis approach.
        myresults[i, 7] = sum(output[seq[i]+2,] <= Rxeff & Rxeff <= output[seq[i]+3,] , na.rm = T) / nrep_conv
      }
    
      #Number of replications managed to converge
      myresults[i, 8] = nrep_conv
      
      #Number of models inside the procedures across all replications that managed to run without issues.
      myresults[i, 9] = sum(output[seq[i] + 4,], na.rm = T)
      
      #Number of warnings
      myresults[i, 10] = sum(str_count(output_status[i,], "1")) + sum(str_count(output_status[i,], "2"))
      
      #Number of errors
      myresults[i, 11] = sum(str_count(output_status[i,], "3")) 
      
    }
    
    #Returning both the results of the calculation and also the output.
    return(list(myresults, output))
    
  }


###########################Set-up for the simulation###########################

#Considering only exchangeable and block-exchangeable models.

Tlist <- c(4, 6) 
nbatchlist <- c(2, 5) 
olaplist <- c(1, 999) #999 represent T-2, this is coded so it can be calculated later.
Klist <- c(2, 5) 
mlist <- c(20, 100) 
sharedtimeeffs <- c(0, 1)
effsize <- c(0, 0.2)
icclist <- c(0.01, 0.05, 0.2) 
caclist <- c(1, 0.9, 0.5)
corrstruct <- c(0) #Exclude 1 to exclude discrete-time-decay model.
treateffvar <- c(0, 0.025)

allcombos <- expand.grid(Tlist, nbatchlist, olaplist, Klist, mlist, sharedtimeeffs, effsize, icclist, caclist, corrstruct, treateffvar)

allcombos <- allcombos %>%
  rename(T = Var1,
         nbatch = Var2,
         nolap = Var3,
         K = Var4,
         M = Var5,
         sharedtime = Var6,
         effsize = Var7,
         ICC = Var8,
         CAC = Var9,
         corrstruct = Var10,
         treateffvar = Var11)

#Calculating T-2 in the olap

allcombos <- allcombos %>%
  mutate(nolap = ifelse(nolap == 999, T-2, nolap))

#Filtering out CAC = 1 and corrstruct = 1, so there won't be duplicate.
table(allcombos$CAC, allcombos$corrstruct) #11,664 observations should be removed.

nrow(allcombos) - nrow(allcombos %>%
                         filter(!(CAC == 1 & corrstruct == 1))) #11,664 removed.

allcombos <- allcombos %>%
  filter(!(CAC == 1 & corrstruct == 1))

#Giving ID to each unique combination so it can be join with other data later.

allcombos <- allcombos %>%
  mutate(comb_ID = row_number())


for (i in 1:nrow(allcombos)) {
  nrep = 500 #Changed from 1,000 to 500
  if (i == 1) {
    start_time <- format(Sys.time(), "%a %b %d %X %Y %Z")
    timeStart<-Sys.time()
    
    set.seed(allcombos[i, 12])
    results = BatchedPowerSim_conts_wrap(nrep, allcombos[i,1], allcombos[i,2], allcombos[i,3],
                                         allcombos[i,4], allcombos[i,5], allcombos[i,6], 
                                         allcombos[i,7], allcombos[i,8],allcombos[i,9],
                                         allcombos[i, 10], allcombos[i, 11])
    
    all_contsresults = results[[1]]
    
    all_contsresults = all_contsresults %>%
      mutate(comb_ID = i,
             T = allcombos[i,1],
             nbatch = allcombos[i,2],
             nolap = allcombos[i,3],
             K = allcombos[i,4],
             M = allcombos[i,5],
             sharedtime = allcombos[i,6],
             effsize = allcombos[i,7],
             ICC = allcombos[i,8],
             CAC = allcombos[i,9],
             corrstruct = allcombos[i,10],
             treateffvar = allcombos[i,11],
             comb_ID = allcombos[i,12]) %>%
      mutate(Bias_minusSE = Bias - 2*Bias_MC_SE,
             Bias_plusSE = Bias + 2*Bias_MC_SE,
             Coverage_MC_SE = sqrt((Coverage*(1-Coverage))/nrep_conv),
             Coverage_minusSE = Coverage - 2*Coverage_MC_SE,
             Coverage_plusSE = Coverage + 2*Coverage_MC_SE,
             EmpSE_MC_SE = Emp_SE/sqrt(2*(nrep_conv - 1)),
             EmpSE_minusSE = Emp_SE - 2*EmpSE_MC_SE,
             EmpSE_plusSE = Emp_SE + 2*EmpSE_MC_SE
      )
    
    all_output = results[[2]]
    
    all_output = cbind(all_output, T = allcombos[i,1],
                       nbatch = allcombos[i,2],
                       nolap = allcombos[i,3],
                       K = allcombos[i,4],
                       M = allcombos[i,5],
                       sharedtime = allcombos[i,6],
                       effsize = allcombos[i,7],
                       ICC = allcombos[i,8],
                       CAC = allcombos[i,9],
                       corrstruct = allcombos[i,10],
                       treateffvar = allcombos[i,11],
                       comb_ID = allcombos[i,12],
                       estimates = rownames(all_output))
    
    end_time <- format(Sys.time(), "%a %b %d %X %Y %Z")
    timeEnd<-Sys.time()
    durt <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))
    
    all_tim_df <- cbind(T = allcombos[i,1],
                        nbatch = allcombos[i,2],
                        nolap = allcombos[i,3],
                        K = allcombos[i,4],
                        M = allcombos[i,5],
                        sharedtime = allcombos[i,6],
                        effsize = allcombos[i,7],
                        ICC = allcombos[i,8],
                        CAC = allcombos[i,9],
                        corrstruct = allcombos[i,10],
                        treateffvar = allcombos[i,11],
                        comb_ID = allcombos[i,12],
                        seed = as.character(allcombos[i, 12]), 
                        start_time = start_time,
                        end_tim = end_time,
                        duration = durt)
    
    write.table(all_contsresults, "all_contsresults_BE&E.csv", sep = ",",
                row.names = FALSE, col.names = i== 1, append = i > 1)
    
    write.table(all_output, "all_output_BE&E.csv", sep = ",",
                row.names = FALSE, col.names = i== 1, append = i > 1)
    
    write.table(all_tim_df, "all_tim_df_BE&E.csv", sep = ",",
                row.names = FALSE, col.names = i== 1, append = i > 1)
    
  }
  
  if (i != 1) {
    start_time <- format(Sys.time(), "%a %b %d %X %Y %Z")
    timeStart<-Sys.time()
    
    set.seed(allcombos[i, 12])
    results = BatchedPowerSim_conts_wrap(nrep, allcombos[i,1], allcombos[i,2], allcombos[i,3],
                                         allcombos[i,4], allcombos[i,5], allcombos[i,6], 
                                         allcombos[i,7], allcombos[i,8],allcombos[i,9],
                                         allcombos[i, 10], allcombos[i, 11])
    
    contsresults = results[[1]] %>%
      mutate(comb_ID = i,
             T = allcombos[i,1],
             nbatch = allcombos[i,2],
             nolap = allcombos[i,3],
             K = allcombos[i,4],
             M = allcombos[i,5],
             sharedtime = allcombos[i,6],
             effsize = allcombos[i,7],
             ICC = allcombos[i,8],
             CAC = allcombos[i,9],
             corrstruct = allcombos[i,10],
             treateffvar = allcombos[i,11],
             comb_ID = allcombos[i,12]) %>%
      mutate(Bias_minusSE = Bias - 2*Bias_MC_SE,
             Bias_plusSE = Bias + 2*Bias_MC_SE,
             Coverage_MC_SE = sqrt((Coverage*(1-Coverage))/nrep_conv),
             Coverage_minusSE = Coverage - 2*Coverage_MC_SE,
             Coverage_plusSE = Coverage + 2*Coverage_MC_SE,
             EmpSE_MC_SE = Emp_SE/sqrt(2*(nrep_conv - 1)),
             EmpSE_minusSE = Emp_SE - 2*EmpSE_MC_SE,
             EmpSE_plusSE = Emp_SE + 2*EmpSE_MC_SE
      )
    
    all_contsresults = all_contsresults %>%
      rbind(contsresults)
    
    output = results[[2]]
    
    output = cbind(output, T = allcombos[i,1],
                   nbatch = allcombos[i,2],
                   nolap = allcombos[i,3],
                   K = allcombos[i,4],
                   M = allcombos[i,5],
                   sharedtime = allcombos[i,6],
                   effsize = allcombos[i,7],
                   ICC = allcombos[i,8],
                   CAC = allcombos[i,9],
                   corrstruct = allcombos[i,10],
                   treateffvar = allcombos[i,11],
                   comb_ID = allcombos[i,12],
                   estimates = rownames(output))
    
    all_output = all_output %>%
      rbind(output)
    
    end_time <- format(Sys.time(), "%a %b %d %X %Y %Z")
    timeEnd<-Sys.time()
    durt <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))
    
    tim_df <- cbind(T = allcombos[i,1],
                    nbatch = allcombos[i,2],
                    nolap = allcombos[i,3],
                    K = allcombos[i,4],
                    M = allcombos[i,5],
                    sharedtime = allcombos[i,6],
                    effsize = allcombos[i,7],
                    ICC = allcombos[i,8],
                    CAC = allcombos[i,9],
                    corrstruct = allcombos[i,10],
                    treateffvar = allcombos[i,11],
                    comb_ID = allcombos[i,12],
                    seed = as.character(allcombos[i, 12]), 
                    start_time = start_time,
                    end_tim = end_time,
                    duration = durt)
    
    all_tim_df = all_tim_df %>%
      rbind(tim_df)
    
    write.table(contsresults, "all_contsresults_BE&E.csv", sep = ",",
                row.names = FALSE, col.names = i== 1, append = i > 1)
    
    write.table(output, "all_output_BE&E.csv", sep = ",",
                row.names = FALSE, col.names = i== 1, append = i > 1)
    
    write.table(tim_df, "all_tim_df_BE&E.csv", sep = ",",
                row.names = FALSE, col.names = i== 1, append = i > 1)
  }
}

