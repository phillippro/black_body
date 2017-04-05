library(foreach)
library(doMC)
#library(doParallel)
#library(parallel)
source('func_post.R')
lambda <- c(70,160,250,350,500,870)#micrometer
#tab <- read.table('all_matchcornish_flux_higal5bands.txt',header=TRUE)
tab <- read.table('data/test.txt',header=TRUE)
Nc <- ncol(tab)
#####set up
c <- 299792458#m/s
tau.up <- 0.3
ind.fix <- 5
par.mean <- c(logA=log(1e-11),T=50,beta=2,tau0=0,sJ=0)
par.sd <- c(logA=0.2,T=20,beta=0.3,tau0=0.05,sJ=2)
par.low <- c(logA=log(1e-13),T=0,beta=1,tau0=1e-5,sJ=0)##lower
par.up <- c(logA=log(1e-8),T=200,beta=5,tau0=tau.up,sJ=50)##upper
cov.start <- 1e-6*diag(abs(par.mean))
Niter <- 1e5
Ncores <- 4
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
frac.adapt <- 0.9
n0 = max(round(Niter/Ncores*(1-frac.adapt)+2),10)#where the adaptive cov. used#at least start from the second interation
Npar <- length(par.up)
eps <- 1e-6
Sd <- 2.4^2/Npar#hyper parameter for covariance
Nupdate <-  1#adaptive frequency
verbose <- FALSE
logFit <- FALSE
error.equal <- FALSE
prior <- 'norm'

#####loop
for(kk in 1:nrow(tab)){
#for(kk in 1){
####read data
    id <- as.character(tab[kk,1])
    fs <- as.numeric(tab[kk,c(2,4,6,8,10,12)])
    dfs <- as.numeric(tab[kk,c(2,4,6,8,10,12)+1])
    dlogfs <- dfs/fs
    logfs <- log(fs)
    if(error.equal){
        dfs <- rep(mean(dfs),length(dfs))
        dlogfs <- rep(mean(dlogfs),length(dlogfs))
    }
    if(logFit){
        df <- data.frame(lambda=lambda,logfs=logfs,dlogfs=dlogfs)
    }else{
        df <- data.frame(lambda=lambda,fs=fs,dfs=dfs)
    }
    par.up['sJ'] <- 100*max(df[,3])
#    par.low['sJ'] <- -100*max(df[,3])

#####AM posterior sampling
    for(jj in 1:2){
        Nsamp <- Niter/Ncores
        if(jj==1){
            out <- AMH(data=df,startvalue=par.mean,cov.start=cov.start,nburn=round(Nsamp/2),Ns=Nsamp,tem=0.1)$out
#            out <- AMsamp(data=df,startvalue=par.mean,cov.start=cov.start,iterations=Niter,tem=0.1)
        }else{
            out <- AMH(data=df,startvalue=par.opt,cov.start=cov.start,nburn=round(Nsamp/2),Ns=Nsamp,tem=1)$out
#            out <- AMsamp(data=df,startvalue=par.opt,cov.start=cov.start,iterations=Niter,tem=1)
        }
        mcmc.out <- out[,1:Npar]
        post.out <- out[,Npar+1]
        loglike.out <- out[,Npar+1]
#        Nburn <- round(Niter/2)
#    mcmc.out <- out$chain[-(1:Nburn),]
#    post.out <- out$post[-(1:Nburn)]
#    loglike.out <- out$like[-(1:Nburn)]
###optimal parameters
        ind <- which.max(post.out)
        par.opt <- mcmc.out[ind,]
        names(par.opt) <- names(par.mean)
#        cat('par.opt=',par.opt,'\n')
    }
    if(ind.fix>0){
        mcmc.out <- mcmc.out[,-ind.fix]
        par.opt <- par.opt[-ind.fix]
    }
####test convergence
    Nsub <- 3##number of subchains
    chain.len <- nrow(mcmc.out)
    subchain.len <- floor(chain.len/Nsub)
    var.sub <- array(data=NA,dim=c(Nsub,ncol(mcmc.out)))
    mean.sub <- array(data=NA,dim=c(Nsub,ncol(mcmc.out)))
    meanvar <- rep(NA,ncol(mcmc.out))
    mean.all <- rep(NA,ncol(mcmc.out))
    var.single.est <- rep(NA,ncol(mcmc.out))
    for(j in 1:ncol(mcmc.out)){ 
        for(i in 1:Nsub){
            var.sub[i,j] <- var(mcmc.out[((i-1)*subchain.len+1):(i*subchain.len),j])
            mean.sub[i,j] <- var(mcmc.out[((i-1)*subchain.len+1):(i*subchain.len),j])
        }
        meanvar[j] <- mean(var.sub[,j])
        mean.all[j] <- mean(mean.sub[,j])
        var.single.est[j] <- subchain.len/Nsub*sum((mean.sub[,j]-mean.all[j])^2)
    }
    var.est <- (subchain.len-1)/subchain.len*meanvar+1/subchain.len*var.single.est
    R <- sqrt(var.est/meanvar)
    if(any(is.na(R))) R[is.na(R)] <- 1
    if(any(R>1.1)){
        conv <- FALSE
    }else{
        conv <- TRUE
    }
#    cat('convergence:',conv,'\n')

###parameters for binning/thining
#    lams <- seq(min(lambda),max(lambda),by=0.1)
    Nsep <- max(round(Niter/10000),1)
    fchop <- 0#0.9
    Nbins = max(min(round(nrow(mcmc.out)/1e2),1000),10)
    Nbreaks = 100#for histograms

####binning
    par.bin <- post.bin <- like.bin <- array(data=NA,dim=c(Nbins,ncol(mcmc.out)))
    for(i in 1:ncol(mcmc.out)){
        tmp <- binning.post2(mcmc.out[,i],post.out,loglike.out,Nbins)
        par.bin[,i] <- tmp$pars
        like.bin[,i] <- tmp$likes
        post.bin[,i] <- tmp$posts
        if(length(tmp$ind.na)>0){
            ind.na <- tmp$ind.na
                                        #cat('the',ind.na,'th bin of ',colnames(mcmc.out)[i],' has no samples!\n')
            like.bin[ind.na,i] <- min(like.bin[-ind.na,i])
            post.bin[ind.na,i] <- min(post.bin[-ind.na,i])
        }
    }

####thining
    if(Nsep>1){
        ind.out <- Nsep*(1:floor(nrow(mcmc.out)/Nsep))
        mcmc.out1 <- mcmc.out[ind.out,]
        post.out1 <- post.out[ind.out]
        loglike.out1 <- loglike.out[ind.out]
        index.mcmc <- ind.out
    }else{
        mcmc.out1 <- mcmc.out
        post.out1 <- post.out
        loglike.out1 <- loglike.out
        index.mcmc <- 1:length(post.out)
    }
    source('black_plot.R',local=TRUE)
}
