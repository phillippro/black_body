library(mvtnorm)
library(MASS)
library(e1071)

blackbody <- function(lambda,pars){
    c <- 299792458#m/s
    nu <- c/(lambda*1e-6)
    logA <- pars['logA']
    T <- pars['T']
    beta <- pars['beta']
    tau0 <- pars['tau0']
    h <- 6.62607004e-34
    k <- 1.38064852e-23
    nu0 <- 1e6*c/850
    tau <- tau0*(nu/nu0)^beta
    S <- 1e26*exp(logA)*(1-exp(-tau))*2*h*nu^3/c^2/(exp(h*nu/(k*T))-1)
    if(logFit) S <- log(S)
    return(S)
}
covariance.n0 <- function(mat){
    cov.par <- Sd*cov(mat)+Sd*eps*diag(ncol(mat))
    return(cov.par)
}

covariance.rep <- function(parms,cov1,mu1,n){
    if(n%%Nupdate==0){
        Sd <- 2.4^2/Npar
        mu2 <- (mu1*(n-1)+parms)/n#mu1 and mu2 are the mean vectors of parameters for n-1 and n iterations
        N <- n-1
        cov2 <- (N-1)*cov1/N+Sd*(N*mu1%*%t(mu1)-(N+1)*mu2%*%t(mu2)+parms%*%t(parms)+eps*diag(length(parms)))/N
    }else{
        cov2 <- cov1
    }
    return(cov2)
}

prior.func <- function(pars){
    logA <- pars['logA']
    T <- pars['T']
    beta <- pars['beta']
    tau0 <- pars['tau0']
    sJ <- pars['sJ']
    logprior <- log(1/(par.up['logA']-par.low['logA']))#log uniform
    if(prior=='norm'){
        logprior <- logprior+log(dnorm(T,par.mean['T'],par.sd['T']))
        logprior <- logprior+log(dnorm(beta,par.mean['beta'],par.sd['beta']))
        logprior <- logprior+log(dnorm(tau0,par.mean['tau0'],par.sd['tau0']))
        logprior <- logprior+log(dnorm(sJ,par.mean['sJ'],par.sd['sJ']))
    }else{
        logprior <- logprior+log(1/(par.up['T']-par.low['T']))
        logprior <- logprior+log(1/(par.up['beta']-par.low['beta']))
        logprior <- logprior+log(1/(par.up['tau0']-par.low['tau0']))
        logprior <- logprior+log(1/(par.up['sJ']-par.low['sJ']))
    }
#
    return(logprior)
}

loglikelihood <- function(data,pars){ 
    lambda <- data[,1]
    if(logFit){
        flux <- data$logfs
        dflux <- data$dlogfs
    }else{
        flux <- data$fs
        dflux <- data$dfs
    }
    sJ <- pars['sJ']
    bb <- blackbody(lambda,pars)
    logLike <- sum(log(1/sqrt(2*pi)/dflux))-sum((flux-bb)^2/(2*(dflux^2+sJ^2)))
#    cat('logLike1=',sum(log(1/sqrt(2*pi)/dflux)),'\n')
#    cat('logLike2=',sum((flux-bb)^2/(2*dflux^2)),'\n')
    return(list(loglike = logLike))
}

posterior <- function(data,param,tem){
    llike <- loglikelihood(data,param)
    pr <- prior.func(param)
    post <- llike$loglike*tem+pr
    return(list(loglike=llike$loglike,logprior=pr,post=post,ind.neg=llike$ind.neg))
}

proposalfunction <- function(param,startvalue,cov.adapt){
    Ntt <- 1000
    for(k in 1:Ntt){
        cov.small <- eps*diag(Npar)
        cov.all <- cov.adapt+cov.small
        param.new <- try(mvrnorm(n=1,mu=param,Sigma=cov.all,tol=1e-10,empirical=FALSE))#this could cause some non-zero values for fix value, but it is too small to be accounted for.
#        cat('param.new=',param.new,'\n')
#        cat('param.new[k]=',param.new['K1'],'\n')
        if(class(param.new)=='try-error'){
            param.new <- try(mvrnorm(n=1,mu=param,Sigma=cov.small,tol=1e-10,empirical=FALSE))
        }
        if(all(param.new>par.low & param.new<par.up)) break()
    }
    names(param.new) <- names(startvalue)
    return(param.new)
}

AMsamp <- function(data,startvalue,cov.start,iterations,tem=1){
    Npar <- length(startvalue)		    
    chain = array(dim=c(iterations+1,Npar))
    post = loglike = rep(NA,iterations+1)
    colnames(chain) <- names(startvalue)
    chain[1,]<- startvalue
    mu1 <- chain[1,]
    post.out = posterior(data,chain[1,],tem)
    post[1] <- post.pre <- post.out$post
    loglike[1] <- loglike.pre <- post.out$loglike
    logprior.pre <- post.out$logprior
    ind.neg <- ind.neg.pre <- post.out$ind.neg
    dt0 <- 0
    cov.adapt <- array(data=0,dim=c(Npar,Npar))
    for(i in 1:iterations){
        if(i == n0){
            cov.adapt <- covariance.n0(mat=chain[1:n0,])
        }else if(i > n0){
            cov.adapt <- covariance.rep(chain[i,],cov.adapt,mu1,i)
        }else{
            cov.adapt <- cov.start
        }
        proposal = proposalfunction(chain[i,],startvalue,cov.adapt)
        if(ind.fix>0){
            proposal[ind.fix] <- startvalue[ind.fix]
        }
        proprop <- posterior(data,proposal,tem)
        post.prop <- proprop$post
	logprior.prop <- proprop$logprior
        loglike.prop <- proprop$loglike
        ind.neg.prop <- proprop$ind.neg

        post.cur <- post.pre
        loglike.cur <- loglike.pre
	logprior.cur <- logprior.pre
        ind.neg.cur <- ind.neg.pre
        if(is.na(post.prop)){
            probab = 0
        }else{
            probab = exp(post.prop-post.cur)
        }
        if(runif(1)<probab){
            chain[i+1,]=proposal
###values for next proposal estimation
            post.pre <- post.prop
            loglike.pre <- loglike.prop
	        logprior.pre <- logprior.prop
            ind.neg.pre <- ind.neg.prop
        }else{
            chain[i+1,]=chain[i,]
            post.pre <- post.cur
            loglike.pre <- loglike.cur
            logprior.pre <- logprior.cur
            ind.neg.pre <- ind.neg.cur
        }
##save values 
        post[i+1]<- logprior.pre + loglike.pre
        loglike[i+1] <- loglike.pre
        mu1 <- (mu1*i+chain[i+1,])/(i+1)
        if(i%%round(Niter/100)==0 & verbose){
            cat('i=',i,'\n')
            cat('variables=',names(par.up),'\n')
            cat('chain[i+1,]=',chain[i+1,],'\n')
            cat('post[i+1,]=',post[i+1],'\n')
            cat('loglike[i+1,]=',loglike[i+1],'\n\n')
        }
    }
    return(list(chain=chain,post=post,like=loglike,cov=cov.adapt))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

data.distr <- function(x,xlab,ylab,main='',oneside=FALSE,plotf=TRUE){
    xs <- seq(min(x),max(x),length.out=1e3)
    fitnorm <- fitdistr(x,"normal")
    p <- hist(x,plot=FALSE)
    xfit <- length(x)*mean(diff(p$mids))*dnorm(xs,fitnorm$estimate[1],fitnorm$estimate[2])
    ylim <- range(xfit,p$counts)
    if(plotf){
        plot(p,xlab=xlab,ylab=ylab,main=main,ylim=ylim)
        lines(xs,xfit,col='red')
    }
    x1=Mode(x)
    x2=mean(x)
    x3=sd(x)
    x4=skewness(x)
    x5=kurtosis(x)
    xs = sort(x)
    x1per = max(min(xs),xs[floor(length(xs)*0.01)])
    x99per = min(xs[ceiling(length(xs)*0.99)],max(xs))
#    abline(v=c(x1per,x99per),col='blue')
    if(plotf){
        if(!oneside){
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3))))),bty='n')
            legend('topright',legend=c(as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }else{
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }
    }
    return(c(x1per=x1per,x99per=x99per,mode=x1,mean=x2,sd=x3,skewness=x4,kurtosis=x5))
}
binning.post2 <- function(par.val,post.val,like.val,Nbins){
    ind <- sort(par.val,index.return=TRUE)$ix
    par.sort <- par.val[ind]
    post.sort <- post.val[ind]
    like.sort <- like.val[ind]
    p1 <- hist(par.sort,breaks=seq(min(par.sort),max(par.sort),length.out=Nbins+1),plot=FALSE)
    index <- c(0,cumsum(p1$counts))
    post.max <- rep(NA,length(p1$mids))
    like.max <- rep(NA,length(p1$mids))
    ind.na <- c()
    for(j in 1:(length(index)-1)){
        if(index[j+1]>index[j]){
            post.max[j] <- max(post.sort[(index[j]+1):index[j+1]])
            like.max[j] <- max(like.sort[(index[j]+1):index[j+1]])
        }
    }
    ind.na <- which(is.na(post.max))
    return(list(likes=like.max,posts=post.max,pars=p1$mids,ind.na=ind.na))
}
###tell foreach how to combine output
comb <- function(x, ...) {  
      mapply(rbind,x,...,SIMPLIFY=FALSE)
}
####parallel computing
AMH <- function(data,startvalue,cov.start,nburn,Ns,tem){
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tmp <- AMsamp(data=data,startvalue=startvalue,cov.start=cov.start,iterations=Ns,tem=tem)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]))
#        cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)])
    }
#    return(out.amh)
    return(list(out=out.amh$out,cov=out.amh$cov))
}
