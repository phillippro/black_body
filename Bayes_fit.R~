c <- 299792458#m/s
model <- function(nu,A,T,beta,tau0,logFit=FALSE){
    h <- 6.62607004e-34
    k <- 1.38064852e-23
    nu0 <- 1e6*c/850
    tau <- tau0*(nu/nu0)^beta
    S <- 1e26*A*(1-exp(-tau))*2*h*nu^3/c^2/(exp(h*nu/(k*T))-1)
    if(logFit) S <- log(S)
    return(S)
}
#par.low <- c(A=1e-13,T=20,beta=1.8,nu0=1e6*c/1000)##lower
tau.up <- 0.3
par.low <- c(A=1e-13,T=0,beta=1,tau0=1e-5)##lower
#par.up <- c(A=1e-8,T=200,beta=2.2,nu0=1e6*c/10)##upper
par.up <- c(A=1e-8,T=200,beta=3,tau0=tau.up)##upper
#par.mean <- c(A=1e-12,T=100,beta=1.8,nu0=1e6*c/100)
par.mean <- c(A=1e-11,T=50,beta=1.8,tau0=0.1)
par.sd <- c(0.2,10,0.2,0.05)
Nsamp <- 1e6
#nls.lm(par,lower=par.low,upper=par.up,fn=model)
lambda <- c(70,160,250,350,500)*1e-6#m
###parameters
prior <- 'uniform'
#prior <- 'norm'
logFit <- FALSE
#logFit <- TRUE
##
###read data
#tab <- read.table('all_matchcornish_flux_higal5bands.txt',header=TRUE)
tab <- read.table('data/test.txt',header=TRUE)
for(kk in 1:nrow(tab)){
ind <- kk
Nc <- ncol(tab)
id <- as.character(tab[ind,1])
#cat('id=',id,'\n')
fs <- as.numeric(tab[ind,c(2,4,6,8,10,12)])
dfs <- as.numeric(tab[ind,c(2,4,6,8,10,12)+1])
#df <- data.frame(nu=c/ls,fs=c(323.913,561.07,371.956,186.321,75.884),dfs=c(0.595,1.209,1.324,1.287,1.432))
df <- data.frame(nu=c/lambda,fs=fs,dfs=dfs)
x <- df[,1]
if(logFit){
    y <- log(df[,2])
    dy <- df[,3]/df[,2]
}else{
    y <- df[,2]
    dy <- df[,3]
}
##
loglikes <- rep(NA,Nsamp)
if(prior=='uniform'){
    A <- exp(runif(Nsamp,log(par.low[1]),log(par.up[1])))
    T <- runif(Nsamp,par.low[2],par.up[2])
    beta <- runif(Nsamp,par.low[3],par.up[3])
#    nu0 <- exp(runif(Nsamp,log(par.low[4]),log(par.up[4])))
    tau0 <- runif(Nsamp,par.low[4],par.up[4])
}else{
    Ntry <- 100
    A <- T <- beta <- tau0 <- c()
    A0 <- exp(rnorm(Nsamp,log(par.mean[1]),par.sd[1]))
    T0 <- rnorm(Nsamp*10,par.mean[2],par.sd[2])
    beta0 <- rnorm(Nsamp*10,par.mean[3],par.sd[3])
    tau0 <- rnorm(Nsamp*10,par.mean[4],par.sd[1])
    A <- A0[A0>par.low[1] & A0<par.up[1]][1:Nsamp]
    T <- T0[T0>par.low[2] & T0<par.up[2]][1:Nsamp]
    beta <- beta0[beta0>par.low[3] & beta0<par.up[3]][1:Nsamp]
    tau0 <- tau0[tau0>par.low[4] & tau0<par.up[4]][1:Nsamp]
}
for(j in 1:Nsamp){
    yp <- model(x,A[j],T[j],beta[j],tau0[j],logFit=logFit)
    loglikes[j] <- 1/sqrt(2*pi)*sum(1/dy)-sum((y-yp)^2/(2*dy^2))
}

###binning
llb <- c()
Tb <- Ab <- betab <- tau0b <- c()
Nbin <- max(50,Nsamp/1e5)
llb <- array(data=NA,dim=c(Nbin,4))
for(j in 1:Nbin){
    dlogA <- (max(log(A))-min(log(A)))/Nbin
    inds <- which(A>exp(min(log(A))+(j-1)*dlogA) & A<exp(min(log(A))+j*dlogA))
    llb[j,1] <- max(loglikes[inds])
    Ab <- c(Ab,median(A[inds]))

    dT <- (max(T)-min(T))/Nbin
    inds <- which(T>(min(T)+(j-1)*dT) & T<(min(T)+j*dT))
    llb[j,2] <- max(loglikes[inds])
    Tb <- c(Tb,median(T[inds]))

    dbeta <- (max(beta)-min(beta))/Nbin
    inds <- which(beta>(min(beta)+(j-1)*dbeta) & beta<(min(beta)+j*dbeta))
    llb[j,3] <- max(loglikes[inds])
    betab <- c(betab,median(beta[inds]))

#    dlogtau0 <- (log(max(tau0))-log(min(tau0)))/Nbin
    dtau0 <- (max(tau0)-min(tau0))/Nbin
    inds <- which(tau0>(min(tau0)+(j-1)*dtau0) & tau0<(min(tau0)+j*dtau0))
    llb[j,4] <- max(loglikes[inds])
    tau0b <- c(tau0b,median(tau0[inds]))
}


####plot
fname <- paste0('results/Bayes_',id,'_',prior,'_',logFit,'_',Nsamp,'.pdf')
pdf(fname,9,6)
par(mfrow=c(2,3),cex=1.2,cex.lab=1.2,cex.axis=1.2,mar=c(4,4,1,1))
###like distr.
plot(Ab,llb[,1],type='l',xlab='A',ylab='Log like',log='x')
plot(Tb,llb[,2],type='l',xlab='T',ylab='Log like')
plot(betab,llb[,3],type='l',xlab='beta',ylab='Log like')
plot(tau0b,llb[,4],type='l',xlab='tau0',ylab='Log like')
#####best fit
ind <- which.max(loglikes)
yp <- model(x,A[ind],T[ind],beta[ind],tau0[ind],logFit=logFit)
xs <- seq(min(x),max(x),length.out=100)
ys <- model(xs,A[ind],T[ind],beta[ind],tau0[ind],logFit=logFit)
if(logFit){
    ylab <- 'log(Flux/Jy)'
}else{
    ylab <- 'Flux/Jy'
}
#xlim=1e9*c(1e2,1e5)
#log='x'
lam <- c/x*1e6#micro
lams <- c/xs*1e6
xlim <- range(lam)
plot(lam,y,xlab='wavelength[micro]',ylab=ylab,xlim=xlim,ylim=range(y,ys))
points(lam,yp,col='blue',pch=20,cex=0.5)
lines(lams,ys,col='red')
#plot(A,loglikes,xlab='A',ylab='log Likelihood')
#plot(T,loglikes,xlab='T',ylab='log Likelihood')
#plot(beta,loglikes,xlab='beta',ylab='log Likelihood')
#plot(tau0,loglikes,xlab='tau0',ylab='log Likelihood')
cat('output:\n')
cat(fname,'\n')
cat('A,T,beta,tau0=',A[ind],T[ind],beta[ind],tau0[ind],'\n\n')
#cat('lambda',c*1e6/tau0[ind],'micro\n\n')
dev.off()
}
