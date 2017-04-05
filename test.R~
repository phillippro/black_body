library(minpack.lm)
model <- function(nu,A,T,beta,nu0){
#    A <- par$A
#    nu <- par$nu
#    T <- par$T
#    beta <- par$beta
#    nu0 <- par$nu0
    h <- 6.62607004e-34
    k <- 1.38064852e-23
#    h <- 6.62606885e-27
#    k <- 1.3806504e-16
    tau <- (nu/nu0)^beta
    e <- exp(h*nu/(k*T))
#    cat('e=',e,'\n')
#    tmp <- log(1-exp(-tau))+3*log(nu)-log(e-1)
    tmp <- (1-exp(-tau))*nu^3/(e-1)
#    s <- h*nu/(k*T)
#    cat('s=',s,'\n')
#    tmp <- -tau+3*log(nu)-s
#    cat('tmp=',tmp,'\n')
#    S <- tmp+A
    S <- A*tmp
    return(S)
}

par <- list(A=exp(-50),T=100,beta=1.5,nu0=1e12)
par.low <- c(A=exp(-100),T=0,beta=0,nu0=1e11)
par.up <- c(A=exp(0),T=1e3,beta=3,nu0=1e13)
#nls.lm(par,lower=par.low,upper=par.up,fn=model)
ls <- c(70,160,250,350,500)*1e-6#micro
fs <- c(114.311,14.958,3.456,2.219,0.936)
dfs <- c(0.339,0.246,0.195,0.256,0.194)
dlogfs <- dfs/fs
c <- 299792458#m/s
#df <- data.frame(nu=c/ls,logfs=log(fs))
df <- data.frame(nu=c/ls,fs=fs)
lm <- nlsLM(fs~model(nu,A,T,beta,nu0),data=df,start=par,lower=par.low,upper=par.up,trace=T,control=nls.lm.control(maxiter=500,ptol=1e-26,ftol=1e-26,maxfev=1e3),weights=1/dfs^2)
nu <- df$nu
#logfs <- df$logfs
fs <- df$fs
pdf('test.pdf',4,4)
nus <- seq(min(nu),max(nu),by=min(nu)/10)
A <- coef(lm)['A']
nu0 <- coef(lm)['nu0']
beta <- coef(lm)['beta']
T <- coef(lm)['T']
#plot(nu,logfs,xlab='Hz',ylab='log Flux',ylim=range(logfs,model(nus,A,T,beta,nu0)))
plot(nu,fs,xlab='Hz',ylab='log Flux',ylim=range(logfs,model(nus,A,T,beta,nu0)),log='x')
points(nu,fitted(lm),col='red')
lines(nus,model(nus,A,T,beta,nu0),col='blue')
dev.off()
