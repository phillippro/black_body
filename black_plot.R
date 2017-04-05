library(magicaxis)
fname <- paste0('results/AMsamp_',id,'_',Niter,'_fix',ind.fix,'_logFit',logFit,'_sJ',par.mean['sJ'],'_errorEqual',error.equal,'_prior',prior)
pdf.name <- paste0(fname,'.pdf')
obj.name <- paste0(fname,'_conv',conv,'.Robj')
####plot
pdf(pdf.name,width=12,height=16)
par(mfrow=c(4,3),mar=c(5,5,5,1),cex.axis=1.5,cex.lab=1.5)
####plot data-model
xs <- seq(10,1000,length.out=100)
y <- fs
dy <- dfs
ys <- blackbody(xs,par.opt)
yp <- blackbody(df$lambda,par.opt)
if(logFit){
    ylab='log(Flux [Jy])'
    ys <- exp(ys)
    yp <- exp(yp)
}else{
    ylab='Flux [Jy]'
}
plot(df[,1],y,xlab=expression(lambda*'['*mu*'m]'),ylab=ylab,ylim=range(min(y,ys[xs>200]),max(ys,y)),log='y',xlim=range(xs),yaxt='n')
magaxis(side=2)
arrows(df[,1],y-dy,df[,1],y+dy,length=0.05,angle=90,code=3,col='black')
points(df[,1],yp,pch=20,col='red',cex=1.2)
lines(xs,ys,col='red')
####tracing plot of posterior and likelihood samples
plot(index.mcmc,post.out1,type='l',xlab='Iterations',ylab='Logarithm Posterior')
plot(index.mcmc,loglike.out1,type='l',xlab='Iterations',ylab='Logarithm Likelihood')
####tracing plot for parameters
par.stat <- c()
for(j in 1:length(par.opt)){
    lab <- names(par.up)[j]
    plot(index.mcmc,mcmc.out1[,j],type='l',xlab='Iterations',ylab=lab,main='')
                                        #posterior
    plot(par.bin[,j],post.bin[,j],type='l',ylab='Logarithm Posterior',xlab=lab,main='')
    p2 <- data.distr(x=mcmc.out[,j],xlab=lab,ylab='Frequency',main='')
    par.stat <- cbind(par.stat,p2)
}
dev.off()
dm <- par.opt-par.stat[1,]
dp <- par.stat[2,]-par.opt
par.stat <- rbind(par.opt,dm,dp,par.stat)
colnames(par.stat) <- names(par.opt)
rownames(par.stat) <- c('opt','dminus','dplus','x1per','x99per','mode','mean','sd','skewness','kurtosis')
cat('pdf output:\n',pdf.name,'\n')
###save data
save(list=ls(all=TRUE),file=obj.name)
cat('Robj output:\n',obj.name,'\n\n')
