## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  if(!requireNamespace('devtools')){
#    install.packages('devtools')
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github(repo="ntdyjack/fasthplus",
#                           ref = "main",
#                           build_vignettes = TRUE,
#                           subdir = NULL)

## -----------------------------------------------------------------------------
library(fasthplus)

## ----motivating-figure, echo=FALSE, fig.cap='[**Figure 1**: The $G_{+}$ discordance metric is biased as function of the proportion of within-cluster distances, which is a function of the group balance.]', out.width = '80%', fig.align='left'----
knitr::include_graphics("motivating_figure.jpeg")

## ----performance-plot, echo=FALSE, fig.cap='[**Figure 2**: Demonstration of computation time (in seconds) for several components of estimating H+]', out.width = '80%', fig.align='left'----
knitr::include_graphics("performance_plot.jpeg")

## ---- fig.width=6,fig.height=4------------------------------------------------
set.seed(1234)
n <- 10000
a <- rnorm(n=n,mean=0.5,sd=1)
b <- rnorm(n=n,mean=-0.5,sd=1)

bins <- seq(min(c(a,b)),max(c(a,b)),length.out=20)
hist(x=a,breaks=bins,main='',xlab='',ylab='Frequency',plot=T, border='blue',col='#0000ff64', freq=T)
hist(x=b,breaks=bins,add=T,border='red',col='#ff000064',freq=T)
legend('topright',legend=c("A","B"), pch=c(22,22),
  col= c('blue','red'),cex=1.5, pt.bg=c('#0000ff64','#ff000064'),bty='n')

## -----------------------------------------------------------------------------
hpe(A=a,B=b,p=1001) # A, B formulation

## ----fig.width=6,fig.height=4-------------------------------------------------
n <- 1000
m <- 100
cl1 <- sapply(1:n, function(i) rnorm(n=m,mean=0.5,sd=1))
cl2 <- sapply(1:n, function(i) rnorm(n=m,mean=-0.5,sd=1))
dat <- t(cbind(cl1,cl2))
d <- dist(dat)
dvec <- as.matrix(d)
dvec <- dvec[upper.tri(dvec)]
l <- c(rep(0,n),rep(1,n))
ind <- sapply(l, function(x) x==l)
ind <- ind[upper.tri(ind)]
iw <- which(ind)
ib <- which(!ind)
dw <- dvec[iw]
db <- dvec[ib]

bins <- seq(min(dvec),max(dvec),length.out=20)
hist(x=dw,breaks=bins,main='',xlab='',ylab='Frequency',plot=T, border='blue',col='#0000ff64', freq=T)
hist(x=db,breaks=bins,add=T,border='red',col='#ff000064',freq=T)
legend('topright',legend=c(expression('D'[W]), expression('D'[B])), pch=c(22,22),
  col= c('blue','red'),cex=1.5, pt.bg=c('#0000ff64','#ff000064'),bty='n')

## -----------------------------------------------------------------------------
hpe(D=d,L=l,p=1001) # D, L formulation
hpb(D=dat,L=l,r=100,t=100) # bootstrap

## -----------------------------------------------------------------------------
library(clusterCrit)
intCriteria(traj=dat, part=as.integer(l), crit='G_plus')

## ----fig.width=6,fig.height=4-------------------------------------------------
n <- 1000
m <- 100
dat <- t(sapply(1:n, function(i) rnorm(n=m,mean=0,sd=1)))
d <- dist(dat)
pc <- prcomp(dat)$x[,1:2]
l <- round(runif(n=n))
cols <- ifelse(l==1,'#0000ff64','#ff000064')
plot(x=pc[,1],y=pc[,2],pch=16,col=cols,cex=0.7,xaxs = "i",yaxs = "i",xlab='PC1',ylab='PC2',xaxt='n',yaxt='n')
legend('topleft',legend=c(expression('Cl'[1]), expression('Cl'[2])), pch=c(21,21),
  col= c('blue','red'),cex=1.5, pt.bg=c('#0000ff64','#ff000064'),bty='n')

## -----------------------------------------------------------------------------
hpe(D=d,L=l) # D,L formulation
intCriteria(traj=dat,part=as.integer(l),crit='G_plus') #G+ value
hpb(D=dat,L=l,r=30,t=30) # bootstrap

## ----fig.width=6,fig.height=4-------------------------------------------------
n <- 500
cl1 <- sapply(1:n, function(i) rnorm(n=m,mean=0.5,sd=1))
cl2 <- sapply(1:n, function(i) rnorm(n=m,mean=-0.5,sd=1))
dat <- t(cbind(cl1,cl2))
d <- dist(dat)
l <- c(rep(0,n),rep(1,n))
pc <- prcomp(dat)$x[,1:2]
cols <- ifelse(l==1,'#0000ff64','#ff000064')
plot(x=pc[,1],y=pc[,2],pch=16,col=cols,cex=0.7,xaxs = "i",yaxs = "i",xlab='PC1',ylab='PC2',xaxt='n',yaxt='n')
legend('top',legend=c(expression('Cl'[1]), expression('Cl'[2])), pch=c(21,21),
  col= c('blue','red'),cex=1.5, pt.bg=c('#0000ff64','#ff000064'),bty='n')

## -----------------------------------------------------------------------------
hpe(D=d,L=l,p=1001) # D,L formulation
hpb(D=dat,L=l,r=100,t=100) # bootstrap
intCriteria(traj=dat,part=as.integer(l),crit='G_plus')

