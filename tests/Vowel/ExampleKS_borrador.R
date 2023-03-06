rm(list=ls())
coordins=read.table("coordin.txt",head=T,dec=",")
subject1=read.table("subjet1Epsd.txt",head=T,dec=",")
subjet1_Spatfd=SpatFD(subject1,coords = coordins,basis="Bsplines",nbasis=25,lambda=0,nharm = 1)
summary(subjet1_Spatfd)
scores(subjet1_Spatfd)
nuevos=expand.grid(x=seq(0,80,len=3),y=seq(0,48,len=4))
modelo1=list(vgm(4.2,"Exp",12))
FunctKrig=KS_scores_lambdas(subjet1_Spatfd, nuevos, model=modelo1, method = "lambda", name=NULL,fill.all=NULL)
#
ggplot_KS(FunctKrig, show.varpred = T, main = "Functional Data", main2 = "Functional Data", ylab = "Value", xlab = "Time", ndigits = 3)



sujeto1SCORES=data.frame(scores(sujetoSpatfd)
coordinates(sujeto1SCORES)=c("x","y")   
sujeto1g=gstat(,id="X1", formula = X1~1, data = sujeto1SCORES)

gstat::fit.variogram(v[[i]], model[[i]], fit.method = 6)