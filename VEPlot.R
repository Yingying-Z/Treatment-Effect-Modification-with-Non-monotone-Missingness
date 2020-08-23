# Demo code for manuscript Evaluation of Treatment Effect Modification by Biomarkers Measured Pre- and Post-randomization in the Presence of Non-monotone Missingness
# It illustrates how the code could be used to plot the VE curves in Fig 3 in the manuscript.
# Date      Programmer   	
# August12-2020  Y Zhuang      
##########################################################################
rm(list=ls(all=TRUE))
#set your working directory
load("Results.RData")


pdf(file="VEplot.pdf")

lwd=1.5
myylim=c(0,1)
leg.x<-1.5
leg.cex<-0.8

par(mar=c(3,4,3,1),oma=c(1,0,0,0))
layout(matrix(1:4,byrow=T,nrow=2))

dotLoc=c(1,10,50,100,200,300,400,500)
plot(Su,VE.S.marginal,type='n',main="VE curve", xlab="", ylab="", xlim=range(Su), ylim=myylim, xaxt="n", yaxt="n", lwd=2.5)
lines(Su,VE.S.marginal,col=1,lwd=lwd)
points(Su[dotLoc],VE.S.marginal[dotLoc], pch = 16, col = 1,lty = 2, lwd = 1)
lines(Su,VE.S.Bpos,col="firebrick",lwd=lwd)
dotLoc=c(5,15,60,170,350,480)
points(Su[dotLoc],VE.S.Bpos[dotLoc], pch = 15, col = "firebrick", lty = 2, lwd = 1)
lines(Su,VE.S.Bneg,col="darkcyan",lwd=lwd)
dotLoc=c(2,35,120,250,400,450,500)
points(Su[dotLoc],VE.S.Bneg[dotLoc], pch = 17, col = "darkcyan", lty = 2, lwd = 1)
axis(side=1, at=c(log10(5),1:5), labels=expression("<10     "," 10",100,10^3,10^4,10^5))
axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25)),cex.axis=1,las=1)
mtext(expression(bold(paste('Mon13 Average Titer=',s[1],' in Vaccinees'))),side=1, las=0, line=2, cex=0.9) 
mtext("Vaccine Efficacy (%)", side=2, las=0, line=2.2, cex=1)
legend(leg.x,0.2,c("Estimated Marginal VE curve"),bty='n',col=1, pch= 16,lwd=lwd,cex=leg.cex)
legend(leg.x,0.13,c("Estimated BL seropositive VE curve"),bty='n',col="firebrick", pch= 15,lwd=lwd,cex=leg.cex)
legend(leg.x,0.06,c('Estimated BL seronegative VE curve'),bty='n',col="darkcyan", pch= 17,lwd=lwd,cex=leg.cex)


#Marginal VE plot
type<-"marginal"
plot(Su,VE.S.marginal,type='n',main="Marginal", xlab="", ylab="", xlim=range(Su), ylim=myylim, xaxt="n", yaxt="n", lwd=2.5)
lines(Su,VE.S.marginal,col=1,lwd=lwd)
dotLoc=c(1,10,50,100,200,300,400,500)
points(Su[dotLoc],VE.S.marginal[dotLoc], pch = 16, col = 1, lty = 2, lwd = 1)
lines(Su,low.ci.Mar,lty=6, col="darkolivegreen4", lwd=2.5)
lines(Su,high.ci.Mar,lty=6, col="darkolivegreen4", lwd=2.5)
lines(Su,low.cb.Mar, lty=2, col="dark blue", lwd=3.5)
lines(Su,high.cb.Mar, lty=2, col="dark blue", lwd=3.5)
axis(side=1, at=c(log10(5),1:5), labels=expression("<10     "," 10",100,10^3,10^4,10^5))
axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25)),cex.axis=1,las=1)
mtext(expression(bold(paste('Mon13 Average Titer=',s[1],' in Vaccinees'))),side=1, las=0, line=2, cex=0.9)
mtext("Vaccine Efficacy (%)", side=2, las=0, line=2.2, cex=1)

leg.x=2
legend(leg.x,0.2,c("Estimated Marginal VE curve"),bty='n',col=1, pch= 16,lwd=lwd,cex=leg.cex)
legend(leg.x,0.13,lty=6,c('95% pointwise CI'),bty='n',col="darkolivegreen4", lwd=2,cex=leg.cex)
legend(leg.x,0.06,lty=2,c('95% simultaneous CI'),bty='n',col="dark blue", lwd=2.5,cex=leg.cex)


#Bpos VE plot
type<-"Bpos"
plot(Su,VE.S.Bpos,type='n',main="BL seropositive", xlab="", ylab="", xlim=range(Su), ylim=myylim, xaxt="n", yaxt="n", lwd=2.5)
lines(Su,VE.S.Bpos,col="firebrick",lwd=lwd)
dotLoc=c(5,15,60,170,350,480)
points(Su[dotLoc],VE.S.Bpos[dotLoc], pch = 15, col = "firebrick", lty = 2, lwd = 1)
lines(Su,low.ci.Bpos,lty=6, col="darkolivegreen4", lwd=2.5)
lines(Su,high.ci.Bpos,lty=6, col="darkolivegreen4", lwd=2.5)
lines(Su,low.cb.Bpos, lty=2, col="dark blue", lwd=3.5)
lines(Su,high.cb.Bpos, lty=2, col="dark blue", lwd=3.5)
axis(side=1, at=c(log10(5),1:5), labels=expression("<10     "," 10",100,10^3,10^4,10^5))
axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25)),cex.axis=1,las=1)

mtext(expression(bold(paste('Mon13 Average Titer=',s[1],' in Vaccinees'))),side=1, las=0, line=2, cex=0.9)
mtext("Vaccine Efficacy (%)", side=2, las=0, line=2.2, cex=1)

leg.x=2
leg.cex=0.75
legend(leg.x,0.2,c("Estimated BL seropositive VE curve"),bty='n',col="firebrick", pch= 15,lwd=lwd,cex=leg.cex)
legend(leg.x,0.13,lty=6,c('95% pointwise CI'),bty='n',col="darkolivegreen4", lwd=2,cex=leg.cex)
legend(leg.x,0.06,lty=2,c('95% simultaneous CI'),bty='n',col="dark blue", lwd=2.5,cex=leg.cex)

#Bneg VE plot
type<-"Bneg"
plot(Su,VE.S.Bneg,type='n',main="BL seronegative", xlab="", ylab="", xlim=range(Su), ylim=myylim, xaxt="n", yaxt="n", lwd=2.5)
lines(Su,VE.S.Bneg,col="darkcyan",lwd=lwd)
dotLoc=c(2,35,120,250,400,450,500)
points(Su[dotLoc],VE.S.Bneg[dotLoc], pch = 17, col = "darkcyan", lty = 2, lwd = 1)
lines(Su,low.ci.Bneg,lty=6, col="darkolivegreen4", lwd=2.5)
lines(Su,high.ci.Bneg,lty=6, col="darkolivegreen4", lwd=2.5)
lines(Su,low.cb.Bneg, lty=2, col="dark blue", lwd=3.5)
lines(Su,high.cb.Bneg, lty=2, col="dark blue", lwd=3.5)
axis(side=1, at=c(log10(5),1:5), labels=expression("<10     "," 10",100,10^3,10^4,10^5))
axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25)),cex.axis=1,las=1)

mtext(expression(bold(paste('Mon13 Average Titer=',s[1],' in Vaccinees'))),side=1, las=0, line=2, cex=0.9)
mtext("Vaccine Efficacy (%)", side=2, las=0, line=2.2, cex=1)

leg.x=1.9
leg.cex=0.75
legend(leg.x,0.16,c('Estimated BL seronegative VE curve'),bty='n',col="darkcyan", pch= 17,lwd=lwd,cex=leg.cex)
legend(leg.x,0.1,lty=6,c('95% pointwise CI'),bty='n',col="darkolivegreen4", lwd=2,cex=leg.cex)
legend(leg.x,0.04,lty=2,c('95% simultaneous CI'),bty='n',col="dark blue", lwd=2.5,cex=leg.cex)

dev.off()




