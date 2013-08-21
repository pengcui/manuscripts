
#------------------------------------------------------------
  
# 10-fold cross validation in svr 
svrcrossvalidation=function(datamatrix,featurevector){     # most modify
  library('e1071')
  kk=c()
  actual=c()
  pre=c()
  n=nrow(datamatrix)/10
  for (i in 1:10){
    if (length(kk) > 0){
      se=sample((1:nrow(datamatrix))[-kk],n)
    }else{
      se=sample((1:nrow(datamatrix)),n)
    }
    te=datamatrix[se,]
    tr=datamatrix[-se,]
    mysvr=svm(tr[,featurevector],log(tr[,1],2),probability=T)
    mysvrpre=predict(mysvr,te[,featurevector])
    actual=c(actual,log(te[,1],2))
    pre=c(pre,mysvrpre)
    kk=c(kk,se)
  }
  myresult=list(actual,pre)
  return(myresult)
}

#-----------------------------------------------------------
samplesize=seq(250,6000,250)
sampleresult=c()
# choose the sample size
data=read.table('/home/ckivip/mywork/Dnase1/zlater/max/data_max/zzzzz',header=T)
alldata=as.matrix(data[,4:41])
alldata=alldata[1:84830,]
for (i in samplesize){
  aa=sample(1:nrow(alldata),i)
  mydata=alldata[aa,]
  mylist=svrcrossvalidation(mydata,c(2:6,8:35))
  mycorr = round(cor(mylist[[1]], mylist[[2]],method=c('spearman')),3)  # compute the correlation
  sampleresult=c(sampleresult,mycorr)
}
sampleresult
plot(seq(250,6000,250),sampleresult,type='o',xlim=c(0,6250),ylim=c(0.6,0.85),xlab='samplesize',ylab='')
mtext('SCC',side=2,line=2)
#-------------------------------------------------------------
  
# the linear model
data=read.table('/home/ckivip/mywork/Dnase1/zlater/max/data_max/zzzzzz_max',header=T)
#head(data)
alldata=data[,4:41]
#head(alldata)
#alldata=alldata[1:84830,]
#aa=sample(1:nrow(alldata),3000)
#mydata=alldata[aa,]
# 10-fold cross validation in lm
mydata=alldata
kk=c()
actual=c()
pre=c()
datamatrix=mydata
n=nrow(datamatrix)/10
for (i in 1:10){
  if (length(kk) > 0){
    se=sample((1:nrow(datamatrix))[-kk],n)
  }else{
    se=sample((1:nrow(datamatrix)),n)
  }
  te=datamatrix[se,]
  tr=datamatrix[-se,]
  #the HM
  #mylm=lm(log(dnase1sig,2)~H3k9me3+H3k27ac+H3k4me3+H3k4me1+H3k4me2+H3k27me3+H3k36me3+H3k9ac+H4k20me1+H3k79me2, data=tr)
  #the TF
  mylm=lm(log(dnase1sig,2)~CTCF+BRCA1+BACH1+CEBPB+CHD1+CHD2+CJUN+GTF2F1+CMYC+CTBP2+JUND+MXI1+MAX+MAFK+RFX5+RAD21+NRF1+SUZ12+SIN3A+TBP+ZNF143+USF2+ZNF274,data=tr)
  # all features
  #mylm=lm(log(dnase1sig,2)~H3k9me3+H3k27ac+H3k4me3+H3k4me1+H3k4me2+H3k27me3+H3k36me3+H3k9ac+H4k20me1+H3k79me2+CTCF+BRCA1+BACH1+CEBPB+CHD1+CHD2+CJUN+GTF2F1+CMYC+CTBP2+JUND+MXI1+MAX+MAFK+RFX5+RAD21+NRF1+SUZ12+SIN3A+TBP+ZNF143+USF2+ZNF274,data=tr)
  #histone=c(2:6,8:12)
  tf=c(13:35)
  #allfeature=c(histone,tf)
  mylmpre=predict(mylm,data.frame(te[,tf]))
  actual=c(actual,log(te[,1],2))
  pre=c(pre,mylmpre)
  kk=c(kk,se)
}
mycorr = round(cor(actual,pre,method=c('spearman')),2)  # compute the correlation
mycorr  

  
--------------------------------------------------------------
  
# the svr model and plot the correlation
data=read.table('/home/ckivip/mywork/Dnase1/zlater/max/data_max/zzzzzz_max',header=T)
#head(data)
alldata=as.matrix(data[,4:41])

#head(alldata)
#alldata=alldata[1:84830,]
#aa=sample(1:nrow(alldata),5000)
#mydata=alldata[aa,]
    
    #mydata=data[aa,]
    #write.table(mydata,file='/home/ckivip/mywork/Dnase1/zlater/data_avg/featuresig/zzzzzz',sep='\t',row.names=F,quote=F)
mydata=alldata
histone=c(2:6,8:12)
tf=c(13:35)
allfeature=c(histone,tf)
mylist=svrcrossvalidation(mydata,c(32))
mycorr = round(cor(mylist[[1]], mylist[[2]],method=c('spearman')),2)  # compute the correlation
mycorr

# no meaning ! using the particular feature
allscc=c()
for (i in histone){
  print(i)
  onelist=svrcrossvalidation(mydata,c(i))
  onecorr=round(cor(onelist[[1]], onelist[[2]],method=c('spearman')),2)
  allscc=c(allscc,onecorr)
}

---------------------------------------------------------------------

# plot the figure
dev.off()
par(mfrow=c(1,2),mar=c(5,4,3,1))
# plot the correlation using all the features
plot(mylist[[1]],mylist[[2]], pch=20, col="blue", cex=0.6,xlab='',ylab='')
title(xlab="The signal of Dnase-seq (log2)", ylab="Predicted Signal",mgp=c(2.3,1,0))
mylm = lm(mylist[[1]]~mylist[[2]])
abline(mylm, lwd=2)
text(0.05,-4.5,paste("SCC=", mycorr, sep=""),cex=0.8)
mtext('(a)',side=3,at=-5.5,adj=0,cex=1.5)
# barplot for the only one feature
#names(allscc)=names(data)[5:15]
#barplot(sort(allscc,T),ylim=c(0,0.7),xlab='',ylab='',axisnames=F)
featurepowers=read.table('/home/ckivip/mywork/Dnase1/zlater/max/data_max/parallel/one1',header=F)
histonevec=c(1,2,4:6,8,10,26,29,30)
tfvec=c(3,7,9,11:14,16:22,24,27,28,31:34,36,37)
barplot(rev(featurepowers[tfvec,2]),ylim=c(0,0.65),xlab='',ylab='',axisnames=F)
# HM: at=seq(0.6,12,1.2) ,TF: at=seq(0.7,28,1.2)
#mtext(names(sort((allscc),T)),line=0.5,side=1,at=seq(0.7,31,1.2),las=2,cex=0.8)

mtext(rev(featurepowers[tfvec,1]),line=0.5,side=1,at=seq(0.7,28,1.2),las=2,cex=0.8)
mtext('SCC',side=2,line=2.3)
mtext('(b)',side=3,at=-1,cex=1.5)

---------------------------------------------------------------------
  
# plot the one, two. three models

mydata1=read.table('~/mywork/Dnase1/zlater/max/data_max/parallel/one2')
mydata2=read.table('~/mywork/Dnase1/zlater/max/data_max/parallel/two2')
mydata3=read.table('~/mywork/Dnase1/zlater/max/data_max/parallel/three2')
library('vioplot')
vioplot(mydata1$V2,mydata2$V3,mydata3$V4,c(0.78),names=c(1,2,3,'all'))
mtext('SCC',side=2,line=2)
library('graphics')
clip(x1=3.85,x2=4.15,y1=0,y2=1)
abline(h=0.78)
points(x=4,y=0.78)
usr=par('usr')
do.call('clip',as.list(usr))
lines(c(1,2,3,4),c(max(mydata1$V2),max(mydata2$V3),max(mydata3$V4),c(0.78)),col=2)
lines(c(1,2,3,4),c(median(mydata1$V2),median(mydata2$V3),median(mydata3$V4),c(0.78)),col=3)
lines(c(1,2,3,4),c(min(mydata1$V2),min(mydata2$V3),min(mydata3$V4),c(0.78)),col=9)
mtext('(b)',side=3,at=-0.15,cex=1.5)
