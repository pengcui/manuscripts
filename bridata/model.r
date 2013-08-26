
# We highly recommand to run these codes using Rstudio piece by piece.

####################################################################

# 1. The support vector regression model and plot the correlation.

setwd('/home/ckivip/mywork/Dnase1/zlater/manuscripts/bridata')

# (1) 10-fold cross validation in svr.

svrcrossvalidation=function(datamatrix,featurevector){
  library('e1071')
  kk=c()
  actual=c()
  pre=c()
  n=nrow(datamatrix)/10
  ##  using different test data in every iteration
  for (i in 1:10){
    if (length(kk) > 0){
      se=sample((1:nrow(datamatrix))[-kk],n)
    }else{
      se=sample((1:nrow(datamatrix)),n)
    }
    te=datamatrix[se,]
    tr=datamatrix[-se,]
    mysvr=svm(tr[,featurevector],log(tr[,1],2),probability=T)
    mysvrpre=round(predict(mysvr,te[,featurevector]),3)
    #  combine the results into actual and pre vectors in every iteration.
    actual=c(actual,log(te[,1],2))
    pre=c(pre,mysvrpre)
    kk=c(kk,se)
  }
  myresult=list(actual,pre)
  return(myresult)
}

# (2) commputation process.

data=read.table('./max_sig.txt',header=T)  ## modify this line to choose avg_sig.txt/max_sig.txt
mydata=as.matrix(data[,4:37])
histone=c(2:11)
tf=c(12:34)
allfeature=c(histone,tf)
mylist=svrcrossvalidation(mydata, allfeature)  ## modify this line to choose different features

# (3) plot the figure.

dev.off()
par(mfrow=c(1,2),mar=c(5,4,3,1))
plot(mylist[[1]], mylist[[2]],pch=20, col="blue", cex=0.6,xlab='',ylab='')
title(xlab="The Dnase-seq signals (log2)",ylab="Predicted signals (log2)",mgp=c(2.3,1,0))
mylm = lm(mylist[[1]]~mylist[[2]])
abline(mylm, lwd=2)
summary(mylm)  ## compute the R squared
text(0.05,-4.5,bquote(R^2 ~ '=' ~ 0.66),cex=0.8)
mtext('(a)',line=1,side=3,at=-5.5,adj=0,cex=1.5)


# barplot for the only one feature
featurepowers=read.table('./one_feature',header=F)
histonevec=1:10
tfvec=11:34

#--------------------------------------
# HM barplot
barplot(featurepowers[histonevec,2],ylim=c(0,0.6),xlab='',ylab='',axisnames=F)
mtext(featurepowers[histonevec,1],line=0.5,side=1,at=seq(0.6,12,1.2),las=2,cex=0.8)

# TF barplot
barplot(featurepowers[tfvec,2],ylim=c(0,0.6),xlab='',ylab='',axisnames=F)
mtext(featurepowers[tfvec,1],line=0.5,side=1,at=seq(0.7,28,1.2),las=2,cex=0.8)
# -------------------------------------

mtext('SCC',side=2,line=2.3)
mtext('(b)',line=1,side=3,at=-1,cex=1.5)

################################################################################

# 2. the linear regression model

data=read.table('./avg_sig.txt',header=T)  ## modify this line to choose avg_sig.txt/max_sig.txt

datamatrix=data[,4:37]

# 10-fold cross validation in lm
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
  
  histone=c(2:11)
  tf=c(12:34)
  allfeature=c(histone,tf)
  
  #---------------------------------------
  
  # you can choose different liner model here!
  
  # the HM model
   mylm=lm(log(dnase1sig,2)~H3k9me3+H3k27ac+H3k4me3+H3k4me1+H3k4me2+H3k27me3+H3k36me3+H3k9ac+H4k20me1+H3k79me2, data=tr)
   mylmpre=predict(mylm,data.frame(te[,histone]))
  # the TF model
  mylm=lm(log(dnase1sig,2)~CTCF+BRCA1+BACH1+CEBPB+CHD1+CHD2+CJUN+GTF2F1+CMYC+CTBP2+JUND+MXI1+MAX+MAFK+RFX5+RAD21+NRF1+SUZ12+SIN3A+TBP+ZNF143+USF2+ZNF274,data=tr)
  mylmpre=predict(mylm,data.frame(te[,tf]))
  # all features model
  mylm=lm(log(dnase1sig,2)~H3k9me3+H3k27ac+H3k4me3+H3k4me1+H3k4me2+H3k27me3+H3k36me3+H3k9ac+H4k20me1+H3k79me2+CTCF+BRCA1+BACH1+CEBPB+CHD1+CHD2+CJUN+GTF2F1+CMYC+CTBP2+JUND+MXI1+MAX+MAFK+RFX5+RAD21+NRF1+SUZ12+SIN3A+TBP+ZNF143+USF2+ZNF274,data=tr)
  mylmpre=predict(mylm,data.frame(te[,allfeature]))

  #---------------------------------------
 
  actual=c(actual,log(te[,1],2))
  pre=c(pre,mylmpre)
  kk=c(kk,se)
}

mylm = lm(actual~pre)
summary(mylm)
################################################################################

# 3. compute correlation coefficients of all possibal one-feature, two-feature, 
# three-feature combinations and plot the results of one, two. three feature models.

# This step includes thousands of models and may spend much computation time. 
# So we use 'snawfall' package to parallel this task.

setwd('/home/ckivip/mywork/Dnase1/zlater/manuscripts/bridata')

# (1) A modified version of 10-fold cross validation in svr 
svrcrossvalidation=function(featurevector){   
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
  myresult=round(cor(actual,pre,method=c('spearman')),3)
  return(myresult)
}

#--------------------------------------------------------------------------

# (2) parallel the task
library('snowfall')
data=read.table('./max_sig.txt',header=T)
allnames=names(data)[4:37]
datamatrix=as.matrix(data[,4:37])

sfInit(parallel=T,cpus=10) # set 10 cpus
sfExport('datamatrix')
sfExport('svrcrossvalidation')

combs=combn(33,3)  ## modify this line to choose one, two, three feature combinations
combs=combs+1  ## because the first column is dnasesig 
mylist=list()
for (i in 1:ncol(combs)){
  mylist[[i]]=combs[,i]
}

mycors=sfClusterApplyLB(mylist,fun='svrcrossvalidation')
sfStop()

# write the result
for (j in 1:length(mycors)){
  lastresult=t(c(allnames[combs[,j]],mycors[[j]]))
  write.table(lastresult,file='./three_feature',append=T,quote=F,sep='\t',row.names=F,col.names=F)
}

#---------------------------------------------------------

# (3) violin plot of feature combinatin results

mydata1=read.table('./one_feature')
mydata2=read.table('./two_feature')
mydata3=read.table('./three_feature')
library('vioplot')
vioplot(mydata1$V2,mydata2$V3,mydata3$V4,c(0.85),names=c(1,2,3,'all'))
mtext('SCC',side=2,line=2.2)
library('graphics')
clip(x1=3.85,x2=4.15,y1=0,y2=1)
abline(h=0.81)
points(x=4,y=0.81)
usr=par('usr')
do.call('clip',as.list(usr))
lines(c(1,2,3,4),c(max(mydata1$V2),max(mydata2$V3),max(mydata3$V4),c(0.81)),col=2)
lines(c(1,2,3,4),c(median(mydata1$V2),median(mydata2$V3),median(mydata3$V4),c(0.81)),col=3)
lines(c(1,2,3,4),c(min(mydata1$V2),min(mydata2$V3),min(mydata3$V4),c(0.81)),col=9)
mtext('(b)',line=1,side=3,at=-0.15,cex=1.5)