setwd('~/mywork/Dnase1/zlater/sampleout/')
mydata250=read.table('mycors250')
mydata500=read.table('mycors500')
mydata750=read.table('mycors750')
mydata1000=read.table('mycors1000')
mydata1250=read.table('mycors1250')
mydata1500=read.table('mycors1500')
mydata1750=read.table('mycors1750')
mydata2000=read.table('mycors2000')
mydata2250=read.table('mycors2250')
mydata2500=read.table('mycors2500')
mydata2750=read.table('mycors2750')
mydata3000=read.table('mycors3000')
mydata3250=read.table('mycors3250')
mydata3500=read.table('mycors3500')
mydata3750=read.table('mycors3750')
mydata4000=read.table('mycors4000')
mydata4250=read.table('mycors4250')
mydata4500=read.table('mycors4500')
mydata4750=read.table('mycors4750')
mydata5000=read.table('mycors5000')
mydata5250=read.table('mycors5250')
mydata5500=read.table('mycors5500')
mydata5750=read.table('mycors5750')
mydata6000=read.table('mycors6000')
boxplot(mydata250[,1],mydata500[,1],mydata750[,1],mydata1000[,1],mydata1250[,1],mydata1500[,1],mydata1750[,1],mydata2000[,1],mydata2250[,1],mydata2500[,1],mydata2750[,1],mydata3000[,1],mydata3250[,1],mydata3500[,1],mydata3750[,1],mydata4000[,1],mydata4250[,1],mydata4500[,1],mydata4750[,1],mydata5000[,1],mydata5250[,1],mydata5500[,1],mydata5750[,1],mydata6000[,1],outpch=NA,xlab='Sample size (×1000)',names=seq(0.25,6,0.25))
mtext('SCC',side=2,line=2.5)
mean(mydata5000[,1])
sd(mydata5000[,1])