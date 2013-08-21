# This script is to compute the histone signal in all the dnase1 peaks
import os
import re
def Binarysearch(list1,search):
	low=0
	high=len(list1)-1
	mid=int((low+high)/2)
	while((high-low)>1):
		if list1[mid]==search:
			low=mid
			break
		if list1[mid]>search:
			high=mid
			mid=int((low+high)/2)
		else:
			low=mid
			mid=int((low+high)/2)
	return low


# read modify levels to moddict
#{'mod':[[start],[end],[level]],...}
def readmodify(moddir,modchr):
	moddict={}
	for eachitem in os.listdir(moddir):
		modname=eachitem[0:-9]
		moddict.setdefault(modname,[[],[],[]])
		print modname
		mydir=moddir+r'/'+'%s' % (eachitem)+r'/'+'%s' % (modchr)
		modfile=open(mydir,'r')
		for eachline in modfile:
			m=re.split('\t',eachline)
			start=int(m[0])
			end=int(m[1])
			level=float(m[2])
			moddict[modname][0].append(start)
			moddict[modname][1].append(end)
			moddict[modname][2].append(level)
		modfile.close()
	return moddict

# read the position of dnase1 of particular chr into dnase1dict
def readdnase1(dnase1dir,dnase1chr):
	dnase1dict={}
	newdir=dnase1dir+r'/'+'%s' % (dnase1chr)
	dnase1file=open(newdir,'r')
	for eachline in dnase1file:
		m=re.split('\t',eachline)
		start=int(m[1])
		end=int(m[2])
		dnase1dict.setdefault(start,end)
	dnase1file.close()
	return dnase1dict

chrlist=[]
for i in range(1,23):
	chrlist.append('chr'+ '%s' % (i))
chrlist.append('chrX')

modnamelist=[]
# get the mod names
for eachitem in os.listdir('/home/ckivip/Epigenetics/Data/H1hesc/Histone/Signal4'):
	modname=eachitem[0:-9]
	modnamelist.append(modname)



# do it!
resultfile=open('/home/ckivip/mywork/Dnase1/zlater/data/histoneresult4','w')
resultfile.write('%s\t%s\t%s\t' % ('Chr','start','end'))

for each in range(0,len(modnamelist)):
	resultfile.write('%s\t' % (modnamelist[each]))
resultfile.write('\n')


for j in range(0,23):
	uniquechr=chrlist[j]
	dnase1dict=readdnase1('/home/ckivip/Epigenetics/Data/H1hesc/DHS/peaks',uniquechr)
	moddict=readmodify('/home/ckivip/Epigenetics/Data/H1hesc/Histone/Signal4',uniquechr)
	for eachkey1 in dnase1dict:
		start=eachkey1
		end=dnase1dict[eachkey1]
		resultfile.write('%s\t%s\t%s\t' % (uniquechr,start,end))
		for k in range(0,len(modnamelist)):
			startindex=Binarysearch(moddict[modnamelist[k]][0],start)+1
			endindex=Binarysearch(moddict[modnamelist[k]][0],end)
			if endindex > startindex:
				findlevel=sum(moddict[modnamelist[k]][2][startindex:endindex+1])*1.0/(endindex-startindex+1)
				resultfile.write('%.3f\t' % (findlevel))
			elif endindex == startindex:
				resultfile.write('%.3f\t' % (moddict[modnamelist[k]][2][startindex]))
			else:
				resultfile.write('NA\t')
		resultfile.write('\n')
	del dnase1dict
	del moddict



