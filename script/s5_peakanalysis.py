# This script is to compute the proportion of tfbs peaks in DHS peaks.

# The binarysearch algorithm.
# list1 is the peak of the  Dnase1 site.
# search is the peak of the  tfbs.
def Binarysearch(list1,search):
	low=0
	high=len(list1)-1
	mid=int((low+high)/2)
	while((high-low)>1):
		if int(list1[mid])==int(search):
			low=mid
			break
		if int(list1[mid])>int(search):
			high=mid
			mid=(low+high)/2
		else:
			low=mid
			mid=(low+high)/2
	return low

import os
import re
# f1 is the peak file of DHS. It generate the dnase1dict
# {'chr1':[[start],[end]],'chr2':[[start],[end]],...}
def credhsdict(f1):
	dnase1dict={}
	for i in range(1,23):
		dnase1dict.setdefault('chr'+'%s' % (i),[[],[]])
	dnase1dict.setdefault('chrX',[[],[]])
	dnase1dict.setdefault('chrY',[[],[]])

	file1=open(f1,'r')
	for eachline1 in file1:
		eachline1=eachline1.strip()
		m=re.split('\s',eachline1)
		dnase1dict[m[0]][0].append(int(m[1]))
		dnase1dict[m[0]][1].append(int(m[2]))
	file1.close()
	return dnase1dict


# get the dnase1dict
dnase1dict=credhsdict('Epigenetics/Data/H1hesc/DHS/Peaks/H1_ESDukeDNaseSeq')


# compute the proportion.
dir1='Epigenetics/Data/H1hesc/Histone/Peaks'
resultfile=open('mywork/Dnase1/zlater/data/aa','w')
allnum=0  # the number of all the 22 tfbs peaks
allinnum=0  # the number of all the 22 tfbs peaks in DHS peaks
for eachfile in os.listdir(dir1):
	n=eachfile.split('.')
	tfbsname=n[0]
	dir2=dir1+r'/'+eachfile
	tfbsfile=open(dir2,'r')
	tfbsnum=0   # the number of the one tfbs peaks
	tfbsinnum=0    # the number of the one tfbs peak in DHS peak
	for eachline in tfbsfile:
		m=re.split('\s',eachline)
		chro=m[0]
		if chro in ['chrM','chrY','chrX']:
			continue
		allnum=allnum+1
		tfbsnum=tfbsnum+1
		start=int(m[1])
		end=int(m[2])
		startindex=Binarysearch(dnase1dict[chro][0],start)
		if dnase1dict[chro][1][startindex] < start and dnase1dict[chro][0][startindex+1] > end:    # the peaks may overlap
			continue
		else:
			allinnum=allinnum+1
			tfbsinnum=tfbsinnum+1
	proportion=tfbsinnum*1.0/tfbsnum
	resultfile.write('%s\t%.3f\n' % (tfbsname,proportion))

print allnum
allproportion=allinnum*1.0/allnum
resultfile.write('%s\t%.3f\n' % ('all',allproportion))
resultfile.close()




























