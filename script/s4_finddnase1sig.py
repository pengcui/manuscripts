#This script is to find the Dnase1 signal in every dnase1 sites.


# The binarysearch algorithm.
# list1 is the sequence of dnase1.
# search is the Dnase1 site.
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



# Read the uniq chromsome's dnase signal to ROM.
# f is the file of the uniq chromsome's dnase signal.
def readtfbs(f):
	numastart=[]
	numaend=[]
	level=[]
	file1=open(f,'r')
	for eachline1 in file1:
		m=re.split('\s',eachline1)
		numastart.append(int(m[1]))
		numaend.append(int(m[2]))
		level.append(float(m[3]))
	file1.close()
	return numastart,numaend,level


# Get the name and Dnase1 region.
# eachline is one line of the file in each chr of dnase1 sites.
def getdnase1(eachline):
	m=re.split('\s',eachline)
	if m[0] != 'Chr':
		chrname=m[0]
		start=int(m[1])
		end=int(m[2])
		return chrname,start,end



# find the Dnase1 signal of dnase1 sites
import re
import os
nu=range(1,23)
list1=[]
for eachnu in nu:
	chrname='chr'+'%s' % eachnu
	list1.append(chrname)

list1.append('chrX')


dir1='/home/ckivip/Epigenetics/Data/H1hesc/DHS/peaks'  
dir2='/home/ckivip/Epigenetics/Data/H1hesc/DHS/Signal'
dir3='/home/ckivip/mywork/Dnase1/zlater/data/zzdnase1sig' 
file3=open(dir3,'w')
file3.write('%s\t%s\t%s\t%s\n' % ('Chr','start','end','dnase1sig'))
for eachchr in list1:
	print eachchr,
	numastart,numaend,level=readtfbs(dir2+r'/'+eachchr)
	file1=open(dir1+r'/'+eachchr,'r')
	for eachline1 in file1:
		eachline1=eachline1.strip()
		chrname,start,end=getdnase1(eachline1)
		file3.write('%s\t%s\t%s\t' % (chrname,start,end))
		startindex=Binarysearch(numastart,start)+1
		endindex=Binarysearch(numastart,end)
		if endindex > startindex:
			findlevel=sum(level[startindex:endindex+1])*1.0/(endindex-startindex+1)
			file3.write('%.3f\n' % (findlevel))
		elif endindex == startindex:
			file3.write('%.3f\n' % (level[startindex]))
		else:
			file3.write('NA\n')
	
	del numastart,numaend,level
	file1.close()
file3.close()

