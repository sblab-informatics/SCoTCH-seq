import sys
import os
from get_gapped_cigar import *
import re



#Read Input File(eg: Splice<strand>_tmp.txt)
filename = sys.argv[1]

#Output(strip "tmp" from Input file)
writefile=filename.replace('_tmp.txt','.txt')

#temporary/intermediate files (useful for troubleshooting). Can comment it out of not needed
Intermediate='Intermediate' + writefile
tmpfile='tmp_' + writefile


# Open files to read and write
read_file = open(filename, 'r')
outfile = open(writefile,'w')
tmp = open(tmpfile, 'w')
Interim = open(Intermediate, 'w')


#Start reading file
Lines = read_file.readlines()

#Set count flags to 0
digitcount=0  
count = 0
dot='.'

# Strips the newline character
for line in Lines:
	values = line.split("\t")
	if values[3] == '*':
		continue
	if 'XA:' in values[4]: 
		continue
	#if re.search('z',values[4], re.IGNORECASE):
	#	continue
	test_str=(values[4].removeprefix('XB:Z:'))
	test_str_len = len(test_str)
	i=0
	pos = values[2]
	cigar=values[3]
	pos1 = int(pos)-1
	K = '.'
	bismark_str=''
	prev=0
	flag=0
	wflag=0
	L = "\n" + values[0] +  "\t" + '+' +  "\t" + values[1] + "\t" + values[3] 
	Interim.writelines(L)	
	for ele in test_str:
		i=i+1
		if ele.isdigit() and flag==1:
			pos1=pos1-prev
			digit=(str(prev)) + (str(ele))
			pos1=pos1+int(digit)
			digitcount=int(digit)
			flag=1
			wflag=1		
		if ele.isdigit() and flag==0:
			digitcount=0
			pos1=pos1+int(ele)
			prev=int(ele)
			digitcount=int(ele)
			flag=1
			wflag=1
		if i==test_str_len and wflag==1:
			pos1=pos1+1
			dots=digitcount*'.'
			digitcount=0
			L=dots
			Interim.writelines(L)
			bismark_str=bismark_str + L
			flag=0
			wflag=0

		elif ele == 'x':
			flag=0
			dots=digitcount*'.'
			digitcount=0
			pos1=pos1+1
			L = dots+ele  
			Interim.writelines(L)
			wflag=0
			bismark_str=bismark_str + L
			continue
		elif ele == 'X':
			flag=0
			dots=digitcount*'.'
			digitcount=0
			pos1=pos1+1
			L = dots+ele
			Interim.writelines(L)
			wflag=0
			bismark_str=bismark_str + L
			continue
		elif ele.isalpha() and ele != 'x' and ele != 'X':  #else:
			pos1=pos1+1
			dots=digitcount*'.'
			digitcount=0
			L=dots+dot
			Interim.writelines(L)
			wflag=0
			bismark_str=bismark_str + L	
			flag=0
			continue	
	new_cgar,new_seq=cigarToList(cigar,bismark_str)
	NL = "\n" + values[0] +  "\t" + values[1] + "\t" + values[2]  + "\t" + values[3]  + "\t" 
	tmp.writelines(NL)
	NL=new_seq
	tmp.writelines(NL)			
	
	meth_pos=([pos for pos, char in enumerate(new_seq) if char == 'x'])
	for i in meth_pos:
		mpos=(int(values[2]) + i)
		L=values[0] +  "\t" + '-' +  "\t" + values[1] + "\t" + str(mpos)  + "\t" + 'z' + "\n"
		outfile.writelines(L)

	meth_posX=([pos for pos, char in enumerate(new_seq) if char == 'X'])
	for i in meth_posX:
		mpos=(int(values[2]) + i)
		L=values[0] +  "\t" + '+' +  "\t" + values[1] + "\t" + str(mpos)  + "\t" + 'Z' + "\n"
		outfile.writelines(L)


tmp.close()
Interim.close()
read_file.close()
outfile.close()

#Remove all temporary files (comment out this section for debugging)
os.remove(Intermediate)
os.remove(tmpfile)

