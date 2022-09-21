#python
# -*- coding: utf-8 -*-
#execute.server=zeta(MRC3)

import os
import sys
import argparse


parser =argparse.ArgumentParser()
parser.add_argument ('-i','--dir_in', action='store', dest="dir_in", help="input directory", type=str)
parser.add_argument ('-o','--dir_out', action='store', dest="dir_out", help="output directory", type=str)
parser.add_argument ('-T','--table', action='store', dest="SRAtable", help="Tab separated SRA table", type=str)
parser.add_argument ('-C','--count',action='store', dest="cout", help="Tab separated count table", type=str)

args = parser.parse_args()
dir_in=args.dir_in.rstrip("/")
dir_out=args.dir_out.rstrip("/")
input=args.SRAtable
cout=args.cout

cnt={}
hh={}
chck={}

#dir_in='/data4/gangtl95/BRCA/breast.GSE194434/fastq/'
#dir_out='/data4/gangtl95/BRCA/breast.GSE194434/fastq/'
#input='/data4/gangtl95/BRCA/breast.GSE194434/SraRunTable.txt'
#cout='/data4/gangtl95/BRCA/breast.GSE194434/fastq/cnt.txt'

if not os.path.exists(input): sys.exit('error in filename: '+input)

n=0

with open(input,"r") as inf:
    for line in inf:
        if n==0:pass
        else:
            line=line.rstrip()
            tmp=line.split(",")
            gsm_id= [s for s in tmp if "GSM" in s]
            srr_id= [s for s in tmp if "SRR" in s]

            id=srr_id + gsm_id

            sampleid=id[0]
            gsm=id[1]

            if input in cnt: cnt[input] +=1
            else: cnt[input] = 1

            if gsm not in hh:
                hh.setdefault(gsm, [])
            hh[gsm].append(sampleid)

        n+=1
with open(cout, "r") as ps:
    for line in ps:
        line=line.rstrip()
        line=line.split("\t")
        nm=line[0]
        num=line[1]
        chck[nm]=num



for key in hh.keys():
    for i in range(0,len (hh[key])):
        if int(chck[hh[key][i]])==2:
            print("paired")
            output1=dir_out+"/"+key+"_1.fastq"
            output2=dir_out+"/"+key+"_2.fastq"
            aa=0
            if os.path.isfile(output1): aa += 1
            if os.path.isfile(output2): aa += 1
            if aa==2:
                print("!! already exists: "+key)
                continue
                
            if len(hh[key]) > 1:
                cmd1="cat"
                cmd2="cat"
                for ff in hh[key]:
                    input1=dir_in+"/"+ff+"_1.fastq"
                    input2=dir_in+"/"+ff+"_2.fastq"
                    if not os.path.isfile(input1): sys.exit('file does not exist '+input1)
                    if not os.path.isfile(input2): sys.exit('file does not exist '+input2)         
                    cmd1 += " "+input1
                    cmd2 += " "+input2
                    cmd1 += " > "+ output1
                    cmd2 += " > "+ output2

            elif len(hh[key])==1:
                ff=hh[key][0]
                input1=dir_in+"/"+ff+"_1.fastq" 
                input2=dir_in+"/"+ff+"_2.fastq"
                if not os.path.isfile(input1): sys.exit('file does not exist '+input1)
                if not os.path.isfile(input2): sys.exit('file does not exist '+input2)
                cmd1 = "mv " +input1+" "+output1
                cmd2 = "mv " +input2+" "+output2
            
            else:
                print("error "+key+"\n")
                continue

            print (cmd1+"\n")
            os.system(cmd1)
            print (cmd2+"\n")
            os.system(cmd2)
            
        else:
            print("single")
            output1=dir_out+"/"+key+".fastq"
            aa=0
            if os.path.isfile(output1): aa += 1

            if aa==2:
                print("!! already exists: "+key)
                continue

            if len(hh[key]) > 1:
                cmd1="cat"
                for ff in hh[key]:
                    input1=dir_in+"/"+ff+".1.fastq"
                    if not os.path.isfile(input1): sys.exit('file does not exist '+input1)
                    cmd1 += " "+input1
                    cmd1 += " > "+ output1

            elif len(hh[key]) ==1:
                ff=hh[key][0]
                input1=dir_in+"/"+ff+".1.fastq"
                if not os.path.isfile(input1): sys.exit('file does not exist '+input1)
                    
                cmd1 = "mv " +input1+" "+output1

            else:
                print("error "+key+"\n")
                continue

        print (cmd1+"\n")
        os.system(cmd1)

sys.stdout=open("../SrrGsm.txt","w")

n=0

with open(input,"r") as inf:
    for line in inf:
        if n==0:pass
        else:
            line=line.rstrip()
            tmp=line.split(",")
            gsm_id= [s for s in tmp if "GSM" in s]
            srr_id= [s for s in tmp if "SRR" in s]

            id=srr_id + gsm_id

            sampleid=id[0]
            gsm=id[1]

            print(sampleid, gsm, sep="\t")
        n+=1
    sys.stdout.close()
