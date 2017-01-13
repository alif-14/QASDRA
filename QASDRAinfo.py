#Ali Fotouhi
#Research Assistant - Istanbul Technical University
#DAMGA Lab - UYBHM Building - ITU Ayazaga Campus - 34469 Istanbul-Turkey
#E-mail: fotouhi@itu.edu.tr

import argparse
import sys

#Phred Quality Base Check
def type_checking(l):
    s="!\"#$%&'()*+,-./0123456789:;<=>?"
    sh="JKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
    qsf3=0
    qsf6=0
    qsf=0
    for i in s:
        if i in l:
            qsf3=33 #quality score format
            break
    for j in sh:
        if j in l:
            qsf6=64
            break
            
    if qsf3==33 and qsf6==64:
        qsf=33
    elif qsf3 != 0:
        qsf = qsf3
    else:
        qsf=64
            
    return qsf


parser = argparse.ArgumentParser(description='Fastq Analyser.')
parser.add_argument('filename',metavar='File Name', help="Name of the fastq file should be like this: 'filename.fastq' ")
inputs=vars(parser.parse_args())
fastqfile=inputs['filename']
check=True
minlen=999999999 ; maxlen=0 ; avg=0 ; fqr=1 ; minqs=256 ; maxqs=0 ; avgrl=0;lenequality=False

try:
    with open(fastqfile,'rU') as fstq:
        l=1 #of line
        for line in fstq:
            if line == " " : continue
            if l % 2 == 0 and l % 4 != 0:
                rlen=len(line)              
            if l % 4 == 0:
                fqs=[]
                qlen=len(line)
                if qlen != rlen:
                    lenequality=True
                    sys.exit() 
                line=line.rstrip()
                rl=len(line) #read length

                avgrl+=rl
                if  rl > maxlen:
                    maxlen=rl
                if rl<minlen:
                    minlen=rl
                if check == True: 
                    base=type_checking(line)
                    check = False
                for item in line:
                    fqs.append(ord(item)-base)
                minfqs=min(fqs)
                maxfqs=max(fqs)
                if minfqs<minqs:
                    minqs=minfqs
                if maxfqs>maxqs:
                    maxqs=maxfqs
                fqr+=1
  
            l+=1
    fqr-=1
    averagel=avgrl/fqr
    print "\nHere we will provide you general information about your file, to give you insight about your data, then you can choose your 'v' value according to these information."
    print "\nYour fastq file: %s, consists of %s reads, with phred quality %s-ASCII base." %(fastqfile,fqr,base)
    print "\nIn this file, Minimum read length is %s, Maximum read length is %s, and Average length of reads is  %s." %(minlen,maxlen,averagel)
    print "\nIn this file, Minimum Phred quality score is %s, and Maximum phred quality score is %s." %(minqs,maxqs)

except IOError:
    print 'cannot Open', fastqfile
except:
    if lenequality == True:
        print "*******\nan ERROR occured while processing your file!! read length and associated qualities doesn't macth at read %s and line %s of your file!!! check your file and try again later!\n*******"  % (str(fqr),str(l))
    else:        
        print "Your File Cannot Be Processed!!! ***" , " *** Try Again!!! "


try:
    print "\nHere we will provide you general information about your file, to give you insight about your data, then you can choose your 'v' value according to these information."
    print "\nYour fastq file: %s, consists of %s reads, with phred quality %s-ASCII base." %(fastqfile,fqr,base)
    print "\nIn this file, Minimum read length is %s, Maximum read length is %s, and Average length of reads is  %s." %(minlen,maxlen,averagel)
    print "\nIn this file, Minimum Phred quality score is %s, and Maximum phred quality score is %s." %(minqs,maxqs)
except:
    print 'Entered file has a problem!!'
