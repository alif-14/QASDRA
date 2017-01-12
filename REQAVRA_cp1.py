#Ali Fotouhi
#Research Assistant - Istanbul Technical University
#DAMGA Lab - UYBHM Building - ITU Ayazaga Campus - 34469 Istanbul-Turkey
#E-mail: fotouhi@itu.edu.tr


import io
import time
import math
import random
import argparse
import sys
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from reportlab.lib import colors 
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image,Table
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

#General Info
def GeneralInfo(fq):
    global minqs,maxqs,maxlen,minlen
    rl=len(fq) #read length
    if  rl > maxlen:
        maxlen=rl
    if rl<minlen:
        minlen=rl
    minfqs=min(fq)
    maxfqs=max(fq)
    if minfqs<minqs:
        minqs=minfqs
    if maxfqs>maxqs:
        maxqs=maxfqs

    return

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
    elif qsf6==64:
        qsf=64
            
    return qsf

def MeanQuality(fqs):
    global meanqualities, qualities, sum_of_rqavgs #read Quality Average
    rq=0
    for i in fqs:
        rq+=i # read quality
        qualities[i]=qualities.get(i,0)+1
    rqavg= rq/len(fqs) #read qualities average
    meanqualities[rqavg]=meanqualities.get(rqavg,0)+1 # dict for storing each mean quality
    sum_of_rqavgs+=rqavg

    return

#Input list of quality scores
#Output: maximal Range with kth smallest element
def InverseRangeSelect(X,k,v):
    A=[]
    z=0
    q=[]
    n=len(X)-1
    for i in range(k+1):
        while X[z] > v and z<=n:
            if z==n: 
                z+=1
                break
            z+=1
            
        q.append(z)
        z+=1
        if z>n: break
    if i<k:
        begin=0
        end=n
        A.append((begin,end))
        return A
    begin=0
    end=q[-1]-1
    A.append((begin,end))
    if q[-1]==n:
        begin=q[0]+1
        end=n
        A.append((begin,end))
    while z<=n:
        begin=q[0]+1
        while X[z]>v and z<=n:
            if z==n: 
                z+=1
                break            
            z+=1
                      
        q.pop(0)
        q.append(z)
        z+=1     
        end=q[-1]-1
        A.append((begin,end))
        if q[-1]==n:
            begin=q[0]+1
            end=n
            A.append((begin,end))  
    return A

#to find Length of Maximal Ranges
#input Maximal Ranages of a each read in Fastq file
#output length of Maximal Ranges of each read

def MaximalRangesLengths(MR): 
    MRL=[] #Lengths of Maximal ranges
    global dispersion,aimr 
    simrpr=0 #simrpr is  cubic means
    for i,j in MR:
        l=(j-i)+1
        MRL.append(l)
        dispersion[l]=dispersion.get(l,0)+1
        simrpr+=(l*l*l)/float(len(MR)) #sum for intensifying the ranges by powering three 
    simrpr=simrpr**(1./3.)    
    aimr+=simrpr 

    return MRL

#giving MRL as argument and getting metrics frequencies
def LMetrics(mrl):
    #for each read we compute longest,shortest,average and number of maximal ranges
    global longest, shortest, average, MRnumber,sum_of_max , sum_of_min , sum_of_avg , sum_of_mrno ,longfilter,shortfilter,avgfilter, maxfr,minfr,avgfr ,var,cv,variances,CV, sf, af, lf #shortest average longest filtration 
    ma=max(mrl)
    mi=min(mrl)
    avg=sum(mrl)/len(mrl)
    
    mrno=len(mrl)
    MRnumber[mrno]=MRnumber.get(mrno,0)+1
    if longfilter:
        if ma>lf:
            longest[ma]=longest.get(ma,0)+1
            sum_of_max += ma 
        else:
            maxfr+=1 #of filtered read
    else:
        longest[ma]=longest.get(ma,0)+1
        sum_of_max += ma #computing sum of maxes for calculating the overall average 
    if shortfilter:        
        if mi>sf:
            shortest[mi]=shortest.get(mi,0)+1
            sum_of_min += mi 
        else:
            minfr+=1 #of filtered read 
    else:
        shortest[mi]=shortest.get(mi,0)+1
        sum_of_min +=mi #computing sum of mins for calculating the overall average 
    if avgfilter:
        if avg>af:
            average[avg]=average.get(avg,0)+1
            sum_of_avg+= avg 
        else:
            avgfr+=1 #of filteration read
    else:
        average[avg]=average.get(avg,0)+1
        sum_of_avg+= avg #computing sum of avgs for calculating the overall average
    # Variance , Coefficient Variation
    try:
        var=0; cv=0
        var=sum([(xi - avg)**2 for xi in mrl]) / float(mrno)
        cv=math.sqrt(var)/avg
    except:
        var=0
    CV+=cv
                
    return

#for finding the average of the metrics and printing terminal based info
def metricsavg(fqr): 
    print "\nprocessed FAStQ file: %s, with %s Base Phred Quality Scores, is consisited of %s reads and %s of them is processed"  %(fastqfile,base,nor,fqr)
    print "Processed Reads of this file, have Minimum read length of %s and Maximum read length of %s." %(minlen,maxlen)
    print "Processed Reads of this file, have Minimum Phred quality score of %s and Maximum phred quality score of %s." %(minqs,maxqs)
    if longfilter == True: print 'for Maximum length, %s of reads have been filtered!!' %(maxfr) 
    if shortfilter == True: print 'for Minimum length, %s of reads have been filtered!!' %(minfr)  
    if avgfilter == True: print 'for Average length, %s of reads have been filtered!!' %(avgfr) 
    print 'Average length of Reads is: ' , fqss/fqr
    print 'Average Mean Quality of the data is: ', "{0:.2f}".format(sum_of_rqavgs/float(fqr))
    print 'Average of Longest Maximal Ranges is: ', "{0:.2f}".format(sum_of_max/float((fqr-maxfr)))
    print 'Average of shortest Maximal Ranges is: ', "{0:.2f}".format(sum_of_min/float(fqr-minfr))
    print 'Grand Average of Maximal Ranges is: ', "{0:.2f}".format(sum_of_avg/float(fqr-avgfr))
    print 'Grand Average of Maximal Ranges Cubic mean is : ', "{0:.2f}".format(aimr/fqr)
    print 'Average of Coefficient Variations is: ' , "{0:.2f}".format(CV/fqr)

    return

#for report
def Report(im0,t):

    doc = SimpleDocTemplate(fastqfile+'k'+str(k)+'v'+str(v)+"Report.pdf",pagesize=letter, rightMargin=0.5,leftMargin=0.5, topMargin=0.5,bottomMargin=0.5)
    Story=[]

    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))

    formatted_time = time.ctime()
    header="<font size=14 ><strong>QAVRA - Quality Assessment via Range Analysis </strong></font>"
    ptext = "<font size=11>File Name: </font><font size=14 color='red'><strong><u>%s</u></strong></font><font size=11> / Date: </font><font size=11 color='blue' > <u>%s</u></font>" % (fastqfile,formatted_time)
    p0 = Paragraph(header, styles["Normal"])
    p1 = Paragraph(ptext, styles["Normal"])
 
    reporttitle=[[p0],[p1]]
    rprttitle=Table(reporttitle,1*[5*inch], 2*[0.2*inch],style=[('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE')]) #
    Story.append(rprttitle)
    Story.append(Spacer(1, 10))

    tabletitle1="<font size=12 ><strong>Input Sequencing Data Digest:</strong></font>"
    tabletitle2="<font size=12 ><strong>Computed QAVRA Vector for k= %d , v= %d:</strong></font>" % (k,v)

    p3 = Paragraph(tabletitle1, styles["BodyText"])
    p4 = Paragraph(tabletitle2, styles["BodyText"]) 
    tabletitles=[[p3,p4]]
    tables=Table(tabletitles,2*[3.5*inch], 1*[0.4*inch])
    Story.append(tables)

    info=[['Quality Score Format:', str(base)+' ASCII-based Phred','Average Longest Maximal Range length',str("{0:.2f}".format(sum_of_max/float(fqr-maxfr)))],['Quality Scores (min,max):', '('+str(minqs)+','+str(maxqs)+')','Average Shortest Maximal Range Length',str("{0:.2f}".format(sum_of_min/float(fqr-minfr)))],['Number of Reads: ',str(nor),'Grand Average Maximal Range Length',str("{0:.2f}".format(sum_of_avg/float(fqr-avgfr)))],['Processed Number of Reads: ',str(fqr),'Cubic Average Maximal Range Length', str("{0:.2f}".format(aimr/fqr)) ],['Read Length (min,max):','('+str(minlen)+','+str(maxlen)+')','Average Coefficient of Variation', str("{0:.2f}".format(CV/fqr))]]

    infotable=Table(info,style=[('LINEBEFORE',(2,0),(-2,-1),1,colors.black),('LINEABOVE',(0,0),(4,0),1,colors.black),])
    Story.append(infotable)
    Story.append(Spacer(1, 7))

    Story.append(im0)    
    Story.append(t)
    ptext = "<font size=6 >developed by: </font><font size=7 color='green' > fotouhi@itu.edu.tr</font>"
    Story.append(Paragraph(ptext, styles["Normal"]))
    doc.build(Story)

    return

    
#Plot function        
def lmetricplot(p):
#fontweight='bold', 
#Mean Qualities
    plt.figure(1)
    if p:
        y=[float(i*100)/sum(meanqualities.values()) for i in meanqualities.values()]
        plt.ylabel('Frequency %')
    else:
        y=meanqualities.values()
        plt.ylabel('Frequency')
    x=meanqualities.keys()
    plt.yscale('log')
    plt.bar(x,y,width=0.7,align='center',color='black')
    plt.title('Distribution of Mean Qualities',y=1.05)
    plt.xlabel('Quality')
    buf5 = io.BytesIO()
    plt.savefig(buf5, format='png')
    im5 = Image(buf5,4*inch, 2*inch)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())

#qualities distribution
    plt.figure(2)
    if p:
        y=[float(i*100)/sum(qualities.values()) for i in qualities.values()]
        plt.ylabel('Frequency %')
    else:
        y=qualities.values()
        plt.ylabel('Frequency')
    x=qualities.keys()
    plt.yscale('log')
    plt.bar(x,y,width=0.7,align='center',color='black')
    plt.title('Distribution of Qualities',y=1.05)
    plt.xlabel('Quality')
    buf6 = io.BytesIO()
    plt.savefig(buf6, format='png')
    im6 = Image(buf6,4*inch, 2*inch)
    mng = plt.get_current_fig_manager()
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())    

#Longest
    plt.figure(3)
    if p:
        y=[float(i*100)/sum(longest.values()) for i in longest.values()]
        plt.ylabel('Frequency %')
    else:
        y=longest.values()
        plt.ylabel('Frequency')
    x=longest.keys()
    plt.yscale('log')
    plt.bar(x,y,width=0.7,align='center',color='black')
    plt.title('Distribution of Longest Lengths of Maximal Ranges'+' ; K= '+str(k)+ ', V= '+str(v),y=1.05)
    plt.xlabel('Length of Maximal Range')
    buf1 = io.BytesIO()
    plt.savefig(buf1, format='png')
    im1 = Image(buf1,4*inch, 2*inch)
    buf1.seek(0)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())

#Shortest,'Shortest Lengths of Maximal Ranges
    plt.figure(4)
    if p:
        y=[float(i*100)/sum(shortest.values()) for i in shortest.values()]
        plt.ylabel('Frequency %')
    else:
        y=shortest.values()
        plt.ylabel('Frequency')
    x=shortest.keys()
    plt.yscale('log')
    plt.bar(x,y,width=0.7,align='center',color='black')
    plt.title('Distribution of Shortest Lengths of Maximal Ranges'+' ; K= '+str(k)+ ', V= '+str(v),y=1.05)
    plt.xlabel('Length of Maximal Range')
    buf2 = io.BytesIO()
    plt.savefig(buf2, format='png')
    im2 = Image(buf2,4*inch, 2*inch)
    buf2.seek(0)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    
#average,'Average Lengths of Maximal Ranges
    plt.figure(5)
    if p:
        y=[float(i*100)/sum(average.values()) for i in average.values()]
        plt.ylabel('Frequency %')
    else:
        y=average.values()
        plt.ylabel('Frequency')
    x=average.keys()
    plt.yscale('log')
    plt.bar(x,y,width=0.7,align='center',color='black')
    plt.title('Distribution of Average Lengths of Maximal Ranges'+' ; K= '+str(k)+ ', V= '+str(v),y=1.05)
    plt.xlabel('Length of Maximal Range')
    buf3 = io.BytesIO()
    plt.savefig(buf3, format='png')
    im3 = Image(buf3,4*inch, 2*inch)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    
#MRnumber,'Number of Maximal Ranges per read
    plt.figure(6)
    if p:
        y=[float(i*100)/sum(MRnumber.values()) for i in MRnumber.values()]
        plt.ylabel('Frequency %')
    else:
        y=MRnumber.values()
        plt.ylabel('Frequency')
    x=MRnumber.keys()
    plt.yscale('log')
    plt.bar(x,y,width=0.7,align='center',color='black')
    plt.title('Distribution of Number of Maximal Ranges Per Read'+' ; K= '+str(k)+ ', V= '+str(v),y=1.05)
    plt.xlabel('Number of Maximal Range')
    buf4 = io.BytesIO()
    plt.savefig(buf4, format='png')
    im4 = Image(buf4,4*inch, 2*inch)
    mng = plt.get_current_fig_manager()
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())

#dispersion,'Dispersion of Maximal Ranges lengths
    plt.figure(7)
    if p:
        y=[float(i*100)/sum(dispersion.values()) for i in dispersion.values()]
        plt.ylabel('Frequency %')
    else:    
        y=dispersion.values()
        plt.ylabel('Frequency')
    x=dispersion.keys()
    plt.yscale('log')
    plt.bar(x,y,width=0.7,align='center',color='black')
    plt.title('Dispersion of Maximal Ranges lengths'+' ; K= '+str(k)+ ', V= '+str(v),y=1.05)
    plt.xlabel('Length of Maximal Range')
    buf0 = io.BytesIO()
    plt.savefig(buf0, format='png')
    im0 = Image(buf0,4*inch, 2*inch)
    buf0.seek(0)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    data=[[im1,im2],[im3,im4],[im5,im6]]
    t=Table(data)
    Report(im0,t)
    plt.show()

    return

###################################################### metrics and frequencies for each metric ###################################################################
longest={} ; shortest={} ; average={} ; MRnumber={} ; dispersion={} ; meanqualities={} ; qualities={}
###################################################### variables for computing averages ###################################################################3
sum_of_max=0 ; sum_of_min=0 ; sum_of_avg=0 ; sum_of_mrno=0 ; sum_of_rqavgs=0 ; maxfr=0; minfr=0; avgfr=0;variances=0;CV=0;minlen=999999999;maxlen=0;minqs=256;maxqs=0;lenequality = False
aimr=0 #average intensified maximal ranges associated with length of reads

######################################################giving the arguments using command line###################################################################
parser = argparse.ArgumentParser(prog='REQAVRA_cp1',description=' ',usage='This program will analyze entered fastq file according to the reads phred quality scores,in order to run program file name, k, and v values are necessary; there are also some other optional features if you wish to use\n%(prog)s [optional features are: -r random sampling, -fl for filtering lengths less than user defined longest length, -fs for filtering lengths less than user defined shortest length, -fa for filtering lengths less than user defined average length, -b ASCII-BASE of Phred quality scores, and -p for plotting distributions instead of their percentage]\nData should be entered like this: python REQAVRA_cp1.py filename.fastq k v [-r number -fl number -fs number -fa number -b {33,64} -p]',epilog='sample command line should be like this:'+'\n'+'python REQAVRA_cp1.py 2 20 -r 50 -fl 20 -fs 5 -fa 15 -b 33 -p'+'\n' + 'if *** Remember you can skip optional features if still it is not clear read the readme.txt file')
parser.add_argument('filename',metavar='File Name', help="Name of the fastq file should be like this: 'filename.fastq' ")
parser.add_argument('k-value',type=int,metavar='k Value',help='value for k')
parser.add_argument('v-value',type=int,metavar='v Value',help='value for v')
parser.add_argument('-b',metavar='ASCII-BASE amount', type=int,choices=[33,64],help='ASCII-BASE of phred quality scores',default=0)
parser.add_argument('-r',metavar='random sampling',help='a number in the range of 1...100 for random sampling',type=int,choices=range(1,101))
parser.add_argument('-fl',type=int,metavar='minimum longest length',help='maximum length for Maximal Ranges length filtration')
parser.add_argument('-fs',type=int,metavar='minimum shortest length',help='minimum length for Maximal Ranges length filtration')
parser.add_argument('-fa',type=int,metavar='minimum average length',help='average length for Maximal Ranges length filtration')
parser.add_argument('-p',help='ploting distribution of metrics instead of percentage',action="store_false",default=True)

inputs=vars(parser.parse_args())
fastqfile=inputs['filename']
k=inputs['k-value']
v=inputs['v-value']
if inputs['r'] != None:
    rv=inputs['r']
    randomization = True
else:
    randomization = False

if inputs['fl'] != None: #longest
    lf=inputs['fl']
    longfilter=True
else:
    longfilter = False
if inputs['fs'] != None: #shortest
    sf=inputs['fs']
    shortfilter=True
else:
    shortfilter = False
if inputs['fa'] != None: #average
    af=inputs['fa']
    avgfilter=True
else:
    avgfilter = False

if inputs['b'] == 0:
    check=True
else:
    check =False
    base=inputs['b']
p=inputs['p']

###########################################################################opening file#####################################
#Extracting FredQualityScores for each read FQS, FQR fredQualityRead
MaximalRanges=[]
fqr=1 #of processed reads
nor=1 #of reads
fqss=0 #fqs size
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
                    lenequality = True
                    sys.exit() 
                line=line.rstrip()
                if check == True: 
                    base=type_checking(line)
                    check = False
                if randomization == True and random.randint(0,100) > rv: 
                    nor+=1
                    l+=1                        
                    continue
                for item in line:
                    fqs.append(ord(item)-base)
                fqss+=len(fqs)
                GeneralInfo(fqs)
                MeanQuality(fqs)
                MaximalRanges=InverseRangeSelect(fqs,k,v)
                MRL=MaximalRangesLengths(MaximalRanges)
                LMetrics(MRL)
                fqr+=1
                nor+=1    
            l+=1
    fqr-=1
    nor-=1  #for counting the len of fastq

    metricsavg(fqr) #for printing averages
    lmetricplot(p) #for plotting metrics

except IOError:
    print 'cannot Open', fastqfile
except:
    if lenequality == True:
        print "*******\nan ERROR occured while processing your file!! read length and associated qualities doesn't macth at read %s and line %s of your file!!! check your file and try again later!\n*******" % (str(fqr),str(l))
    else:        
        print "Your File Cannot Be Processed!!! ***" , " *** Try Again!!! "

