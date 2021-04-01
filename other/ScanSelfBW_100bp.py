import optparse
import subprocess
import twobitreader
import sys

def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac

def sperr(cmd):
    a=subprocess.Popen(cmd,stderr=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac

def isFloat(x):
    try:
        float(x)
        if str(x) in ['inf', 'infinity', 'INF', 'INFINITY', 'True', 'NAN', 'nan', 'False', '-inf', '-INF', '-INFINITY', '-infinity', 'NaN', 'Nan']:
            return False
        else:
            return True
    except:
        return False

#cmd1 = "bigWigSummary a.bw chr1 100 200 100 "
#res = sp(cmd1)
#signal = str(res[0], encoding = "utf-8")
#signal.strip("\n").split("\t")
#bigWigSummary a.bw chr1 3073253 3074322 100

def get_Command(Ch,S,E,bw,interval,extend):
    center=round((S+E)/2)
    cmd1 = ["bigWigSummary",bw,Ch,center-extend,center+extend,interval]
    cmd1 = ' '.join([ str(i) for i in cmd1 ])
    return cmd1

def get_Signal(inputfile,outputfile,bw,interval,extend):
    infile = open(inputfile)
    outfile = open(outputfile,"w")
    line = infile.readline()
    count = 0
    while line:
        array = line.strip("\n").split()
        Ch = array[0]
        S = int(array[1])
        E = int(array[2])
        cmd1 = get_Command(Ch,S,E,bw,interval,extend)

        res = sp(cmd1)
        signal = str(res[0], encoding = "utf-8")
        signal = signal.strip("\n").split("\t")
        if signal[0]=="":
            signal = [0]*interval

        for i in range(0,len(signal)):
            if isFloat(signal[i]):
                signal[i] = float(signal[i])
            else:
                signal[i] = 0

        array.extend(signal)
        outfile.write("\t".join(map(str,array))+"\n")
        
        count = count+1
        if count % 500 == 0:
            print(count)

        line = infile.readline()

    infile.close()
    outfile.close()

#------------------------main---------------------------
def main():
    parser = optparse.OptionParser()
#========OPTIONS============
    parser.add_option("-n","--name",dest="name",type="str",help="")
    parser.add_option("-b","--bed",dest="bed",type="str",help="")
    parser.add_option("-w","--bw",dest="bw",type="str",help="")
    parser.add_option("-o","--outdir",dest="outdir",type="str",help="")

#========GET OPTIONS========
    (options,args) = parser.parse_args()
    name=options.name
    bed=options.bed
    bw=options.bw
    outdir=options.outdir

    interval = 40
    extend = 2000
    outfile=outdir+name+".peak2kbNearby_100bp.bed"
    get_Signal(bed,outfile,bw,interval,extend)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)
