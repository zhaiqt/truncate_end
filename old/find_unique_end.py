#import sys
import argparse

parser = argparse.ArgumentParser( prog='find_unique_end',description="find the either 5' or 3' unique seqs", epilog='python find_unique_end -i inputfile -o outputfile')
parser.add_argument ('-i','--input',help='Input File Name', default="mouse_VH.txt")
#parser.add_argument('-5','--head', help="trim 5' region", action='store_true')
#parser.add_argument('-o','--output',help='Output File Name',required="true")
parser.add_argument('-3','--direction', help="trim 3' region", action='store_false', default='true')
parser.add_argument('-n','--primerlen', type=int, help='length of nucleotide', default='22')
parser.print_help()
args=parser.parse_args()

##show values##
#print (parser.print_help())
print ("input file: %s" % args.input)
#print ("output file: %s" % args.output)
print ("start with 5prime?")
print args.direction
print ("primer length: %s" % args.primerlen)

def truncate(oligo,is5prime,length):
	if is5prime == true:
		return oligo[:length-1]
	else:
		return oligo[length+1:]


Infilename1 = args.input
Infile1 = open(Infilename1, 'r')

#Outfilename1= args.output
if args.direction==True:
	Outfilename2=rgs.input.rstrip(".txt")+"_sense_primerlen_"+str(args.primerlen)+".txt"
else:
	Outfilename2=args.input.rstrip(".txt")+"_AntiSense_primerlen_"+str(args.primerlen)+".txt"

Outfilename1=args.input.rstrip(".txt")+"_fulllength_noN"+str(args.primerlen)+".txt"
Outfile1 = open(Outfilename1, 'w')
Outfile2 = open(Outfilename2, 'w')



####parse the input file, generate key value pair of each sequences
germName=''
germSeq=''
allDict={}
for line in Infile1:
	line = line.strip()
	if line and line.startswith(">"):
		line=line.lstrip('>')
		#process previous entry
		if not germSeq.upper().startswith("N"):
			allDict[germName]= germSeq
		#begin a new entry
		germName=line
		germSeq=''
	else: 
		germSeq = germSeq + line

for key in allDict:
	Outfile1.write(key + "\t" + allDict[key] + '\n')
print str(len(allDict))+' sequences have been writen into ' + Outfilename1 + '\n'


Infile1.close()
Outfile1.close()
Outfile2.close()