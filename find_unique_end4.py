#import sys
import argparse

parser = argparse.ArgumentParser( prog='find_unique_end',description="find the either 5' or 3' unique seqs", epilog='python find_unique_end.py -i inputfile -o outputfile')
parser.add_argument ('-i','--input',help='Input File Name', default="./data/mouse_VH.txt")
#parser.add_argument('-5','--head', help="trim 5' region", action='store_true')
#parser.add_argument('-o','--output',help='Output File Name',required="true")
parser.add_argument('-3','--direction', help="trim 3' region", action='store_false', default=True)
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
	if is5prime == True:
		#print oligo
		return oligo[:length-1]
	else:
		return oligo[length+1:]

####parse the input file, generate key value pair of each sequences
def read_dnaFastas(fastaFile):
	Infile1 = open(fastaFile, 'r')
	germName=''
	germSeq=''
	allDict={}
	nucleotides=set('ATGCatgcNn')
	for line in Infile1:
		line = line.strip().strip('\n')
		if line and line.startswith(">"):
			line=line.lstrip('>')
			#process previous entry
			allDict[germName]= germSeq
			#begin a new entry
			germName=line
			germSeq=''
		elif line =='':  #/check whether it is the last line
				allDict[germName]= germSeq
		elif  set(line) <= nucleotides and len(line)>10: 
			germSeq = germSeq + line
	Infile1.close()
	return allDict;

########## prepare output files ###############################
#Infilename1 = args.input
#Infile1 = open(Infilename1, 'r')
#Outfilename1= args.output
	
Outfilename1=args.input.rstrip(".txt")+"_fulllength_noN"+str(args.primerlen)+".txt"
Outfile1 = open(Outfilename1, 'w')

if args.direction == True:
	Outfilename2=args.input.rstrip('.txt')+"Sense_primerlen_"+str(args.primerlen)
	Outfilename3=args.input.rstrip('.txt')+"Unique_5prime_lenth_"+str(args.primerlen)
else:
	Outfilename2=args.input.rstrip('.txt')+"AntiSense_primerlen_"+str(args.primerlen)
	Outfilename3=args.input('.txt')+"Unique_3prime_lenth_"+str(args.primerlen)
Outfile2 = open(Outfilename2, 'w')
Outfile3 = open(Outfilename3, 'w')

############# main ###############
allFasta={}
allFasta = read_dnaFastas(args.input)
trunFasta={}
unique_trunFasta={}
unique_Oligo=[]
all_Oligo=[]

for key in allFasta:
	Outfile1.write('>'+str(key) + "\t" + allFasta[key] + '\n')
	if not allFasta[key].upper().startswith('N'):
		
		currentOligo= truncate(allFasta[key],args.direction,args.primerlen) 
		# generate trimmed fasta file and dictinary
		all_Oligo.append(currentOligo)
		trunFasta[key] = currentOligo 
		Outfile2.write(str(key) + "\t" + trunFasta[key] + '\n')
print str(len(allFasta))+' sequences have been writen into ' + Outfilename1 + '\n'
print str(len(trunFasta)) + ' truncated sequences have been writen into ' + Outfilename2 + '\n'

unique_Oligo=set(all_Oligo)

for seq in unique_Oligo:
	count=0
	for key in trunFasta:
		if trunFasta[key] == seq:
			count = count +1
		unique_trunFasta[seq]= count;
	Outfile3.write(">"+str(count) +'_counts_'+ seq+"\n" + seq + '\n')


print str(len(unique_Oligo)) +" unique oligos whose lenth is "+str(args.primerlen)+ 'have been writen into ' + Outfilename3 + '\n'
Outfile1.close()
Outfile2.close()
Outfile3.close()