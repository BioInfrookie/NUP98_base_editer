# !/usr/bin/env python
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--infile", dest="infile", help="give a fasta file to me", metavar="FILE")
parser.add_option("--outfile",  dest="outfile", help="the name of oufput file [fastq]", metavar="FILE")
parser.add_option("--outfile_2",  dest="outfile_2", help="the name of oufput file [fastq]", metavar="FILE")
parser.add_option("--scale",  dest="scale", help="the name of oufput file [fastq]", metavar="FILE")
parser.add_option("--fname",  dest="fname", help="the name of oufput file [Sample]", metavar="FILE")

(options, args) = parser.parse_args()

pileup = open(options.infile, "r")
output = open(options.outfile,'w')
output_2 = open(options.outfile_2,'w')
scale_factor = float(options.scale)
Sample_fname = options.fname

output.write("chr"+"\t"+"position"+"\t"+"ref"+"\t"+"rate"+"\t"+"num"+"\t"+"mutate"+"\t"+"type"+"\t"+"Sample"+"\n")
output_2.write("chr"+"\t"+"position"+"\t"+"ref"+"\t"+"A"+"\t"+"T"+"\t"+"C"+"\t"+"G"+"\n")

for line in pileup:
	linesplit = line.split("\t")
	chrom = linesplit[0]
	position = linesplit[1]
	ref = linesplit[2]
	A = "A:" + str(linesplit[4].count("A") + linesplit[4].count("a"))
	T = "T:" + str(linesplit[4].count("T") + linesplit[4].count("t"))
	C = "C:" + str(linesplit[4].count("C") + linesplit[4].count("c"))
	G = "G:" + str(linesplit[4].count("G") + linesplit[4].count("g"))

	if ref == "A":
		A = "A:" + str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "T":
		T = "T:" + str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "C":
		C = "C:" + str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "G":
		G = "G:" + str(linesplit[4].count(".") + linesplit[4].count(","))
	else:
		continue

	A_num = float(int(A.split(":")[1])*scale_factor)
	T_num = float(int(T.split(":")[1])*scale_factor)
	C_num = float(int(C.split(":")[1])*scale_factor)
	G_num = float(int(G.split(":")[1])*scale_factor)

	all_reads = float(A_num) + float(T_num) + float(C_num) + float(G_num)
	if all_reads != 0:
		A_rate = format(float(A_num)/all_reads, '.10')
		G_rate = format(float(G_num)/all_reads, '.10')
		C_rate = format(float(C_num)/all_reads, '.10')
		T_rate = format(float(T_num)/all_reads, '.10')
	else:
		continue
	# output = "\t".join([chrom, position, A, T, C, G])
	output.write(chrom+"\t"+position+"\t"+ref+"\t"+str(A_rate)+"\t"+str(A_num)+"\t"+str(ref)+"->"+"A"+"\t"+"A"+ "\t" + Sample_fname+"\n"+
	chrom+"\t"+position+"\t"+ref+"\t"+str(T_rate)+"\t"+str(T_num)+"\t"+str(ref)+"->"+"T"+"\t"+"T"+ "\t"+Sample_fname+"\n"+
	chrom+"\t"+position+"\t"+ref+"\t"+str(C_rate)+"\t"+str(C_num)+"\t"+str(ref)+"->"+"C"+"\t"+"C"+ "\t"+Sample_fname+"\n"+
	chrom+"\t"+position+"\t"+ref+"\t"+str(G_rate)+"\t"+str(G_num)+"\t"+str(ref)+"->"+"G"+"\t"+"G"+ "\t"+Sample_fname+"\n")
	output_2.write(chrom+"\t"+position+"\t"+ref+"\t"+str(A_num)+"\t"+str(T_num)+"\t"+str(C_num)+"\t"+str(G_num)+"\n")

pileup.close()
output.close()