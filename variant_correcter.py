from Bio import SeqIO, Entrez, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq

import matplotlib.pyplot as plt
import os, sys

file_name = 'MATN3_DES136-7_2-1RD072018-02-14-16-04-15_copy.ab1'
Entrez.email = "fsysy@naver.com"
reference_gene_NM_number = "NM_002381.4"
record = SeqIO.read(file_name, 'abi')
blue = record.annotations['abif_raw']['DATA9']
red = record.annotations['abif_raw']['DATA10']
green = record.annotations['abif_raw']['DATA11']
yellow = record.annotations['abif_raw']['DATA12']
call_time = record.annotations['abif_raw']['PLOC2']
sequence = record.annotations['abif_raw']['PBAS1']
color = ['blue','red','green','yellow']
#if trigger is True, keydown makes a nucleotide sequence changing.
nucleotide_change_trigger = False



def id_search(NM_number):
    handle = Entrez.esearch(db='nucleotide', term=NM_number)
    record = Entrez.read(handle)
    idlist = record["IdList"]
    if len(idlist) == 1:
        print "The id is %s"%idlist[0]
        return idlist[0]
    else:
        if len(idlist) == 0:
            print "No search NM_number in genbank. please search other names"
            return False
        else:
            print "Multple ids were found! %s"%s(",".join(idlist))
            return False

def get_sequence(Id):
    handle = Entrez.efetch(db="nucleotide", id=Id,rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    return record

def abi_trim(seq_record):
	start = False   # flag for starting position of trimmed sequence
	segment = 20    # minimum sequence length
	trim_start = 0  # init start index
	cutoff = 0.05   # default cutoff value for calculating base score
	if len(seq_record) <= segment:
		return seq_record
	else:
		score_list = [cutoff - (10 ** (qual / -10.0)) for qual in
			seq_record.letter_annotations['phred_quality']]
		cummul_score = [0]
		for i in range(1, len(score_list)):
			score = cummul_score[-1] + score_list[i]
			if score < 0:
				cummul_score.append(0)
			else:
				cummul_score.append(score)
				if not start:
					trim_start = i
					start = True
	trim_finish = cummul_score.index(max(cummul_score))
	return [trim_start,trim_finish]

def sequence_writer(ax):
    trim = abi_trim(record)
    for i in range(trim[0],trim[1]):
        ax.text(call_time[i],-15,sequence[i],ha='center', va='center')

def count_gap(seq1, refseq):
        alignments = pairwise2.align.localms(seq1,refseq,1,-1,-2,-2)
        formatted_align1 = alignments[0][0]
        return formatted_align1

def nucleotide_click(xpoint,ypoint):
        if ypoint<-5 and ypoint >-20:
                print xpoint, ypoint
                with open("result.ab1","w") as output_handle:
                        SeqIO.write(record, output_handle, "abi")
                        print "check the save data"
        
def press(event):
    print('press', event.key)
    sys.stdout.flush()

        
class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        #print('click', event)
        #print(event.xdata, event.ydata)
        nucleotide_click(event.xdata,event.ydata)




#find sequence
gene_id = id_search("NM_002381.4")
if gene_id:
    reference_gb = get_sequence(gene_id)


#sequence alignment
trimmed_position = abi_trim(record)
trimmed_sequence = sequence[trimmed_position[0]:trimmed_position[1]]
trimmed_sequence_seq = Seq(trimmed_sequence)

print count_gap(trimmed_sequence_seq, reference_gb)


#figure making!
fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot([0], [0])
fig.canvas.mpl_connect('key_press_event', press)
linebuilder = LineBuilder(line)
sequence_writer(ax)
plt.plot(blue, color='black')
plt.plot(red, color='green')
plt.plot(green, color='red')
plt.plot(yellow, color='blue')

plt.grid(True)

plt.show()
