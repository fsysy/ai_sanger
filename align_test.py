from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
file_name = 'MATN3_DES136-7_2-1RD072018-02-14-16-04-15_copy.ab1'
record = SeqIO.read(file_name, 'abi')
sequence = record.annotations['abif_raw']['PBAS1']
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
trimmed_position = abi_trim(record)
trimmed_sequence = sequence[trimmed_position[0]:trimmed_position[1]]

alignments = pairwise2.align.globalxx(trimmed_sequence,sequence)
print(format_alignment(*alignments[0]))

