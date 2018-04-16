from Bio import SeqIO
import os.path
import csv
import matplotlib.pyplot as plt
from collections import defaultdict
channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12','PLOC2']
GATC = {"G":(1,0,0,0),"A":(0,1,0,0),"T":(0,0,1,0),"C":(0,0,0,1),"R":(1,1,0,0),"Y":(0,0,1,1),"M":(0,1,0,1),"K":(1,0,1,0),"S":(1,0,0,1),"W":(0,1,1,0),"H":(0,1,1,1),"B":(1,0,1,1),"V":(1,1,0,1),"D":(1,1,1,0),"N":(1,1,1,1)}
f = open('data.csv', 'w')
wr = csv.writer(f)
 
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
 
def make_AI_training_set_form(seq_record,view_port=50):
	trimmed_position = abi_trim(seq_record)
	raw = record.annotations['abif_raw']
	peak_time = raw["PLOC2"]
	start_time = peak_time[trimmed_position[0]]
	end_time = peak_time[trimmed_position[1]]
	q = trimmed_position[0]
	while q < trimmed_position[1]:
		current_time_position = peak_time[q]
		if current_time_position > start_time + view_port and current_time_position < end_time - view_port:
			
			plt.plot(raw["DATA9"][current_time_position-view_port: current_time_position+view_port], color='blue')
			plt.plot(raw["DATA10"][current_time_position-view_port: current_time_position+view_port], color='red')
			plt.plot(raw["DATA11"][current_time_position-view_port: current_time_position+view_port], color='green')
			plt.plot(raw["DATA12"][current_time_position-view_port: current_time_position+view_port], color='yellow')
			
			
			d9 = raw["DATA9"][current_time_position-view_port:current_time_position+view_port]
			d10 = raw["DATA10"][current_time_position-view_port:current_time_position+view_port]
			d11 = raw["DATA11"][current_time_position-view_port:current_time_position+view_port]
			d12 = raw["DATA12"][current_time_position-view_port:current_time_position+view_port]
			#ai_array = d9+d10+d11+d12+GATC[raw["PBAS2"][q]]
			wr.writerow(d9)
			wr.writerow(d10)
			wr.writerow(d11)
			wr.writerow(d12)
			wr.writerow(GATC[raw["PBAS2"][q]])
			wr.writerow('@')
			plt.show()
			break
		q = q + 1
	return [start_time, peak_time[trimmed_position[0]+1]]
 
	
folder = os.getcwd()
for filename in os.listdir(folder):
	if filename[-3:] == "ab1":
		record = SeqIO.read(filename, 'abi')
		print(make_AI_training_set_form(record,100))
f.close()
exit(0)