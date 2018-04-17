from Bio import SeqIO
import matplotlib.pyplot as plt
import csv
record = SeqIO.read('MATN3_DES136-7_2-1RD072018-02-14-16-04-15_copy.ab1', 'abi')
blue = record.annotations['abif_raw']['DATA9']
red = record.annotations['abif_raw']['DATA10']
green = record.annotations['abif_raw']['DATA11']
yellow = record.annotations['abif_raw']['DATA12']
call_time = record.annotations['abif_raw']['PLOC2']
sequence = record.annotations['abif_raw']['PBAS1']
color = ['blue','red','green','yellow']
sampling_half_length = 100

#trim_abi_file:return trim start point, end point
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

def peak_calling_choice(trim_point):
        trim_start = trim_point[0] #20
        trim_end = trim_point[1] #577
        
        trim_start_time = call_time[trim_start] #233
        trim_end_time = call_time[trim_end] #6964

        q = trim_start #20

        while q<trim_end: #20<577
                q_time = call_time[q] #233
                if q_time - trim_start_time > sampling_half_length:
                        break
                q=q+1
        peak_call_start_q = q

        q = trim_end
        while q>trim_start:
                q_time = call_time[q]
                if trim_end_time - q_time > sampling_half_length:
                        break
                q=q-1
        peak_call_end_q = q
        return [peak_call_start_q, peak_call_end_q]
        

trim_point = abi_trim(record)

peak_cc = peak_calling_choice(trim_point)
q=0
plot_q = peak_cc[0]
for i in range(peak_cc[0]:peak_cc[1]):
        if q>5:
                break
        print i
        peak_cc_time = [call_time[peak_cc[i]],call_time[peak_cc[i]]]
        plot_range = [peak_cc_time[0]-sampling_half_length,peak_cc_time[0]+sampling_half_length]
        plt.plot(list(blue)[plot_range[0]:plot_range[1]], color='blue')
        plt.plot(list(red)[plot_range[0]:plot_range[1]], color='red')
        plt.plot(list(green)[plot_range[0]:plot_range[1]], color='green')
        plt.plot(list(yellow)[plot_range[0]:plot_range[1]], color='yellow')	
        plt.show()
        q = q +1




#data = []
#with open('data.csv','rb') as f:
#	reader = csv.reader(f)
#	color_code = 0
#	for line in reader:
#		data.append(line)
		#print line
#print(data[0])
