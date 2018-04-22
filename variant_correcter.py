from Bio import SeqIO
import matplotlib.pyplot as plt
import os
file_name = 'MATN3_DES136-7_2-1RD072018-02-14-16-04-15_copy.ab1'
record = SeqIO.read(file_name, 'abi')
blue = record.annotations['abif_raw']['DATA9']
red = record.annotations['abif_raw']['DATA10']
green = record.annotations['abif_raw']['DATA11']
yellow = record.annotations['abif_raw']['DATA12']
call_time = record.annotations['abif_raw']['PLOC2']
sequence = record.annotations['abif_raw']['PBAS1']
color = ['blue','red','green','yellow']

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
        ax.text(call_time[i],-15,sequence[i])
    
class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        print('click', event)
        print(event.xdata, event.ydata)
        if event.dblclick:
            print "double"
        else:
            print'single'

fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot([0], [0])
linebuilder = LineBuilder(line)
sequence_writer(ax)
plt.plot(blue, color='black')
plt.plot(red, color='green')
plt.plot(green, color='red')
plt.plot(yellow, color='blue')

plt.grid(True)

plt.show()
