from Bio import SeqIO
import matplotlib.pyplot as plt
import csv
record = SeqIO.read('MATN3_DES136-7_2-1RD072018-02-14-16-04-15_copy.ab1', 'abi')
blue = record.annotations['abif_raw']['DATA9']
red = record.annotations['abif_raw']['DATA10']
green = record.annotations['abif_raw']['DATA11']
yellow = record.annotations['abif_raw']['DATA12']
call_time = record.annotations['abif_raw']['PLOC2']
color = ['blue','red','green','yellow']
#print call_time
data = []
with open('data.csv','rb') as f:
	reader = csv.reader(f)
	color_code = 0
	for line in reader:
		data.append(line)
		print line
print(data[0])
plt.plot(tuple(data[0]), color='blue')
plt.plot(tuple(data[1]), color='red')
plt.plot(tuple(data[2]), color='green')
plt.plot(tuple(data[3]), color='yellow')	
	
	
#plt.plot(data[0], color='blue')
#plt.plot(data[1], color='red')
#plt.plot(data[2], color='green')
#plt.plot(data[3], color='yellow')
#for x in call_time:
#	plt.text(x,-10,'C')
#print data[0]
#print data[1]
#print data[2]
#print data[3]
#print data[4]	


plt.show()