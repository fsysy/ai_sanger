from Bio import SeqIO
import matplotlib.pyplot as plt

record = SeqIO.read('MATN3_DES136-7_2-1RD072018-02-14-16-04-15_copy.ab1', 'abi')
blue = record.annotations['abif_raw']['DATA9']
red = record.annotations['abif_raw']['DATA10']
green = record.annotations['abif_raw']['DATA11']
yellow = record.annotations['abif_raw']['DATA12']
call_time = record.annotations['abif_raw']['PLOC2']
#print call_time

plt.plot(blue, color='blue')
plt.plot(red, color='red')
plt.plot(green, color='green')
plt.plot(yellow, color='yellow')
for x in call_time:
	plt.text(x,-10,'C')

plt.show()