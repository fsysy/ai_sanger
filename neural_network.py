from sklearn.neural_network import MLPClassifier
clf = MLPClassifier(solver = 'lbfgs', alpha=1e-5, hidden_layer_sizes=(100,2), random_state=1)

file_name = "training_set/text/MATN3_DES136-7_2-1RD072018-02-14-16-04-15_copy.ab1_training_set.txt"
f = open(file_name,"r")
line_data = f.readline()
i = 1
X=[]
y=[]

def comma_parser(line):
    parsed_data =  line.split(',')
    parsed_data[-1].replace("\n","")
    return parsed_data

while line_data:
    i += 1
    if (i<400):
        parsed_line = comma_parser(line_data)
        answer = parsed_line[0]
        question = [float(text) for text in parsed_line[3:]]
        X.append(question)
        y.append(answer)
        line_data = f.readline()
    else:
        break
clf.fit(X,y)
#print([question])

total = 0
correct_number = 0
error_number = 0
while line_data:
    parsed_line = comma_parser(line_data)
    answer = parsed_line[0]
    question = [float(text) for text in parsed_line[3:]]
    prediction = clf.predict([question])
    total += 1
    if answer == prediction[0]:
        correct_number += 1
    else:
        error_number +=1
    #print "correct: %s predict: %s"% (answer, clf.predict([question]))
    line_data = f.readline()
print "total: %d, correct: %d, error: %d, error_rate: %f%%"%(total, correct_number, error_number, float(error_number)/float(total)*100)
f.close()
