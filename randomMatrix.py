import random

rows = 1000
cols = 1000
fout = open('matrix.txt','w')

for row in range(rows):
    for col in range(cols):
        fout.write('%.2f ' % random.uniform(0,10))
    fout.write('\n')
