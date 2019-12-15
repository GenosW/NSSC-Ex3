import numpy as np
import matplotlib.pyplot as plt
import csv
import os 


def applyFilter(x,a):
    y = np.array([])
    for i in range(len(x)):
        if(i%a == 0):
            y = np.append(y,x[i])
        
    return y


## Code starts here ##

filtervalue = 50
rnorm = np.array([])
enorm = np.array([])
scriptdir = os.path.dirname(__file__)
data_file = os.path.join(scriptdir,'data.txt')
with open(data_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter =';')
    line_count = 0
    for row in csv_reader:
        if(line_count > 0):
            rnorm = np.append(rnorm,float(row[0]))
            enorm = np.append(enorm,float(row[1]))
        line_count += 1
       
rnorm = applyFilter(rnorm,filtervalue)
enorm = applyFilter(enorm,filtervalue)

fig, ax = plt.subplots()
plt.semilogy(np.arange(len(rnorm))*filtervalue,rnorm,label='|rk|2 / |r0|2')
plt.title('|rk|2 / |r0|2 in respect to iteration number',fontsize=15)
plt.xlabel('Number of iterations')
plt.ylabel('|rk|2 / |r0|2')
ax.tick_params(labelsize=10)
ax.yaxis.label.set_size(10)
ax.xaxis.label.set_size(10)
plt.grid()
plt.legend()
ax.legend(prop={'size': 10})
output = os.path.join(scriptdir,'R-norm.png')
plt.savefig(output)
fig, ax = plt.subplots()
plt.semilogy(np.arange(len(enorm))*filtervalue,enorm,label='|ek|A')
plt.title('Error on the A-norm in respect to iteration number',fontsize=15)
plt.xlabel('Number of iterations')
plt.ylabel('|ek|A')
ax.tick_params(labelsize=10)
ax.yaxis.label.set_size(10)
ax.xaxis.label.set_size(10)
plt.grid()
plt.legend()
ax.legend(prop={'size': 10})
output = os.path.join(scriptdir,'A-norm_of_ek.png')
plt.savefig(output)
            
