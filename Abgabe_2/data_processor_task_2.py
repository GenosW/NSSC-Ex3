import numpy as np
import matplotlib.pyplot as plt
import csv
import os


## Code starts here ##

errorOwn = np.array([])
errorDiagonal = np.array([])
errorCholesky = np.array([])
scriptdir = os.path.dirname(__file__)
data_file1 = os.path.join(scriptdir,'data.txt')
data_file2 = os.path.join(scriptdir,'dataEigen.txt')
data_file3 = os.path.join(scriptdir,'dataEigenCholesky.txt')
with open(data_file1) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter =';')
    line_count = 0
    for row in csv_reader:
        if(line_count > 0):
            errorOwn = np.append(errorOwn,float(row[0]))
        line_count += 1
        
        
with open(data_file2) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter =';')
    line_count = 0
    for row in csv_reader:
        if(line_count > 0):
            errorDiagonal = np.append(errorDiagonal,float(row[0]))
        line_count += 1
        
        
with open(data_file3) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter =';')
    line_count = 0
    for row in csv_reader:
        if(line_count > 0):
            errorCholesky = np.append(errorCholesky,float(row[0]))
        line_count += 1
       
fig, ax = plt.subplots()
plt.semilogy(np.arange(len(errorOwn))*20, errorOwn, label='Our implementation')
plt.semilogy(np.arange(len(errorOwn))*20, errorDiagonal, label='Eigen: Diagonal')
plt.semilogy(np.arange(len(errorOwn))*20, errorCholesky, label='Eigen: Incomplete Cholseky')
plt.title('Performance of different implementations in respect to iteration number',fontsize=10)
plt.xlabel('Number of iterations')
plt.ylabel('|rk|2 / |r0|2')
ax.tick_params(labelsize=10)
ax.yaxis.label.set_size(10)
ax.xaxis.label.set_size(10)
plt.grid()
plt.legend()
ax.legend(prop={'size': 10})
output = os.path.join(scriptdir,'Errors.png')
plt.savefig(output)

