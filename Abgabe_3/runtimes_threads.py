import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv
name_of_txt_file = 'runtimes_Exa3.csv'

threats1 = []
time_sin1 = []
time_sin2inv1 = []
xcube1 = []

threats2 = []
time_sin2 = []
time_sin2inv2 = []
xcube2 = []


with open(name_of_txt_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=';')
    line_count = 0

    for row in csv_reader:
        if (line_count > 0 and line_count < 9):
            threats1.append(int(row[0]))
            time_sin1.append(float(row[1]))
            time_sin2inv1.append(float(row[2]))
            xcube1.append(float(row[3]))
        if(line_count >= 9):
            threats2.append(int(row[0]))
            time_sin2.append(float(row[1]))
            time_sin2inv2.append(float(row[2]))
            xcube2.append(float(row[3]))
        line_count += 1


x = np.arange(len(threats1))  # the label locations
plt.figure()
matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)


plt.title("1 million samples per thread", fontsize=15)
plt.plot(threats1, time_sin1, '-o', label='sin')
plt.plot(threats1, time_sin2inv1, '-o', label='sininv')
plt.plot(threats1, xcube1, '-o', label='xcube')
plt.legend(fontsize=15)
plt.xlabel("number of threads", fontsize=15)
plt.ylabel("runtime[s]", fontsize=15)
plt.grid()
plt.savefig("runtimes_1m_abs.png")

plt.figure()
plt.title("100 million samples per thread", fontsize=15)
plt.xlabel("threads", fontsize=15)
plt.ylabel("runtime", fontsize=15)
matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)
plt.plot(threats2, time_sin2, '-o', label='sin')
plt.plot(threats2, time_sin2inv2, '-o', label='sininv')
plt.plot(threats2, xcube2, '-o', label='xcube')
plt.legend(fontsize=15)
plt.xlabel("number of threads")
plt.ylabel("runtime[s]")
plt.grid()
plt.savefig("runtimes_100m_abs.png")