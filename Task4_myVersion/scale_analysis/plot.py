import numpy
import csv
import os
import matplotlib.pyplot as plt

if __name__ == "__main__":
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir,"slurm-794.csv")
    i = -1
    runtime = []
    iterations = []
    threads = []
    tempth = []
    temp = []
    temp2 = []
    with open(filename) as file:
        reader = csv.reader(file,delimiter=";")
        for row in reader:
            if row[0][0] == "#": #comment char
                pass
            else:
                tempth.append(int(row[0]))
                temp.append(int(row[1]))
                temp2.append(float(row[2]))
                if int(row[0]) == 40:
                    threads.append(tempth)
                    iterations.append(temp)
                    runtime.append(temp2)
                    temp = []
                    temp2 = []
                    tempth = []

    itPerc = []
    for i in iterations:
        temp = []
        for ele in i:
            temp.append(float(ele)/float(i[-1]))
        itPerc.append(temp)

    labels = ["256","512","1024","2048"]
    for x,y,lab in zip(threads,itPerc,labels):
        plt.plot(x,y,marker="o",label=lab)
    plt.suptitle("Iterations done (scaled)\nschedule=dynamic")
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(dir,"dynamic_percent"))

    for x,y,lab in zip(threads,iterations,labels):
        plt.plot(x,y,marker="o",label=lab)
    plt.suptitle("Iterations done\nschedule=dynamic")
    plt.savefig(os.path.join(dir,"dynamic_iterations"))

    for x,y,lab in zip(threads,runtime,labels):
        plt.plot(x,y,marker="o",label=lab)
    plt.suptitle("Average runtime per iteration\nschedule=dynamic")
    plt.savefig(os.path.join(dir,"dynamic_runtime"))
                    
                