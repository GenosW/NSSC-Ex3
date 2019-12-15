import numpy
import csv
import os
import matplotlib.pyplot as plt

if __name__ == "__main__":
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir,"slurm-797.csv")
    policy = "static,1"
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
    should = []
    for i,thr in zip(iterations[:1],threads):
        temp = []
        temp2 = []
        for ele in i:
            temp.append(float(i[-1])/float(ele))
        for ind in range(len(i)):
            if ind == 0:
                temp2.append(temp[0])
            else:
                print(ind," : ",(thr[ind]/thr[ind-1]))
                val = temp[ind-1]/(thr[ind]/thr[ind-1])
                print(val)
                temp2.append(val)
        itPerc.append(temp)
        should.append(temp2)

    labels = ["256","512","1024","2048"]
    for x,y,z,lab in zip(threads,itPerc,should,labels):
        plt.plot(x,y,marker="o",label=lab)
        plt.plot(x,z,marker="x",label=lab)
    plt.suptitle("Iterations done (scaled)\nschedule="+policy)
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(dir,policy+"_percent"))

    for x,y,lab in zip(threads,iterations,labels):
        plt.plot(x,y,marker="o",label=lab)
    plt.suptitle("Iterations done\nschedule="+policy)
    plt.savefig(os.path.join(dir,policy+"_iterations"))

    for x,y,lab in zip(threads,runtime,labels):
        plt.plot(x,y,marker="o",label=lab)
    plt.suptitle("Average runtime per iteration\nschedule="+policy)
    plt.savefig(os.path.join(dir,policy+"_runtime"))
                    
                