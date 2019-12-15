import numpy
import csv
import os
import matplotlib.pyplot as plt

if __name__ == "__main__":
    dir = os.path.dirname(__file__)
    resolutions = ["256","512","1024","2048"]
    for reso in resolutions:
        dataName = "data_"+reso+".csv"
        filename = os.path.join(dir,dataName)
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

        labels = ["dynamic","static","static,1","dynamic,innerLen"]
        for x,y,lab in zip(threads,runtime,labels):
            if lab=="static" or lab=="dynamic,innerLen": lineStyle = "--"
            else: lineStyle="-"
            plt.semilogy(x,y,marker="x",label=lab,ls=lineStyle,lw=2,basey=2)
        plt.title("Average runtime per iterationn\n resolution="+reso)
        plt.xlabel("# of threads")
        plt.ylabel("runtime [s]")
        plt.grid()
        plt.legend(title="schedule(,chunk size)")
        plt.savefig(os.path.join(dir,"resolution"+reso+".pdf"))
        plt.close()
                        