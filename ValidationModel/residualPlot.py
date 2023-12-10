import numpy as np
import csv
import json
import matplotlib.pyplot as plt

'''
Run 
pyFoamPlotWatcher.py --write-files --single-data-files-only  log.mpirun
in consol before this script.
'''
path = 'log.foamRun.analyzed/linear_'

def csvMaker(param):
    
    with open('%s%s'%(path,param), 'r') as input_file:
        lines = input_file.readlines()
        newLines = []
        for line in lines:
            newLine = line.strip('|').split()
            newLines.append(newLine)

    with open('%s%s.csv'%(path,param), 'w') as output_file:
        file_writer = csv.writer(output_file)
        file_writer.writerows(newLines)

def resPlotter(paramList):


    for param in paramList:
        csvMaker(param)
        with open('%s%s.csv'%(path,param),'r',newline="") as fp:
            data = csv.reader(fp)

            time = []
            residuals = []
            counter = 0
            for dat in data:
                if (counter >= 1) and not(counter % 2):
                    time.append(float(dat[0]))
                    residuals.append(float(dat[1]))
                counter += 1
        plt.plot(time,residuals, label = param)
    plt.yscale('log')
    plt.ylabel('Residual')
    plt.xlabel('Time (s)')
    plt.title('Change in initial residuals during simulation')
    plt.legend()
    plt.grid()
    plt.show()
    
    fp.close()


parameters = ['Ux', 'Uy', 'p_rgh']
resPlotter(parameters)