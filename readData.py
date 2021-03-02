import numpy as np
import pandas as pd
import glob

relativePathList = {}
relativePathList["ageDist_time/before"] = "/logBin/"
relativePathList["ageDist_time/during"] = "/logBin/"
relativePathList["ageDist_time/after"] = "/logBin/"
relativePathList["ageDist_op/before"] = "/logBin/"
relativePathList["ageDist_op/during"] = "/logBin/"
relativePathList["ageDist_op/after"] = "/logBin/"
relativePathList["clusterSizeDist"] = "/logBin/"
relativePathList["clusterSizeDist_exact"] = "/logBin/"
relativePathList["clusterSizeDist_time"] = "/logBin/"
relativePathList["deltaUpperBoundDist_time/before"] = "/logBin/"
relativePathList["deltaUpperBoundDist_time/during"] = "/logBin/"
relativePathList["deltaUpperBoundDist_time/after"] = "/logBin/"
relativePathList["deltaUpperBoundDist_op/before"] = "/logBin/"
relativePathList["deltaUpperBoundDist_op/during"] = "/logBin/"
relativePathList["deltaUpperBoundDist_op/after"] = "/logBin/"
relativePathList["deltaUpperBoundDist_tot"] = "/logBin/"
relativePathList["interEventTime"] = "/average/"
relativePathList["interEventTimeDist_time/before"] = "/logBin/"
relativePathList["interEventTimeDist_time/during"] = "/logBin/"
relativePathList["interEventTimeDist_time/after"] = "/logBin/"
relativePathList["interEventTimeDist_op/before"] = "/logBin/"
relativePathList["interEventTimeDist_op/during"] = "/logBin/"
relativePathList["interEventTimeDist_op/after"] = "/logBin/"
relativePathList["interEventTimeDist_tot"] = "/logBin/"
relativePathList["interEventTime_deltaUpperBound/before"] = "/logBin/"
relativePathList["interEventTime_deltaUpperBound/during"] = "/logBin/"
relativePathList["interEventTime_deltaUpperBound/after"] = "/logBin/"
relativePathList["meanClusterSize"] = "/average/"
relativePathList["meanClusterSize_trial"] = "/average/"
relativePathList["orderParameter"] = "/average/"
relativePathList["orderParameter_trial"] = "/average/"
relativePathList["orderParameterDist"] = "/linBin/"
relativePathList["orderParameterVariance"] = "/average/"
relativePathList["orderParameterVariance_trial"] = "/average/"
relativePathList["sampled_deltaUpperBound_interEventTime"] = "/logBin/"
relativePathList["deltaUpperBound_interEventTime_tot"] = "/logBin/"
relativePathList["interEventTime_deltaUpperBound_tot"] = "/logBin/"
relativePathList["sampled_upperBound_interEventTime"] = "/logBin/"
relativePathList["interEventTime_upperBound_tot"] = "/"
relativePathList["upperBound_interEventTime_tot"] = "/"
relativePathList["sampled_time_interEventTime"] = "/logBin/"
relativePathList["time_interEventTime_tot"] = "/"
relativePathList["interEventTime_time_tot"] = "/"

relativePathList["points"] = "/"
relativePathList["dynamics"] = "/"


baseDirectory = "../data/mBFW/"
absolutePathList = {}
for observable, relativePath in relativePathList.items():
    absolutePathList[observable] = baseDirectory + observable + relativePath

#* CSV Reader
def readCSV(t_fileName):
    data = pd.read_csv(t_fileName, sep=',', header=None)
    data = data.values.transpose()
    if (len(data) == 1):
        return data[0]
    else:
        return tuple([row for row in data])


#* File Name Convections
def NG(t_networkSize, t_acceptanceThreshold):
    return "N{:.1e},G{:.1f}*".format(t_networkSize, t_acceptanceThreshold)

def NGT(t_networkSize, t_acceptanceThreshold, t_time):
    return "N{:.1e},G{:.1f}*,T{:.4f}*".format(t_networkSize, t_acceptanceThreshold, t_time)

def NGOP(t_networkSize, t_acceptanceThreshold, t_orderParameter):
    return "N{:.1e},G{:.1f}*,OP{:.4f}*".format(t_networkSize, t_acceptanceThreshold, t_orderParameter)


#* Get the order parameter/time value in directory
def extractRepeater(t_type, t_networkSize, t_acceptanceThreshold):
    repeaterList = set()
    if (t_type == "clusterSizeDist_time" or t_type == "orderParameterDist"):
        target = "T"
    elif (t_type == "clusterSizeDist" or t_type == "clusterSizeDist_exact"):
        target = "OP"
    for file in glob.glob(absolutePathList[t_type] + NG(t_networkSize, t_acceptanceThreshold)):
        repeaterList.add(float(file[file.find(target)+len(target) : file.find(".txt")]))
    return repeaterList

#* Read various points
def readPoint(t_type, t_networkSize, t_acceptanceThreshold):
    with open(absolutePathList["points"] + NG(t_networkSize, t_acceptanceThreshold)[:-1] + ".txt") as file:
        point = 0
        for line in file.readlines():
            if t_type in line:
                point = float(line[line.find(": ")+2 : ])
    return point

def readPoints(t_networkSize, t_acceptanceThreshold):
    points = {}
    with open(absolutePathList["points"] + NG(t_networkSize, t_acceptanceThreshold)[:-1] + ".txt") as file:
        content = file.readlines()
        for line in content:
            if "t_a" in line:
                points["t_a"] = float(line[line.find(": ")+2 : ])
            elif "m_a" in line:
                points["m_a"] = float(line[line.find(": ")+2 : ])
            elif "t_inflection" in line:
                points["t_inflection"] = float(line[line.find(": ")+2 : ])
            elif "m_inflection" in line:
                points["m_inflection"] = float(line[line.find(": ")+2 : ])
            elif "t_c_var" in line:
                points["t_c_var"] = float(line[line.find(": ")+2 : ])
            elif "m_c_var" in line:
                points["m_c_var"] = float(line[line.find(": ")+2 : ])
            elif "t_c_mcs" in line:
                points["t_c_mcs"] = float(line[line.find(": ")+2 : ])
            elif "m_c_mcs" in line:
                points["m_c_mcs"] = float(line[line.find(": ")+2 : ])
            elif "t_c_csd" in line:
                points["t_c_csd"] = float(line[line.find(": ")+2 : ])
            elif "m_c_csd" in line:
                points["m_c_csd"] = float(line[line.find(": ")+2 : ])
    return points

#* Read Observables
def read(t_type, t_networkSize, t_acceptanceThreshold, t_reapeater=None):
    #* Read time-accumulated distributions
    if (t_type == "clusterSizeDist_time" or t_type == "orderParameterDist"):
        file = glob.glob(absolutePathList[t_type] + NGT(t_networkSize, t_acceptanceThreshold, t_reapeater))

    #* Read orderparameter-accumulated distributions
    elif (t_type == "clusterSizeDist" or t_type == "clusterSizeDist_exact"):
        file = glob.glob(absolutePathList[t_type] + NGOP(t_networkSize, t_acceptanceThreshold, t_reapeater))

    #* Read general data
    else:
        file = glob.glob(absolutePathList[t_type] + NG(t_networkSize, t_acceptanceThreshold))

    #* Check found files
    if (len(file) != 1):
        print("There is problem at reading " + t_type + " at N={:.1e}".format(t_networkSize) + ", G={:.1f}".format(t_acceptanceThreshold))
        return
    return readCSV(file[0])

if __name__=="__main__":
    print("This is a module readData.py")