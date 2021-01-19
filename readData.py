import numpy as np
import pandas as pd
import glob

relativePathList = {}
relativePathList["ageDist_time/before"] = "/logBin/"
relativePathList["ageDist_time/during"] = "/logBin/"
relativePathList["ageDist_time/after"] = "/logBin/"
relativePathList["clusterSizeDist"] = "/logBin/"
relativePathList["clusterSizeDist_exact"] = "/logBin/"
relativePathList["clusterSizeDist_time"] = "/logBin/"
relativePathList["deltaUpperBoundDist_time/before"] = "/logBin/"
relativePathList["deltaUpperBoundDist_time/during"] = "/logBin/"
relativePathList["deltaUpperBoundDist_time/after"] = "/logBin/"
relativePathList["interEventTime"] = "/average/"
relativePathList["interEventTimeDist_time/before"] = "/logBin/"
relativePathList["interEventTimeDist_time/during"] = "/logBin/"
relativePathList["interEventTimeDist_time/after"] = "/logBin/"
relativePathList["meanClusterSize"] = "/average/"
relativePathList["meanClusterSize_trial"] = "/average/"
relativePathList["orderParameter"] = "/average/"
relativePathList["orderParameter_trial"] = "/average/"
relativePathList["orderParameterDist"] = "/linBin/"
relativePathList["orderParameterVariance"] = "/average/"
relativePathList["orderParameterVariance_trial"] = "/average/"



rootPath = "../data/mBFW_hybrid/"
absolutePathList = {}
for observable, relativePath in relativePathList.items():
    absolutePathList[observable] = rootPath + observable + relativePath


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


#* Read t_a
def readt_a(t_networkSize, t_acceptanceThreshold):
    fileName = rootPath + "t_a.txt"
    target = "N{:.1e},G{:.1f}".format(t_networkSize, t_acceptanceThreshold)
    with open(fileName) as file:
        content = file.readlines()
        for line in content:
            if target in line:
                line = line[line.find("inflection: ")+12 : ]
                inflectionTime = line[:line.find(",")]
                inflectionOP = line[line.find(",")+1 : line.find("\t")]
                line = line[line.find("\t") : ]
                t_a = line[line.find("t_a:")+4 : line.find(",")]
                m_a = line[line.find("m_a:")+4 :]
    return float(inflectionTime), float(inflectionOP), float(t_a), float(m_a)

#* Read t_c
def readt_c(t_networkSize, t_acceptanceThreshold):
    fileName = rootPath + "t_c.txt"
    target = "N{:.1e},G{:.1f}".format(t_networkSize, t_acceptanceThreshold)
    t_c_var = 0
    m_c_var = 0
    t_c_mcs = 0
    m_c_mcs = 0

    with open(fileName) as file:
        content = file.readlines()
        for line in content:
            if target in line:
                line = line[line.find("\t"):]
                if "var" in line:
                    t_c_var = line[line.find("t_c_var:")+8 : line.find(",")]
                    m_c_var = line[line.find("m_c:")+4 : -1]
                elif "mcs" in line:
                    t_c_mcs = line[line.find("t_c_mcs:")+8 : line.find(",")]
                    m_c_mcs = line[line.find("m_c:")+4 : -1]
    return float(t_c_var), float(m_c_var), float(t_c_mcs), float(m_c_mcs)

if __name__=="__main__":
    print("This is a module readData.py")