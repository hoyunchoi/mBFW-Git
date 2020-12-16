import numpy as np
import pandas as pd
import glob

relativePathList = {}
relativePathList["clusterSizeDist"] = "/logBin/"
relativePathList["clusterSizeDist_exact"] = "/logBin/"
relativePathList["clusterSizeDist_time"] = "/logBin/"
relativePathList["meanClusterSize"] = "/average/"
relativePathList["meanClusterSize_trial"] = "/average/"
relativePathList["orderParameter"] = "/average/"
relativePathList["orderParameter_trial"] = "/average/"
relativePathList["orderParameterVariance"] = "/average/"
relativePathList["orderParameterVariance_trial"] = "/average/"


rootPath = "../data/mBFW_hybrid/"
absolutePathList = {}
for observable, relativePath in relativePathList.items():
    absolutePathList[observable] = rootPath + observable + relativePath


#* CSV Reader
def readCSV(t_filename):
    data = pd.read_csv(t_filename, sep=',', header=None)
    data = data.values.transpose()
    if (len(data) == 1):
        return data[0]
    else:
        return tuple([row for row in data])

#* File Name Convections
def nameGlobFile(t_networkSize, t_acceptanceThreshold):
    return "N{:.1e},G{:.1f}*".format(t_networkSize, t_acceptanceThreshold)

def nameGlobFile_T(t_networkSize, t_acceptanceThreshold, t_time):
    return "N{:.1e},G{:.1f}*,T{:.4f}*".format(t_networkSize, t_acceptanceThreshold, t_time)

def nameGlobFile_OP(t_networkSize, t_acceptanceThreshold, t_orderParameter):
    return "N{:.1e},G{:.1f}*,OP{:.4f}*".format(t_networkSize, t_acceptanceThreshold, t_orderParameter)


#* Read Observables
def read(t_observable, t_networkSize, t_acceptanceThreshold, t_reapeater=None):
    #* Read time-accumulated distributions
    if (t_observable == "clusterSizeDist_time" or t_observable == "orderParameterDist"):
        file = glob.glob(absolutePathList[t_observable] + nameGlobFile_T(t_networkSize, t_acceptanceThreshold, t_reapeater))

    #* Read orderparameter-accumulated distributions
    elif (t_observable == "clusterSizeDist" or t_observable == "clusterSizeDist_exact"):
        file = glob.glob(absolutePathList[t_observable] + nameGlobFile_OP(t_networkSize, t_acceptanceThreshold, t_reapeater))

    #* Read general data
    else:
        file = glob.glob(absolutePathList[t_observable] + nameGlobFile(t_networkSize, t_acceptanceThreshold))

    #* Check found files
    if (len(file) != 1):
        print("There is problem at reading " + t_observable + " at N={:.1e}".format(t_networkSize) + ", G={:.1f}".format(t_acceptanceThreshold))
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