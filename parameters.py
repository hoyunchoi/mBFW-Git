import numpy as np
from os import walk

networkSizeList = [1e4, 2e4, 4e4, 8e4, 1.6e5, 3.2e5, 6.4e5, 1.28e6, 2.56e6, 5.12e6, 1.024e7]
N_strToInt = {"1.0e+04":10000, "2.0e+04":20000, "4.0e+04":40000, "8.0e+04":80000, "1.6e+05":160000, "3.2e+05":320000, "6.4e+05":640000, "1.3e+06":1280000, "2.6e+06":2560000, "5.1e+06":5120000, "1.0e+07":10240000}
acceptanceThresholdList = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
time = {}
for networkSize in networkSizeList:
    deltaT = 1/networkSize
    time[networkSize] = np.arange(0, 1, 1/networkSize)

orderParameter_clusterSizeDist = {}
directory = "../data/mBFW_hybrid/clusterSizeDist/logBin/"
for _, _, fileNameList in walk(directory):
    for fileName in fileNameList:
        networkSize = N_strToInt[fileName[1:8]]
        acceptanceThreshold = float(fileName[10:13])
        op = float(fileName[fileName.find("OP")+2 : -4])

        if (networkSize, acceptanceThreshold) in orderParameter_clusterSizeDist:
            orderParameter_clusterSizeDist[networkSize, acceptanceThreshold].add(op)
        else:
            orderParameter_clusterSizeDist[networkSize, acceptanceThreshold] = set([op])

time_clusterSizeDist = {}
directory = "../data/mBFW_hybrid/clusterSizeDist_time/logBin/"
for _, _, fileNameList in walk(directory):
    for fileName in fileNameList:
        networkSize = N_strToInt[fileName[1:8]]
        acceptanceThreshold = float(fileName[10:13])
        t = float(fileName[fileName.find(",T")+2 : -4])

        if (networkSize, acceptanceThreshold) in time_clusterSizeDist:
            time_clusterSizeDist[networkSize, acceptanceThreshold].add(t)
        else:
            time_clusterSizeDist[networkSize, acceptanceThreshold] = set([t])

nearCritical = {}
nearCritical[0.2] = [0.0, 0.0]
nearCritical[0.3] = [0.0, 0.0]
nearCritical[0.4] = [0.0, 0.0]
nearCritical[0.5] = [0.9, 0.95]
nearCritical[0.6] = [0.0, 0.0]
nearCritical[0.7] = [0.0, 0.0]
nearCritical[0.8] = [0.0, 0.0]
nearCritical[0.9] = [0.0, 0.0]
nearCritical[1.0] = [0.0, 0.0]

clusterSizeDist_fitRange = {}
clusterSizeDist_fitRange[10000, 0.5] = [10, 23]
clusterSizeDist_fitRange[20000, 0.5] = [10, 25]
clusterSizeDist_fitRange[40000, 0.5] = [10, 26]
clusterSizeDist_fitRange[80000, 0.5] = [10, 28]
clusterSizeDist_fitRange[160000, 0.5] = [10, 30]
clusterSizeDist_fitRange[320000, 0.5] = [10, 31]
clusterSizeDist_fitRange[640000, 0.5] = [10, 33]
clusterSizeDist_fitRange[1280000, 0.5] = [10, 35]
clusterSizeDist_fitRange[2560000, 0.5] = [10, 36]
clusterSizeDist_fitRange[5120000, 0.5] = [10, 37]
clusterSizeDist_fitRange[10240000, 0.5] = [10, 38]





