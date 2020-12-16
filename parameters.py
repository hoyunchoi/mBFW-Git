networkSizeList = [1e4, 2e4, 4e4, 8e4, 1.6e5, 3.2e5, 6.4e5, 1.28e6, 2.56e6, 5.12e6, 1.024e7]
acceptanceThresholdList = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

plot_lim = {}
fit_lim = {}
t_c_inf = {}



for networkSize in networkSizeList:
    plot_lim[0.2, networkSize] = [0.99, 0.999]
    plot_lim[0.3, networkSize] = [0.96, 0.99]
    plot_lim[0.4, networkSize] = [0.945, 0.975]
    plot_lim[0.5, networkSize] = [0.92, 0.95]
    plot_lim[0.6, networkSize] = [0.885, 0.915]
    plot_lim[0.7, networkSize] = [0.84, 0.87]
    plot_lim[0.8, networkSize] = [0.78, 0.81]
    plot_lim[0.9, networkSize] = [0.7, 0.73]


#* g=0.2
for networkSize in networkSizeList[-2:]:
    fit_lim[0.2, networkSize] = [0.9953,0.9958]
for networkSize in networkSizeList[5:-2]:
    fit_lim[0.2, networkSize] = [0.9952,0.9958]
for networkSize in networkSizeList[2:5]:
    fit_lim[0.2, networkSize] = [0.995,0.9957]
#! For g=0.2, N=1e4,2e4 critical point from mean cluster size is inaccurate
for networkSize in networkSizeList[:2]:
    fit_lim[0.2, networkSize] = [0.9945,0.9955]

#* g=0.3
for networkSize in networkSizeList[3:]:
    fit_lim[0.3, networkSize] = [0.983,0.985]
fit_lim[0.3, 4e4] = [0.982,0.984]
for networkSize in networkSizeList[:2]:
    fit_lim[0.3, networkSize] = [0.981,0.984]

#* g=0.4
for networkSize in networkSizeList[3:]:
    fit_lim[0.4, networkSize] = [0.9625,0.9655]
for networkSize in networkSizeList[:3]:
    fit_lim[0.4, networkSize]= [0.9605,0.965]

#* g=0.5
t_c_inf[0.5] = 0.93712
for networkSize in networkSizeList[3:]:
    fit_lim[0.5, networkSize] = [0.935,0.939]
for networkSize in networkSizeList[:3]:
    fit_lim[0.5, networkSize] = [0.933,0.937]

#* g=0.6
for networkSize in networkSizeList[3:]:
    fit_lim[0.6, networkSize] = [0.899,0.903]
for networkSize in networkSizeList[1:3]:
    fit_lim[0.6, networkSize] = [0.897,0.902]
fit_lim[0.6, 1e4] = [0.8955,0.9005]

#* g=0.7
for networkSize in networkSizeList[6:]:
    fit_lim[0.7, networkSize] = [0.855, 0.858]
for networkSize in networkSizeList[3:6]:
    fit_lim[0.7, networkSize] = [0.8535, 0.857]
for networkSize in networkSizeList[1:3]:
    fit_lim[0.7, networkSize] = [0.851, 0.857]
fit_lim[0.7, 1e4] = [0.845, 0.855]

#* g= 0.8
for networkSize in networkSizeList[5:]:
    fit_lim[0.8, networkSize] = [0.797, 0.801]
for networkSize in networkSizeList[2:5]:
    fit_lim[0.8, networkSize] = [0.793, 0.799]
for networkSize in networkSizeList[:2]:
    fit_lim[0.8, networkSize] = [0.788, 0.796]

#* g=0.9