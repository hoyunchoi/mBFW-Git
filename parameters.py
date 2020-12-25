networkSizeList = [1e4, 2e4, 4e4, 8e4, 1.6e5, 3.2e5, 6.4e5, 1.28e6, 2.56e6, 5.12e6, 1.024e7]
acceptanceThresholdList = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

t_c_var_inf = {}
t_c_mcs_inf = {}
nu_bar_var = {}
nu_bar_mcs = {}

#* x limit of plotting near critical point
chi_plotRange = {}
chi_plotRange[0.2] = [0.99, 0.999]
chi_plotRange[0.3] = [0.96, 0.99]
chi_plotRange[0.4] = [0.945, 0.975]
chi_plotRange[0.5] = [0.92, 0.95]
chi_plotRange[0.6] = [0.885, 0.915]
chi_plotRange[0.7] = [0.84, 0.87]
chi_plotRange[0.8] = [0.78, 0.81]
chi_plotRange[0.9] = [0.69, 0.73]

#* x limit of fitting Susceptibility(order parameter variance/ mean cluster size) by order 3 polynomial
chi_fitRange = {}
#* g=0.2
t_c_var_inf[0.2] = 0.99549    #? Excluding 2560000, 5120000, 10240000
nu_bar_var[0.2] = 1/0.49527   #? Excluding 2560000, 5120000, 10240000, inaccurate
t_c_mcs_inf[0.2] = 0.99563  #? Excluding 10000, 20000, 560000, 5120000, 10240000
nu_bar_mcs[0.2] = 1/0.77473   #? Excluding 10000, 20000, 560000, 5120000, 10240000
for networkSize in networkSizeList:
    #! For g=0.2, N=1e4,2e4 critical point from mean cluster size is inaccurate
    if networkSize <= 2e4:
        chi_fitRange[networkSize, 0.2] = [0.9945,0.9955]
    elif networkSize <= 1.6e5:
        chi_fitRange[networkSize, 0.2] = [0.995,0.9957]
    elif networkSize <= 2.56e6:
        chi_fitRange[networkSize, 0.2] = [0.9952,0.9958]
    else:
        chi_fitRange[networkSize, 0.2] = [0.9953,0.9958]

#* g=0.3
t_c_var_inf[0.3] = 0.983683     #? Excluding 2560000, 5120000, 10240000
nu_bar_var[0.3] = 1/0.76412   #? Excluding 2560000, 5120000, 10240000
t_c_mcs_inf[0.3] = 0.98369  #? Excluding 5120000, 10240000
nu_bar_mcs[0.3] = 1/0.79163   #? Excluding 5120000, 10240000
for networkSize in networkSizeList:
    if networkSize <= 2e4:
        chi_fitRange[networkSize, 0.3] = [0.981,0.984]
    elif networkSize <= 4e4:
        chi_fitRange[networkSize, 0.3] = [0.982,0.984]
    else:
        chi_fitRange[networkSize, 0.3] = [0.983,0.985]

#* g=0.4
t_c_var_inf[0.4] = 0.964274     #? Excluding 5120000, 10240000
nu_bar_var[0.4] = 1/0.68074   #? Excluding 5120000, 10240000
t_c_mcs_inf[0.4] = 0.96416  #? Excluding 10240000
nu_bar_mcs[0.4] = 1/0.57756   #? Excluding 10240000
for networkSize in networkSizeList:
    if networkSize <= 4e4:
        chi_fitRange[networkSize, 0.4]= [0.9605,0.965]
    else:
        chi_fitRange[networkSize, 0.4] = [0.9625,0.9655]

#* g=0.5
t_c_var_inf[0.5] = 0.93706
nu_bar_var[0.5] = 1/0.68494
t_c_mcs_inf[0.5] = 0.93697  #?
nu_bar_mcs[0.5] = 1/0.48079   #?
for networkSize in networkSizeList:
    if networkSize <= 4e4:
        chi_fitRange[networkSize, 0.5] = [0.933,0.937]
    else:
        chi_fitRange[networkSize, 0.5] = [0.935,0.939]

#* g=0.6
t_c_var_inf[0.6] = 0.90164
nu_bar_var[0.6] = 1/0.64606
t_c_mcs_inf[0.6] = 0.90137  #?
nu_bar_mcs[0.6] = 1/0.55483   #?
for networkSize in networkSizeList:
    if networkSize <= 1e4:
        chi_fitRange[networkSize, 0.6] = [0.8955,0.9005]
    elif networkSize <= 4e4:
        chi_fitRange[networkSize, 0.6] = [0.897,0.902]
    else:
        chi_fitRange[networkSize, 0.6] = [0.899,0.903]

#* g=0.7
t_c_var_inf[0.7] = 0.85666
nu_bar_var[0.7] = 1/0.69739
t_c_mcs_inf[0.7] = 0.85631  #?
nu_bar_mcs[0.7] = 1/0.55880   #?
for networkSize in networkSizeList:
    if networkSize <= 1e4:
        chi_fitRange[networkSize, 0.7] = [0.845, 0.855]
    elif networkSize <= 4e4:
        chi_fitRange[networkSize, 0.7] = [0.851, 0.857]
    elif networkSize <= 3.2e5:
        chi_fitRange[networkSize, 0.7] = [0.8535, 0.857]
    else:
        chi_fitRange[networkSize, 0.7] = [0.855, 0.858]


#* g= 0.8
t_c_var_inf[0.8] = 0.79885
nu_bar_var[0.8] = 1/0.65410
t_c_mcs_inf[0.8] = 0.79837
nu_bar_mcs[0.8] = 1/0.62671
for networkSize in networkSizeList:
    if networkSize <= 2e4:
        chi_fitRange[networkSize, 0.8] = [0.788, 0.796]
    elif networkSize <= 1.6e5:
        chi_fitRange[networkSize, 0.8] = [0.793, 0.799]
    else:
        chi_fitRange[networkSize, 0.8] = [0.797, 0.801]

#* g=0.9
t_c_var_inf[0.9] = 0.71853
nu_bar_var[0.9] = 1/0.58615
t_c_mcs_inf[0.9] = 0.71827
nu_bar_mcs[0.9] = 1/0.54507
for networkSize in networkSizeList:
    if networkSize <= 2e4:
        chi_fitRange[networkSize, 0.9] = [0.701, 0.714]
    elif networkSize <= 8e4:
        chi_fitRange[networkSize, 0.9] =  [0.7095, 0.7165]
    elif networkSize <= 6.4e5:
        chi_fitRange[networkSize, 0.9] = [0.714, 0.718]
    else:
        chi_fitRange[networkSize, 0.9] = [0.7165, 0.719]


#* Cluster Size Distribution
csd_fitRange = {}
m_c_csd = {}

#* g=0.2
csd_fitRange[10000, 0.2] = [11, 20] ; m_c_csd[10000, 0.2] = 0.96
csd_fitRange[20000, 0.2] = [11, 22]; m_c_csd[20000, 0.2] = 0.96
csd_fitRange[40000, 0.2] = [11, 25]; m_c_csd[40000, 0.2] = 0.96
csd_fitRange[80000, 0.2] = [11, 27]; m_c_csd[80000, 0.2] = 0.97
csd_fitRange[160000, 0.2] = [11, 30]; m_c_csd[160000, 0.2] = 0.97
csd_fitRange[320000, 0.2] = [11, 33]; m_c_csd[320000, 0.2] = 0.97
csd_fitRange[640000, 0.2] = [11, 35]; m_c_csd[640000, 0.2] = 0.97
csd_fitRange[1280000, 0.2] = [11, 38]; m_c_csd[1280000, 0.2] = 0.97
csd_fitRange[2560000, 0.2] = [11, 41]; m_c_csd[2560000, 0.2] = 0.97
csd_fitRange[5120000, 0.2] = [11, 44]; m_c_csd[5120000, 0.2] = 0.97
csd_fitRange[10240000, 0.2] = [11, 46]; m_c_csd[10240000, 0.2] = 0.97

# #* g=0.3
# csd_fitRange[10000, 0.3] = []
# csd_fitRange[20000, 0.3] = []
# csd_fitRange[40000, 0.3] = []
# csd_fitRange[80000, 0.3] = []
# csd_fitRange[160000, 0.3] = []
# csd_fitRange[320000, 0.3] = []
# csd_fitRange[640000, 0.3] = []
# csd_fitRange[1280000, 0.3] = []
# csd_fitRange[2560000, 0.3] = []
# csd_fitRange[5120000, 0.3] = []
# csd_fitRange[10240000, 0.3] = []


# #* g=0.4
# csd_fitRange[10000, 0.4] = []
# csd_fitRange[20000, 0.4] = []
# csd_fitRange[40000, 0.4] = []
# csd_fitRange[80000, 0.4] = []
# csd_fitRange[160000, 0.4] = []
# csd_fitRange[320000, 0.4] = []
# csd_fitRange[640000, 0.4] = []
# csd_fitRange[1280000, 0.4] = []
# csd_fitRange[2560000, 0.4] = []
# csd_fitRange[5120000, 0.4] = []
# csd_fitRange[10240000, 0.4] = []

#* g=0.5
csd_fitRange[10000, 0.5] = [11, 21]; m_c_csd[10000, 0.5] = 0.80
csd_fitRange[20000, 0.5] = [11, 23]; m_c_csd[20000, 0.5] = 0.80
csd_fitRange[40000, 0.5] = [11, 25]; m_c_csd[40000, 0.5] = 0.81
csd_fitRange[80000, 0.5] = [11, 27]; m_c_csd[80000, 0.5] = 0.81
csd_fitRange[160000, 0.5] = [11, 29]; m_c_csd[160000, 0.5] = 0.82
csd_fitRange[320000, 0.5] = [11, 31]; m_c_csd[320000, 0.5] = 0.82
csd_fitRange[640000, 0.5] = [11, 33]; m_c_csd[640000, 0.5] = 0.82
csd_fitRange[1280000, 0.5] = [11, 35]; m_c_csd[1280000, 0.5] = 0.82
csd_fitRange[2560000, 0.5] = [11, 36]; m_c_csd[2560000, 0.5] = 0.82
csd_fitRange[5120000, 0.5] = [11, 37]; m_c_csd[5120000, 0.5] = 0.82
csd_fitRange[10240000, 0.5] = [11, 38]; m_c_csd[10240000, 0.5] = 0.82


# #* g=0.6
# csd_fitRange[10000, 0.6] = []
# csd_fitRange[20000, 0.6] = []
# csd_fitRange[40000, 0.6] = []
# csd_fitRange[80000, 0.6] = []
# csd_fitRange[160000, 0.6] = []
# csd_fitRange[320000, 0.6] = []
# csd_fitRange[640000, 0.6] = []
# csd_fitRange[1280000, 0.6] = []
# csd_fitRange[2560000, 0.6] = []
# csd_fitRange[5120000, 0.6] = []
# csd_fitRange[10240000, 0.6] = []

# #* g=0.7
# csd_fitRange[10000, 0.7] = []
# csd_fitRange[20000, 0.7] = []
# csd_fitRange[40000, 0.7] = []
# csd_fitRange[80000, 0.7] = []
# csd_fitRange[160000, 0.7] = []
# csd_fitRange[320000, 0.7] = []
# csd_fitRange[640000, 0.7] = []
# csd_fitRange[1280000, 0.7] = []
# csd_fitRange[2560000, 0.7] = []
# csd_fitRange[5120000, 0.7] = []
# csd_fitRange[10240000, 0.7] = []

# #* g=0.8
# csd_fitRange[10000, 0.8] = []
# csd_fitRange[20000, 0.8] = []
# csd_fitRange[40000, 0.8] = []
# csd_fitRange[80000, 0.8] = []
# csd_fitRange[160000, 0.8] = []
# csd_fitRange[320000, 0.8] = []
# csd_fitRange[640000, 0.8] = []
# csd_fitRange[1280000, 0.8] = []
# csd_fitRange[2560000, 0.8] = []
# csd_fitRange[5120000, 0.8] = []
# csd_fitRange[10240000, 0.8] = []

# #* g=0.9
# csd_fitRange[10000, 0.9] = []
# csd_fitRange[20000, 0.9] = []
# csd_fitRange[40000, 0.9] = []
# csd_fitRange[80000, 0.9] = []
# csd_fitRange[160000, 0.9] = []
# csd_fitRange[320000, 0.9] = []
# csd_fitRange[640000, 0.9] = []
# csd_fitRange[1280000, 0.9] = []
# csd_fitRange[2560000, 0.9] = []
# csd_fitRange[5120000, 0.9] = []
# csd_fitRange[10240000, 0.9] = []
