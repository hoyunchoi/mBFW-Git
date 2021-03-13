from readData import read
import numpy as np

networkSizeList = [1e4, 2e4, 4e4, 8e4, 1.6e5, 3.2e5, 6.4e5, 1.28e6, 2.56e6, 5.12e6, 1.024e7]
acceptanceThresholdList = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# * Plot range for near critical point
chi_plotRange = {}
chi_plotRange[0.2] = [0.99, 0.999]
chi_plotRange[0.3] = [0.97, 0.99]
chi_plotRange[0.4] = [0.945, 0.975]
chi_plotRange[0.5] = [0.92, 0.95]
chi_plotRange[0.6] = [0.885, 0.915]
chi_plotRange[0.7] = [0.84, 0.87]
chi_plotRange[0.8] = [0.78, 0.81]
chi_plotRange[0.9] = [0.69, 0.73]

# *---------------------------------------------------------------------------------------------------
# * Fit range for chi
chi_fitRange = {}

#! G=0.2
chi_fitRange[10000, 0.2] = [0.9935, 0.9955]
chi_fitRange[20000, 0.2] = [0.994, 0.9955]
chi_fitRange[40000, 0.2] = [0.995, 0.9957]
chi_fitRange[80000, 0.2] = [0.995, 0.9957]
chi_fitRange[160000, 0.2] = [0.995, 0.9957]
chi_fitRange[320000, 0.2] = [0.9952, 0.9958]
chi_fitRange[640000, 0.2] = [0.9952, 0.9958]
chi_fitRange[1280000, 0.2] = [0.9952, 0.9958]
chi_fitRange[2560000, 0.2] = [0.9952, 0.9958]
chi_fitRange[5120000, 0.2] = [0.9953, 0.9958]
chi_fitRange[10240000, 0.2] = [0.9953, 0.9958]

#! G=0.3
chi_fitRange[10000, 0.3] = [0.981, 0.984]
chi_fitRange[20000, 0.3] = [0.981, 0.984]
chi_fitRange[40000, 0.3] = [0.982, 0.984]
chi_fitRange[80000, 0.3] = [0.983, 0.9845]
chi_fitRange[160000, 0.3] = [0.983, 0.9845]
chi_fitRange[320000, 0.3] = [0.983, 0.9845]
chi_fitRange[640000, 0.3] = [0.983, 0.9845]
chi_fitRange[1280000, 0.3] = [0.983, 0.9845]
chi_fitRange[2560000, 0.3] = [0.9833, 0.9845]
chi_fitRange[5120000, 0.3] = [0.9833, 0.9845]
chi_fitRange[10240000, 0.3] = [0.9833, 0.9845]

#! G=0.4
chi_fitRange[10000, 0.4] = [0.9605, 0.965]
chi_fitRange[20000, 0.4] = [0.9605, 0.965]
chi_fitRange[40000, 0.4] = [0.9605, 0.965]
chi_fitRange[80000, 0.4] = [0.9625, 0.9655]
chi_fitRange[160000, 0.4] = [0.9625, 0.9655]
chi_fitRange[320000, 0.4] = [0.9630, 0.9655]
chi_fitRange[640000, 0.4] = [0.9630, 0.9655]
chi_fitRange[1280000, 0.4] = [0.9630, 0.9655]
chi_fitRange[2560000, 0.4] = [0.9630, 0.9655]
chi_fitRange[5120000, 0.4] = [0.9630, 0.9655]
chi_fitRange[10240000, 0.4] = [0.9633, 0.9655]

#! G=0.5
chi_fitRange[10000, 0.5] = [0.933, 0.937]
chi_fitRange[20000, 0.5] = [0.933, 0.937]
chi_fitRange[40000, 0.5] = [0.933, 0.937]
chi_fitRange[80000, 0.5] = [0.935, 0.938]
chi_fitRange[160000, 0.5] = [0.935, 0.938]
chi_fitRange[320000, 0.5] = [0.935, 0.938]
chi_fitRange[640000, 0.5] = [0.9355, 0.938]
chi_fitRange[1280000, 0.5] = [0.9355, 0.938]
chi_fitRange[2560000, 0.5] = [0.9355, 0.938]
chi_fitRange[5120000, 0.5] = [0.9355, 0.938]
chi_fitRange[10240000, 0.5] = [0.9355, 0.938]

#! G=0.6
chi_fitRange[10000, 0.6] = [0.8955, 0.9005]
chi_fitRange[20000, 0.6] = [0.897, 0.902]
chi_fitRange[40000, 0.6] = [0.897, 0.902]
chi_fitRange[80000, 0.6] = [0.899, 0.903]
chi_fitRange[160000, 0.6] = [0.899, 0.903]
chi_fitRange[320000, 0.6] = [0.899, 0.903]
chi_fitRange[640000, 0.6] = [0.899, 0.903]
chi_fitRange[1280000, 0.6] = [0.899, 0.903]
chi_fitRange[2560000, 0.6] = [0.899, 0.903]
chi_fitRange[5120000, 0.6] = [0.899, 0.903]
chi_fitRange[10240000, 0.6] = [0.899, 0.903]

#! G=0.7
chi_fitRange[10000, 0.7] = [0.845, 0.855]
chi_fitRange[20000, 0.7] = [0.851, 0.857]
chi_fitRange[40000, 0.7] = [0.851, 0.857]
chi_fitRange[80000, 0.7] = [0.8535, 0.857]
chi_fitRange[160000, 0.7] = [0.8535, 0.857]
chi_fitRange[320000, 0.7] = [0.8535, 0.857]
chi_fitRange[640000, 0.7] = [0.855, 0.858]
chi_fitRange[1280000, 0.7] = [0.855, 0.858]
chi_fitRange[2560000, 0.7] = [0.855, 0.858]
chi_fitRange[5120000, 0.7] = [0.855, 0.858]
chi_fitRange[10240000, 0.7] = [0.855, 0.858]

#! G=0.8
chi_fitRange[10000, 0.8] = [0.788, 0.796]
chi_fitRange[20000, 0.8] = [0.788, 0.796]
chi_fitRange[40000, 0.8] = [0.793, 0.799]
chi_fitRange[80000, 0.8] = [0.793, 0.799]
chi_fitRange[160000, 0.8] = [0.793, 0.799]
chi_fitRange[320000, 0.8] = [0.797, 0.801]
chi_fitRange[640000, 0.8] = [0.797, 0.801]
chi_fitRange[1280000, 0.8] = [0.797, 0.801]
chi_fitRange[2560000, 0.8] = [0.797, 0.801]
chi_fitRange[5120000, 0.8] = [0.797, 0.801]
chi_fitRange[10240000, 0.8] = [0.797, 0.801]

#! G=0.9
chi_fitRange[10000, 0.9] = [0.701, 0.714]
chi_fitRange[20000, 0.9] = [0.701, 0.714]
chi_fitRange[40000, 0.9] = [0.7095, 0.7165]
chi_fitRange[80000, 0.9] = [0.7095, 0.7165]
chi_fitRange[160000, 0.9] = [0.714, 0.718]
chi_fitRange[320000, 0.9] = [0.714, 0.718]
chi_fitRange[640000, 0.9] = [0.714, 0.718]
chi_fitRange[1280000, 0.9] = [0.7165, 0.719]
chi_fitRange[2560000, 0.9] = [0.7165, 0.719]
chi_fitRange[5120000, 0.9] = [0.7165, 0.719]
chi_fitRange[10240000, 0.9] = [0.7165, 0.719]

# *---------------------------------------------------------------------------------------------------
# * Critical point in thermodynamic limit and nu_bar
t_c_var_inf = {}
nu_bar_var = {}
t_c_mcs_inf = {}
nu_bar_mcs = {}

#! G=0.2
t_c_var_inf[0.2] = 0.995474
nu_bar_var[0.2] = 1/0.52155
t_c_mcs_inf[0.2] = 0.995635
nu_bar_mcs[0.2] = 1/0.68175  # ? Excluding 10000, 20000, 2560000, 5120000, 10240000

#! G=0.3
t_c_var_inf[0.3] = 0.983737
nu_bar_var[0.3] = 1/0.70173  # ? Excluding 2560000, 5120000, 10240000
t_c_mcs_inf[0.3] = 0.983689
nu_bar_mcs[0.3] = 1/0.78967  # ? Excluding 2560000, 5120000, 10240000

#! G=0.4
t_c_var_inf[0.4] = 0.964284
nu_bar_var[0.4] = 1/0.65864  # ? Excluding 5120000, 10240000
t_c_mcs_inf[0.4] = 0.964058
nu_bar_mcs[0.4] = 1/0.62966

#! G=0.5
t_c_var_inf[0.5] = 0.937061
nu_bar_var[0.5] = 1/0.68135  # ? Excluding 10240000
t_c_mcs_inf[0.5] = 0.936897
nu_bar_mcs[0.5] = 1/0.59502  # ? Excluding 10240000
# ------------------------------------

#! G=0.6
t_c_var_inf[0.6] = 0.901616
nu_bar_var[0.6] = 1/0.69724
t_c_mcs_inf[0.6] = 0.901384
nu_bar_mcs[0.6] = 1/0.54742

#! G=0.7
t_c_var_inf[0.7] = 0.856656
nu_bar_var[0.7] = 1/0.71512
t_c_mcs_inf[0.7] = 0.856363
nu_bar_mcs[0.7] = 1/0.55401  # ? excluding 2560000

#! G=0.8
t_c_var_inf[0.8] = 0.798838
nu_bar_var[0.8] = 1/0.66728
t_c_mcs_inf[0.8] = 0.798378
nu_bar_mcs[0.8] = 1/0.62711

#! G=0.9
t_c_var_inf[0.9] = 0.718534
nu_bar_var[0.9] = 1/0.58429
t_c_mcs_inf[0.9] = 0.718161
nu_bar_mcs[0.9] = 1/0.56890

# *---------------------------------------------------------------------------------------------------
# * Fit range for cluster size distribution
csd_fitRange = {}
m_c_csd = {}

#! G=0.2
csd_fitRange[10000, 0.2] = [11, 20]
m_c_csd[10000, 0.2] = 0.962
csd_fitRange[20000, 0.2] = [11, 22]
m_c_csd[20000, 0.2] = 0.966
csd_fitRange[40000, 0.2] = [11, 25]
m_c_csd[40000, 0.2] = 0.967
csd_fitRange[80000, 0.2] = [11, 27]
m_c_csd[80000, 0.2] = 0.969
csd_fitRange[160000, 0.2] = [11, 30]
m_c_csd[160000, 0.2] = 0.969
csd_fitRange[320000, 0.2] = [11, 32]
m_c_csd[320000, 0.2] = 0.971
csd_fitRange[640000, 0.2] = [11, 34]
m_c_csd[640000, 0.2] = 0.972
csd_fitRange[1280000, 0.2] = [11, 36]
m_c_csd[1280000, 0.2] = 0.973
csd_fitRange[2560000, 0.2] = [11, 39]
m_c_csd[2560000, 0.2] = 0.972
csd_fitRange[5120000, 0.2] = [11, 41]
m_c_csd[5120000, 0.2] = 0.973
csd_fitRange[10240000, 0.2] = [11, 43]
m_c_csd[10240000, 0.2] = 0.972

#! G=0.3
csd_fitRange[10000, 0.3] = [11, 20]
m_c_csd[10000, 0.3] = 0.92
csd_fitRange[20000, 0.3] = [11, 22]
m_c_csd[20000, 0.3] = 0.93
csd_fitRange[40000, 0.3] = [11, 24]
m_c_csd[40000, 0.3] = 0.93
csd_fitRange[80000, 0.3] = [11, 26]
m_c_csd[80000, 0.3] = 0.93
csd_fitRange[160000, 0.3] = [11, 28]
m_c_csd[160000, 0.3] = 0.93
csd_fitRange[320000, 0.3] = [11, 31]
m_c_csd[320000, 0.3] = 0.93
csd_fitRange[640000, 0.3] = [11, 34]
m_c_csd[640000, 0.3] = 0.93
csd_fitRange[1280000, 0.3] = [11, 36]
m_c_csd[1280000, 0.3] = 0.93
csd_fitRange[2560000, 0.3] = [11, 38]
m_c_csd[2560000, 0.3] = 0.93
csd_fitRange[5120000, 0.3] = [11, 40]
m_c_csd[5120000, 0.3] = 0.93
csd_fitRange[10240000, 0.3] = [11, 42]
m_c_csd[10240000, 0.3] = 0.93

#! G=0.4
csd_fitRange[10000, 0.4] = [11, 21]
m_c_csd[10000, 0.4] = 0.850
csd_fitRange[20000, 0.4] = [11, 23]
m_c_csd[20000, 0.4] = 0.865
csd_fitRange[40000, 0.4] = [11, 25]
m_c_csd[40000, 0.4] = 0.875
csd_fitRange[80000, 0.4] = [11, 27]
m_c_csd[80000, 0.4] = 0.880
csd_fitRange[160000, 0.4] = [11, 28]
m_c_csd[160000, 0.4] = 0.880
csd_fitRange[320000, 0.4] = [11, 30]
m_c_csd[320000, 0.4] = 0.880
csd_fitRange[640000, 0.4] = [11, 31]
m_c_csd[640000, 0.4] = 0.885
csd_fitRange[1280000, 0.4] = [11, 32]
m_c_csd[1280000, 0.4] = 0.885
csd_fitRange[2560000, 0.4] = [11, 33]
m_c_csd[2560000, 0.4] = 0.885
csd_fitRange[5120000, 0.4] = [11, 34]
m_c_csd[5120000, 0.4] = 0.885
csd_fitRange[10240000, 0.4] = [11, 35]
m_c_csd[10240000, 0.4] = 0.885

#! G=0.5
csd_fitRange[10000, 0.5] = [11, 22]
m_c_csd[10000, 0.5] = 0.780
csd_fitRange[20000, 0.5] = [11, 24]
m_c_csd[20000, 0.5] = 0.796
csd_fitRange[40000, 0.5] = [11, 26]
m_c_csd[40000, 0.5] = 0.806
csd_fitRange[80000, 0.5] = [11, 27]
m_c_csd[80000, 0.5] = 0.815
csd_fitRange[160000, 0.5] = [11, 29]
m_c_csd[160000, 0.5] = 0.818
csd_fitRange[320000, 0.5] = [11, 31]
m_c_csd[320000, 0.5] = 0.818
csd_fitRange[640000, 0.5] = [11, 33]
m_c_csd[640000, 0.5] = 0.818
csd_fitRange[1280000, 0.5] = [11, 35]
m_c_csd[1280000, 0.5] = 0.817
csd_fitRange[2560000, 0.5] = [11, 36]
m_c_csd[2560000, 0.5] = 0.817
csd_fitRange[5120000, 0.5] = [11, 37]
m_c_csd[5120000, 0.5] = 0.817
csd_fitRange[10240000, 0.5] = [11, 38]
m_c_csd[10240000, 0.5] = 0.817

#! G=0.6
csd_fitRange[10000, 0.6] = [11, 22]
m_c_csd[10000, 0.6] = 0.690
csd_fitRange[20000, 0.6] = [11, 24]
m_c_csd[20000, 0.6] = 0.710
csd_fitRange[40000, 0.6] = [11, 26]
m_c_csd[40000, 0.6] = 0.720
csd_fitRange[80000, 0.6] = [11, 27]
m_c_csd[80000, 0.6] = 0.740
csd_fitRange[160000, 0.6] = [11, 29]
m_c_csd[160000, 0.6] = 0.740
csd_fitRange[320000, 0.6] = [11, 31]
m_c_csd[320000, 0.6] = 0.740
csd_fitRange[640000, 0.6] = [11, 33]
m_c_csd[640000, 0.6] = 0.740
csd_fitRange[1280000, 0.6] = [11, 35]
m_c_csd[1280000, 0.6] = 0.740
csd_fitRange[2560000, 0.6] = [11, 37]
m_c_csd[2560000, 0.6] = 0.740
csd_fitRange[5120000, 0.6] = [11, 39]
m_c_csd[5120000, 0.6] = 0.740
csd_fitRange[10240000, 0.6] = [11, 41]
m_c_csd[10240000, 0.6] = 0.740

#! G=0.7
csd_fitRange[10000, 0.7] = [11, 22]
m_c_csd[10000, 0.7] = 0.580
csd_fitRange[20000, 0.7] = [11, 24]
m_c_csd[20000, 0.7] = 0.610
csd_fitRange[40000, 0.7] = [11, 26]
m_c_csd[40000, 0.7] = 0.630
csd_fitRange[80000, 0.7] = [11, 28]
m_c_csd[80000, 0.7] = 0.640
csd_fitRange[160000, 0.7] = [11, 30]
m_c_csd[160000, 0.7] = 0.650
csd_fitRange[320000, 0.7] = [11, 31]
m_c_csd[320000, 0.7] = 0.650
csd_fitRange[640000, 0.7] = [11, 33]
m_c_csd[640000, 0.7] = 0.650
csd_fitRange[1280000, 0.7] = [11, 35]
m_c_csd[1280000, 0.7] = 0.650
csd_fitRange[2560000, 0.7] = [11, 37]
m_c_csd[2560000, 0.7] = 0.650
csd_fitRange[5120000, 0.7] = [11, 39]
m_c_csd[5120000, 0.7] = 0.650
csd_fitRange[10240000, 0.7] = [11, 41]
m_c_csd[10240000, 0.7] = 0.650

#! G=0.8
csd_fitRange[10000, 0.8] = [11, 22]
m_c_csd[10000, 0.8] = 0.430
csd_fitRange[20000, 0.8] = [11, 24]
m_c_csd[20000, 0.8] = 0.450
csd_fitRange[40000, 0.8] = [11, 26]
m_c_csd[40000, 0.8] = 0.500
csd_fitRange[80000, 0.8] = [11, 28]
m_c_csd[80000, 0.8] = 0.520
csd_fitRange[160000, 0.8] = [11, 30]
m_c_csd[160000, 0.8] = 0.520
csd_fitRange[320000, 0.8] = [11, 31]
m_c_csd[320000, 0.8] = 0.540
csd_fitRange[640000, 0.8] = [11, 33]
m_c_csd[640000, 0.8] = 0.540
csd_fitRange[1280000, 0.8] = [11, 35]
m_c_csd[1280000, 0.8] = 0.540
csd_fitRange[2560000, 0.8] = [11, 37]
m_c_csd[2560000, 0.8] = 0.540
csd_fitRange[5120000, 0.8] = [11, 39]
m_c_csd[5120000, 0.8] = 0.540
csd_fitRange[10240000, 0.8] = [11, 41]
m_c_csd[10240000, 0.8] = 0.540

#! G=0.9
csd_fitRange[10000, 0.9] = [11, 22]
m_c_csd[10000, 0.9] = 0.230
csd_fitRange[20000, 0.9] = [11, 23]
m_c_csd[20000, 0.9] = 0.300
csd_fitRange[40000, 0.9] = [11, 25]
m_c_csd[40000, 0.9] = 0.300
csd_fitRange[80000, 0.9] = [11, 27]
m_c_csd[80000, 0.9] = 0.360
csd_fitRange[160000, 0.9] = [11, 29]
m_c_csd[160000, 0.9] = 0.380
csd_fitRange[320000, 0.9] = [11, 31]
m_c_csd[320000, 0.9] = 0.380
csd_fitRange[640000, 0.9] = [11, 33]
m_c_csd[640000, 0.9] = 0.380
csd_fitRange[1280000, 0.9] = [11, 35]
m_c_csd[1280000, 0.9] = 0.380
csd_fitRange[2560000, 0.9] = [11, 36]
m_c_csd[2560000, 0.9] = 0.380
csd_fitRange[5120000, 0.9] = [11, 38]
m_c_csd[5120000, 0.9] = 0.380
csd_fitRange[10240000, 0.9] = [11, 40]
m_c_csd[10240000, 0.9] = 0.380


def findTa(networkSize, acceptanceThreshold):
    op = read("orderParameter", networkSize, acceptanceThreshold)
    slope = op[1:] - op[:-1]
    maxIndex = np.argmax(slope)
    inflection_t = maxIndex / networkSize
    inflection_m = op[maxIndex]
    maxSlope = slope[maxIndex] * networkSize
    t_a = inflection_t - maxSlope/inflection_m
    m_a = op[int(t_a*networkSize)]
    with open("../data/mBFW/points/" + "N{:.1e},G{:.1f}".format(networkSize, acceptanceThreshold) + ".txt", 'a') as file:
        file.write("t_a: {:.15f}".format(t_a))
        file.write("m_a: {:.15f}".format(m_a))
        file.write("t_inflection: {:.15f}".format(inflection_t))
        file.write("m_inflection: {:.15f}".format(inflection_m))


# *---------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    # print("This is a module parameters.py")
    import readData

    for acceptanceThreshold in acceptanceThresholdList:
        for networkSize in networkSizeList:
            orderParameter = readData.read("orderParameter", networkSize, acceptanceThreshold)
            t_c_csd = 0.0
            for index, op in enumerate(orderParameter):
                if (op > m_c_csd[networkSize, acceptanceThreshold]):
                    t_c_csd = index/networkSize
                    break

            with open("../data/mBFW/points/" + "N{:.1e},G{:.1f}".format(networkSize, acceptanceThreshold) + ".txt", 'a') as file:
                file.write("t_c_csd: " + str(t_c_csd) + "\n")
                file.write("m_c_csd: " + str(m_c_csd[networkSize, acceptanceThreshold]) + "\n")
