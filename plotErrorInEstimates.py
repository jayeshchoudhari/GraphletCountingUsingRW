from __future__ import division
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

# graphType = "g32"
graphType = "g46"

graphName = "orkut"
# graphName = "sinaweibo"

exact = {"g32" : {"orkut": 524643952, "sinaweibo": 212977684}, "g46" : {"orkut": 2427699127, "sinaweibo": 662717336}}

graphTypeToName = {"g32": "Triangle", "g46": "4-Clique"}

staticExact = exact[graphType][graphName]

allEstimates = defaultdict(list)

percInd = 0

# f = open("./allOutput/4-Clique-sinaweibo-singleRW.out")
# f = open("allOutput/4-Clique-sinaweibo-high-newStart.out")
# f = open("allOutput/4-Clique-sinaweibo-high.out")

# f = open("./allOutput/3-Clique-orkut-singleRW.out")
f = open("./allOutput/4-Clique-orkut-singleRW.out")
# f = open("allOutput/4-Clique-orkut.out")

xarr = []
yarr = []
fivePercentLineXarr = []
fivePercentLineYarr = []
fivePercentLineYarrNeg = []
exactCurve = []


varXarr = []
varYarr = []
valuesYarr = []
valuesStd = []

prevFractionEdges = 0

for line in f:

    if line[0] == "%":
        flds = line.strip().split()
        fractionEdges = float(flds[-1])
        if percInd > 0:
            # print(xarr)
            # print(yarr)
            if xarr[0] > 0.0:
                plt.plot(xarr, yarr, 'o', label=str(prevFractionEdges)+"%")
                # plt.plot(xarr[0], np.std(valuesYarr), 'o', label=str(prevFractionEdges)+"%")
                fivePercentLineXarr.append(xarr[0])
                fivePercentLineYarr.append(5.0)
                fivePercentLineYarrNeg.append(-5.0)
                exactCurve.append(0.0)
                varXarr.append(xarr[0])
                varYarr.append(np.std(yarr))
                valuesStd.append(np.std(valuesYarr))
            # allEstimates[xarr[0]] = yarr
        percInd += 1
        xarr = []
        yarr = []
        valuesYarr = []
        prevFractionEdges = fractionEdges
        

    else:
        flds = line.strip().split()
        estVal = float(flds[2])
        # if estVal > 0:
        xarr.append(fractionEdges)
        err = 100 * (estVal - staticExact)/staticExact
        yarr.append(err)
        valuesYarr.append(estVal)

        if err > 100:
            print("Here ---- ", estVal)
    

plt.plot(xarr, yarr, 'o', label=str(prevFractionEdges)+"%")
# plt.plot(xarr[0], np.std(valuesYarr), 'o', label=str(prevFractionEdges)+"%")
fivePercentLineXarr.append(xarr[0])
fivePercentLineYarr.append(5.0)
fivePercentLineYarrNeg.append(-5.0)
exactCurve.append(0.0)
varXarr.append(xarr[0])
varYarr.append(np.std(yarr))
valuesStd.append(np.std(valuesYarr))

plt.plot(fivePercentLineXarr, fivePercentLineYarr, 'r--', label="5% Error")
plt.plot(fivePercentLineXarr, fivePercentLineYarrNeg, 'r--', label="-5% Error")
plt.plot(fivePercentLineXarr, exactCurve, 'r--*', label="Exact")
# plt.plot(varXarr, varYarr, 'kD-', label="Std. Dev.")
# plt.plot(varXarr, valuesStd, 'k*-', label="Std. Dev.(Estimate)")

f.close()

plt.title("soc " + graphName + "-- Error in Estimate for " + graphTypeToName[graphType])
plt.ylabel("% Error in Estimate")
plt.xlabel("% Edges Observed")
plt.legend(loc = "upper right")

plt.show()