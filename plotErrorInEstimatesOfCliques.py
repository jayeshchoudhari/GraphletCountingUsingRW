from __future__ import division
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np


maxK = int(sys.argv[1])
graphTypeArr = ["g32", "g46", "g510"]
cliqueType = {3: '3-Clique', 4: '4-Clique', 5: '5-Clique', 6: '6-Clique'}
allKVals = range(3, maxK+1)


print(allKVals)

# graphType = "g32"
# graphType = "g46"
# graphType = "g45"
graphType = "g510"
# graphType = "g615"

# graphName = "orkut"
graphName = "sinaweibo"
# graphName = "twitter-higgs"

exact = {"g32" : {"orkut": 524643952, "sinaweibo": 212977684, "twitter-higgs": 83023455}, "g46" : {"orkut": 2427699127, "sinaweibo": 662717336, "twitter-higgs": 429733066}, "g45" : {"orkut": 33254726607, "sinaweibo": 27186599063, "twitter-higgs": 23813747782}, "g510" : {"orkut": 10781372075, "sinaweibo": 3279184432, "twitter-higgs":2170177171}, "g59" : {"orkut": 84290005191, "sinaweibo": 46555727590, "twitter-higgs":37220621763}, "g615": {"orkut": 0, "sinaweibo": 19324925499, "twitter-higgs":0}}

graphTypeToName = {"g32": "Triangle", "g45": "4-Chord-Cycle", "g46": "4-Clique", "g510": "5-Clique", "g59": "Almost-5-Clique", "g615": "6-Clique"}

staticExact = exact[graphType][graphName]

allEstimates = defaultdict(list)

percInd = 0


# f = open("./allOutput/demet-twitter_higgs_g510_IntCliqueEstimates_RW.out")
f = open("./allOutput/demet-sinaweibo_g510_IntCliqueEstimates_RW.out")
# f = open("./allOutput/demet-orkut_g510_IntCliqueEstimates_RW.out")

# f = open("./allOutput/demet-sinaweibo_g615_80_100_100_RW.out")
# f = open("./allOutput/demet-sinaweibo_g510_80_100_100_RW.out")
# f = open("./allOutput/demet-sinaweibo_g510_80_100_100_UAR.out")
# f = open("./allOutput/demet-sinaweibo_g510_80_100_100.out")
# f = open("./allOutput/demet-sinaweibo_g45_80_100.out")
# f = open("./allOutput/demet-sinaweibo_g46_70_70.out")
# f = open("./allOutput/demet-sinaweibo_g46_60_60.out")
# f = open("./allOutput/4-Clique-sinaweibo-moreEdges.out")
# f = open("./allOutput/3-Clique-sinaweibo-DiffStartPoint.out")
# f = open("./allOutput/4-Clique-sinaweibo-singleRW.out")
# f = open("allOutput/4-Clique-sinaweibo-high-newStart.out")
# f = open("allOutput/4-Clique-sinaweibo-high.out")

# f = open("./allOutput/demet-orkut_g510_80_100_100_UAR.out")
# f = open("./allOutput/demet-orkut_g510_80_100_100.out")
# f = open("./allOutput/demet-orkut_g45_80_100.out")
# f = open("./allOutput/outputFiles/3-Clique-orkut-DiffStartPoint.out")
# f = open("./allOutput/3-Clique-orkut-moreEdges.out")
# f = open("./allOutput/4-Clique-orkut-moreEdges.out")
# f = open("./allOutput/3-Clique-orkut-singleRW.out")
# f = open("./allOutput/4-Clique-orkut-singleRW.out")
# f = open("allOutput/4-Clique-orkut.out")


percentageWiseRelErrXYArr = defaultdict(lambda: defaultdict(list))


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


medXarr = []
medYarr = []

for line in f:

    if line[0] == "%":
        flds = line.strip().split()
        fractionEdges = float(flds[-1])
        if percInd > 0:
            # print(xarr)
            # print(yarr)
            if xarr[0] > 0.0:
                # plt.plot(xarr, yarr, 'o', label=str(prevFractionEdges)+"%")

                for kval in allKVals:
                    # print(len(xarr), )
                    percentageWiseRelErrXYArr[kval]['xarr'].append(xarr)
                    percentageWiseRelErrXYArr[kval]['yarr'].append(yarr[kval])

                # plt.plot(xarr[0], np.std(yarr), 'ko', label="std " + str(prevFractionEdges)+"%")
                
                # medXarr.append(xarr[0])
                # medYarr.append(np.median(yarr))

                # fivePercentLineXarr.append(xarr[0])
                # fivePercentLineYarr.append(5.0)
                # fivePercentLineYarrNeg.append(-5.0)
                # exactCurve.append(0.0)
                # varXarr.append(xarr[0])
                # varYarr.append(np.std(yarr))
                # valuesStd.append(np.std(valuesYarr))
            # allEstimates[xarr[0]] = yarr
        percInd += 1
        xarr = []
        yarr = defaultdict(list)
        valuesYarr = []
        prevFractionEdges = fractionEdges
        

    else:
        flds = line.strip().split()
        xarr.append(fractionEdges)

        allEstVals = [float(x) for x in flds[2:-1]]

        for indVal in range(len(allEstVals)):
            # print(indVal) 
            # print(graphTypeArr[indVal], graphName)
            staticExact = exact[graphTypeArr[indVal]][graphName]
            err = 100 * (allEstVals[indVal] - staticExact)/staticExact
            yarr[allKVals[indVal]].append(err)
            # valuesYarr.append(estVal)

        if err > 100:
            print("Here ---- ", estVal)
    

# plt.plot(xarr, yarr, 'o', label=str(prevFractionEdges)+"%")
print(len(xarr), )
for kval in allKVals:
    percentageWiseRelErrXYArr[kval]['xarr'].append(xarr)
    percentageWiseRelErrXYArr[kval]['yarr'].append(yarr[kval])


f.close()

print("got the values...")


numRows = 1
fig = plt.figure(figsize=(20, 6))
fig.subplots_adjust(hspace=0.35, wspace=0.15)
plt.subplots_adjust(left=0.07, right=0.93)


# print(percentageWiseRelErrXYArr[4])

did = 1
for cliqueVal in percentageWiseRelErrXYArr:
    ax = fig.add_subplot(numRows, maxK-2, did)
    fivePercentLineXarr = []
    fivePercentLineYarr = []
    fivePercentLineYarrNeg = []
    exactCurve = []
    varXarr = []
    varYarr = []

    for eachXArrInd in range(len(percentageWiseRelErrXYArr[cliqueVal]['xarr'])):
        pxarr = percentageWiseRelErrXYArr[cliqueVal]['xarr'][eachXArrInd]
        pyarr = percentageWiseRelErrXYArr[cliqueVal]['yarr'][eachXArrInd]
        # print(pxarr, pyarr)
        ax.plot(pxarr, pyarr, 'o', label=str(pxarr[0])+"%")

        fivePercentLineXarr.append(pxarr[0])
        fivePercentLineYarr.append(5.0)
        fivePercentLineYarrNeg.append(-5.0)
        exactCurve.append(0.0)
        varXarr.append(pxarr[0])
        varYarr.append(np.std(pyarr))

    ax.set_ylim([-10,10])
    ax.plot(fivePercentLineXarr, fivePercentLineYarr, 'r--', label="5% Error")
    ax.plot(fivePercentLineXarr, fivePercentLineYarrNeg, 'r--', label="-5% Error")
    ax.plot(fivePercentLineXarr, exactCurve, 'r--*', label="Exact", markersize=10)
    ax.plot(varXarr, varYarr, 'kD-', label="Std. Dev.")

    # ax.legend(loc = "upper right", ncol=2, fontsize=14)
    ax.set_title("soc " + graphName + " " + cliqueType[cliqueVal], fontsize=20)
    ax.set_ylabel("% Error in Estimate", fontsize=20)
    ax.set_xlabel("% Edges Observed", fontsize=20)
    ax.tick_params(axis='x', which='major', labelsize=18)
    ax.tick_params(axis='y', which='major', labelsize=18)
    
    did += 1

# plt.plot(xarr[0], np.std(yarr), 'ko', label=str(prevFractionEdges)+"%")
# medXarr.append(xarr[0])
# medYarr.append(np.median(yarr))

# plt.plot(medXarr, medYarr, 'k^-', label="Median Rel Err", markersize = 12)

# fivePercentLineXarr.append(xarr[0])
# fivePercentLineYarr.append(5.0)
# fivePercentLineYarrNeg.append(-5.0)
# exactCurve.append(0.0)
# varXarr.append(xarr[0])
# varYarr.append(np.std(yarr))
# valuesStd.append(np.std(valuesYarr))

# plt.plot(fivePercentLineXarr, fivePercentLineYarr, 'r--', label="5% Error")
# plt.plot(fivePercentLineXarr, fivePercentLineYarrNeg, 'r--', label="-5% Error")
# plt.plot(fivePercentLineXarr, exactCurve, 'r--*', label="Exact", markersize=10)
# plt.plot(varXarr, varYarr, 'kD-', label="Std. Dev.")
# # plt.plot(varXarr, valuesStd, 'k*-', label="Std. Dev.(Estimate)")


# plt.title("soc " + graphName + "-- Error in Estimate for " + graphTypeToName[graphType] + "-- (Single RW)")
# plt.title("soc " + graphName + "-- Error in Estimate for " + graphTypeToName[graphType] + "-- (more edges)")
# plt.title("soc " + graphName + "-- Error in Estimate for " + graphTypeToName[graphType] + "-- (Random start node for each run)")
# plt.ylabel("% Error in Estimate")
# plt.xlabel("% Edges Observed")
# plt.legend(loc = "upper right", ncol=3, fontsize=14)

newHandles, newLabels = ax.get_legend_handles_labels()
fig.legend(newHandles, newLabels, loc= 'upper center', ncol=8, fontsize=16)
# fig.suptitle('This is a somewhat long figure title', fontsize=20)
plt.show()