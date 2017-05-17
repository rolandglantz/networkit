# extension imports
from _NetworKit import Coverage, Modularity, CommunityDetector, PLP, LPDegreeOrdered, PLM
from networkit import *
from networkit.community import Partition
from networkit.community import detectCommunities
from networkit.correspondences import Corres
from networkit.correspondences import Correspondences
from random import shuffle

import time
import sys
import subprocess
import os

setLogLevel('TRACE')

# external imports
import math
import numpy as np
try:
	import tabulate
except ImportError:
	print(""" WARNING: module 'tabulate' not found, please install it to use the full functionality of NetworKit """)

setNumberOfThreads(1);	
#setNumberOfThreads(8);	
#setNumberOfThreads(16);	

def _createPairToyPartitions():
        partition1 = Partition(30)
        partition1.setUpperBound(6)
        partition1.addToSubset(0, 0)
        partition1.addToSubset(0, 1)
        partition1.addToSubset(0, 2)
        partition1.addToSubset(0, 3)
        partition1.addToSubset(1, 4)
        partition1.addToSubset(1, 5)
        partition1.addToSubset(1, 6)
        partition1.addToSubset(1, 7)
        partition1.addToSubset(2, 8)
        partition1.addToSubset(2, 9)
        partition1.addToSubset(2, 10)
        partition1.addToSubset(2, 11)
        partition1.addToSubset(2, 12)
        partition1.addToSubset(2, 13)
        partition1.addToSubset(3, 14)
        partition1.addToSubset(3, 15)
        partition1.addToSubset(3, 16)
        partition1.addToSubset(3, 17)
        partition1.addToSubset(3, 18)
        partition1.addToSubset(3, 19)
        partition1.addToSubset(4, 20)
        partition1.addToSubset(4, 21)
        partition1.addToSubset(4, 22)
        partition1.addToSubset(4, 23)
        partition1.addToSubset(4, 24)
        partition1.addToSubset(5, 25)
        partition1.addToSubset(5, 26)
        partition1.addToSubset(5, 27)
        partition1.addToSubset(5, 28)
        partition1.addToSubset(5, 29)

        partition2 = Partition(30)
        partition2.setUpperBound(6)
        partition2.addToSubset(0, 0)
        partition2.addToSubset(0, 1)
        partition2.addToSubset(0, 2)
        partition2.addToSubset(0, 3)
        partition2.addToSubset(0, 4)
        partition2.addToSubset(0, 5)
        partition2.addToSubset(0, 6)
        partition2.addToSubset(1, 7)
        partition2.addToSubset(1, 8)
        partition2.addToSubset(1, 9)
        partition2.addToSubset(1, 10)
        partition2.addToSubset(1, 11)
        partition2.addToSubset(1, 12)
        partition2.addToSubset(2, 13)
        partition2.addToSubset(2, 14)
        partition2.addToSubset(2, 15)
        partition2.addToSubset(2, 16)
        partition2.addToSubset(3, 17)
        partition2.addToSubset(3, 18)
        partition2.addToSubset(3, 19)
        partition2.addToSubset(4, 20)
        partition2.addToSubset(4, 21)
        partition2.addToSubset(4, 25)
        partition2.addToSubset(4, 26)
        partition2.addToSubset(4, 27)
        partition2.addToSubset(5, 22)
        partition2.addToSubset(5, 23)
        partition2.addToSubset(5, 24)
        partition2.addToSubset(5, 28)
        partition2.addToSubset(5, 29)

        partition11 = Partition(8)
        partition11.setUpperBound(8)
        partition11.addToSubset(0, 0)
        partition11.addToSubset(1, 1)
        partition11.addToSubset(2, 2)
        partition11.addToSubset(3, 3)
        partition11.addToSubset(4, 4)
        partition11.addToSubset(5, 5)
        partition11.addToSubset(6, 6)
        partition11.addToSubset(7, 7)

        partition22 = Partition(8)
        partition22.setUpperBound(2)
        partition22.addToSubset(0, 0)
        partition22.addToSubset(0, 1)
        partition22.addToSubset(0, 2)
        partition22.addToSubset(0, 3)
        partition22.addToSubset(1, 4)
        partition22.addToSubset(1, 5)
        partition22.addToSubset(1, 6)
        partition22.addToSubset(1, 7)

        #partition11 = Partition(10)
        #partition11.setUpperBound(10)
        #partition11.addToSubset(0, 0)
        #partition11.addToSubset(1, 1)
        #partition11.addToSubset(2, 2)
        #partition11.addToSubset(3, 3)
        #partition11.addToSubset(4, 4)
        #partition11.addToSubset(5, 5)
        #partition11.addToSubset(6, 6)
        #partition11.addToSubset(7, 7)
        #partition11.addToSubset(8, 8)
        #partition11.addToSubset(9, 9)

        #partition22 = Partition(10)
        #partition22.setUpperBound(5)
        #partition22.addToSubset(0, 0)
        #partition22.addToSubset(0, 1)
        #partition22.addToSubset(1, 2)
        #partition22.addToSubset(1, 3)
        #partition22.addToSubset(2, 4)
        #partition22.addToSubset(2, 5)
        #partition22.addToSubset(3, 6)
        #partition22.addToSubset(3, 7)
        #partition22.addToSubset(4, 8)
        #partition22.addToSubset(4, 9)

        partition222 = Partition(30)
        partition222.setUpperBound(4)

        partition222.addToSubset(0, 0)
        partition222.addToSubset(0, 1)
        partition222.addToSubset(0, 2)
        partition222.addToSubset(0, 3)
        partition222.addToSubset(0, 4)
        partition222.addToSubset(0, 5)
        partition222.addToSubset(0, 6)
        partition222.addToSubset(0, 8)

        partition222.addToSubset(1, 7)
        partition222.addToSubset(1, 9)
        partition222.addToSubset(1, 10)
        partition222.addToSubset(1, 11)
        partition222.addToSubset(1, 14)
        partition222.addToSubset(1, 15)
        partition222.addToSubset(1, 16)

        partition222.addToSubset(2, 12)
        partition222.addToSubset(2, 17)
        partition222.addToSubset(2, 18)
        partition222.addToSubset(2, 19)

        partition222.addToSubset(3, 13)
        partition222.addToSubset(3, 20)
        partition222.addToSubset(3, 21)
        partition222.addToSubset(3, 22)
        partition222.addToSubset(3, 23)
        partition222.addToSubset(3, 24)
        partition222.addToSubset(3, 25)
        partition222.addToSubset(3, 26)
        partition222.addToSubset(3, 27)
        partition222.addToSubset(3, 28)
        partition222.addToSubset(3, 29)

        partitionA = Partition(31)
        partitionA.setUpperBound(6)
        partitionA.addToSubset(0, 0)
        partitionA.addToSubset(0, 1)
        partitionA.addToSubset(0, 2)
        partitionA.addToSubset(0, 3)
        partitionA.addToSubset(1, 4)
        partitionA.addToSubset(1, 5)
        partitionA.addToSubset(1, 6)
        partitionA.addToSubset(1, 7)
        partitionA.addToSubset(2, 8)
        partitionA.addToSubset(2, 9)
        partitionA.addToSubset(2, 10)
        partitionA.addToSubset(2, 11)
        partitionA.addToSubset(2, 12)
        partitionA.addToSubset(2, 13)
        partitionA.addToSubset(3, 14)
        partitionA.addToSubset(3, 15)
        partitionA.addToSubset(3, 16)
        partitionA.addToSubset(3, 17)
        partitionA.addToSubset(3, 18)
        partitionA.addToSubset(3, 19)
        partitionA.addToSubset(3, 20)
        partitionA.addToSubset(4, 21)
        partitionA.addToSubset(4, 22)
        partitionA.addToSubset(4, 23)
        partitionA.addToSubset(4, 24)
        partitionA.addToSubset(4, 25)
        partitionA.addToSubset(5, 26)
        partitionA.addToSubset(5, 27)
        partitionA.addToSubset(5, 28)
        partitionA.addToSubset(5, 29)
        partitionA.addToSubset(5, 30)

        partitionB = Partition(31)
        partitionB.setUpperBound(5)
        partitionB.addToSubset(0, 0)
        partitionB.addToSubset(0, 1)
        partitionB.addToSubset(0, 2)
        partitionB.addToSubset(0, 3)
        partitionB.addToSubset(0, 4)
        partitionB.addToSubset(0, 5)
        partitionB.addToSubset(0, 6)
        partitionB.addToSubset(0, 8)
        partitionB.addToSubset(1, 7)
        partitionB.addToSubset(1, 9)
        partitionB.addToSubset(1, 10)
        partitionB.addToSubset(1, 11)
        partitionB.addToSubset(1, 14)
        partitionB.addToSubset(1, 15)
        partitionB.addToSubset(1, 16)
        partitionB.addToSubset(2, 12)
        partitionB.addToSubset(2, 17)
        partitionB.addToSubset(2, 18)
        partitionB.addToSubset(2, 19)
        partitionB.addToSubset(3, 13)
        partitionB.addToSubset(3, 21)
        partitionB.addToSubset(3, 22)
        partitionB.addToSubset(3, 23)
        partitionB.addToSubset(3, 26)
        partitionB.addToSubset(3, 27)
        partitionB.addToSubset(4, 20)
        partitionB.addToSubset(4, 24)
        partitionB.addToSubset(4, 25)
        partitionB.addToSubset(4, 28)
        partitionB.addToSubset(4, 29)
        partitionB.addToSubset(4, 30)
        return (partitionB, partitionA)

def _createPairZeroPartitions():
        nrClusters = 200 # must be even
        nrClustersHalf = round(nrClusters/2)
        nrElements = 1000000 # must be multiple of nrClusters
        nrElementsPerCluster = round(nrElements / nrClusters)

        ########## initialize partitionZ1
        partitionZ1 = Partition(nrElements)
        partitionZ1.setUpperBound(nrClusters)

        ########## set partitionZ1
        index = 0
        for clusterID1 in range(nrClusters):
                for j in range(nrElementsPerCluster):
                        partitionZ1.addToSubset(clusterID1, index)
                        index = index + 1

        ########## shuffle cluster IDs of partitionZ1
        shuffledClusterIDs1 = list(range(nrClusters))
        shuffle(shuffledClusterIDs1)#\mathcal{B} given by shuffledClusterIDs1[0], ...,  shuffledClusterIDs1[nrClusters/2-1]

        ########## get shuffled list of elements in U_B
        shuffledElementsFromB = [0] * ((nrClustersHalf) * nrElementsPerCluster) 
        index = 0
        for i in range(nrClustersHalf):
                for j in range(nrElementsPerCluster):
                        shuffledElementsFromB[index] = j + nrElementsPerCluster * shuffledClusterIDs1[i] 
                        index = index + 1
        shuffle(shuffledElementsFromB)

        ########## get shuffled list of elements in U_B'
        shuffledElementsFromB_comp = [0] * ((nrClustersHalf) * nrElementsPerCluster) 
        index = 0
        for i in range(nrClustersHalf,nrClusters):
                for j in range(nrElementsPerCluster):
                        shuffledElementsFromB_comp[index] = j + (nrElementsPerCluster * shuffledClusterIDs1[i]) 
                        index = index + 1
        shuffle(shuffledElementsFromB_comp)

        ########## initialize partitionZ2
        partitionZ2 = Partition(nrElements)
        partitionZ2.setUpperBound(nrClusters)

        ########## shuffle cluster IDs of partitionZ2
        shuffledClusterIDs2 = list(range(nrClusters))
        shuffle(shuffledClusterIDs2)#\mathcal{B}' given by shuffledClusterIDs2[0], ...,  shuffledClusterIDs2[nrClusters/2-1]

        ########## set partitionZ2
        index = 0
        for i in range(nrClustersHalf):
                for j in range(nrElementsPerCluster):
                        partitionZ2.addToSubset(shuffledClusterIDs2[i], shuffledElementsFromB[index])
                        index = index + 1
        index = 0
        for i in range(nrClustersHalf, nrClusters):
                for j in range(nrElementsPerCluster):
                        partitionZ2.addToSubset(shuffledClusterIDs2[i], shuffledElementsFromB_comp[index])
                        index = index + 1

        return (partitionZ1, partitionZ2)

        
(partition1, partition2) = _createPairToyPartitions()
#(partition1, partition2) = _createPairZeroPartitions()

#check whether partition1, resp. partition11 and partition2, resp. partition22, are compact
#if not, exit
nrElements1Before = partition1.numberOfElements() 
nrElements2Before = partition2.numberOfElements() 
nrClusters1Before = partition1.upperBound()
nrClusters2Before = partition2.upperBound()
partition1.compact()
partition2.compact()
nrElements1After = partition1.numberOfElements() 
nrElements2After = partition2.numberOfElements() 
nrClusters1After = partition1.upperBound()
nrClusters2After = partition2.upperBound()
partitionsAreCompact = 1
if nrElements1Before > nrElements1After:
        print("Elements of ground set of first partition are not compact: possible reduction from ", nrElements1Before, "to ", nrElements1After)
        partitionsAreCompact = 0
if nrClusters1Before > nrClusters1After:
        print("Clusters of first partition are not compact: possible reduction from ", nrClusters1Before, "to ", nrClusters1After)
        partitionsAreCompact = 0
if nrElements2Before > nrElements2After:
        print("Elements of ground set of second partition are not compact: possible reduction from ", nrElements2Before, "to ", nrElements2After)
        partitionsAreCompact = 0
if nrClusters2Before > nrClusters2After:
        print("Clusters of second partition are not compact: possible reduction from ", nrClusters2Before, "to ", nrClusters2After)
        partitionsAreCompact = 0
if partitionsAreCompact == 0:
        print("Both partitions must be compact. EXITING!")
        raise SystemExit

graphDir = os.path.expanduser('~/c+/KarMa/data/graphs/extra_social/')

Gset = ['as-22july06',#index 0
        #'as-skitter',#index 1
        #'citationCiteseer',#index 2
        #'coAuthorsCiteseer',#index 3
        #'coAuthorsDBLP',#index 4
        #'coPapersCiteseer',#index 5
        #'coPapersDBLP',#index 6
        #'email-EuAll',#index 7
        #'loc-brightkite_edges',#index 8
        #'loc-gowalla_edges',#index 9
        #'p2p-Gnutella04',#index 10
        #'PGPgiantcompo',#index 11
        #'soc-Slashdot0902',#index 12
        #'web-Google',#index 13
        #'wiki-Talk'#index 14
]

GsetNIX = ['as-22july06',#index 7
           'as-skitter',#index 8
           'email-EuAll',#index 9
           'loc-brightkite_edges',#index 10
           'loc-gowalla_edges',#index 11
           'p2p-Gnutella04',#index 12
           'soc-Slashdot0902',#index 13
           'wiki-Talk'#index 14
]

nrSeeds = 1

#for graphName in Gset:
#        print ("graphName is ", graphName)
#        graphFile = graphDir + graphName + '.graph'
#        G = readGraph(graphFile, Format.METIS)
#        print ("graphFile is ", graphFile)
#        for gammaValue in range (1081, 1092):
#                anzPLM = 0
#                anzPLP = 0
#                for seed in range(1,nrSeeds + 1):
#                        #print("seed is ", seed)
#                        #algo3 = PLM(G, refine=False)
#                        #partition3 = detectCommunities(G, algo3, True)
#                        #algo4 = PLM(G, refine=False,gamma=gammaValue)
#                        #partition4 = detectCommunities(G, algo4, True)
#                        algo3 = PLM(G, refine=False, gamma=gammaValue)
#                        partition3 = detectCommunities(G, algo3, True)
#                        algo4 = PLP(G)
#                        partition4 = detectCommunities(G, algo4, True)
#                        partition3.compact()
#                        partition4.compact()
#                        anzPLM = anzPLM + partition3.upperBound() 
#                        anzPLP = anzPLP + partition4.upperBound() 
#                        print ("Anzahl Cluster in PLM_", gammaValue, " is ", anzPLM)
#                        print ("Anzahl Cluster in PLP is ", anzPLP)
#                        
#                        #minCut = Corres().run(partition3, partition4)
#                        
#                        #print("minCut is", minCut)
#                #meanAnzPLM = anzPLM / nrSeeds
#                #meanAnzPLP = anzPLP / nrSeeds
#                #print("gammaValue, meanAnzPLM und meanAnzPLP sind ", gammaValue, ", ", meanAnzPLM, " und ", meanAnzPLP)

gammaValues = [246,221,183,69,103,76,1089]
#gammaValues = [5,15,290,100,100,100,100,100,100,100,100,79,100,100,5]

graphID = 0;
for graphName in Gset:
        print ("graphName ist ", graphName)
        graphFile = graphDir + graphName + '.graph'
        print ("graphFile ist ", graphFile)
        start = time.time()
        G = readGraph(graphFile, Format.METIS)
        finish = time.time()
        print ("Einlesen von ", graphFile, " dauerte ", finish - start, " Sekunden")
        anzPLM3 = 0
        anzPLM4 = 0
        for seed in range(1,nrSeeds + 1):
                #print("seed ist ", seed)
                start = time.time()
                algo3 = PLM(G, refine=False)
                partition3 = detectCommunities(G, algo3, True)
                partition3.compact()
                finish = time.time()
                print ("Ausrechnen von partition3 dauerte ", finish - start, " Sekunden")
                start = time.time()
                #algo4 = PLM(G, refine=True, 100)
                algo4 = PLM(G, refine=True)
                partition4 = detectCommunities(G, algo4, True)
                partition4.compact()
                finish = time.time()
                print ("Ausrechnen von partition4 dauerte ", finish - start, " Sekunden")
                anzPLM3 = anzPLM3 + partition3.upperBound() 
                anzPLM4 = anzPLM4 + partition4.upperBound() 
                print ("Anzahl Cluster in partition3 ist ", anzPLM3, " und Anzahl Cluster in partition4 ist ", anzPLM4)                
                start = time.time()
                minCut = Corres().run(partition3, partition4)
                finish = time.time()
                print ("Ausrechnen der Korrespondenzen dauerte ", finish - start, " Sekunden")
                
                #print("minCut is", minCut)
        graphID = graphID + 1;
#        meanAnzPLM3 = anzPLM3 / nrSeeds
#        meanAnzPLP4 = anzPLP4 / nrSeeds
