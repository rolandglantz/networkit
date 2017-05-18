/*
 * Corres.cpp
 *
 *  Created on September 20, 2016
 *  Author: Roland Glantz (roland.glantz@kit.edu)
 */

#include <iterator>
#include <functional>
#include <numeric>
#include <iostream>
#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <fstream>
#include <sstream>

#include "Corres.h"

#define NR_CLUSTERS 32
#define NR_CLUSTERS_HALF 16
#define NR_ELEMENTS 1000000
#define NR_ELEMENTS_HALF 500000
#define NR_ELEMENTS_PER_CLUSTER 31250
#define WRITE_PARTITION_A 0
#define WRITE_PARTITION_B 0
#define WRITE_SEGMENTATION_A 0
#define WRITE_SEGMENTATION_B 0
#define WRITE_SHUFFLED_PARTITIONS 0
#define SEED1 1
#define SEED2 2
#define SEED3 3

#define CONFLICT 1
#define VALID 2

#define GREEDY 1
#define CHECK 100000

namespace NetworKit {

/**********************************************************************/
/*                         createToyPartitions                        */
/**********************************************************************/
void Corres::createToyPartitions(Partition& partitionA, Partition& partitionB) {
    partitionA = Partition(31);
    partitionA.setUpperBound(6);
    partitionA.addToSubset(0, 0);
    partitionA.addToSubset(0, 1);
    partitionA.addToSubset(0, 2);
    partitionA.addToSubset(0, 3);
    partitionA.addToSubset(1, 4);
    partitionA.addToSubset(1, 5);
    partitionA.addToSubset(1, 6);
    partitionA.addToSubset(1, 7);
    partitionA.addToSubset(2, 8);
    partitionA.addToSubset(2, 9);
    partitionA.addToSubset(2, 10);
    partitionA.addToSubset(2, 11);
    partitionA.addToSubset(2, 12);
    partitionA.addToSubset(2, 13);
    partitionA.addToSubset(3, 14);
    partitionA.addToSubset(3, 15);
    partitionA.addToSubset(3, 16);
    partitionA.addToSubset(3, 17);
    partitionA.addToSubset(3, 18);
    partitionA.addToSubset(3, 19);
    partitionA.addToSubset(3, 20);
    partitionA.addToSubset(4, 21);
    partitionA.addToSubset(4, 22);
    partitionA.addToSubset(4, 23);
    partitionA.addToSubset(4, 24);
    partitionA.addToSubset(4, 25);
    partitionA.addToSubset(5, 26);
    partitionA.addToSubset(5, 27);
    partitionA.addToSubset(5, 28);
    partitionA.addToSubset(5, 29);
    partitionA.addToSubset(5, 30);
    
    partitionB = Partition(31);
    partitionB.setUpperBound(5);
    partitionB.addToSubset(0, 0);
    partitionB.addToSubset(0, 1);
    partitionB.addToSubset(0, 2);
    partitionB.addToSubset(0, 3);
    partitionB.addToSubset(0, 4);
    partitionB.addToSubset(0, 5);
    partitionB.addToSubset(0, 6);
    partitionB.addToSubset(0, 8);
    partitionB.addToSubset(1, 7);
    partitionB.addToSubset(1, 9);
    partitionB.addToSubset(1, 10);
    partitionB.addToSubset(1, 11);
    partitionB.addToSubset(1, 14);
    partitionB.addToSubset(1, 15);
    partitionB.addToSubset(1, 16);
    partitionB.addToSubset(2, 12);
    partitionB.addToSubset(2, 17);
    partitionB.addToSubset(2, 18);
    partitionB.addToSubset(2, 19);
    partitionB.addToSubset(3, 13);
    partitionB.addToSubset(3, 21);
    partitionB.addToSubset(3, 22);
    partitionB.addToSubset(3, 23);
    partitionB.addToSubset(3, 26);
    partitionB.addToSubset(3, 27);
    partitionB.addToSubset(4, 20);
    partitionB.addToSubset(4, 24);
    partitionB.addToSubset(4, 25);
    partitionB.addToSubset(4, 28);
    partitionB.addToSubset(4, 29);
    partitionB.addToSubset(4, 30);
}

/**********************************************************************/
/*                      createPairZeroPartitions                      */
/**********************************************************************/
void Corres::createPairZeroPartitions(Partition& partitionA, Partition& partitionB) {
    //get seed for random partitions partitionA, partitionB
    unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
    
    //initialize partitionA
    partitionA = Partition(NR_ELEMENTS);
    partitionA.setUpperBound(NR_CLUSTERS);
    
    //set partitionA
    index indx = 0;
    for(index clusterID1 = 0; clusterID1 < NR_CLUSTERS; clusterID1++) {
        for(index j = 0; j < NR_ELEMENTS_PER_CLUSTER;j++) {
            partitionA.addToSubset(clusterID1, indx);
            indx++;
        }
    }
    
    //shuffle cluster IDs of partitionA
    std::array<index,NR_CLUSTERS> shuffledClusterIDs1;
    for(index id1 = 0;id1 < NR_CLUSTERS; id1++) {
        shuffledClusterIDs1[id1] = id1;
    }
    shuffle(shuffledClusterIDs1.begin(), shuffledClusterIDs1.end(), std::default_random_engine(seed));
    //\mathcal{B} given by shuffledClusterIDs1[0], ...,  shuffledClusterIDs1[NR_CLUSTERS/2-1]
    
    //get shuffled list of elements in U_B
    std::array<index,NR_ELEMENTS_HALF> shuffledElementsFromB;
    indx = 0;
    for(index i = 0; i < NR_CLUSTERS_HALF; i++) {
        for(index j = 0; j < NR_ELEMENTS_PER_CLUSTER; j++) {
            shuffledElementsFromB[indx] = j + NR_ELEMENTS_PER_CLUSTER * shuffledClusterIDs1[i];
            indx++;
        }
    }
    shuffle(shuffledElementsFromB.begin(), shuffledElementsFromB.end(), std::default_random_engine(seed));
    
    
    //get shuffled list of elements in U_B'
    std::array<index,NR_ELEMENTS_HALF> shuffledElementsFromB_comp;
    indx = 0;
    for(index i =  NR_CLUSTERS_HALF; i < NR_CLUSTERS; i++) {
        for(index j = 0; j < NR_ELEMENTS_PER_CLUSTER; j++) {
            shuffledElementsFromB_comp[indx] = j + (NR_ELEMENTS_PER_CLUSTER * shuffledClusterIDs1[i]);
            indx++;
        }
    }
    shuffle(shuffledElementsFromB_comp.begin(), shuffledElementsFromB_comp.end(), std::default_random_engine(seed));
    
    //initialize partitionB
    partitionB = Partition(NR_ELEMENTS);
    partitionB.setUpperBound(NR_CLUSTERS);
    
    //shuffle cluster IDs of partitionB
    std::array<index,NR_CLUSTERS> shuffledClusterIDs2;
    for(index id2 = 0;id2 < NR_CLUSTERS; id2++) {
        shuffledClusterIDs2[id2] = id2;
    }
    shuffle(shuffledClusterIDs2.begin(), shuffledClusterIDs2.end(), std::default_random_engine(seed));
    //\mathcal{B}' given by shuffledClusterIDs2[NR_CLUSTERS/2], ...,  shuffledClusterIDs2[NR_CLUSTERS-1]
    
    //set partitionB
    indx = 0;
    for(index i = 0; i <  NR_CLUSTERS_HALF; i++) {
        for(index j = 0; j <  NR_ELEMENTS_PER_CLUSTER; j++) {
            //            std::cout << "shuffledClusterIDs2[i] ist " << shuffledClusterIDs2[i] << " und shuffledElementsFromB[indx] ist " << shuffledElementsFromB[indx] << std::endl;
            partitionB.addToSubset(shuffledClusterIDs2[i], shuffledElementsFromB[indx]);
            indx++;
        }
    }
    indx = 0;
    for(index i = NR_CLUSTERS_HALF; i < NR_CLUSTERS; i++) {
        for(index j = 0; j < NR_ELEMENTS_PER_CLUSTER; j++) {
            //            std::cout << "shuffledClusterIDs2[i] ist " << shuffledClusterIDs2[i] << " und shuffledElementsFromB_comp[indx] ist " << shuffledElementsFromB_comp[indx] << std::endl;
            partitionB.addToSubset(shuffledClusterIDs2[i], shuffledElementsFromB_comp[indx]);
            indx++;
        }
    }
}

/**********************************************************************/
/*                        kickOutIdenticalTwins                       */
/**********************************************************************/
void Corres::kickOutIdenticalTwins(const Partition& tmpNormalPartitionA, const Partition& tmpNormalPartitionB,
                                            Partition& normalPartitionA, Partition& normalPartitionB) {
    
    
    count cardTmpPartitionA = tmpNormalPartitionA.upperBound();
    count cardTmpPartitionB = tmpNormalPartitionB.upperBound();
    std::map<index, count> cardinalityOfTmpClusterA = tmpNormalPartitionA.subsetSizeMap();
    std::map<index, count> cardinalityOfTmpClusterB = tmpNormalPartitionB.subsetSizeMap();
    
    //represent clusters by sets, partitions by vectors of sets
    std::vector<std::set<index>> clustersOfTmpPartitionA(cardTmpPartitionA);
    std::vector<std::set<index>> clustersOfTmpPartitionB(cardTmpPartitionB);
    for(index i1 = 0; i1 < cardTmpPartitionA; i1++) {
        clustersOfTmpPartitionA[i1] = tmpNormalPartitionA.getMembers(i1);
    }
    for(index i2 = 0; i2 < cardTmpPartitionB; i2++) {
        clustersOfTmpPartitionB[i2] = tmpNormalPartitionB.getMembers(i2);
    }
    
    //find identical twins
    std::vector<bool> identicalTwinA;
    identicalTwinA.resize(tmpNormalPartitionA.upperBound());
    std::fill(identicalTwinA.begin(), identicalTwinA.end(), false);
    std::vector<bool> identicalTwinB;
    identicalTwinB.resize(tmpNormalPartitionB.upperBound());
    std::fill(identicalTwinB.begin(), identicalTwinB.end(), false);
    count nrIdenticalTwins = 0;
    for(index i1 = 0; i1 < cardTmpPartitionA; i1++) {
        for(index i2 = 0; i2 < cardTmpPartitionB; i2++) {
            std::vector<index> interset;
            std::set_intersection((clustersOfTmpPartitionA[i1]).begin(), (clustersOfTmpPartitionA[i1]).end(),
                                  (clustersOfTmpPartitionB[i2]).begin(), (clustersOfTmpPartitionB[i2]).end(),
                                  std::back_inserter(interset));
            if(((count) interset.size() == cardinalityOfTmpClusterA[i1]) && ((count) interset.size() == cardinalityOfTmpClusterB[i2])) {
                nrIdenticalTwins++;
                identicalTwinA[i1] = true;
                identicalTwinB[i2] = true;
            }
        }
    }
    TRACE("XXX:Number of identical twins is ", nrIdenticalTwins);
    
    //get relation between old and new cluster IDs
    std::vector<index> newClusterID_A(cardTmpPartitionA, cardTmpPartitionA);//no such cluster
    count nrElementsInTwins = 0;
    count newClusterID = 0;
    for(index i1 = 0; i1 < cardTmpPartitionA; i1++) {
        if(identicalTwinA[i1] == false) {
            newClusterID_A[i1] = newClusterID;
            newClusterID++;
        } else {
            nrElementsInTwins += cardinalityOfTmpClusterA[i1];
        }
        
    }
    std::vector<index> newClusterID_B(cardTmpPartitionB, cardTmpPartitionB);//no such cluster
    newClusterID = 0;
    for(index i2 = 0; i2 < cardTmpPartitionB; i2++) {
        if(identicalTwinB[i2] == false) {
            newClusterID_B[i2] = newClusterID;
            newClusterID++;
        }
    }
    
    //write normalPartitionA
    normalPartitionA = Partition(tmpNormalPartitionA.numberOfElements() - nrElementsInTwins);
    normalPartitionA.setUpperBound(tmpNormalPartitionA.numberOfSubsets() - nrIdenticalTwins);
    count newElementID = 0;
    for(count oldElementID = 0;  oldElementID < tmpNormalPartitionA.numberOfElements();  oldElementID++) {
        index oldClusterID = tmpNormalPartitionA.subsetOf(oldElementID);
        if(identicalTwinA[oldClusterID] == false) {
            normalPartitionA.addToSubset(newClusterID_A[oldClusterID], newElementID);
            newElementID++;
        }
    }
    
    //write normalPartitionB
    normalPartitionB = Partition(tmpNormalPartitionB.numberOfElements() - nrElementsInTwins);
    normalPartitionB.setUpperBound(tmpNormalPartitionB.numberOfSubsets() - nrIdenticalTwins);
    newElementID = 0;
    for(count oldElementID = 0;  oldElementID < tmpNormalPartitionB.numberOfElements();  oldElementID++) {
        index oldClusterID = tmpNormalPartitionB.subsetOf(oldElementID);
        if(identicalTwinB[oldClusterID] == false) {
            normalPartitionB.addToSubset(newClusterID_B[oldClusterID], newElementID);
            newElementID++;
        }
    }
}

/**********************************************************************/
/*                         normalizeElements                          */
/**********************************************************************/
void Corres::normalizeElements(const Partition& partitionA, const Partition& partitionB,
                                        Partition& normalPartitionA, Partition& normalPartitionB,
                                        std::vector<index>& old2newElement) {
    
    //TRACE("partitionA is ", partitionA.getVector());
    //TRACE("partitionB is ", partitionB.getVector());
    
    //consecutively renumber the elements of partitionA such that
    //the new number of any element in cluster i is smaller than
    //the new number of any element in cluster j whenever i < j
    std::map<index, count> cardinalityOfClusterA = partitionA.subsetSizeMap();
    std::vector<count> currSlotInCluster(partitionA.numberOfSubsets(), 0);
    for(count cluster = 1; cluster < partitionA.numberOfSubsets(); cluster++) {
        currSlotInCluster[cluster] = currSlotInCluster[cluster-1] + cardinalityOfClusterA[cluster-1];
    }
    for(count oldElement = 0; oldElement < partitionA.numberOfElements(); oldElement++) {
        old2newElement[oldElement] = currSlotInCluster[partitionA.subsetOf(oldElement)];
        //std::cout << "oldElement is " << oldElement << " and old2newElement[oldElement] is " << old2newElement[oldElement] << std::endl;
        (currSlotInCluster[partitionA.subsetOf(oldElement)])++;
    }
    
    //normalization of partitionA w.r.t newElementNumber
    Partition tmpNormalPartitionA = Partition(partitionA.numberOfElements());
    tmpNormalPartitionA.setUpperBound(partitionA.numberOfSubsets());
    for(count oldElement = 0; oldElement < partitionA.numberOfElements(); oldElement++) {
        tmpNormalPartitionA.addToSubset(partitionA.subsetOf(oldElement), old2newElement[oldElement]);
    }
    
    //normalization of partitionB w.r.t newElementNumber
    Partition  tmpNormalPartitionB = Partition(partitionB.numberOfElements());
    tmpNormalPartitionB.setUpperBound(partitionB.numberOfSubsets());
    for(count oldElement = 0; oldElement < partitionB.numberOfElements(); oldElement++) {
        tmpNormalPartitionB.addToSubset(partitionB.subsetOf(oldElement), old2newElement[oldElement]);
    }
    
    kickOutIdenticalTwins(tmpNormalPartitionA, tmpNormalPartitionB, normalPartitionA, normalPartitionB);
}

/**********************************************************************/
/*                           getDistributions                         */
/**********************************************************************/
void Corres::getDistributions(Partition& partition1, Partition& partition2) {
    
    cardPartition1 = partition1.upperBound();
    cardPartition2 = partition2.upperBound();
    cardPartition2Dec = cardPartition2 - 1;
    cardinalityOfCluster1 = partition1.subsetSizeMap();
    cardinalityOfCluster2 = partition2.subsetSizeMap();
    for(index i1 = 0; i1 < cardPartition1; i1++) {
        cardinalityOfCluster1HalfFloor[i1] = (cardinalityOfCluster1[i1])/2;
        cardinalityOfCluster1HalfCeil[i1] = cardinalityOfCluster1HalfFloor[i1];
        if(cardinalityOfCluster1[i1] % 2 == 1) {
            (cardinalityOfCluster1HalfCeil[i1])++;
        }
    }
    for(index i2 = 0; i2 < cardPartition2; i2++) {
        cardinalityOfCluster2HalfFloor[i2] = (cardinalityOfCluster2[i2])/2;
        cardinalityOfCluster2HalfCeil[i2] = cardinalityOfCluster2HalfFloor[i2];
        if(cardinalityOfCluster2[i2] % 2 == 1) {
            (cardinalityOfCluster2HalfCeil[i2])++;
        }
    }
    
    //represent clusters by sets, partitions by vectors of sets
    std::vector<std::set<index>> clustersOfPartition1(cardPartition1);
    std::vector<std::set<index>> clustersOfPartition2(cardPartition2);
    for(index i1 = 0; i1 < cardPartition1; i1++) {
        clustersOfPartition1[i1] = partition1.getMembers(i1);
    }
    for(index i2 = 0; i2 < cardPartition2; i2++) {
        clustersOfPartition2[i2] = partition2.getMembers(i2);
    }
    
    //Get distributions from vectors of sets
    distributions.resize(cardPartition1, std::vector<count> (cardPartition2));
    distributionsPrime.resize(cardPartition2, std::vector<count> (cardPartition1));
    
    for(index i1 = 0; i1 < cardPartition1; i1++) {
        for(index i2 = 0; i2 < cardPartition2; i2++) {
            std::vector<index> interset;
            std::set_intersection((clustersOfPartition1[i1]).begin(), (clustersOfPartition1[i1]).end(),
                                  (clustersOfPartition2[i2]).begin(), (clustersOfPartition2[i2]).end(),
                                  std::back_inserter(interset));
            distributions[i1][i2] = (count) interset.size();
            //std::cout << "distributions[i1][i2] is " << distributions[i1][i2] << std::endl;
            distributionsPrime[i2][i1] = distributions[i1][i2];
        }
    }
    //TRACE("cardinalityOfCluster1 is ", cardinalityOfCluster1);
    //TRACE("distributions is ", distributions);
}

/**********************************************************************/
/*                         makeBipartiteGraph                         */
/**********************************************************************/
void Corres::makeBipartiteGraph(Partition& partition1, Partition& partition2) {
    
    //for each cluster i1 in partition1 find the neighboring (intersecting) clusters in partition 2
    bipartForwardNeigh.resize(cardPartition1);
    for(index i1 = 0; i1 < cardPartition1; i1++) {
        count numberNeighs = 0;
        for(index i2 = 0; i2 < cardPartition2; i2++) {
            if(distributions[i1][i2] > 0) {
                numberNeighs++;
            }
        }
        //std::cout << "1:numberNeighs is " << numberNeighs << std::endl;
        (bipartForwardNeigh[i1]).resize(numberNeighs+1);
        count counter = 0;
        for(index i2 = 0; i2 < cardPartition2; i2++) {
            if(distributions[i1][i2] > 0) {
                (bipartForwardNeigh[i1])[counter] = i2;
                counter++;
            }
        }
        (bipartForwardNeigh[i1])[counter] = cardPartition2; //no valid neighbor
    }
    
    //for each cluster i2 in partition2 find the neighboring (intersecting) clusters in partition 1
    bipartBackwardNeigh.resize(cardPartition2);
    for(index i2 = 0; i2 < cardPartition2; i2++) {
        count numberNeighs = 0;
        for(index i1 = 0; i1 < cardPartition1; i1++) {
            if(distributionsPrime[i2][i1] > 0) {
                numberNeighs++;
            }
        }
        //std::cout << "2:numberNeighs is " << numberNeighs << std::endl;
        (bipartBackwardNeigh[i2]).resize(numberNeighs+1);
        count counter = 0;
        for(index i1 = 0; i1 < cardPartition1; i1++) {
            if(distributionsPrime[i2][i1] > 0) {
                (bipartBackwardNeigh[i2])[counter] = i1;
                counter++;
            }
        }
        (bipartBackwardNeigh[i2])[counter] = cardPartition1; //no valid neighbor
    }
}

/**********************************************************************/
/*                                 peak                               */
/**********************************************************************/
count Corres::peak(count cardCluster2, count overlap) {
    return(std::min(overlap, cardCluster2 - overlap));
}

/**********************************************************************/
/*                              getBoundPeak                          */
/**********************************************************************/
count Corres::getBoundPeak(void) {
    count ret = 0;
    for(count i2 = 0; i2 < distriB_s.size(); i2++) {
        ret+= std::min(distriB_s[i2], distriB_t[i2]);
    }
    for(count i1 = 0; i1 < cardPartition1; i1++) {
        if(belongs[i1] == 0){
            ret+= std::min(incBound_s[i1], incBound_t[i1]);
        }
    }
    return(ret);
}

/**********************************************************************/
/*                            incrementS2                             */
/**********************************************************************/
count Corres::incrementS2(index newCluster) {
    count status = VALID;
    belongs[newCluster] = 1;
    index partner2 = (bipartForwardNeigh[newCluster])[0];
    count counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_s
        count oldDistrib_s = distriB_s[partner2];
        (distriB_s[partner2])+= distributions[newCluster][partner2];
        
        //update incBound_t
        if((distriB_s[partner2] >= cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_s < cardinalityOfCluster2HalfCeil[partner2])) {
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_t[partner1] += distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
	counter2++;
        partner2 = (bipartForwardNeigh[newCluster])[counter2];
    }
    
    //compute fragmentation of correspondence
    fracInter = 0;
    for(index i2 = 0; i2 < cardPartition2; i2++) {
        fracInter+= std::min(distriB_s[i2], cardinalityOfCluster2[i2] - distriB_s[i2]);
    }
    return(status);
}

/**********************************************************************/
/*                            incrementT2                             */
/**********************************************************************/
count Corres::incrementT2(index newCluster) {
    count status = VALID;
    belongs[newCluster] = 2;
    index partner2 = (bipartForwardNeigh[newCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_t
        count oldDistrib_t = distriB_t[partner2];
        (distriB_t[partner2])+= distributions[newCluster][partner2];
        
        //update incBound_s
        if((distriB_t[partner2] >= cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_t < cardinalityOfCluster2HalfCeil[partner2])) {
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_s[partner1] += distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[newCluster])[counter2];
    }
    //compute fragmentation of correspondence
    fracInter = 0;
    for(index i2 = 0; i2 < cardPartition2; i2++) {
        fracInter+= std::min(distriB_t[i2], cardinalityOfCluster2[i2] - distriB_t[i2]);
    }
    return(status);
}

/**********************************************************************/
/*                            decrementS2                             */
/**********************************************************************/
void Corres::decrementS2(index oldCluster) {
    belongs[oldCluster] = 0;
    index partner2 = (bipartForwardNeigh[oldCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_s
        count oldDistrib_s = distriB_s[partner2];
        (distriB_s[partner2])-= distributions[oldCluster][partner2];
        
        //update incBound_t
        if((distriB_s[partner2] < cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_s >= cardinalityOfCluster2HalfCeil[partner2])) {
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_t[partner1] -= distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[oldCluster])[counter2];
    }
}

/**********************************************************************/
/*                             decrementT2                            */
/**********************************************************************/
void Corres::decrementT2(index oldCluster) {
    belongs[oldCluster] = 0;
    index partner2 = (bipartForwardNeigh[oldCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        
        //update distriB_t
        count oldDistrib_t = distriB_t[partner2];
        (distriB_t[partner2])-= distributions[oldCluster][partner2];
        
        //update incBound_s
        if((distriB_t[partner2] < cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_t >= cardinalityOfCluster2HalfCeil[partner2])) {
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_s[partner1] -= distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[oldCluster])[counter2];
    }
}

/**********************************************************************/
/*                            incrementS3                             */
/**********************************************************************/
count Corres::incrementS3(index newCluster) {
    belongs[newCluster] = 1;
    index partner2 = (bipartForwardNeigh[newCluster])[0];
    count counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_s
        count oldDistrib_s = distriB_s[partner2];
         (distriB_s[partner2])+= distributions[newCluster][partner2];
        
        //update incBound_t and numberPartnersOfS
        if((distriB_s[partner2] >= cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_s < cardinalityOfCluster2HalfCeil[partner2])) {
            numberPartnersOfS++;
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_t[partner1] += distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[newCluster])[counter2];
    }

    count status = 0;
    if((numberPartnersOfS > 0) && (numberPartnersOfS < cardPartition2Dec)) {
        //we have a valid correspondence cutting off $\mathcal{S}$ or $\mathcal{T}$
        status = VALID;
        //compute fragmentation of correspondence
        fracInter = 0;
        for(index i2 = 0; i2 < cardPartition2; i2++) {
            fracInter+= std::min(distriB_s[i2], cardinalityOfCluster2[i2] - distriB_s[i2]);
        }
    }
    return(status);
}

/**********************************************************************/
/*                            incrementT3                             */
/**********************************************************************/
count Corres::incrementT3(index newCluster) {
    belongs[newCluster] = 2;
    index partner2 = (bipartForwardNeigh[newCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_t
        count oldDistrib_t = distriB_t[partner2];
        (distriB_t[partner2])+= distributions[newCluster][partner2];
        
        //update incBound_s and numberPartnersOfT
        if((distriB_t[partner2] >= cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_t < cardinalityOfCluster2HalfCeil[partner2])) {
            numberPartnersOfT++;
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_s[partner1] += distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[newCluster])[counter2];
    }
    int status = 0;
    if((numberPartnersOfT > 0) && (numberPartnersOfT < cardPartition2Dec)) {
        //we have a valid correspondence cutting off $\mathcal{S}$ or $\mathcal{T}$
        status = VALID;
        //compute fragmentation of correspondence
        fracInter = 0;
        for(index i2 = 0; i2 < cardPartition2; i2++) {
            fracInter+= std::min(distriB_t[i2], cardinalityOfCluster2[i2] - distriB_t[i2]);
        }
    }
    return(status);
}

/**********************************************************************/
/*                            decrementS3                             */
/**********************************************************************/
void Corres::decrementS3(index oldCluster) {
    belongs[oldCluster] = 0;
    index partner2 = (bipartForwardNeigh[oldCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_s
        count oldDistrib_s = distriB_s[partner2];
        (distriB_s[partner2])-= distributions[oldCluster][partner2];
        
        //update incBound_t and numberPartnersOfS
        if((distriB_s[partner2] < cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_s >= cardinalityOfCluster2HalfCeil[partner2])) {
            numberPartnersOfS--;
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_t[partner1] -= distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[oldCluster])[counter2];
    }
}

/**********************************************************************/
/*                             decrementT3                            */
/**********************************************************************/
void Corres::decrementT3(index oldCluster) {
    belongs[oldCluster] = 0;
    index partner2 = (bipartForwardNeigh[oldCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        
        //update distriB_t
        count oldDistrib_t = distriB_t[partner2];
        (distriB_t[partner2])-= distributions[oldCluster][partner2];
        
        //update incBound_s and numberPartnersOfT
        if((distriB_t[partner2] < cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_t >= cardinalityOfCluster2HalfCeil[partner2])) {
            numberPartnersOfT--;
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_s[partner1] -= distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[oldCluster])[counter2];
    }
}

/**********************************************************************/
/*                            incrementS4                             */
/**********************************************************************/
count Corres::incrementS4(index newCluster) {
    belongs[newCluster] = 1;
    numberBelongs2s++;
    count status= 0;
    index partner2 = (bipartForwardNeigh[newCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_s
        count oldDistrib_s = distriB_s[partner2];
        (distriB_s[partner2])+= distributions[newCluster][partner2];
        
        //update incBound_t, distriB_sPrime
        if((distriB_s[partner2] >= cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_s < cardinalityOfCluster2HalfCeil[partner2])) {
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_t[partner1] += distributionsPrime[partner2][partner1];
                distriB_sPrime[partner1]+= distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        
        //update belongsPrime and status
        if((distriB_s[partner2] > cardinalityOfCluster2HalfFloor[partner2]) && (belongsPrime[partner2] == 0)) {
            belongsPrime[partner2] = 1;
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                if((distriB_sPrime[partner1] > cardinalityOfCluster1HalfFloor[partner1]) && belongs[partner1] == 2) {
                    status= CONFLICT;
                }
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[newCluster])[counter2];
    }
    
    count numberBelongs2sCovered = 0;
    for(count i1 = 0; i1 < cardPartition1; i1++) {
        if((belongs[i1] == 1) && (distriB_sPrime[i1] >= cardinalityOfCluster1HalfCeil[i1])) {
            numberBelongs2sCovered++;
        }
    }
    status = CONFLICT;
    if(numberBelongs2s == numberBelongs2sCovered) {//we have a mutual correspondence cutting off $\mathcal{S}$
        status = VALID;
        //compute fragmentation of mutual correspondence
        fracInter = 0;
        for(index i2 = 0; i2 < cardPartition2; i2++) {
            fracInter+= std::min(distriB_s[i2], cardinalityOfCluster2[i2] - distriB_s[i2]);
        }
    }
    return(status);
}

/**********************************************************************/
/*                            incrementT4                             */
/**********************************************************************/
count Corres::incrementT4(index newCluster) {
    belongs[newCluster] = 2;
    //std::cout << "Adding " << newCluster << " to " << "the t-side" << std::endl;
    numberBelongs2t++;
    int status = 0;
    index partner2 = (bipartForwardNeigh[newCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_t
        count oldDistrib_t = distriB_t[partner2];
        (distriB_t[partner2])+= distributions[newCluster][partner2];
        
        //update incBound_s, distriB_tPrime and numberBelongs2tCovered
        if((distriB_t[partner2] >= cardinalityOfCluster2HalfFloor[partner2]) && (oldDistrib_t < cardinalityOfCluster2HalfFloor[partner2])) {
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_s[partner1] += distributionsPrime[partner2][partner1];
                //count oldDistriB_tPrime = distriB_tPrime[partner1];
                distriB_tPrime[partner1]+= distributionsPrime[partner2][partner1];
//                if((distriB_tPrime[partner1] >= cardinalityOfCluster1HalfFloor[partner1]) &&
//                   (oldDistriB_tPrime < cardinalityOfCluster1HalfFloor[partner1]) && (belongs[partner1] == 2)) {
//                    numberBelongs2tCovered++;
//                }
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        
        //update belongsPrime and status
        if((distriB_t[partner2] > cardinalityOfCluster2HalfFloor[partner2]) && (belongsPrime[partner2] == 0)) {
            belongsPrime[partner2] = 2;
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                if((distriB_tPrime[partner1] > cardinalityOfCluster1HalfFloor[partner1]) && belongs[partner1] == 1) {
                    status = CONFLICT;
                }
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        counter2++;
        partner2 = (bipartForwardNeigh[newCluster])[counter2];
    }

    count numberBelongs2tCovered = 0;
    for(count i1 = 0; i1 < cardPartition1; i1++) {
        if((belongs[i1] == 2) && (distriB_tPrime[i1] >= cardinalityOfCluster1HalfCeil[i1])) {
            numberBelongs2tCovered++;
        }
    }

    status = CONFLICT;
    if(numberBelongs2t == numberBelongs2tCovered) {
        status = VALID;
        //compute fragmentation of mutual correspondence
        fracInter = 0;
        for(index i2 = 0; i2 < cardPartition2; i2++) {
            fracInter+= std::min(distriB_t[i2], cardinalityOfCluster2[i2] - distriB_t[i2]);
        }
    }
    return(status);
}

/**********************************************************************/
/*                            decrementS4                             */
/**********************************************************************/
void Corres::decrementS4(index oldCluster) {
    //std::cout << "This is decrementS4" << std::endl;
    belongs[oldCluster] = 0;
    numberBelongs2s--;
    index partner2 = (bipartForwardNeigh[oldCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_s
        count oldDistrib_s = distriB_s[partner2];
        (distriB_s[partner2])-= distributions[oldCluster][partner2];
        
        //update incBound_t, distriB_sPrime
        if((distriB_s[partner2] < cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_s >= cardinalityOfCluster2HalfCeil[partner2])) {
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_t[partner1] -= distributionsPrime[partner2][partner1];
                distriB_sPrime[partner1]-= distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        
        //update belongsPrime
        if((distriB_s[partner2] <= cardinalityOfCluster2HalfFloor[partner2]) &&
           (belongsPrime[partner2] == 1)) {
            belongsPrime[partner2] = 0;
        }
        counter2++;
        partner2 = (bipartForwardNeigh[oldCluster])[counter2];
    }
}

/**********************************************************************/
/*                             decrementT4                            */
/**********************************************************************/
void Corres::decrementT4(index oldCluster) {
    //std::cout << "This is decrementT4" << std::endl;
    belongs[oldCluster] = 0;
    numberBelongs2t--;
    index partner2 = (bipartForwardNeigh[oldCluster])[0];
    index counter2 = 0;
    while(partner2 != cardPartition2) {
        //update distriB_t
        count oldDistrib_t = distriB_t[partner2];
        (distriB_t[partner2])-= distributions[oldCluster][partner2];
        
        //update incBound_s, distriB_tPrime
        if((distriB_t[partner2] < cardinalityOfCluster2HalfCeil[partner2]) && (oldDistrib_t >= cardinalityOfCluster2HalfCeil[partner2])) {
            index partner1 = (bipartBackwardNeigh[partner2])[0];
            index counter1 = 0;
            while(partner1 != cardPartition1) {
                incBound_s[partner1] -= distributionsPrime[partner2][partner1];
                distriB_tPrime[partner1]-= distributionsPrime[partner2][partner1];
                counter1++;
                partner1 = (bipartBackwardNeigh[partner2])[counter1];
            }
        }
        //update belongsPrime
        if((distriB_t[partner2] <= cardinalityOfCluster2HalfFloor[partner2]) &&
           (belongsPrime[partner2] == 2)) {
            belongsPrime[partner2] = 0;
        }
        counter2++;
        partner2 = (bipartForwardNeigh[oldCluster])[counter2];
    }
}

/**********************************************************************/
/*                               potFracs                             */
/**********************************************************************/
count Corres::potFracs(count& sFracPotNew, count& tFracPotNew) {
  //srand (time(NULL));  
  std::vector<count> sFracPot(cardPartition1);
  std::vector<count> tFracPot(cardPartition1);

  std::fill(sFracPot.begin(), sFracPot.end(), 0);
  std::fill(tFracPot.begin(), tFracPot.end(), 0);
    
    for(count i2 = 0; i2 < cardPartition2; i2++) {
        index counter = 0;
        index partner1 = (bipartBackwardNeigh[i2])[0];
        while(partner1 != cardPartition1) {
            if(belongs[partner1] == 0) {
                sFracPot[partner1]+= std::min(distriB_s[i2] + distributionsPrime[i2][partner1], distriB_t[i2]);
                tFracPot[partner1]+= std::min(distriB_s[i2], distriB_t[i2] + distributionsPrime[i2][partner1]);
            }
            counter++;
            partner1 = bipartBackwardNeigh[i2][counter];
        }
    }

    index bestCluster1=cardPartition1;//no such cluster
    count fittestness = 0;
    for(index i1 = 0; i1 < cardPartition1; i1++) {
        //for(index ii1 = 0; ii1 < cardPartition1; ii1++) {
        //count i1 = cardPartition1 - 1 - ii1;
        if(belongs[i1] == 0) {//i1 has not been inserted yet
            count fitness = (signed long int)(labs(((signed long int)(sFracPot[i1])) - ((signed long int)(tFracPot[i1]))));
            //fitness += std::max((signed long int)(sFracPot[i1]), (signed long int)(tFracPot[i1]));
            fitness += rand() % (1 + std::min((signed long int)(sFracPot[i1]), (signed long int)(tFracPot[i1])));
            if(fitness >= fittestness) {
                fittestness = fitness;
                bestCluster1 = i1;
                sFracPotNew = sFracPot[i1];
                tFracPotNew = tFracPot[i1];
            }
        }
    }
    
    
       // index bestCluster1=cardPartition1;//no such cluster
       // count maxDiff = 0;//absolute difference between sFracPot and tFracPot
       // for(index i1 = 0; i1 < cardPartition1; i1++) {
       //     if(belongs[i1] == 0) {//i1 has not been inserted yet
       //         count diff = labs(sFracPot[i1] - tFracPot[i1]);
       //         if(diff >= maxDiff) {
       //             maxDiff = diff;
       //             bestCluster1 = i1;
       //             sFracPotNew = sFracPot[i1];
       //             tFracPotNew = tFracPot[i1];
       //         }
       //     }
       // }
    //std::cout << "   sFracPotNew and tFracPotNew are " << sFracPotNew << " and " << tFracPotNew << std::endl;
    return(bestCluster1);
}

/**********************************************************************/
/*                 resetBestBelongsMutual(bestBelongs                 */
/**********************************************************************/
void Corres::resetBestBelongsMutual(std::vector<int>& bestBelongs, bool sWon) {

    if(sWon == true) {
        for(index i1 = 0; i1 < cardPartition1; i1++) {
            if(belongs[i1] == 1) {
                bestBelongs[i1] = 1;
            } else {
                bestBelongs[i1] = 2;
            }
        }
    } else {//ie, tWon == true
        for(index i1 = 0; i1 < cardPartition1; i1++) {
            if(belongs[i1] == 2) {
                bestBelongs[i1] = 2;
            } else {
                bestBelongs[i1] = 1;
            }
        }
    }
}

/**********************************************************************/
/*                           greedyDescent3                           */
/**********************************************************************/
index Corres::greedyDescent3(index s, index t, count& bestFrac, count& currentFrac,
                                     std::vector<index>& position2cluster, count& position,
                                     std::vector<int>& bestBelongs,
                                     count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                     std::vector<bool>& doneWith) {
    count status = 0;
    count newCluster = cardPartition1;//no such cluster
    while(currentFrac < bestFrac) {
        count sFracPotNew = 0;
        count tFracPotNew = 0;
        newCluster = potFracs(sFracPotNew, tFracPotNew);
        //std::cout << "s, t and newCluster are " << s << " " << t << " and " << newCluster << std::endl;
        if(newCluster < cardPartition1) {//newCluster is valid cluster ID for greedy descent
            //            if((sFracPotNew == 0) || (tFracPotNew == 0)) {
            //                doneWith[position] = true;
            //            }
            position2cluster[position] = newCluster;
            //cluster with ID newCluster goes to U_{\mathcal{B}_s} or to U_{\mathcal{B}_t}
            bool sWon = false;
            if(sFracPotNew < tFracPotNew) {
                status = (*this.*incrementS)(newCluster);
                sWon = true;
            } else {
                status = (*this.*incrementT)(newCluster);
            }
            
            currentFrac = getBoundPeak();
            if((status == VALID) && (fracInter < bestFrac)) {
                bestFrac = fracInter;
                resetBestBelongsMutual(bestBelongs, sWon);
            }
            position++;
            
        } else {//newCluster == cardPartition1, i.e., the s-t cut is complete
            
            //check if S has canonical partners
            bool lonely_s = true;
            index partner2 = 0;
            while((lonely_s == true) && (partner2 < cardPartition2)) {
                if(distriB_s[partner2] >= cardinalityOfCluster2HalfCeil[partner2]) {
                    lonely_s = false;
                }
                partner2++;
            }
            
            //if not, find best partner for s
            if(lonely_s == true) {
                //std::cout << "lonely s" << std::endl;
                count minDamage2s = numberOfElements;
                //index bestPartner_s = cardPartition2;
                for(index i2 = 0; i2 < cardPartition2; i2++) {
                    index potDamage = cardinalityOfCluster2[i2] - 2 * distriB_s[i2];
                    if(potDamage < minDamage2s) {
                        minDamage2s = potDamage;
                        //bestPartner_s = i2;
                    }
                }
                //std::cout << "1: currentFrac is " << currentFrac << std::endl;
                currentFrac += minDamage2s;
                //std::cout << "2: currentFrac is " << currentFrac << std::endl;
            } else {
                //check if T has canonical partners
                bool lonely_t = true;
                index partner2 = 0;
                while((lonely_t == true) && (partner2 < cardPartition2)) {
                    if(distriB_t[partner2] >= cardinalityOfCluster2HalfCeil[partner2]) {
                        lonely_t = false;
                    }
                    partner2++;
                }
                
                //is not find best partner for t
                if(lonely_t == true) {
                    //std::cout << "lonely t" << std::endl;
                    count minDamage2t = numberOfElements;
                    //index bestPartner_t = cardPartition2;
                    for(index i2 = 0; i2 < cardPartition2; i2++) {
                        index potDamage = cardinalityOfCluster2[i2] - 2 * distriB_t[i2];
                        if(potDamage < minDamage2t) {
                            minDamage2t = potDamage;
                            //bestPartner_t = i2;
                        }
                    }
                    //std::cout << "1: currentFrac is " << currentFrac << std::endl;
                    currentFrac += minDamage2t;
                    //std::cout << "2: currentFrac is " << currentFrac << std::endl;
                }
            }
            if(currentFrac <  bestFrac) {
                bestFrac = currentFrac;
                std::copy(belongs.begin(), belongs.end(), bestBelongs.begin());
                
                //                for(count i1 = 0; i1 < cardPartition1; i1++) {
                //                    std::cout << "   belongs[i1] is " << belongs[i1] << std::endl;
                //                }
                //                count minSum = 0;
                //                for(count i2 = 0; i2 < cardPartition2; i2++) {
                //                    std::cout << "   distriB_s[i2] is " << distriB_s[i2] << std::endl;
                //                    std::cout << "   distriB_t[i2] is " << distriB_t[i2] << std::endl;
                //                    minSum+= std::min(distriB_s[i2], distriB_t[i2]);
                //                }
                //                std::cout << "   minSum is " << minSum << std::endl;
            }
        }
    }
    //std::cout << "returning bestFrac = " << bestFrac << std::endl;
    return(position - 1);//position >= 1
}

/**********************************************************************/
/*                           greedyDescent4                           */
/**********************************************************************/
index Corres::greedyDescent4(index s, index t, count& bestFrac, count& currentFrac,
                                     std::vector<index>& position2cluster, count& position,
                                     std::vector<int>& bestBelongs,
                                     count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                     std::vector<bool>& doneWith) {
    //std::cout << "greedyDescent4: position is " << position + 1 << std::endl;
    count status = 0;
    count newCluster = cardPartition1;//no such cluster
    while((currentFrac < bestFrac) && (status != CONFLICT)) {//initially the case
        count sFracPotNew = 0;
        count tFracPotNew = 0;
        newCluster = potFracs(sFracPotNew, tFracPotNew);
        //std::cout << std::endl;
//        if((newCluster == cardPartition1)&&(position < cardPartition1 - 2)) std::cout << "Was soll denn das? position is " << position + 1 << std::endl;
        //std::cout << "s, t and newCluster are " << s << " " << t << " and " << newCluster << std::endl;
        if(newCluster < cardPartition1) {//newCluster is valid cluster ID for greedy descent, initially the case
            //            if((sFracPotNew == 0) || (tFracPotNew == 0)) {
            //                doneWith[position] = true;
            //            }
            position2cluster[position] = newCluster;
            //cluster with ID newCluster goes to U_{\mathcal{B}_s} or to U_{\mathcal{B}_t}
            bool sWon = false;
            if(sFracPotNew < tFracPotNew) {
                status = (*this.*incrementS)(newCluster);
                sWon = true;
            } else {
                belongs[newCluster] = 2;
                status = (*this.*incrementT)(newCluster);
            }
            //if(conflict > 0) std::cout << "CONFLICT at position " << position + 1 << std::endl;
            //std::cout << "Last cluster added at position " << position + 1 << " is " << newCluster + 1 << std::endl;
            //            for(index i1 = 0; i1 < cardPartition1; i1++) {
            //                if(belongs[i1] == 1) std::cout << i1  + 1 << " belongs to s" << std::endl;
            //                if(belongs[i1] == 2) std::cout << i1  + 1 << " belongs to t" << std::endl;
            //            }
            //            for(index i2 = 0; i2 < cardPartition1; i2++) {
            //                if(belongsPrime[i2] == 1) std::cout << i2 + 1 << "' belongs to s" << std::endl;
            //                if(belongsPrime[i2] == 2) std::cout << i2 + 1 << "' belongs to t" << std::endl;
            //            }
            

            if(status != CONFLICT) {
                currentFrac = getBoundPeak();
                if((status == VALID) && (fracInter < bestFrac)) {
                    //std::cout << "greedyDescent4::Early exit? bestFrac goes from " << bestFrac << " to " << fracInter << std::endl;
                    //std::cout << "greedyDescent4::currentFrac is " << currentFrac << std::endl;
                    bestFrac = fracInter;
                    resetBestBelongsMutual(bestBelongs, sWon);
                }
            }
            position++;
        } else {//newCluster == cardPartition1, i.e., the s-t cut is complete
            //std::cout << "GOAL: position, currentFrac and status are " << position << " and " << currentFrac << " and " << status << std::endl;
            if(currentFrac <  bestFrac) {
                if(status == VALID) {
                    bestFrac = currentFrac;
                    //std::cout << "New bestFrac is " << bestFrac << std::endl;
                    std::copy(belongs.begin(), belongs.end(), bestBelongs.begin());
                } else {
                    status = CONFLICT;
                }
            }
        }
    }
    //std::cout << std::endl;
    return(position - 1);//position >= 1
}

/**********************************************************************/
/*                             greedyBB                               */
/**********************************************************************/
count Corres::greedyBB(index s, index t, count bestFrac, count currentFrac,
                                std::vector<int>& bestBelongs,
                                count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                void (Corres::*decrementS)(index newCluster), void (Corres::*decrementT)(index newCluster),
                                index (Corres::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
                                                                       std::vector<index>& insertedAt, count& position,
                                                                       std::vector<int>& bestBelongs,
                                                                       count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                                                       std::vector<bool>& doneWith))
{
    //std::cout << "greedyBB: bestFrac and current Frac are " << bestFrac << " and " << currentFrac << std::endl;
    count position = 0;//initialization
    count maxPosition = cardPartition1 - 2;//s and t are never inserted
    //next vector is for keeping track of when clusters were inserted
    std::vector<index> position2cluster(maxPosition, cardPartition1);//cardPartition1 means no cluster
    std::vector<bool> doneWith(maxPosition, false);//needed to avoid infinite loops
    
    count numberForwardGreedyDescents = 0;
    bool done = false;
    if(currentFrac == bestFrac) {
        done = true;
    }
    while(done == false) {
//        if((numberForwardGreedyDescents >= 1) && (bestFrac > currentFrac)) {
//            bestFrac = currentFrac;
//        }
        count position1 = position;
        // if(numberForwardGreedyDescents == CHECK) {
	//   TRACE("Initial values of bestFrac and currentFrac are ", bestFrac, " and ", currentFrac);
	// }
        position = (*this.*greedyDescent)(s, t, bestFrac, currentFrac, position2cluster, position, bestBelongs,
                                          incrementS, incrementT, doneWith);
        
        if(position >= position1) {
            numberForwardGreedyDescents ++;
        }
//        if(position < position1) {
//            std::cout << "Was soll denn das? position < position1" << std::endl;
//        }
        // if((numberForwardGreedyDescents >= CHECK)&&(numberForwardGreedyDescents < CHECK + 10)) {
	//   TRACE("Went from ", position1 + 1, " to ", position + 1, ": bestFrac and currentFrac are ", bestFrac, " and ", currentFrac);
	// }
	// if(numberForwardGreedyDescents >= CHECK + 10) {
	//   TRACE("taking too long");
	//    exit(0);
	//  }        
//        if(currentFrac >= bestFrac) {
//            std::cout << "   backtracking because of high currentFrac" << std::endl;
//        } else {
//            std::cout << "   backtracking because of conflict" << std::endl;
//        }
//        std::cout << std::endl;
        
        if(bestFrac == 0) {
            return(0);
        }

        
        if((GREEDY == 1)&&(bestFrac < numberOfElements)) {
            done = true;
        } else {
            //backtracking
            bool blameHighCurrentFrac = false;
            if(currentFrac >= bestFrac) {//i.e., we have been kicked out of greedyDescent because currentFrac was too high
                blameHighCurrentFrac = true;
            }
            bool keepGoingBack = true;
            while(keepGoingBack == true) {//position > 0
                if(blameHighCurrentFrac == true) {//i.e., backtracking needs to continue because currentFrac is too high
                    doneWith[position] = true;
                }
                while((doneWith[position] == true) && (position > 0)) {//undo latest cluster assignments
                    //std::cout << "going back: position is " << position + 1 << std::endl;
                    index cluster = position2cluster[position];
                    if(belongs[cluster] == 1) {
                        (*this.*decrementS)(cluster);
                    }
                    if(belongs[cluster] == 2) {
                        (*this.*decrementT)(cluster);
                    }
                    //std::cout << "Resetting belongs of cluster " << cluster + 1 << " at position " << position + 1 << std::endl;
                    doneWith[position] = false;
                    position--;
                }
                
                if(doneWith[position] == false) {//new first step of extension
                    count status = 0;
                    //std::cout << "bin drin" << std::endl;
                    index cluster = position2cluster[position];
                    bool sWon = false;
                    if(belongs[cluster] == 1) {
                        sWon = true;
                        //std::cout << "Cluster " << cluster + 1 << " at position " << position + 1 << " switches from s to t" << std::endl;
                        (*this.*decrementS)(cluster);
                        status = (*this.*incrementT)(cluster);
                    } else {
                        if(belongs[cluster] == 2) {
                            //std::cout << "Cluster " << cluster + 1 << " at position " << position + 1 << " switches from t to s" << std::endl;
                            (*this.*decrementT)(cluster);
                            status = (*this.*incrementS)(cluster);
                        }
                    }
                    //std::cout << "Should I go back?";
                    doneWith[position] = true;
                    keepGoingBack = false;
                    blameHighCurrentFrac = false;
                    if(status == CONFLICT) {
                        //std::cout << " Going back because of conflict" << std::endl;
                        keepGoingBack = true;
                    } else {
                        currentFrac = getBoundPeak();
                        if((status == VALID) && (fracInter < bestFrac)) {
                            //std::cout << "greedyBB::Early exit? bestFrac goes from " << bestFrac << " to " << fracInter << std::endl;
                            //std::cout << "greedyBB::currentFrac is " << currentFrac << std::endl;
                            bestFrac = fracInter;
                            resetBestBelongsMutual(bestBelongs, sWon);
                        }
                        if(currentFrac >= bestFrac) {
                            //std::cout << " Going back because of high currentFrac" << std::endl;
                            blameHighCurrentFrac = true;
                            keepGoingBack = true;
                        }                }
                    if(position == 0) {
                        if(keepGoingBack == true) {//conflict or currentFrac too high
                            done = true;
                        }
                        keepGoingBack = false;
                        //std::cout << "Changed my mind because I am back to square 1" << std::endl;
                    }
                } else {//doneWith[position] == true AND position == 0 (see do-loop above)
                    //std::cout << "bin draussen" << std::endl;
                    keepGoingBack = false;
                    done = true;
                }
            }
            position++;
        }
        //std::cout << "Number of forwardGreedyDescents is " << numberForwardGreedyDescents << std::endl;
        
    }
    if(numberForwardGreedyDescents > maxNumberForwardGreedyDescents) {
        maxNumberForwardGreedyDescents = numberForwardGreedyDescents;
    }
    return(bestFrac);
}

/**********************************************************************/
/*                    getQualityTrivialSolutions2                     */
/**********************************************************************/
count Corres::getQualityTrivialSolutions2(index s, index t, bool& tWins) {
    count ret = getQualityTrivialSolutions4(s, t, tWins);
    if(cardinalityOfCluster1[s] < ret) {
        ret = cardinalityOfCluster1[s];
        tWins = false;
    }
    if(cardinalityOfCluster1[t] < ret) {
        ret = cardinalityOfCluster1[t];
        tWins = true;
    }
    return(ret);
}

/**********************************************************************/
/*                    getQualityTrivialSolutions3                     */
/**********************************************************************/
count Corres::getQualityTrivialSolutions3(index s, index t, bool& tWins) {
    count symDiff_s = numberOfElements;
    count neigh = (bipartForwardNeigh[s])[0];
    count counter = 0;
    count sCard = cardinalityOfCluster1[s];
    while(neigh != cardPartition2) {
        count symDiff = sCard + cardinalityOfCluster2[neigh] - (2 * distributions[s][neigh]);
        if(symDiff < symDiff_s) {
            symDiff_s = symDiff;
            //std::cout << "   s found a partner, we have a trivial solution!" << std::endl;
        }
        counter++;
        neigh = (bipartForwardNeigh[s])[counter];
    }
    tWins = false;
    count symDiff_t = numberOfElements;
    count tCard = cardinalityOfCluster1[t];
    neigh = (bipartForwardNeigh[t])[0];
    counter = 0;
    while(neigh != cardPartition2) {
        count symDiff = tCard + cardinalityOfCluster2[neigh] - (2 * distributions[t][neigh]);
        if(symDiff < symDiff_t) {
            symDiff_t = symDiff;
            //std::cout << "   t found a partner, we have a trivial solution!" << std::endl;
        }
        counter++;
        neigh = (bipartForwardNeigh[t])[counter];
    }
    if(symDiff_s < symDiff_t) {
        return(symDiff_s);
    } else {
        tWins = true;
        return(symDiff_t);
    }
}

/**********************************************************************/
/*                    getQualityTrivialSolutions4                     */
/**********************************************************************/
count Corres::getQualityTrivialSolutions4(index s, index t, bool& tWins) {
    
    count qualityTrivialSolutions = numberOfElements;
    count neigh = (bipartForwardNeigh[s])[0];
    count counter = 0;
    count sCard = cardinalityOfCluster1[s];
    while(neigh != cardPartition2) {
        count twiceIntersect = 2 * distributions[s][neigh];
        if((twiceIntersect >= cardinalityOfCluster1[s]) && (twiceIntersect >= cardinalityOfCluster2[neigh])) {
            count symDiff = sCard + cardinalityOfCluster2[neigh] - twiceIntersect;
            if(symDiff < qualityTrivialSolutions) {
                qualityTrivialSolutions = symDiff;
                //std::cout << "   s found a partner, we have a trivial solution!" << std::endl;
            }
        }
        counter++;
        neigh = (bipartForwardNeigh[s])[counter];
    }
    tWins = false;
    count tCard = cardinalityOfCluster1[t];
    neigh = (bipartForwardNeigh[t])[0];
    counter = 0;
    while(neigh != cardPartition2) {
        count twiceIntersect = 2 * distributions[t][neigh];
        if((twiceIntersect >= cardinalityOfCluster1[t]) && (twiceIntersect >= cardinalityOfCluster2[neigh])) {
            count symDiff = tCard + cardinalityOfCluster2[neigh] - twiceIntersect;
            if(symDiff < qualityTrivialSolutions) {
                qualityTrivialSolutions = symDiff;
                //std::cout << "   t found a partner, we have a trivial solution!" << std::endl;
                tWins = true;
            }
        }
        counter++;
        neigh = (bipartForwardNeigh[t])[counter];
    }
    return(qualityTrivialSolutions);
}

/**********************************************************************/
/*                              minCut4                               */
/**********************************************************************/
count Corres::minCut4(index s, index t, std::vector<int>& bestBelongs,
                               count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                               void (Corres::*decrementS)(index newCluster), void (Corres::*decrementT)(index newCluster),
                               count (Corres::*getQualityTrivialSolutions)(index s, index t, bool& tWins),
                               index (Corres::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
                                                                      std::vector<index>& insertedAt, count& position,
                                                                      std::vector<int>& bestBelongs,
                                                                      count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                                                      std::vector<bool>& doneWith))
{
    //reset belongs, belongsPrime, distriB_s, distriB_t, distriB_sPrime, distriB_tPrime
    std::fill(belongs.begin(), belongs.end(), 0);
    std::fill(belongsPrime.begin(), belongsPrime.end(), 0);
    std::fill(distriB_s.begin(), distriB_s.end(), 0);
    std::fill(distriB_t.begin(), distriB_t.end(), 0);
    std::fill(distriB_sPrime.begin(), distriB_sPrime.end(), 0);
    std::fill(distriB_tPrime.begin(), distriB_tPrime.end(), 0);
    
    std::fill(incBound_s.begin(), incBound_s.end(), 0);
    std::fill(incBound_t.begin(), incBound_t.end(), 0);
//    std::fill(  totalDamage_s.begin(),   totalDamage_s.end(), 0);
//    std::fill(  totalDamage_t.begin(),   totalDamage_t.end(), 0);

//    for (count i2 = 0; i2 < cardPartition2; i2++) {
//        std::fill(damage_s[i2].begin(), damage_s[i2].end(), 0);
//    }
//    for (count i2 = 0; i2 < cardPartition2; i2++) {
//        std::fill(damage_t[i2].begin(), damage_t[i2].end(), 0);
//    }

    numberBelongs2s = 0;
    numberBelongs2t = 0;
    numberPartnersOfS = 0;
    numberPartnersOfT = 0;
    
    (*this.*incrementS)(s);
    (*this.*incrementT)(t);
    if(cardPartition1 > 2) {
        bool tWins = false;
        count qualityTrivialSolutions = (*this.*getQualityTrivialSolutions)(s, t, tWins);
        count minNontrivCut = greedyBB(s, t, qualityTrivialSolutions, getBoundPeak(), bestBelongs,
                                       incrementS, incrementT, decrementS, decrementT, greedyDescent);
        if(minNontrivCut < numberOfElements) {
            //std::cout << "Weight of cut and qualityTrivialSolutions is " << minNontrivCut << " and " << qualityTrivialSolutions << std::endl;
            std::cout << "Weight of cut is " << minNontrivCut << std::endl;
        } else {
            std::cout << "No correspondence found" << std::endl;
        }
        
        if(minNontrivCut < qualityTrivialSolutions) {
            return(minNontrivCut);
        } else {
            if(tWins == false) {
                for(count i = 0; i < cardPartition1; i++) {
                    bestBelongs[i] = 2;
                }
                bestBelongs[s] = 1;
            } else {
                for(count i = 0; i < cardPartition1; i++) {
                    bestBelongs[i] = 1;
                }
                bestBelongs[t] = 2;
            }
            return(qualityTrivialSolutions);
        }
    } else {
        std::copy(belongs.begin(), belongs.end(), bestBelongs.begin());
        return(getBoundPeak());
    }
}

/**********************************************************************/
/*                            bucketSort                              */
/**********************************************************************/
void Corres::bucketSort(count cMax,
                                 std::vector<count>& gomoryHuParent, std::vector<count>& cutWithGomoryHuParent,
                                 std::vector<count>& sortedGomoryHuParent, std::vector<count>& sortedCutWithGomoryHuParent) {
    
    std::vector<std::vector<index> > bucket(cMax+1);
    for(count i = 0; i < cardPartition1; i++) {
        if(cutWithGomoryHuParent[i] <= cMax) {
            (bucket[cutWithGomoryHuParent[i]]).push_back(i);
        }
    }
    count counter = 0;
    for(count j = 0; j <= cMax ; j++) {
        for(index k = 0; k < (bucket[j]).size(); k++) {
            sortedGomoryHuParent[counter] = (bucket[j])[k];
            sortedCutWithGomoryHuParent[counter] = j;
            counter++;
        }
    }
}


/**********************************************************************/
/*                              gusfield                              */
/**********************************************************************/

count Corres::gusfield(std::vector<index>& gomoryHuParent, std::vector<count>& cutWithGomoryHuParent,
                                index& bestS, index& bestT, std::vector<int>& bestBestBelongs,
                                count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                void (Corres::*decrementS)(index newCluster), void (Corres::*decrementT)(index newCluster),
                                count (Corres::*getQualityTrivialSolutions)(index s, index t, bool& tWins),
                                index (Corres::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
                                                                       std::vector<index>& insertedAt, count& position,
                                                                       std::vector<int>& bestBelongs,
                                                                       count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                                                       std::vector<bool>& doneWith),
                                count (Corres::*minCut)(index s, index t, std::vector<int>& bestBelongs,
                                                                 count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                                                 void (Corres::*decrementS)(index newCluster), void (Corres::*decrementT)(index newCluster),
                                                                 count (Corres::*getQualityTrivialSolutions)(index s, index t, bool& tWins),
                                                                 index (Corres::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
                                                                                                        std::vector<index>& insertedAt, count& position,
                                                                                                        std::vector<int>& bestBelongs,
                                                                                                        count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
                                                                                                        std::vector<bool>& doneWith))) {
 
 //belongs[] expresses membership of a cluster $c$ from partition1 as follows.
 //belongs[c] = 1: $c \in $\mathcal{B}_s$
 //belongs[c] = 2: $c \in $\mathcal{B}_t$
 //belongs[c] = 0: $c \notin $\mathcal{B}_s$, $c \notin $\mathcal{B}_t$
 belongs.resize(cardPartition1);
 
 //belongsPrime[] expresses membership of a cluster $c'$ from $\mathcal{B}' (partition2) .
 //belongsPrime[c] = 1: $c \in $\mathcal{B}'_s$
 //belongsPrime[c] = 2: $c \in $\mathcal{B}'_t$
 //belongsPrime[c] = 0: $c \notin $\mathcal{B}'_s$, $c \notin $\mathcal{B}'_t$
 belongsPrime.resize(cardPartition2);
 
 distriB_s.resize(cardPartition2); //distribution of U_{\mathcal{B}_s} over the clusters of partition2
 distriB_t.resize(cardPartition2); //distribution of U_{\mathcal{B}_t} over the clusters of partition2
 distriB_sPrime.resize(cardPartition1); //distribution of U_{\mathcal{B}'_s} over the clusters of partition1
 distriB_tPrime.resize(cardPartition1); //distribution of U_{\mathcal{B}'_t} over the clusters of partition1
 
 //sFracPot.resize(cardPartition1, 0);
 //tFracPot.resize(cardPartition1, 0);
 
 incBound_s.resize(cardPartition1, 0);
 incBound_t.resize(cardPartition1, 0);
// totalDamage_s.resize(cardPartition1, 0);
// totalDamage_t.resize(cardPartition1, 0);
// damage_s.resize(cardPartition2, std::vector<count> (cardPartition1));
// damage_t.resize(cardPartition2, std::vector<count> (cardPartition1));

 maxNumberForwardGreedyDescents = 0;
 
 //build first Gomory-Hu tree (star with vertex 0 in the center
 for(index s = 1; s < cardPartition1; s++) {
     gomoryHuParent[s] = 0;
 }
 
 //rebuild the Gomory-Hu tree
 bestBestBelongs.resize(cardPartition1, 0);
 count cMin = numberOfElements;
 count cMax = 0;
 count totalCut = 0;
 for(index s = 1; s < cardPartition1; s++) {
     index t = gomoryHuParent[s];
     //TRACE(" s und t sind ", s, " und  ", t);
     std::vector<int> bestBelongs(cardPartition1, 0);//best belongs per s-t cut
     //std::cout << "s and t are " << s << " and " << t << std::endl;
     
     //std::cout << "cardinalityOfCluster1[3] is " << cardinalityOfCluster1[3] << std::endl;
     //std::cout << "cardinalityOfCluster1[5] is " << cardinalityOfCluster1[3] << std::endl;
     count c = minCut4(s, t, bestBelongs, incrementS, incrementT, decrementS, decrementT, getQualityTrivialSolutions, greedyDescent);
     //std::cout << std::endl;
     totalCut += c;
     //TRACE("   minCut zwischen ", s, " und  ", t, " ist ", c);
     
     //TRACE("   bestBelongs is ", bestBelongs);
     cutWithGomoryHuParent[s] = c;//label edge of Gomory-Hu tree
     if(c < cMin) {//update value of minimal Cut
         cMin = c;
         bestS = s;
         bestT = t;
         std::copy(bestBelongs.begin(), bestBelongs.end(), bestBestBelongs.begin());
     }
     if(c > cMax) {//update value of minimal Cut
         cMax = c;
     }
     //relink edges of Gomory-Hu tree
     //        for(index i = s + 1; i < cardPartition1; i++) {
     //            if(gomoryHuParent[i] == t) {//i has same parent as s, that is t
     //                if(bestBelongs[i] == 1) {//i is on the same side of the cut as s
     //                    //TRACE("       gomoryHuParent of ", i, " is set to ", s);
     //                    gomoryHuParent[i] = s;
     //                }
     //            }
     //        }
     for(index i = 0; i < cardPartition1; i++) {
         if((i != s)&&(gomoryHuParent[i] == t)) {//i has same parent as s, that is t
             if(bestBelongs[i] == 1) {//i is on the same side of the cut as s
                 //TRACE("       gomoryHuParent of ", i, " is set to ", s);
                 gomoryHuParent[i] = s;
             }
         }
     }
     if(bestBelongs[gomoryHuParent[t]] == 1) {
         gomoryHuParent[s] = gomoryHuParent[t];
         gomoryHuParent[t] = s;
         cutWithGomoryHuParent[s] = cutWithGomoryHuParent[t];
         cutWithGomoryHuParent[t] = c;
     }
     
 }
 
 //TRACE("gomoryHuParent is ", gomoryHuParent);
 //TRACE("cutWithGomoryHuParent is ", cutWithGomoryHuParent);
 //TRACE("cardPartition1, cardPartition2, cMin, cMax, totalCut and numberOfElements are ", cardPartition1, ", ", cardPartition2, ", ", cMin, ", ", cMax, ", ", totalCut, " and ", numberOfElements);
 
 //sort gomoryHuParent, cutWithGomoryHuParent
 //bucketSort(cMax, gomoryHuParent, cutWithGomoryHuParent, sortedGomoryHuParent, sortedCutWithGomoryHuParent);
  
std::cout << "Total cut has weight " << totalCut << std::endl;
return(totalCut);
}

//0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5
//0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5

//0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 4, 4, 4, 5, 5
//0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 4, 4, 4, 5, 5


/**********************************************************************/
/*                        getBestBelongsPrime                         */
/**********************************************************************/
void Corres::getBestBelongsPrime(std::vector<int>& bestBelongs,  std::vector<int>& bestBelongsPrime) {
    for(index i2 = 0; i2 < cardPartition2; i2++) {
        count thresh = (cardinalityOfCluster2[i2]) / 2;
        count fit = 0;
        for(index i1 = 0; i1 < cardPartition1; i1++) {
            if(bestBelongs[i1] == 1) {
                fit+=distributions[i1][i2];
            }
        }
        if(fit >= thresh) {
            bestBelongsPrime[i2] = 1;
        } else {
            bestBelongsPrime[i2] = 2;
        }
    }
}

/**********************************************************************/
/*                       evaluateCorrespondence                       */
/**********************************************************************/
void Corres::evaluateCorrespondence(std::vector<int>& bestBelongs, std::vector<int>& bestBelongsPrime,
                                             count& clusterCardinality, count& clusterCardinalityPrime,
                                             count& elementCardinality, count& elementCardinalityPrime,
                                             double symDiff, double& size, double& quality) {
    //cardinalities w.r.t. bestBelongs
    count clusterCardinalityA = 0;
    count elementCardinalityA = 0;
    for(index i1 = 0; i1 < cardPartition1; i1++) {
        if(bestBelongs[i1] == 1) {
            clusterCardinalityA++;
            elementCardinalityA+= cardinalityOfCluster1[i1];
        }
    }
    
    //cardinalities w.r.t. bestBelongsPrime
    count clusterCardinalityAPrime = 0;
    count elementCardinalityAPrime = 0;
    for(index i2 = 0; i2 < cardPartition2; i2++) {
        if(bestBelongsPrime[i2] == 1) {
            clusterCardinalityAPrime++;
            elementCardinalityAPrime+= cardinalityOfCluster2[i2];
        }
    }
    
    //go with $\mathcal{B}$, as opposed to $\mathcal{C} \setminus \mathcal{B}$, if $U_{\mathcal{B}}$ has fewer elements
    count elementCardinalityB = numberOfElements - elementCardinalityA;
    if(elementCardinalityA <= elementCardinalityB) {
        clusterCardinality = clusterCardinalityA;
        elementCardinality = elementCardinalityA;
        clusterCardinalityPrime = clusterCardinalityAPrime;
        elementCardinalityPrime = elementCardinalityAPrime;
    } else {
        clusterCardinality = cardPartition1 - clusterCardinalityA;
        elementCardinality = elementCardinalityB;
        clusterCardinalityPrime = cardPartition2 - clusterCardinalityAPrime;
        elementCardinalityPrime = numberOfElements - elementCardinalityAPrime;
    }
    double sizeIntersection = (elementCardinality + elementCardinalityPrime - symDiff) / 2.0;
    size = elementCardinality + elementCardinalityPrime - sizeIntersection;
    quality = sizeIntersection / size;

}

/**********************************************************************/
/*                     evaluateAllCorrespondences                     */
/**********************************************************************/
double Corres::evaluateAllCorrespondences(std::vector<index>& gomoryHuParent,
					  std::vector<count>& cutWithGomoryHuParent) {
    
  //2D vector to store the number of clusters in partition1 vs number of clusters in partition2
  // for each correspondence
  //entry at (i1, i2) stands for number of such correspondences
  std::vector<std::vector<count> > pairs(cardPartition1, std::vector<index>(cardPartition2,0));
    
  //2D vector to store total dissimilarities corresponding to 2D vector pairs
  std::vector<std::vector<count> > dissim(cardPartition1, std::vector<index>(cardPartition2));
    
  count counter = 0;
  count totalDissimilarity = 0;
  //std::cout << std::endl;
  for(count i1 = 0; i1 < cardPartition1; i1++) {
    if(gomoryHuParent[i1] < cardPartition1) {//{i1, gomoryHuParent[i1]} is an edge of the Gomory-Hu tree
      counter++;
      //find bestBelongs, i.e., clusters on the i1-side of Gomory-Hu tree
      std::vector<int> bestBelongs(cardPartition1, 0);
            
      //make a copy of gomoryHuParent
      //std::vector<count> copyGomoryHuParent(gomoryHuParent);
            
      std::vector<count> copyGomoryHuParent(cardPartition1, 0);
      for(count i1 = 0; i1 < cardPartition1; i1++) {
	copyGomoryHuParent[i1] = gomoryHuParent[i1];
      }
            
      //... first, turn copyGomoryHuParent into transitive closure of gomoryHuParent
      copyGomoryHuParent[i1] = i1;
      for(count j1 = 0; j1 < cardPartition1; j1++) {
	count k1 = j1;
	while((k1 != i1) && (k1 != cardPartition1)) {
	  k1 = copyGomoryHuParent[k1];
	}
	index forefather = k1;
	k1 = j1;
	while(k1 != forefather) {
	  copyGomoryHuParent[k1] = forefather;
	  k1 = copyGomoryHuParent[k1];
	}
      }
            
      //... then derive bestBelongs from transitive closure of copyGomoryHuParent
      for(count j1 = 0; j1 < cardPartition1; j1++) {
	if(copyGomoryHuParent[j1] == i1) {
	  bestBelongs[j1] = 1;
	}
      }
            
      //find bestBelongsPrime
      std::vector<int> bestBelongsPrime(cardPartition2, 0);
      getBestBelongsPrime(bestBelongs, bestBelongsPrime);
            
      //evaluate correspondence (bestBelongs, bestBelongsPrime)
      count clusterCardinality = 0;
      count clusterCardinalityPrime = 0;
      count elementCardinality = 0;
      count elementCardinalityPrime = 0;
      count symDiff = cutWithGomoryHuParent[i1];
      totalDissimilarity += symDiff;
      double size = 0.0;
      double quality = 0.0;
      evaluateCorrespondence(bestBelongs, bestBelongsPrime,
			     clusterCardinality, clusterCardinalityPrime,
			     elementCardinality, elementCardinalityPrime,
			     symDiff, size, quality);
      (pairs[clusterCardinality][clusterCardinalityPrime])++;
      (dissim[clusterCardinality][clusterCardinalityPrime])+= cutWithGomoryHuParent[i1];
    }
  }
    
  count diagDissim = 0;
  count offDiagDissim = 0;
  for(count i1 = 0; i1 < cardPartition1; i1++) {
    for(count i2 = 0; i2 < cardPartition2; i2++) {
      if(i1 == i2) {
	diagDissim+= dissim[i1][i2];
      } else {
	offDiagDissim+= dissim[i1][i2];
      }
    }
  }

  for(count i1 = 0; i1 < cardPartition1; i1++) {
    for(count i2 = 0; i2 < cardPartition2; i2++) {
      if(pairs[i1][i2] > 0) {
	TRACE("XXX:Number of ", i1, " : ", i2, " pairs is ", pairs[i1][i2]);
      }
    }
  }
  //TRACE("XXX:cardPartition1 is ", cardPartition1);

  //double imbalance = ((double) offDiagDissim) / ((double) diagDissim);
  double imbalance = ((double) offDiagDissim) / ((double) totalDissimilarity);
    
  TRACE("Normalized total dissimilarity and imbalance are ",
	((double) totalDissimilarity) / ((double)numberOfElements), " and ", imbalance);
  return(imbalance);
}

/**********************************************************************/
/*                   readSegmentationNormalizeClusters                */
/**********************************************************************/
bool Corres::readSegmentationNormalizeClusters(std::string filename, Partition& partition,
                                                        std::map<index, index>& clusterID2segmentID) {
    
    // open file for reading
    std::ifstream infile(filename.c_str());
    if (!infile) {
        std::cerr << "Error opening file" << filename << std::endl;
        return(false);
    }
    
    //get numberOfLines in file = number of elements partitioned
    count nrElements = 0;
    std::string line;
    while(std::getline(infile, line)) {nrElements++;}
    
    
    //return to first line
    infile.clear();
    infile.seekg(0);
    
    //read file, put content in vector named segmentation
    std::vector<index> segmentation(nrElements);
    index elementID = 0;
    while(std::getline(infile, line)) {
        std::stringstream ss(line);
        ss >> segmentation[elementID];
        elementID++;
    }
    //    for(elementID = 0; elementID < nrElements; elementID++) {
    //        std::cout << "elementID und segmentation[elementID] sind " << elementID << " und " << segmentation[elementID] << std::endl;
    //    }
    
    //map segmentIDs, i.e., the entries of the vector named segmentation,
    //onto (consectutive) clusterIDs and vice versa (cluster normalization)
    count nrClusters = 0;
    std::map<index, index> segmentID2clusterID;
    for(elementID = 0; elementID < nrElements; elementID++) {
        index segID = segmentation[elementID];
        //std::cout << "segID und nrClusters sind " << segID << " und " << nrClusters << std::endl;
        if(segmentID2clusterID.find(segID) == segmentID2clusterID.end()) {//i.e., new clusterID
            segmentID2clusterID[segID] = nrClusters;
            clusterID2segmentID[nrClusters] = segID;
            //std::cout << "IF::segID und nrClusters sind " << segID << " und " << nrClusters << std::endl;
            nrClusters++;
        }
    }
    infile.close();
    
    //specfy partition
    partition = Partition(nrElements);
    partition.setUpperBound(nrClusters);
    for(elementID = 0; elementID < nrElements; elementID++) {
      partition.addToSubset(segmentID2clusterID[segmentation[elementID]],elementID);
    }

    // count nrClustersDec = nrClusters - 1;
    // for(elementID = 0; elementID < nrElements; elementID++) {
    //     partition.addToSubset(nrClustersDec - segmentID2clusterID[segmentation[elementID]],elementID);
    // }

    
    return(true);
}

  /**********************************************************************/
  /*                                run                                 */
  /**********************************************************************/
  double Corres::run(const Partition& partitionA, const Partition& partitionB) {
    int level = 3;
    //level == 1: cuts through bipartite graph, not yet implemented
    //level == 2: nontrivial (one-sided) correspondences
    //level == 3: non-degenerate correspondences
    //level == 4: mutual correspondences
    if(level == 1) {
      //std::cout << "Level == 1 not yet implemented. Sorry!" << std::endl;
      TRACE("Level == 1 not yet implemented. Sorry!");
      return(0);
    }
    //depending on level, choose the right functions for incrementing, decrementing S or T
    count (Corres::*incrementS)(index newCluster) = &Corres::incrementS4;
    count (Corres::*incrementT)(index newCluster) = &Corres::incrementT4;
    void (Corres::*decrementS)(index newCluster) = &Corres::decrementS4;
    void (Corres::*decrementT)(index newCluster) = &Corres::decrementT4;
    count (Corres::*getQualityTrivialSolutions)(index s, index t, bool& tWins) = &Corres::getQualityTrivialSolutions4;
    index (Corres::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
				   std::vector<index>& insertedAt, count& position,
				   std::vector<int>& bestBelongs,
				   count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
				   std::vector<bool>& doneWith) = &Corres::greedyDescent4;
        
    count (Corres::*minCut)(index s, index t, std::vector<int>& bestBelongs,
			    count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
			    void (Corres::*decrementS)(index newCluster), void (Corres::*decrementT)(index newCluster),
			    count (Corres::*getQualityTrivialSolutions)(index s, index t, bool& tWins),
			    index (Corres::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
							   std::vector<index>& insertedAt, count& position,
							   std::vector<int>& bestBelongs,
							   count (Corres::*incrementS)(index newCluster), count (Corres::*incrementT)(index newCluster),
							   std::vector<bool>& doneWith)) = &Corres::minCut4;
    if(level == 2) {
      incrementS = &Corres::incrementS2;
      incrementT = &Corres::incrementT2;
      decrementS = &Corres::decrementS2;
      decrementT = &Corres::decrementT2;
      getQualityTrivialSolutions = &Corres::getQualityTrivialSolutions2;
      greedyDescent = &Corres::greedyDescent4;
      minCut = &Corres::minCut4;
    }
        
    if(level == 3) {
      incrementS = &Corres::incrementS3;
      incrementT = &Corres::incrementT3;
      decrementS = &Corres::decrementS3;
      decrementT = &Corres::decrementT3;
      getQualityTrivialSolutions = &Corres::getQualityTrivialSolutions3;
      greedyDescent = &Corres::greedyDescent3;
      minCut = &Corres::minCut4;
    }
        
    if(level == 4) {
      incrementS = &Corres::incrementS4;
      incrementT = &Corres::incrementT4;
      decrementS = &Corres::decrementS4;
      decrementT = &Corres::decrementT4;
      getQualityTrivialSolutions = &Corres::getQualityTrivialSolutions4;
      greedyDescent = &Corres::greedyDescent4;
      minCut = &Corres::minCut4;
    }
        
    index bestS, bestT;
    std::vector<int> bestBestBelongs;
    Partition partition1, partition2;//element normalizations of partitionA and partitionB
    std::vector<index> old2newElement(partitionA.numberOfElements());//needed to undo normalization of elelments
    normalizeElements(partitionA, partitionB, partition1, partition2, old2newElement);
    numberOfElements = partition1.numberOfElements();
        
    //std::cout << "Number of elements: " <<  numberOfElements << std::endl;
    //std::cout << "Number of clusters in partition 1: " << partition1.upperBound() << std::endl;
    //std::cout << "Number of clusters in partition 2: " << partition2.upperBound() << std::endl;
    //TRACE("Number of elements is ", numberOfElements);
    //TRACE("Number of clusters in partition 1: ", partition1.upperBound());
    //TRACE("Number of clusters in partition 2: ", partition2.upperBound());
        
        
    getDistributions(partition1, partition2);
    //TRACE("partition1 is ", partition1.getVector());
    //TRACE("partition2 is ", partition2.getVector());
        
    //make bipartite graph
    makeBipartiteGraph(partition1, partition2);
        
    //Gomory-Hu tree
    std::vector<index> gomoryHuParent(cardPartition1, cardPartition1);
    std::vector<count> cutWithGomoryHuParent(cardPartition1, numberOfElements);
    //std::vector<count> sortedGomoryHuParent(cardPartition1, cardPartition1);
    //std::vector<count> sortedCutWithGomoryHuParent(cardPartition1, numberOfElements);
    //std::cout << "Gusfield's algorithm ..." << std::endl;
    //TRACE("Gusfield's algorithm ...");
 
    time_t seconds0 = time(NULL);       
    count totalCut = gusfield(gomoryHuParent, cutWithGomoryHuParent,
				//sortedGomoryHuParent, sortedCutWithGomoryHuParent,
				bestS, bestT, bestBestBelongs,
				incrementS, incrementT, decrementS, decrementT, getQualityTrivialSolutions,
				greedyDescent, minCut);
    time_t seconds1 = time(NULL);
    std::cout << "Gusfield took " << seconds1 - seconds0 << " seconds." << std::endl;

        
    //std::cout << "Minimum cut has weight " << minimumCut << std::endl;
    //(std::cout << "Maximum number of greedyDescents is " <<  maxNumberForwardGreedyDescents << std::endl;
    //TRACE("Minimum cut has weight ", minimumCut);
    //TRACE("Maximum number of greedyDescents ", maxNumberForwardGreedyDescents);
        
    evaluateAllCorrespondences(gomoryHuParent, cutWithGomoryHuParent);
        
    //TRACE("bestS is ", bestS);
    //TRACE("bestT is ", bestT);
    // std::vector<int> bestBestBelongsPrime(cardPartition2, 0);
    // getBestBelongsPrime(bestBestBelongs, bestBestBelongsPrime);
    // TRACE("bestBestBelongs is ", bestBestBelongs);
    // TRACE("bestBestBelongsPrime is ", bestBestBelongsPrime);
        
    return(((double)totalCut)/((double)numberOfElements));
  }

} /* namespace NetworKit */

