/*
 * File:   Correspondences.h
 * Author: Roland Glantz (roland.glantz@kit.edu)
 *
 */

#ifndef CORRESPONDENCES_H
#define	CORRESPONDENCES_H

#include "../structures/Partition.h"

namespace NetworKit {

class Correspondences {

public:
    count numberOfElements;

    count cardPartition1, cardPartition2, cardPartition2Dec;

    std::vector<count> incBound_s;
    std::vector<count> incBound_t;
//    std::vector<count> totalDamage_s;
//    std::vector<count> totalDamage_t;
//    std::vector<std::vector<count>> damage_s;
//    std::vector<std::vector<count>> damage_t;

    count numberBelongs2s = 0;
    count numberBelongs2t = 0;
    count numberPartnersOfS = 0;
    count numberPartnersOfT = 0;

    count fracInter = 0;

    std::map<index, count> cardinalityOfCluster1;
    std::map<index, count> cardinalityOfCluster1HalfFloor;
    std::map<index, count> cardinalityOfCluster1HalfCeil;

    std::map<index, count> cardinalityOfCluster2;
    std::map<index, count> cardinalityOfCluster2HalfFloor;
    std::map<index, count> cardinalityOfCluster2HalfCeil;

    std::vector<std::vector<count>> distributions;
    std::vector<std::vector<count>> distributionsPrime;

    std::vector<std::vector<index>> bipartForwardNeigh;
    std::vector<std::vector<index>> bipartBackwardNeigh;

    count maxNumberForwardGreedyDescents = 0;


    //belongs[] expresses membership of a cluster $c$ from $\mathcal{B} (partition1) as follows.
    //belongs[c] = 1: $c \in $\mathcal{B}_s$
    //belongs[c] = 2: $c \in $\mathcal{B}_t$
    //belongs[c] = 0: $c \notin $\mathcal{B}_s$, $c \notin $\mathcal{B}_t$
    std::vector<int> belongs;
    //belongsPrime[] expresses membership of a cluster $c'$ from $\mathcal{B}' (partition2) .
    //belongsPrime[c] = 1: $c \in $\mathcal{B}'_s$
    //belongsPrime[c] = 2: $c \in $\mathcal{B}'_t$
    //belongsPrime[c] = 0: $c \notin $\mathcal{B}'_s$, $c \notin $\mathcal{B}'_t$
    std::vector<int> belongsPrime;

    std::vector<count> distriB_s;//distribution of U_{\mathcal{B}_s} over the clusters of partition2
    std::vector<count> distriB_t;//distribution of U_{\mathcal{B}_t} over the clusters of partition2
    std::vector<count> distriB_sPrime;//distribution of U_{\mathcal{B}'_s} over the clusters of partition1
    std::vector<count> distriB_tPrime;//distribution of U_{\mathcal{B}'_t} over the clusters of partition1

    std::vector<index> gomoryHuParent;
    std::vector<count> cutWithGomoryHuParent;

    void createToyPartitions(Partition& PartitionA, Partition& PartitionB);

    void createPairZeroPartitions(Partition& partitionA, Partition& partitionB);

    void normalizeElements(const Partition& partitionA, const Partition& partitionB,
                           Partition& normalPartitionA, Partition& normalPartitionB,
                           std::vector<index>& old2newElement);

    void kickOutIdenticalTwins(const Partition& partitionA, const Partition& partitionB,
                               Partition& normalPartitionA, Partition& normalPartitionB);

    void getDistributions(Partition& partition1, Partition& partition2);

    void makeBipartiteGraph(Partition& partition1, Partition& partition2);

    count peak(count cardCluster2, count overlap);

    count getBoundPeak(void);

    count potFracs(count& sFracPotNew, count& tFracPotNew);

    void resetBestBelongsMutual(std::vector<int>& bestBelongs, bool sWon);

    index greedyDescent3(index s, index t, count& bestFrac, count& currentFrac,
                        std::vector<index>& position2cluster, count& position,
                        std::vector<int>& bestBelongs,
                        count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                        std::vector<bool>& doneWith);

    index greedyDescent4(index s, index t, count& bestFrac, count& currentFrac,
                        std::vector<index>& insertedAt, count& position,
                        std::vector<int>& bestBelongs,
                        count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                        std::vector<bool>& doneWith);

    count greedyBB(index s, index t, count bestFrac, count currentFrac, std::vector<int>& bestBelongs,
                   count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                   void (Correspondences::*decrementS)(index newCluster), void (Correspondences::*decrementT)(index newCluster),
                   index (Correspondences::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
                                                          std::vector<index>& insertedAt, count& position,
                                                          std::vector<int>& bestBelongs,
                                                          count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                                                          std::vector<bool>& doneWith));

    void getBestBelongsPrime(std::vector<int>& bestBelongs,  std::vector<int>& bestBelongsPrime);

    void evaluateCorrespondence(std::vector<int>& bestBelongs, std::vector<int>& bestBelongsPrime,
                                count& clusterCardinality, count& clusterCardinalityPrime,
                                count& elementCardinality, count& elementCardinalityPrime,
                                double symDiff, double& size, double& quality);

    count minCut4(index s, index t, std::vector<int>& bestBelongs,
                  count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                  void (Correspondences::*decrementS)(index newCluster), void (Correspondences::*decrementT)(index newCluster),
                  count (Correspondences::*getQualityTrivialSolutions)(index s, index t, bool& tWins),
                  index (Correspondences::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
                                                         std::vector<index>& insertedAt, count& position,
                                                         std::vector<int>& bestBelongs,
                                                         count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                                                         std::vector<bool>& doneWith));

    void bucketSort(count cMax,
                    std::vector<count>& gomoryHuParent, std::vector<count>& cutWithGomoryHuParent,
                    std::vector<count>& sortedGomoryHuParent, std::vector<count>& sortedCutWithGomoryHuParent);

    count incrementS2(index newCluster);
    count incrementT2(index newCluster);
    void decrementS2(index newCluster);
    void decrementT2(index newCluster);
    count getQualityTrivialSolutions2(index s, index t, bool& tWins);

    count incrementS3(index newCluster);
    count incrementT3(index newCluster);
    void decrementS3(index newCluster);
    void decrementT3(index newCluster);
    count getQualityTrivialSolutions3(index s, index t, bool& tWins);

    count incrementS4(index newCluster);
    count incrementT4(index newCluster);
    void decrementS4(index newCluster);
    void decrementT4(index newCluster);
    count getQualityTrivialSolutions4(index s, index t, bool& tWins);

    count gusfield(std::vector<index>& gomoryHuParent, std::vector<count>& cutWithGomoryHuParent,
                   index& bestS, index& bestT, std::vector<int>& bestBestBelongs,
                   count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                   void (Correspondences::*decrementS)(index newCluster), void (Correspondences::*decrementT)(index newCluster),
                   count (Correspondences::*getQualityTrivialSolutions)(index s, index t, bool& tWins),
                   index (Correspondences::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
                                                          std::vector<index>& insertedAt, count& position,
                                                          std::vector<int>& bestBelongs,
                                                          count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                                                          std::vector<bool>& doneWith),
                   count (Correspondences::*minCut)(index s, index t, std::vector<int>& bestBelongs,
                                                    count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                                                    void (Correspondences::*decrementS)(index newCluster), void (Correspondences::*decrementT)(index newCluster),
                                                    count (Correspondences::*getQualityTrivialSolutions)(index s, index t, bool& tWins),
                                                    index (Correspondences::*greedyDescent)(index s, index t, count& bestFrac, count& currentFrac,
                                                                                           std::vector<index>& insertedAt, count& position,
                                                                                           std::vector<int>& bestBelongs,
                                                                                           count (Correspondences::*incrementS)(index newCluster), count (Correspondences::*incrementT)(index newCluster),
                                                                                           std::vector<bool>& doneWith)));

    double evaluateAllCorrespondences(std::vector<index>& gomoryHuParent, std::vector<count>& cutWithGomoryHuParent);


    bool readSegmentationNormalizeClusters(std::string filename, Partition& partition,
                                           std::map<index, index>& clusterID2segmentationID);

    bool writePartition(std::string filename, Partition& partition);

    bool writeSegmentation(std::string filename, Partition& partition,
                           std::map<index, index>& clusterID2segemntationID,
                           std::vector<index>& old2newElemet);

    bool writeShuffledPartition(std::string filename, Partition& partition,
                                std::vector<index>& elementPermutation,
                                std::vector<index>& clusterPermutation);

    void instruction();

    double run(const Partition& partition1, const Partition& partition2);

    /**
     * Detect runs the correspondences detection algorithm for the passed in partitions with the specified level.
     * It's results are stored in the instance's attributes.
     *
     * Available Levels:
     *
     * level == 1: cuts through bipartite graph, not yet implemented
     * level == 2: nontrivial (one-sided) correspondences
     * level == 3: non-degenerate correspondences
     * level == 4: mutual correspondences
     *
     */
    void detect(unsigned int level, const Partition& partitionA, const Partition& partitionB);
};
} /* namespace NetworKit */
#endif	/* CORRESPONDENCES_H */

