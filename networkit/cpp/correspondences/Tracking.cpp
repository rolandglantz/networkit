/*
 * Tracking.cpp
 *
 *  Created on: July 8, 2017
 *      Author: Paul Skopnik
 */

#include "Tracking.h"

#include <cmath>
#include <deque>
#include <fstream>
#include <initializer_list>
#include <numeric>
#include <set>
#include <string>

namespace NetworKit {

std::vector<count> partitionSubsetSizes(const Partition& p) {
	std::vector<count> sizes(p.upperBound(), 0);

	p.forEntries([&](index, index s){
		if (s != none)
			sizes[s] += 1;
	});

	return sizes;
}

GHGraph::GHGraph(
	const std::vector<index>& parents,
	const std::vector<count>& weights,
	std::vector<NodeData> nodes,
	std::vector<index> children
	) : size(parents.size()),
		parents(parents),
		weights(weights),
		nodes(nodes),
		children(children) {
	}


GHGraph GHGraph::build(const std::vector<index>& parents, const std::vector<count>& weights) {
	std::vector<GHGraph::NodeData> nodes(parents.size(), {0, 0, 1, 0});
	std::vector<index> children(parents.size() - 1);

	for (auto it = parents.cbegin(); it != parents.cend(); ++it) {
		if (*it < parents.size())
			++nodes[*it].totalChildren;
	}

	std::deque<index> queue;

	index cumChildrenIndex = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it) {
		it->childrenIndex = cumChildrenIndex;
		cumChildrenIndex += it->totalChildren;

		if (it->totalChildren == 0)
			queue.push_back(it - nodes.begin());
	}

	while (!queue.empty()) {
		index i = queue.front();
		queue.pop_front();

		index p = parents[i];

		if (p < parents.size()) {
			nodes[p].subtreeSize += nodes[i].subtreeSize;

			children[nodes[p].childrenIndex + nodes[p].foundChildren] = i;

			if (++nodes[p].foundChildren == nodes[p].totalChildren)
				queue.push_back(p);
		}
	}

	return GHGraph(
		parents,
		weights,
		std::move(nodes),
		std::move(children)
	);
}

GHGraph::Neighbours GHGraph::neighbours(index node) const {
	return GHGraph::Neighbours(*this, node);
}


std::vector<count> calculateSubtreeSizes(const std::vector<index>& parents) {
	struct NodeData {
		count totalChildren, foundChildren;
		count subtreeSize;
	};

	std::vector<NodeData> nodes(parents.size(), {0, 0, 1});
	std::vector<count> subtreeSizes(parents.size(), 0);

	for (auto it = parents.cbegin(); it != parents.cend(); ++it) {
		if (*it < parents.size())
			++nodes[*it].totalChildren;
	}

	std::deque<index> queue;

	for (auto it = nodes.cbegin(); it != nodes.cend(); ++it) {
		if (it->totalChildren == 0)
			queue.push_back(it - nodes.cbegin());
	}

	while (!queue.empty()) {
		index i = queue.front();
		queue.pop_front();

		index p = parents[i];

		if (p < parents.size()) {
			nodes[p].subtreeSize += nodes[i].subtreeSize;

			if (++nodes[p].foundChildren == nodes[p].totalChildren)
				queue.push_back(p);
		}

		subtreeSizes[i] = nodes[i].subtreeSize;
	}

	return subtreeSizes;
}

CorrespondencesExtractor::CorrespondencesExtractor(Correspondences& correspondences)
	: correspondences(correspondences) {}

TimestepData::Correspondence  CorrespondencesExtractor::extract(const std::vector<index>& parts) const {
	TimestepData::Correspondence correspondence = {
		parts,
		std::vector<index>(0),

		std::vector<index>(parts.size(), 0),
		std::vector<index>(0),

		std::vector<std::vector<count>>(0)
	};

	count intersection = 0;
	count sizePPrime = 0;
	count sizeP = 0;

	for (index i = 0; i < this->correspondences.cardPartition2; ++i) {
		count sum = 0;

		for (index part : parts) {
			sum += this->correspondences.distributions[part][i];
		}

		if (2 * sum > this->correspondences.cardinalityOfCluster2[i]) {
			correspondence.pPrime.push_back(i);
			correspondence.pPrimeSizes.push_back(sum);

			for (auto it = parts.cbegin(); it != parts.cend(); ++it) {
				correspondence.pSizes[it - parts.cbegin()]
					+= this->correspondences.distributions[*it][i];

				correspondence.pToPPrimeSizes[it - parts.cbegin()].push_back(
					this->correspondences.distributions[*it][i]
				);
			}

			intersection += sum;
			sizePPrime += this->correspondences.cardinalityOfCluster2[i];
		}
	}

	for (index part : parts) {
		sizeP += this->correspondences.cardinalityOfCluster1[part];
	}

	correspondence.intersectionSize = intersection;
	correspondence.unionSize = sizePPrime + sizeP - intersection;

	return correspondence;
}

TimestepData::Timestep LeafExpansion::analyseFirst(const Partition& partition) const {
	// DEBUG
	// Write out the partition
	std::ofstream partFile("1_parts.tsv");
	for (index i = 0; i < partition.upperBound(); ++i) {
		auto s = partition.getMembers(i);
		for (auto it = s.cbegin(); it != s.cend(); ++it) {
			partFile << *it << " ";
		}
		partFile << std::endl;
	}
	partFile.close();

	return {
		partitionSubsetSizes(partition),
		std::vector<TimestepData::Correspondence>()
	};
}

TimestepData::Timestep LeafExpansion::analyseStep(
	const Partition& partition1,
	const Partition& partition2
) const {
	TimestepData::Timestep timestep;
	timestep.partitionSizes = partitionSubsetSizes(partition2);

	std::vector<TimestepData::Correspondence>& correspondences = timestep.correspondences;

	Correspondences c;
	c.detect(2, partition1, partition2);
	CorrespondencesExtractor corrExtractor(c);

	std::vector<count> subtreeSizes = calculateSubtreeSizes(c.gomoryHuParent);


	// DEBUG
	// Count timesteps
	static count timestepIndex = 1;
	++timestepIndex;

	// DEBUG
	// Write out the partition
	std::ofstream partFile(std::to_string(timestepIndex).append("_parts.tsv").c_str());
	for (index i = 0; i < partition2.upperBound(); ++i) {
		auto s = partition2.getMembers(i);
		for (auto it = s.cbegin(); it != s.cend(); ++it) {
			partFile << *it << " ";
		}
		partFile << std::endl;
	}
	partFile.close();

	// DEBUG
	// Write out gomory hu tree
	std::ofstream file(std::to_string(timestepIndex).append("_ght.dot").c_str());
	file << Strings::digraph << Strings::opCurly;
	for (auto it = c.gomoryHuParent.cbegin(); it != c.gomoryHuParent.cend(); ++it) {
		if (*it < c.gomoryHuParent.size()) {
			file << "p" << Strings::underscore << it - c.gomoryHuParent.cbegin();
			file << Strings::arrow << "p" << Strings::underscore << *it;
			file << Strings::opBracket << "label=" << c.cutWithGomoryHuParent[it - c.gomoryHuParent.cbegin()] << Strings::clBracket;
			file << Strings::semicolon;
		} else {
			file << "p" << Strings::underscore << it - c.gomoryHuParent.cbegin() << Strings::semicolon;
		}
	}
	file << Strings::clCurly;
	file.close();


	// Analyse using leaf expansion

	for (auto it = subtreeSizes.cbegin(); it != subtreeSizes.cend(); ++it) {
		if (*it == 1) {
			// Leaf node
			index parent = c.gomoryHuParent[it - subtreeSizes.cbegin()];
			if (subtreeSizes[parent] == 2) {
				// Parent node only has one child
				if (c.cutWithGomoryHuParent[parent] <= c.cutWithGomoryHuParent[it - subtreeSizes.cbegin()]) {
					// The "cut" is at least as cheap when including the parent
					//  -> include the parent node in the correspondence
					std::vector<count> parts;
					parts.push_back(it - subtreeSizes.cbegin());
					parts.push_back(parent);
					correspondences.push_back(corrExtractor.extract(parts));
				} else {
					// The "cut" is cheaper when excluding the parent
					//  -> only include the leaf in the correspondence
					correspondences.push_back(corrExtractor.extract(std::vector<index>(1, it - subtreeSizes.cbegin())));
				}
			} else {
				// Parent node has more than one child -> ignore parent, only include leaf in the correspondence
				correspondences.push_back(corrExtractor.extract(std::vector<index>(1, it - subtreeSizes.cbegin())));
			}
		} else if (*it == partition2.upperBound()) {
			// Root node
			std::vector<count> childrenNo = this->calculateChildrenNo(c.gomoryHuParent);
			if (childrenNo[it - subtreeSizes.cbegin()] == 0) {
				// The root node is the only node (only one partition)
				//  -> search for correspondence with this one node
				correspondences.push_back(corrExtractor.extract(std::vector<index>(1, it - subtreeSizes.cbegin())));
			} else if (childrenNo[it - subtreeSizes.cbegin()] == 1) {
				// Root node only has one child, there is a distinct cut separating it from all other nodes
				index child;
				for (auto childIt = c.gomoryHuParent.cbegin(); childIt != c.gomoryHuParent.cend(); ++childIt) {
					// Search for the root node's child
					if (*childIt == static_cast<unsigned int>(it - subtreeSizes.cbegin()))
						child = childIt - c.gomoryHuParent.cbegin();
				}
				if (childrenNo[child] == 1) {
					// The root node's child only has one child, there is a distinct cut separating the two from other nodes
					index childsChild;
					for (auto childIt = c.gomoryHuParent.cbegin(); childIt != c.gomoryHuParent.cend(); ++childIt) {
						// Search for the root node's child's child
						if (*childIt == static_cast<unsigned int>(it - subtreeSizes.cbegin()))
							childsChild = childIt - c.gomoryHuParent.cbegin();
					}
					if (c.cutWithGomoryHuParent[childsChild] <= c.cutWithGomoryHuParent[child]) {
						// The "cut" is at least as cheap when including the child
						//  -> include the child node in the correspondence
						std::vector<count> parts;
						parts.push_back(it - subtreeSizes.cbegin());
						parts.push_back(child);
						correspondences.push_back(corrExtractor.extract(parts));
					} else {
						// The "cut" is cheaper excluding the child
						//  -> only include root node in correspondence
						correspondences.push_back(corrExtractor.extract(std::vector<index>(1, it - subtreeSizes.cbegin())));
					}
				} else {
					// Child has no more than one child, there is no single "cut"
					//  -> only include root node in correspondence
					correspondences.push_back(corrExtractor.extract(std::vector<index>(1, it - subtreeSizes.cbegin())));
				}
			// } else {
			// 	// Ignore as there is no single "cut" using the root
			}
		}
	}

	return timestep;
}

std::vector<count> LeafExpansion::calculateChildrenNo(const std::vector<index>& parents) const {
	std::vector<count> childrenNo(parents.size(), 0);

	for (auto it = parents.cbegin(); it != parents.cend(); ++it) {
		if (*it < parents.size())
			++childrenNo[*it];
	}

	return childrenNo;
}

DCGOwnershipExtractor::DCGOwnershipExtractor(
	const DynamicCommunitiesGenerator& generator,
	const std::vector<std::vector<index>>& parts
) : generator(GeneratorState(generator)), parts(parts) {
	// Calculate subcluster to part mapping
	this->subclusterToPart = this->invertParts(parts);
	// this->noOfParts = parts.size();
}

std::vector<DCGTimestepData::Ownership> DCGOwnershipExtractor::extract(
	const Partition& partition
) const {
	return this->extract(partition, partitionSubsetSizes(partition));
}

std::vector<DCGTimestepData::Ownership> DCGOwnershipExtractor::extract(
	const Partition& partition,
	const std::vector<count>& subsetSizes
) const {
	std::vector<DCGTimestepData::Ownership> ownership(partition.upperBound());

	// Build ownership info
	std::vector<count> partCounters;
	for (index i = 0; i < partition.upperBound(); ++i) {
		partCounters.assign(this->parts.size(), 0);

		auto memberSet = partition.getMembers(i);
		for (index member : memberSet) {
			index subcluster = this->generator.getIndividualSubluster(member);
			++partCounters[this->subclusterToPart[subcluster]];
		}


		auto largestPartIt = std::max_element(partCounters.begin(), partCounters.end());
		count largestPartSize = *largestPartIt;

		*largestPartIt = 0;
		auto scndLargestPartIt = std::max_element(partCounters.cbegin(), partCounters.cend());

		double ownershipMargin = static_cast<double>(largestPartSize - *scndLargestPartIt)
			/ static_cast<double>(2 * subsetSizes[i]);

		ownership[i] = {
			static_cast<index>(largestPartIt - partCounters.cbegin()),
			static_cast<double>(largestPartSize) / subsetSizes[i],
			ownershipMargin
		};
	}

	return ownership;
}

std::vector<index> DCGOwnershipExtractor::invertParts(const std::vector<std::vector<index>>& parts) {
	std::vector<index> inverted(1, 0);

	index partsInd;

	for (auto partsIt = parts.cbegin(); partsIt != parts.cend(); ++partsIt) {
		partsInd = partsIt - parts.cbegin();
		for (auto subcluster : *partsIt) {
			if (inverted.size() <= subcluster)
				inverted.resize(subcluster + 1);

			inverted[subcluster] = partsInd;
		}
	}

	return inverted;
}

DCGLeafExpansion::DCGLeafExpansion(
	const DynamicCommunitiesGenerator& generator,
	const std::vector<std::vector<index>>& parts
) : ownershipExtractor(generator, parts) {
}

DCGTimestepData::Timestep DCGLeafExpansion::analyseFirst(const Partition& partition) const {
	TimestepData::Timestep simpleTimestep = this->leafExpansion.analyseFirst(partition);

	return {
		simpleTimestep.partitionSizes,
		this->ownershipExtractor.extract(partition),
		simpleTimestep.correspondences
	};
}

DCGTimestepData::Timestep DCGLeafExpansion::analyseStep(
	const Partition& partition1,
	const Partition& partition2
) const {
	TimestepData::Timestep simpleTimestep = this->leafExpansion.analyseStep(
		partition1, partition2);

	return {
		simpleTimestep.partitionSizes,
		this->ownershipExtractor.extract(partition2),
		simpleTimestep.correspondences
	};
}

CheapestSetsGenerator::CheapestSetsGenerator(Correspondences& c, bool)
 : reachedEnd(true),
	c(c),
	parents(c.gomoryHuParent),
	weights(c.cutWithGomoryHuParent),
	graph(c.gomoryHuParent, c.cutWithGomoryHuParent) {}

CheapestSetsGenerator::CheapestSetsGenerator(Correspondences& c)
	: CheapestSetsGenerator(c, std::vector<index>(c.gomoryHuParent.size(), 0), 0, 0) {}

CheapestSetsGenerator::CheapestSetsGenerator(
	Correspondences& c, std::vector<index> partSets, index rootSet, count cost
) : c(c),
	parents(c.gomoryHuParent),
	weights(c.cutWithGomoryHuParent),
	partSets(partSets),
	graph(GHGraph::build(c.gomoryHuParent, c.cutWithGomoryHuParent))
{
	this->weightIndices = this->indexSortWeight();
	this->weightIt = this->weightIndices.cbegin();
	this->weightItEnd = this->weightIndices.cend();

	this->minSetIndex = *std::max_element(
		this->partSets.cbegin(),
		this->partSets.cend()
	);

	this->refResultSet(0) = this->createRootResult(rootSet, cost);
}

std::vector<index> CheapestSetsGenerator::indexSortWeight() const {
	std::vector<index> weightIndices(this->graph.getSize(), 0);

	std::iota(weightIndices.begin(), weightIndices.end(), 0);

	std::sort(weightIndices.begin(), weightIndices.end(),
		[this](int i, int j) {
			return this->weights[i] < this->weights[j];
		}
	);

	return weightIndices;
}

HT::Result CheapestSetsGenerator::createRootResult(index rootSet, count cost) const {
	std::vector<index> nodes;

	for (auto it = this->partSets.cbegin(); it != this->partSets.cend(); ++it) {
		if (*it == rootSet)
			nodes.push_back(it - this->partSets.cbegin());
	}

	return {
		nodes,
		std::vector<index>(0),
		cost
	};
}

void CheapestSetsGenerator::advance() {
	while (this->weightIt != this->weightItEnd) {
		index parentSetIndex = this->partSets[*this->weightIt];

		this->partSets[*this->weightIt] = this->minSetIndex + this->sets.size();
		this->partSets[this->parents[*this->weightIt]] = this->minSetIndex + this->sets.size() + 1;

		this->propagateResult(*this->weightIt, parentSetIndex);
		this->propagateResult(this->parents[*this->weightIt], parentSetIndex);

		++this->weightIt;

		if (this->weightIt != this->weightItEnd
			&& this->resultQueue.top()->cost < this->weights[*this->weightIt]
		)
			// Early exit possible: All not yet explored results have a cost of at least the
			//  upcoming weight. The cost of the result on top of the result queue is already
			//  cheaper and can thus be yielded next.
			return;
	}

	this->reachedEnd = true;
}

HT::Result& CheapestSetsGenerator::propagateResult(index l, index parentSetIndex) {
	index setIndex = this->partSets[l];

	HT::Result& result = this->refResultSet(setIndex);
	result.cost = this->bfsExpand(l, parentSetIndex, result.p);

	this->resultQueue.push(&result);

	return result;
}

count CheapestSetsGenerator::bfsExpand(index from, index expandInto, std::vector<index>& nodes) {
	count outerWeight = 0;
	index ownSet = this->partSets[from];

	nodes.clear();

	std::deque<index> queue;
	queue.push_back(from);

	while (!queue.empty()) {
		auto neighbours = this->graph.neighbours(queue.front());

		for (GHGraph::Edge edge : neighbours) {
			if (partSets[edge.b] == expandInto) {
				queue.push_back(edge.b);
				this->partSets[edge.b] = ownSet;
			} else if (partSets[edge.b] != ownSet) {
				outerWeight += edge.weight;
			}
		}

		nodes.push_back(queue.front());
		queue.pop_front();
	}

	return outerWeight;
}

HT::Result& CheapestSetsGenerator::refResultSet(int set) {
	index vectorIndex = set - this->minSetIndex;

	if (vectorIndex >= this->sets.size())
		this->sets.resize(vectorIndex + 1, {});

	return this->sets[vectorIndex];
}


SmallestMutual::SmallestMutual(Correspondences& c)
: c(c),
	parents(c.gomoryHuParent),
	weights(c.cutWithGomoryHuParent),
	graph(GHGraph::build(c.gomoryHuParent, c.cutWithGomoryHuParent)),
	partSets(c.gomoryHuParent.size(), 0)
{
}


std::vector<HT::Result> SmallestMutual::run() {
	this->weightIndices = this->indexSortWeight();

	HT::Result rootResult = this->createRootResult();

	HT::ResultSet<count> resultSet = this->processTree(0, rootResult);
	// HT::ResultSet<count> resultSet = this->processTree(0, HT::Result());

	if (resultSet.results.size() < 1
		|| resultSet.results.front().p.size() < this->graph.getSize()
	)
		// Don't return rootResult
		return resultSet.results;
	else
		return {};
}

std::vector<index> SmallestMutual::indexSortWeight() const {
	std::vector<index> weightIndices(this->graph.getSize(), 0);

	std::iota(weightIndices.begin(), weightIndices.end(), 0);

	std::sort(weightIndices.begin(), weightIndices.end(),
		[this](int i, int j) {
			return this->weights[i] < this->weights[j];
		}
	);

	return weightIndices;
}

HT::Result SmallestMutual::createRootResult() const {
	std::vector<index> p(this->graph.getSize());

	HT::Result result{
		std::vector<index>(this->graph.getSize())
	};

	std::iota(result.p.begin(), result.p.end(), 0);

	this->calculatePPrime(result);

	return result;
}

count SmallestMutual::bfsExpand(index from, index expandInto, std::vector<index>& nodes) {
	count outerWeight = 0;
	index ownSet = this->partSets[from];

	nodes.clear();

	std::deque<index> queue;
	queue.push_back(from);

	while (!queue.empty()) {
		auto neighbours = this->graph.neighbours(queue.front());

		for (GHGraph::Edge edge : neighbours) {
			if (partSets[edge.b] == expandInto) {
				queue.push_back(edge.b);
				this->partSets[edge.b] = ownSet;
			} else if (partSets[edge.b] != ownSet) {
				outerWeight += edge.weight;
			}
		}

		nodes.push_back(queue.front());
		queue.pop_front();
	}

	return outerWeight;
}

HT::ResultSet<count> SmallestMutual::processTree(index setIndex, const HT::Result parent) {
	for (auto it = this->weightIndices.cbegin(); it != this->weightIndices.cend() - 1; ++it) {
		if (this->partSets[*it] != setIndex || this->partSets[this->parents[*it]] != setIndex)
			continue;

		// Copy partSets to be able to revert the split
		std::vector<index> previousPartSets(this->partSets);

		index parentSetIndex = this->partSets[*it];

		// Assign parts a and b of edge (a -> b) new set indices
		this->partSets[*it] = ++this->maxSetIndex;
		this->partSets[this->parents[*it]] = ++this->maxSetIndex;

		HT::Result resultA = this->buildResult(*it, parentSetIndex);
		HT::Result resultB = this->buildInverseResult(this->parents[*it], parentSetIndex,
			resultA, parent.pPrime);

		bool isMutualA = this->isMutual(resultA);
		bool isMutualB = this->isMutual(resultB);

		HT::ResultSet<count> treeA, treeB;

		if (isMutualA)
			treeA = this->processTree(
				this->partSets[*it],
				resultA
			);

		if (isMutualB)
			treeB = this->processTree(
				this->partSets[this->parents[*it]],
				resultB
			);

		if (isMutualA && isMutualB) {
			HT::ResultSet<count> combined = {
				0,
				std::vector<HT::Result>(treeA.results)
			};
			combined.results.insert(
				combined.results.end(),
				// Breaking from c++11 to c++14:
				// combined.results.cend(),
				treeB.results.cbegin(),
				treeB.results.cend()
			);

			return combined;
		} else if (isMutualA)
			return treeA;
		else if (isMutualB)
			return treeB;
		else
			this->partSets.swap(previousPartSets);
	}

	return HT::ResultSet<count>{
		0,
		{parent}
	};
}

HT::Result SmallestMutual::buildResult(index l, index parentSetIndex) {
	HT::Result result;

	this->bfsExpand(l, parentSetIndex, result.p);

	this->calculatePPrime(result);

	return result;
}

HT::Result SmallestMutual::buildInverseResult(
	index l,
	index parentSetIndex,
	const HT::Result& other,
	const std::vector<index>& superSet
) {
	HT::Result result;

	this->bfsExpand(l, parentSetIndex, result.p);

	this->calculateInversePPrime(result, superSet, other.pPrime);

	return result;
}


void SmallestMutual::calculatePPrime(HT::Result& result) const {
	for (index i = 0; i < this->c.cardPartition2; ++i) {
		count sum = 0;

		for (index part : result.p) {
			sum += this->c.distributions[part][i];
		}

		if (2 * sum > this->c.cardinalityOfCluster2[i])
			result.pPrime.push_back(i);
		// Maybe in pPrime
		// else if (2 * sum == this->c.cardinalityOfCluster2[i])
	}
}

void SmallestMutual::calculateInversePPrime(
	HT::Result& result,
	const std::vector<index>& superSet,
	const std::vector<index>& other
) const {
	std::set<index> all(superSet.cbegin(), superSet.cend());

	for (auto e : other) {
		all.erase(all.find(e));
	}

	std::copy(all.cbegin(), all.cend(), std::back_inserter(result.pPrime));
}

bool SmallestMutual::isMutual(const HT::Result& result) const {
	std::set<index> pSet;
	std::set<index> pSetMaybe;

	for (index i = 0; i < this->c.cardPartition1; ++i) {
		count sum = 0;

		for (index part : result.pPrime) {
			sum += this->c.distributions[i][part];
		}

		if (2 * sum > this->c.cardinalityOfCluster1[i])
			pSet.insert(i);
		else if (2 * sum == this->c.cardinalityOfCluster1[i])
			pSetMaybe.insert(i);
	}

	if (result.p.size() < pSet.size()
		|| result.p.size() > pSet.size() + pSetMaybe.size()
	)
		return false;

	count maybeCount = result.p.size() - pSet.size();

	for (index part : result.p) {
		if (pSet.count(part) != 1) {
			if (pSetMaybe.count(part) == 1)
				--maybeCount;
			else
				return false;
		}
	}

	if (maybeCount != 0)
		return false;

	return true;
}

} /* namespace NetworKit */
