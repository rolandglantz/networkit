/*
 * DynamicCommunitiesGenerator.cpp
 *
 *  Created on: May 26, 2017
 *      Author: Paul Skopnik
 */

#include "DynamicCommunitiesGenerator.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/Enforce.h"
#include "../auxiliary/Log.h"

#include <utility>
#include <string>
#include <random>

namespace NetworKit {

void infoVector(const std::vector<index> v, const char* prefix) {
	std::string s = prefix;

	s += "[";

	for (auto val : v) {
		s += std::to_string(val);
		s += " ";
	}

	s += "]";

	INFO(s.c_str());
}

void infoPartialTrees(const std::vector<std::vector<index>>& partialTrees,
		const std::vector<DynamicCommunitiesGenerator::Subcluster>& subclusters,
		const char* prefix) {
	std::string s = prefix;


	for (auto partialTree : partialTrees) {
		s += "[";
		for (auto node : partialTree) {
			s += "(";
			s += std::to_string(node);
			s += " @ ";
			s += std::to_string(subclusters[node].postOrder);
			s += " p=";
			s += std::to_string(subclusters[node].parent);
			s += ")";
			s += " ";
		}
		s += "] ";
	}


	INFO(s.c_str());
}

void infoCluster(const std::vector<index>& cluster,
		const std::vector<DynamicCommunitiesGenerator::Subcluster>& subclusters,
		const char* prefix) {
	std::string s = prefix;


	s += "[";
	for (auto node : cluster) {
		s += "(";
		s += std::to_string(node);
		s += " @ ";
		s += std::to_string(subclusters[node].postOrder);
		s += " p=";
		s += std::to_string(subclusters[node].parent);
		s += ")";
		s += " ";
	}
	s += "] ";


	INFO(s.c_str());
}

SubclusterEdges::SubclusterEdges(SymmetricMatrix<double> affinities)
	: affinities(affinities),
		availableEdges(affinities.getN()),
		availableEdgesWeight(0),
		selectedEdgesWeight(0),
		availableEdgesNumber(0),
		selectedEdgesNumber(0) {
	for (count i = 1; i <= affinities.getN(); i++) {
		for (count j = i + 1; j <= affinities.getN(); j++) {
			this->availableEdges.set(i, j, 1);
			this->availableEdgesWeight += affinities.get(i, j);
		}
	}
}

count SubclusterEdges::getN() const {
	return this->affinities.getN();
}

void SubclusterEdges::selectEdge(index i, index j) {
	if (this->availableEdges.get(i, j) == 1) {
		this->availableEdges.set(i, j, 0);

		double edgeWeight = this->affinities.get(i, j);

		this->availableEdgesWeight -= edgeWeight;
		this->selectedEdgesWeight += edgeWeight;

		--this->availableEdgesNumber;
		++this->selectedEdgesNumber;
	}
	// else {
	// 	// Already selected
	// }
}

void SubclusterEdges::unselectEdge(index i, index j) {
	if (this->availableEdges.get(i, j) == 0) {
		this->availableEdges.set(i, j, 1);

		double edgeWeight = this->affinities.get(i, j);

		this->availableEdgesWeight += edgeWeight;
		this->selectedEdgesWeight -= edgeWeight;

		++this->availableEdgesNumber;
		--this->selectedEdgesNumber;
	}
	// else {
	// 	// Already unselected
	// }
}

double SubclusterEdges::getAvailableEdgesWeight() const {
	return this->availableEdgesWeight;
}

double SubclusterEdges::getSelectedEdgesWeight() const {
	return this->selectedEdgesWeight;
}

count SubclusterEdges::getAvailableEdgesNumber() const {
	return this->availableEdgesNumber;
}

count SubclusterEdges::getSelectedEdgesNumber() const {
	return this->selectedEdgesNumber;
}

SubclusterEdges::Index SubclusterEdges::availableEdgeAtDistribution(double affinity) const {
	if (affinity >= this->availableEdgesWeight || affinity < 0) {
		throw("SubclusterEdges::availableEdgeAtDistribution: Invalid affinity value");
	}

	double cumulativeAffinity = 0;

	auto availableIt = this->availableEdges.cbegin();
	auto affinitiesIt = this->affinities.cbegin();

	for (; availableIt != this->availableEdges.cend(); ++availableIt, ++affinitiesIt) {
		if (*availableIt == 0)
			continue;

		cumulativeAffinity += *affinitiesIt;
		if (cumulativeAffinity > affinity) {
			break;
		}
	}
	// TODO
	// If there is no break, default to last edge
	// Current approach may not work when availableIt is already at the post-end (::end()) stage.
	// (Due to floating point calculation errors)

	Aux::enforce(availableIt != this->availableEdges.cend(),
		"SubclusterEdges::availableEdgeAtDistribution: Invalid iterator (implementation is wip)");

	return this->availableEdges.indexFromIterator(availableIt);
}

AffinityGenerator::Affinities AffinityGenerator::halfHalf(count k, double w) const {
	AffinityGenerator::Affinities affinities(k);

	double intraEdgeWeight = 1;
	double interEdgeWeight = w;

	for (index i = 1; i <= k / 2; ++i) {
		for (index j = i + 1; j <= k / 2; ++j) {
			affinities.set(i, j, intraEdgeWeight);
		}
	}

	for (index i = k / 2 + 1; i <= k; ++i) {
		for (index j = i + 1; j <= k; ++j) {
			affinities.set(i, j, intraEdgeWeight);
		}
	}

	for (index i = 1; i <= k / 2; ++i) {
		for (index j = k / 2 + 1; j <= k; ++j) {
			affinities.set(i, j, interEdgeWeight);
		}
	}

	return affinities;
}

DynamicCommunitiesGenerator::DynamicCommunitiesGenerator(const Parameters parameters)
	: parameters(parameters),
		timestep(0),
		preGenerated(false),
		individuals(parameters.n + 1, {0, 0}),
		clusters(parameters.affinities.getN() + 1, std::vector<index>()),
		subclusterEdges(parameters.affinities),
		subclusters(parameters.affinities.getN() + 1, {0, 0, 0}) {
	// Assign each subclusters to its own clusters
	for (index i = 1; i < this->subclusters.size(); ++i) {
		this->subclusters[i].cluster = i;
		this->clusters[i].push_back(i);
	}

	count perSubcluster = parameters.n / (this->subclusters.size() - 1);
	count left = parameters.n % (this->subclusters.size() - 1);

	index i = 1;
	for (index q = 1; q < this->subclusters.size(); ++q) {
		for (index c = 0; c < perSubcluster; ++c) {
			this->individuals[i++] = {q, q};
		}

		if (left > 0) {
			this->individuals[i++] = {q, q};
			--left;
		}
	}

	#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATE
	(GeneratorValidator(*this)).validate();
	#endif
}

Partition DynamicCommunitiesGenerator::next() {
	if (!this->preGenerated)
		this->preGenerate();

	#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATE
	(GeneratorValidator(*this)).validateClusterNum();
	#endif

	this->performMerge();
	this->performSplit();
	this->performIndividualMoves();

	Partition partition(this->individuals.size() - 1);
	partition.setUpperBound(this->subclusters.size());

	for (auto it = this->individuals.cbegin() + 1; it != this->individuals.cend(); ++it) {
		partition.addToSubset((*it).subcluster - 1, it - this->individuals.cbegin() - 1);
	}

	return partition;
}

void DynamicCommunitiesGenerator::preGenerate() {
	// Perform k - l merge operations, then the designated number of clusters is reached.
	while (this->clusters.size() - 1 > this->parameters.l)
		this->performMerge();

	// for (int i = 0; i < 100; ++i) {
	// 	this->performMerge();
	// 	this->performSplit();
	// 	this->performIndividualMoves();
	// }

	this->preGenerated = true;
}

void DynamicCommunitiesGenerator::performMerge() {
	double availableEdgesWeight = this->subclusterEdges.getAvailableEdgesWeight();

	// Choose r \in [0, availableEdgesWeight)
	double r = Aux::Random::real(availableEdgesWeight);

	SubclusterEdges::Index index = this->subclusterEdges.availableEdgeAtDistribution(r);

	this->addSubclusterEdge(index.i, index.j);

	#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATE
	(GeneratorValidator(*this)).validate();
	#endif
}

void DynamicCommunitiesGenerator::performSplit() {
	count numberOfEdges = this->subclusters.size() - this->clusters.size();

	count sum = 0;
	index i = 0;

	// Choose r \in [0, numberOfEdges - 1]
	index r = Aux::Random::integer(numberOfEdges - 1);

	// Skip subclusters[0], as it is not used
	for (auto it = this->subclusters.cbegin() + 1; it != this->subclusters.cend(); ++it) {
		if ((*it).parent == 0)
			continue;

		if (++sum > r) {
			// *it was chosen
			i = it - this->subclusters.cbegin();
			break;
		}
	}

	// Development artefact
	// Aux::enforce(i != 0, "Found i = 0");
	// Aux::enforce(this->subclusters[i].parent != 0, "Found subcluster edge with parent == 0!");

	this->removeSubclusterEdge(i, this->subclusters[i].parent);

	#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATE
	(GeneratorValidator(*this)).validate();
	#endif
}

void DynamicCommunitiesGenerator::performIndividualMoves() {
	count no_of_moves = this->parameters.p_move_v * this->parameters.n;

	std::vector<index> individuals, subcluster;
	individuals.reserve(no_of_moves);
	subcluster.reserve(no_of_moves);

	std::mt19937_64& urng = Aux::Random::getURNG();
	std::uniform_int_distribution<uint64_t> dist{1, this->parameters.n};

	index i;

	for (; individuals.size() < no_of_moves;) {
		i = dist(urng);

		// Ignore, if the individual is already marked to be moved
		if (this->individuals[i].subcluster == 0)
			continue;

		individuals.push_back(i);
		subcluster.push_back(this->individuals[i].subcluster);
		this->individuals[i].subcluster = 0;
	}

	for (i = 0; i < no_of_moves; ++i) {
		this->moveIndividual(individuals[i], subcluster[i]);
	}

	#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATE
	(GeneratorValidator(*this)).validate();
	#endif
}

void DynamicCommunitiesGenerator::moveIndividual(index v, index subcluster) {
	index q;

	if (subcluster == this->individuals[v].homeSubcluster) {
		// In home subcluster
		q = Aux::Random::integer(1, this->subclusters.size() - 1);

		if (q >= this->individuals[v].homeSubcluster)
			++q;
	} else {
		// Not in home subcluster
		double p = Aux::Random::probability();

		if (p < this->parameters.alpha) {
			// Decided to move to home subcluster
			q = this->individuals[v].homeSubcluster;
		} else {
			// Decided to move to arbitrary subcluster, but not the home subcluster
			q = Aux::Random::integer(1, this->subclusters.size() - 2);

			if (q >= this->individuals[v].homeSubcluster)
				// Skip the home subcluster, if necessary
				++q;
			if (q >= subcluster)
				// Skip the current subcluster, if necessary
				++q;
		}
	}

	this->individuals[v].subcluster = q;
}

void DynamicCommunitiesGenerator::addSubclusterEdge(index i, index j) {
	// TODO
	//  Figure out which cluster to make parent, child
	this->makeClusterRoot(j);

	// TODO? Is this a copy? Use a reference instead?
	std::vector<index> iCluster = this->clusters[this->subclusters[i].cluster];
	std::vector<index> jCluster = this->clusters[this->subclusters[j].cluster];

	std::vector<index> mergedCluster;
	mergedCluster.reserve(iCluster.size() + jCluster.size());

	auto iClusterIt = iCluster.cbegin();
	auto jClusterIt = jCluster.cbegin();
	index iPostOrder = this->subclusters[i].postOrder;
	count jClusterSize = jCluster.size();

	// Indices of the clusters, newClusterIndex < unusedClusterIndex
	index newClusterIndex = this->subclusters[i].cluster,
		unusedClusterIndex = this->subclusters[j].cluster;
	if (unusedClusterIndex < newClusterIndex) {
		index temp = newClusterIndex;
		newClusterIndex = unusedClusterIndex;
		unusedClusterIndex = temp;
	}

	// TODO
	// Line too long
	for(; *iClusterIt != i; ++iClusterIt) {
		mergedCluster.push_back(*iClusterIt);
		this->subclusters[*iClusterIt].cluster = newClusterIndex;
	}
	for (; iClusterIt != iCluster.cend() && this->subclusters[*iClusterIt].postOrder <= iPostOrder; ++iClusterIt) {
		mergedCluster.push_back(*iClusterIt);
		this->subclusters[*iClusterIt].cluster = newClusterIndex;
	}
	for (; jClusterIt != jCluster.cend(); ++jClusterIt) {
		mergedCluster.push_back(*jClusterIt);
		this->subclusters[*jClusterIt].postOrder += iPostOrder;
		this->subclusters[*jClusterIt].cluster = newClusterIndex;
	}
	for(; iClusterIt != iCluster.cend(); ++iClusterIt) {
		mergedCluster.push_back(*iClusterIt);
		this->subclusters[*iClusterIt].postOrder += jClusterSize;
		this->subclusters[*iClusterIt].cluster = newClusterIndex;
	}

	this->subclusters[j].parent = i;

	index z = j;
	while (this->subclusters[z].parent != 0) {
		z = this->subclusters[z].parent;
		this->subclusters[z].postOrder += jClusterSize;
	}

	// Remove all edges between the two new clusters to the available edge pool
	for (iClusterIt = iCluster.cbegin(); iClusterIt != iCluster.cend(); ++iClusterIt) {
		for (jClusterIt = jCluster.cbegin(); jClusterIt != jCluster.cend(); ++jClusterIt) {
			this->subclusterEdges.selectEdge(*iClusterIt, *jClusterIt);
		}
	}

	// TODO move instead of copy
	this->clusters[newClusterIndex] = mergedCluster;
	if (unusedClusterIndex != this->clusters.size() - 1) {
		for (auto it = this->clusters.back().cbegin(); it != this->clusters.back().cend(); ++it)
			this->subclusters[*it].cluster = unusedClusterIndex;
		// TODO move instead of copy
		this->clusters[unusedClusterIndex] = this->clusters.back();
	}
	this->clusters.pop_back();
}

void DynamicCommunitiesGenerator::removeSubclusterEdge(index i, index j) {
	index oldClusterIndex = this->subclusters[i].cluster;

	// Development artefact
	// Aux::enforce(oldClusterIndex > 0, "oldClusterIndex invalid (>0)");
	// Aux::enforce(oldClusterIndex < this->clusters.size(),
	// 	"oldClusterIndex invalid (>=clusters.size)");
	// Aux::enforce(oldClusterIndex == this->subclusters[j].cluster,
	// 	"Edge invalid (cluster[i] != cluster[j]");

	std::vector<index> oldCluster = this->clusters[oldClusterIndex];

	std::vector<index> newParentCluster, newChildCluster;

	index postOrderBound = 0;
	index newChildClusterIndex = this->clusters.size();
	index iPostOrder = this->subclusters[i].postOrder;

	// oldCluster is preOrder sorted...

	auto it = oldCluster.cbegin();

	// Iterate over items before i
	for(; *it != i; ++it) {
		newParentCluster.push_back(*it);
		if (this->subclusters[*it].postOrder < iPostOrder
				&& this->subclusters[*it].postOrder >= postOrderBound)
			postOrderBound = this->subclusters[*it].postOrder + 1;
	}

	// TODO more elegantly
	// Iterate over items in the i subtree
	for (; it != oldCluster.cend() && this->subclusters[*it].postOrder <= iPostOrder; ++it) {
		newChildCluster.push_back(*it);
		this->subclusters[*it].postOrder -= postOrderBound;
		this->subclusters[*it].cluster = newChildClusterIndex;
	}

	// Iterate over items behind the i subtree (in terms of preOrder)
	for (; it != oldCluster.cend(); ++it) {
		newParentCluster.push_back(*it);
		this->subclusters[*it].postOrder -= newChildCluster.size();
	}

	// Reduce postOrder of items "above" i (transitive parents)
	index z = i;
	while(this->subclusters[z].parent != 0) {
		z = this->subclusters[z].parent;
		this->subclusters[z].postOrder -= newChildCluster.size();
	}

	// Re-add all edges between the two new clusters to the available edge pool
	for (auto pIt = newParentCluster.cbegin(); pIt != newParentCluster.cend(); ++pIt) {
		for (auto cIt = newChildCluster.cbegin(); cIt != newChildCluster.cend(); ++cIt) {
			this->subclusterEdges.unselectEdge(*pIt, *cIt);
		}
	}

	this->subclusters[i].parent = 0;
	// TODO move instead of copy
	this->clusters[oldClusterIndex] = newParentCluster;
	// TODO move instead of copy
	this->clusters.push_back(newChildCluster);
}

void DynamicCommunitiesGenerator::makeClusterRoot(index i) {
	this->rejoinPartialTrees(this->extractPartialTrees(i));

	#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATE
	(GeneratorValidator(*this)).validate();
	#endif
}

std::vector<std::vector<index>> DynamicCommunitiesGenerator::extractPartialTrees(index i) {
	// List of trees, return value
	std::vector<std::vector<index>> partialTrees;

	// Current state
	//  postOrderBound is an upper, outer bound for encountered postOrder values.
	//  postOrderBound = 0, reflecting that postOrder 0 has not been reached yet.
	//  removedSize stores the size of the extracted tree and is used to amend postOrderBound
	//  correctly after a partial tree has been comitted.
	std::vector<index> partialTree;
	index postOrderBound = 0, postOrderDelta, removedSize = 0;

	// Stacks used during processing
	std::vector<std::vector<index>> partialTreeStack;
	std::vector<index> postOrderDeltaStack;

	// TODO? Is this a copy? Use a reference instead?
	std::vector<index> cluster = this->clusters[this->subclusters[i].cluster];

	// Stack of root nodes of partial trees
	std::vector<index> rootStack;

	// Fill rootStack with all nodes on the path from cluster root to node i
	//  These will be the roots of the partial trees.
	index z = i;
	while (z != 0) {
		rootStack.push_back(z);
		z = this->subclusters[z].parent;
	}

	for (auto it = cluster.cbegin(); it != cluster.cend(); ++it) {
		if (!rootStack.empty() && *it == rootStack.back()) {
			// Encountered root node.

			// Only if there are elements in the partialTree, i.e. don't push the empty partialTree
			//  when encountering the cluster root in the first iteration
			if (partialTree.size() > 0) {
				// Push the partial tree and the postOrderDelta onto their stack

				// TODO move instead of copy
				partialTreeStack.push_back(partialTree);
				partialTree.clear();

				postOrderDeltaStack.push_back(postOrderDelta);
			}

			// Remove the encountered root node
			rootStack.pop_back();

			// Add the root to the partial tree.
			//  The postOrder is not changed yet, as it is required to recognise the end of the
			//  partial tree.
			partialTree.push_back(*it);
			postOrderDelta = postOrderBound;
		} else if (this->subclusters[*it].postOrder > this->subclusters[partialTree.front()].postOrder) {
			// Reached end of current partial tree.

			while (this->subclusters[*it].postOrder > this->subclusters[partialTree.front()].postOrder) {

				// Amend the postOrder of the partial tree's root node.
				//  Simply set to the tree size - 1 , i.e. last post order occuring in the tree instead
				//  of using postOrderDelta
				this->subclusters[partialTree.front()].postOrder = partialTree.size() - 1;

				// Restore postOrderDelta, add the size of the just extracted partial tree.
				removedSize += partialTree.size();
				postOrderDelta = postOrderDeltaStack.back() + removedSize;
				postOrderDeltaStack.pop_back();

				// Commit partial tree
				// TODO move instead of copy
				partialTrees.push_back(partialTree);

				// Restore partialTree
				// TODO move instead of copy
				partialTree = partialTreeStack.back();
				partialTreeStack.pop_back();
			}

			// Process encountered node.
			partialTree.push_back(*it);
			// Note: The post order of the encountered node will be greater, as the last processed
			//  subtree was on the same level but with a lower pre order (post order completeness
			//  of subtrees).
			postOrderBound = this->subclusters[*it].postOrder + 1;
			this->subclusters[*it].postOrder -= postOrderDelta;
		} else {
			// Encountered arbitrary node.

			// Push onto partial tree node list
			partialTree.push_back(*it);

			// Increase postOrderBound, if necessary
			if (this->subclusters[*it].postOrder >= postOrderBound)
				postOrderBound = this->subclusters[*it].postOrder + 1;

			// Amend postOrder to reflect (absolute) position in partial tree
			this->subclusters[*it].postOrder -= postOrderDelta;
		}
	}

	// Commit currenly open partial tree.
	this->subclusters[partialTree.front()].postOrder = partialTree.size() - 1;
	// TODO move instead of copy
	partialTrees.push_back(partialTree);

	while (!partialTreeStack.empty()) {
		// Restore state
		// TODO move instead of copy
		partialTree = partialTreeStack.back();
		partialTreeStack.pop_back();

		postOrderDelta = postOrderDeltaStack.back();
		postOrderDeltaStack.pop_back();

		this->subclusters[partialTree.front()].postOrder = partialTree.size() - 1;

		// TODO move instead of copy
		partialTrees.push_back(partialTree);
	}


	#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATE
	(GeneratorValidator(*this)).validatePartialTrees(partialTrees);
	#endif

	return partialTrees;
}

void DynamicCommunitiesGenerator::rejoinPartialTrees(std::vector<std::vector<index>> partialTrees) {
	std::vector<index> cluster;
	index postOrderDelta = 0;


	for (auto treeIt = partialTrees.cbegin(); treeIt != partialTrees.cend(); ++treeIt) {
		for (auto nodeIt = (*treeIt).cbegin(); nodeIt != (*treeIt).cend(); ++nodeIt) {
			// If the node *nodeIt is not the root of the partial tree
			if (nodeIt != (*treeIt).cbegin())
			if (nodeIt != (*treeIt).cbegin())
				this->subclusters[*nodeIt].postOrder += postOrderDelta;
			cluster.push_back(*nodeIt);
		}

		// Increase postOrderDelta
		postOrderDelta += (*treeIt).size() - 1;
	}

	for (auto treeIt = partialTrees.crbegin(); treeIt != partialTrees.crend(); ++treeIt) {
		this->subclusters[(*treeIt).front()].postOrder = postOrderDelta;
		++postOrderDelta;

		// Set the parent pointer to the partial tree's root node
		//  - If it is the first partial tree, set to 0. The root node is the cluster's new root
		//  - Otherwise set to the root node of the previous partial tree
		this->subclusters[(*treeIt).front()].parent =
			(treeIt + 1 != partialTrees.crend())? (*(treeIt + 1)).front() : 0;
	}

	if (cluster.size() > 0) {
		// If there are any elements in the cluster...

		// Set the cluster at the index to the cluster created from the partial trees.
		index clusterIndex = this->subclusters[cluster.front()].cluster;
		this->clusters[clusterIndex] = cluster;
	}
}

#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATOR

void GeneratorValidator::validate0Cluster() {
	Aux::enforce(this->state.getCluster(0).size() == 0,
		"cluster[0] has size > 0");
}

void GeneratorValidator::validateClusterNum() {
	if (this->state.getPreGenerated())
		Aux::enforce(this->state.getClusters().size() - 1 == this->state.getParameters().l,
			"clusters has size != l after pre-generation");
}

void GeneratorValidator::validateSubclustersLength() {
	const std::vector<DynamicCommunitiesGenerator::Subcluster> subclusters = this->state.getSubclusters();

	Aux::enforce(subclusters.size() - 1 == this->state.getParameters().affinities.getN(),
		"Invalid length of subclusters");
}

void GeneratorValidator::validateSubclusterCluster() {
	count subclusterNum = 0;
	bool found;

	const std::vector<std::vector<index>> clusters = this->state.getClusters();

	for (auto clustersIt = clusters.cbegin() + 1; clustersIt != clusters.cend(); ++clustersIt) {
		for (auto subclusterIt = (*clustersIt).cbegin(); subclusterIt != (*clustersIt).cend(); ++subclusterIt) {
			Aux::enforce(
				this->state.getSubclusterCluster(*subclusterIt)
					== static_cast<index>(clustersIt - clusters.cbegin()),
				"Found inconsistent state: subcluster not in cluster subcluster.cluster"
			);

			found = false;
			for (auto findIt = (*clustersIt).cbegin(); findIt != subclusterIt; ++findIt) {
				if (*findIt == *subclusterIt)
					found = true;
			}

			Aux::enforce(!found, "Found duplicate subcluster in cluster");

			++subclusterNum;
		}
	}


	Aux::enforce(subclusterNum == this->state.getParameters().affinities.getN(),
		"Invalid number of subclusters found when iterating over clusters");
}

void GeneratorValidator::validateSubclusterParent() {
	index subclusterParent, subclusterCluster;

	for (index i = 1; i <= this->state.getParameters().affinities.getN(); ++i) {
		subclusterParent = this->state.getSubclusterParent(i);
		subclusterCluster = this->state.getSubclusterCluster(i);

		Aux::enforce(subclusterParent == 0
			|| this->state.getSubclusterCluster(subclusterParent) == subclusterCluster,
			"Subcluster parent in different cluster than subcluster");
	}
}

GeneratorValidator::Forest GeneratorValidator::buildSubclusterForest(
		const std::vector<std::vector<index>>& clusters,
		index offset
	) {
	const std::vector<DynamicCommunitiesGenerator::Subcluster>& subclusters
		= this->state.getSubclusters();

	Forest forest = {
		// List of root nodes
		std::vector<index>(),
		// All nodes contained in the forest
		std::vector<TreeNode>(subclusters.size())
	};

	// Manually set nodes[0] to a sentry
	forest.nodes[0] = {
		0,
		0,
		std::vector<index>()
	};

	// Stacks for tree traversal
	std::vector<index> parentStack;
	std::vector<std::vector<index>> childrenStack;

	// Use sentry as bottom element on stacks (collects roots)
	parentStack.push_back(0);
	childrenStack.emplace_back();

	for (auto clusterIt = clusters.cbegin() + offset; clusterIt != clusters.cend(); ++clusterIt) {
		for (auto it = (*clusterIt).cbegin(); it != (*clusterIt).cend(); ++it) {
			while (parentStack.size() > 1
					&& subclusters[parentStack.back()].postOrder < subclusters[*it].postOrder) {
				// Commit subtree rooted by parentStack.back().
				index parentInd = parentStack.back();
				parentStack.pop_back();

				forest.nodes[parentInd] = {
					parentInd,
					parentStack.back(),
					std::move(childrenStack.back())
				};

				childrenStack.pop_back();

				// Add parentInd to its parent's children list
				childrenStack.back().push_back(parentInd);
			}

			// Push *it on stacks
			parentStack.push_back(*it);
			childrenStack.emplace_back();
		}

		// Empty stack
		while (parentStack.size() > 1) {
			// Commit subtree
			index parentInd = parentStack.back();

			parentStack.pop_back();

			forest.nodes[parentInd] = {
				parentInd,
				parentStack.back(),
				std::move(childrenStack.back())
			};

			childrenStack.pop_back();

			// Add parentInd to its parent's children list
			childrenStack.back().push_back(parentInd);
		}

		Aux::enforce(parentStack.size() == 1,
			"Elements remaining on stack after cluster traversal");
		Aux::enforce(parentStack.back() == 0,
			"Wrong element on stack after cluster traversal");
	}

	Aux::enforce(parentStack.size() == 1,
		"Elements remaining on stack after cluster traversal");
	Aux::enforce(parentStack.back() == 0,
		"Wrong element on stack after cluster traversal");

	// Swap the virtual root's (at nodes[0]) children with the roots list
	forest.roots.swap(childrenStack.back());

	return forest;
}

void GeneratorValidator::validateClusters(Forest& forest) {
	count clustersEncountered = 0;

	const std::vector<DynamicCommunitiesGenerator::Subcluster>& subclusters
		= this->state.getSubclusters();

	Aux::enforce(this->state.getClusters().size() - 1 == forest.roots.size(),
		"Wrong number of forest roots discovered");

	std::vector<index> stack;
	index clusterIndex;
	index nodeIndex;

	for (auto it = forest.roots.cbegin(); it != forest.roots.cend(); ++it) {
		stack.push_back(*it);
		clusterIndex = subclusters[*it].cluster;

		++clustersEncountered;

		while (!stack.empty()) {
			nodeIndex = stack.back();
			stack.pop_back();

			Aux::enforce(subclusters[nodeIndex].cluster == clusterIndex,
				"Wrong cluster index of subcluster in (cluster) tree");

			for (auto cIt = forest.nodes[nodeIndex].children.cbegin();
					cIt != forest.nodes[nodeIndex].children.cend(); ++cIt) {
				stack.push_back(*cIt);
			}
		}
	}

	Aux::enforce(this->state.getClusters().size() - 1 == clustersEncountered,
		"Wrong number of clusters encountered");
}

void GeneratorValidator::validatePartialTreesParents(
		const std::vector<std::vector<index>>& partialTrees,
		Forest& forest
	) {
	const std::vector<DynamicCommunitiesGenerator::Subcluster>& subclusters
		= this->state.getSubclusters();

	// Validate parents.
	for (auto treeIt = partialTrees.cbegin(); treeIt != partialTrees.cend(); ++treeIt) {
		if (treeIt + 1 == partialTrees.cend())
			Aux::enforce(subclusters[(*treeIt).at(0)].parent == 0,
				"Wrong subcluster parent (of partial tree root)");
		else
			Aux::enforce(subclusters[(*treeIt).at(0)].parent == (*(treeIt + 1)).at(0),
				"Wrong subcluster parent (of partial tree root)");

		for (auto it = (*treeIt).cbegin() + 1; it != (*treeIt).cend(); ++it) {
			Aux::enforce(subclusters[*it].parent == forest.nodes[*it].parent,
				"Wrong subcluster parent (of partial tree node)");
		}
	}
}

void GeneratorValidator::validateParents(Forest& forest) {
	const std::vector<DynamicCommunitiesGenerator::Subcluster>& subclusters
		= this->state.getSubclusters();

	for (auto it = subclusters.cbegin(); it != subclusters.cend(); ++it) {
		Aux::enforce((*it).parent == forest.nodes[it - subclusters.cbegin()].parent,
			"Wrong subcluster parent");
	}
}

void GeneratorValidator::validateOrder(
		Forest& forest,
		const std::vector<std::vector<index>>& clusters
	) {
	// Calculate preOrder for each subcluster
	index preOrder;
	std::vector<index> preOrders(this->state.getSubclusters().size());
	for (auto clusterIt = clusters.cbegin(); clusterIt != clusters.cend(); ++clusterIt) {
		preOrder = -1;
		for (auto it = (*clusterIt).cbegin(); it != (*clusterIt).cend(); ++it) {
			preOrders[*it] = ++preOrder;
		}
	}

	DFS dfs(forest, this->state.getSubclusters(), preOrders);

	for (auto it = forest.roots.cbegin(); it != forest.roots.cend(); ++it) {
		dfs.startAt(*it);
	}
}

void GeneratorValidator::DFS::startAt(index node) {
	this->postOrder = -1;
	this->preOrder = -1;

	this->process(node);
}

void GeneratorValidator::DFS::process(index node) {
	Aux::enforce(++this->preOrder == this->preOrders[node],
		"Wrong preOrder encountered");

	for (auto cIt = this->forest.nodes[node].children.cbegin();
			cIt != this->forest.nodes[node].children.cend(); ++cIt) {
		this->process(*cIt);
	}

	Aux::enforce(++this->postOrder == this->subclusters[node].postOrder,
		"Wrong postOrder encountered");
}

void GeneratorValidator::validateIndividuals() {
	Aux::enforce(this->state.getIndividuals().size() - 1 == this->state.getParameters().n,
		"Wrong number of individuals");

	const std::vector<DynamicCommunitiesGenerator::Individual>& individuals
		= this->state.getIndividuals();

	for (auto it = individuals.cbegin() + 1; it != individuals.cend(); ++it) {
		Aux::enforce((*it).subcluster != 0,
			"Individual with invalid subcluster 0 encountered");
		Aux::enforce((*it).homeSubcluster != 0,
			"Individual with invalid home subcluster 0 encountered");
	}
}

#endif /* DYNAMICCOMMUNITIESGENERATOR_VALIDATOR */

} /* namespace NetworKit */
