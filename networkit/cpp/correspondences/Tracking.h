/*
 * Tracking.h
 *
 *  Created on: July 8, 2017
 *      Author: Paul Skopnik
 */

#ifndef TRACKING_H_
#define TRACKING_H_

// #define DEBUG_TOPDOWN_BACKTRACK
// #define DEBUG_WRITE_GHT
// #define DEBUG_WRITE_HDT
// #define DEBUG_WRITE_PARTITION

#include <algorithm>
#include <deque>
#include <fstream>
#include <iterator>
#include <ostream>
#include <queue>
#include <vector>

#include "../Globals.h"
#include "../structures/Partition.h"
#include "../generators/DynamicCommunitiesGenerator.h"
#include "Correspondences.h"

namespace NetworKit {

std::vector<count> partitionSubsetSizes(const Partition& p);

class GHGraph {
public:
	struct Edge {
		index a, b;
		double weight;
	};

	class EdgeIterator {
	// class EdgeIterator : std::iterator<std::random_access_iterator_tag, Edge> {
	public:
		EdgeIterator() = default;
		EdgeIterator(const GHGraph& graph, index node) : graph(&graph), node(node) {
			this->initialise();
		}
		EdgeIterator(const GHGraph& graph, index node, bool) : graph(&graph), node(node) {
			this->state = after;
		}

		inline Edge operator*() {
			switch (this->state) {
				case unintialised:
					return Edge();
				case parent:
					return {
						this->node,
						this->getParent(this->node),
						this->getWeight(this->node)
					};
				case children:
					return {
						this->node,
						*this->iterator,
						this->getWeight(*this->iterator)
					};
				default:
					return Edge();
			}
		}

		// inline Edge operator->() {
		// 	return **this;
		// }

		inline EdgeIterator& operator++() {
			const GHGraph::NodeData& node = this->graph->nodes[this->node];

			switch (this->state) {
				case parent:
					if (node.totalChildren == 0) {
						this->state = after;
						break;
					}

					this->state = children;
					this->iterator = this->graph->children.cbegin() + node.childrenIndex;

					break;
				case children:
					++this->iterator;

					if (this->iterator
							== this->graph->children.cbegin()
								+ node.childrenIndex + node.totalChildren)
					{
						this->state = after;
					}

					break;
				case after:
				default:
					break;
			}

			return *this;
		}

		inline EdgeIterator operator++(int) {
			EdgeIterator tmp(*this);
			++*this;
			return tmp;
		}

		friend inline bool operator==(const EdgeIterator& lhs, const EdgeIterator& rhs) {
			if (lhs.state != rhs.state)
				return false;

			if (lhs.state == unintialised)
				return true;

			if (lhs.graph != rhs.graph)
				return false;

			if (lhs.node != rhs.node)
				return false;

			switch (lhs.state) {
				case parent:
					return true;
				case children:
					return lhs.iterator == rhs.iterator;
				case after:
					return true;
				default:
					return false;
			}
		}

		friend inline bool operator!=(const EdgeIterator& lhs, const EdgeIterator& rhs) {
			return !(lhs == rhs);
		}

	protected:
		enum State
		{
			unintialised,
			parent,
			children,
			after
		};

		State state = unintialised;

		const GHGraph* graph;
		index node;
		std::vector<index>::const_iterator iterator;

		inline void initialise() {
			if (this->hasParent(this->node)) {
				this->state = parent;
			} else if (this->hasChildren(this->node)) {
				this->state = children;
				const GHGraph::NodeData& node = this->graph->nodes[this->node];
				this->iterator = this->graph->children.cbegin() + node.childrenIndex;
			} else {
				this->state = after;
			}
		}

		inline index getParent(index node) {
			return this->graph->parents[node];
		}

		inline bool hasParent(index node) {
			return this->graph->parents[node] < this->graph->size;
		}

		inline bool hasChildren(index node) {
			return this->graph->nodes[node].totalChildren > 0;
		}

		inline double getWeight(index node) {
			return this->graph->weights[node];
		}
	};

	class Neighbours {
	public:
		Neighbours(const GHGraph& graph, index node) : graph(graph), node(node) {}
		~Neighbours() = default;

		inline EdgeIterator begin() const {
			return EdgeIterator(this->graph, this->node);
		}
		inline EdgeIterator end() const {
			return EdgeIterator(this->graph, this->node, false);
		}
		inline EdgeIterator cbegin() const {
			return EdgeIterator(this->graph, this->node);
		}
		inline EdgeIterator cend() const {
			return EdgeIterator(this->graph, this->node, false);
		}
	protected:
		const GHGraph& graph;
		index node;
	};

	/**
	 * Constructs an empty GHGraph.
	 *
	 * The static build() function is needed to properly build the GHGraph structure.
	 */
	GHGraph(const std::vector<index>& parents, const std::vector<count>& weights)
		: size(0), parents(parents), weights(weights) {}
	~GHGraph() = default;

	static GHGraph build(const std::vector<index>& parents, const std::vector<count>& weights);

	Neighbours neighbours(index node) const;
	inline count getSize() const {
		return this->size;
	}
protected:
	struct NodeData {
		count totalChildren, foundChildren;
		count subtreeSize;
		index childrenIndex;
	};

	const count size;
	const std::vector<index>& parents;
	const std::vector<count>& weights;
	const std::vector<NodeData> nodes;
	const std::vector<index> children;

	GHGraph(
		const std::vector<index>& parents,
		const std::vector<count>& weights,
		std::vector<NodeData> nodes,
		std::vector<index> children
	);
};

/**
 * Tracking is the interface for timestep snapshot based trackers of dynamic networks.
 */
class Tracking {
public:
	Tracking() = default;
	virtual ~Tracking() = default;

	/**
	 * Add adds a partition to the tracker.
	 * Tracking expects it to be called once for each timestep, in order.
	 */
	virtual void add(const Partition& partition) = 0;
};

template<class Data, class StepAnalyser, class Outputer>
class StepByStep : public Tracking {
public:
	StepByStep(Data data = Data(), StepAnalyser analyser = StepAnalyser(), Outputer outputer = Outputer())
		: data(data), analyser(analyser), outputer(outputer) {
			this->outputer.setData(this->data);
		}

	virtual ~StepByStep() = default;

	virtual void add(const Partition& partition) {
		if (first) {
			this->data.addTimestep(
				this->analyser.analyseFirst(partition)
			);
			this->first = false;
		} else {
			this->data.addTimestep(
				this->analyser.analyseStep(this->lastPartition, partition)
			);
		}

		// Copy partition over so it can be used in the next add() call.
		this->lastPartition = partition;
	}

	virtual Data& getData() {
		return this->data;
	}

	virtual Outputer& getOutputer() {
		return this->outputer;
	}

protected:
	Partition lastPartition;
	bool first = true;

	Data data;
	StepAnalyser analyser;
	Outputer outputer;
};

template<class Data, class Analyser, class Outputer>
class Simple : public Tracking {
public:
	Simple(Data data = Data(), Analyser analyser = Analyser(), Outputer outputer = Outputer())
		: data(data), analyser(analyser), outputer(outputer) {
			this->analyser.setData(this->data);
			this->outputer.setData(this->data);
		}
	virtual ~Simple() = default;

	virtual void add(const Partition& partition) {
		this->analyser->add(partition);
	}

protected:
	Data data;
	Analyser analyser;
	Outputer outputer;
};

class TimestepData {
public:
	/**
	 * Represents a correspondence between the partition of a timestep t and t+1.
	 * p comprises the involved parts' indices for t, pPrime those for t+1.
	 * pSizes contains for all parts in p the size of the union of the part with the union of
	 * pPrime, pPrimeSizes contains the size of the union of the part with the union of p for all
	 * parts in pPrime. The indices correspond to the indices of parts in p and pPrime.
	 * pToPPrimeSizes contains the number of elements in the intersection of each part in p and
	 * each part in pPrime: pToPPrimeSizes[i][j] = p[i] \cap pPrime[j].
	 * intersectionSize is the number of elements in the intersection of the unions of p and
	 * pPrime, while unionSize is the number of elements in the union of the unions of p and
	 * pPrime.
	 */
	struct Correspondence {
		std::vector<index> p;
		std::vector<index> pPrime;

		std::vector<count> pSizes;
		std::vector<count> pPrimeSizes;

		std::vector<std::vector<count>> pToPPrimeSizes;

		count intersectionSize;
		count unionSize;
	};

	/**
	 * Timestep comprises all required information about a timestep.
	 * When the Timestep object describes timestep t, the contained correspondences are between
	 * timestep t-1 and t.
	 */
	struct Timestep {
		std::vector<count> partitionSizes;
		std::vector<Correspondence> correspondences;
	};

	inline void addTimestep(Timestep timestep) {
		this->timesteps.push_back(timestep);
	}

	inline const std::vector<Timestep>& getTimesteps() const {
		return this->timesteps;
	}

	inline std::vector<Timestep>::iterator begin() {
		return this->timesteps.begin();
	}

	inline std::vector<Timestep>::iterator end() {
		return this->timesteps.end();
	}

	inline std::vector<Timestep>::const_iterator cbegin() const {
		return this->timesteps.cbegin();
	}

	inline std::vector<Timestep>::const_iterator cend() const {
		return this->timesteps.cend();
	}
protected:
	std::vector<Timestep> timesteps;
};

class DCGTimestepData {
public:
	/**
	 * Represents a correspondence between the partition of a timestep t and t+1.
	 * p comprises the involved parts' indices for t, pPrime those for t+1.
	 * pSizes contains for all parts in p the size of the union of the part with the union of
	 * pPrime, pPrimeSizes contains the size of the union of the part with the union of p for all
	 * parts in pPrime. The indices correspond to the indices of parts in p and pPrime.
	 * pToPPrimeSizes contains the number of elements in the intersection of each part in p and
	 * each part in pPrime: pToPPrimeSizes[i][j] = p[i] \cap pPrime[j].
	 * intersectionSize is the number of elements in the intersection of the unions of p and
	 * pPrime, while unionSize is the number of elements in the union of the unions of p and
	 * pPrime.
	 */
	typedef TimestepData::Correspondence Correspondence;
	// struct Correspondence {
	// 	std::vector<index> p;
	// 	std::vector<index> pPrime;

	// 	std::vector<count> pSizes;
	// 	std::vector<count> pPrimeSizes;

	// std::vector<std::vector<count>> pToPPrimeSizes;

	// 	count intersectionSize;
	// 	count unionSize;
	// };

	/**
	 * Comprises ownership information about a cluster.
	 * largestPart is the index of the cluster holding the largest stake of the cluster, while
	 * stake expresses the relative size of that stake (fraction of cluster owned).
	 * ownershipMargin expresses the minium fraction of the cluster the owning part must lose to
	 * be replaced as the largest part. That is: Half the difference between the sizes (in terms
	 * of individuals) of the largest and second largest parts divided by the cluster size.
	 */
	struct Ownership {
		index largestPart;
		double stake;
		double ownershipMargin;
	};

	/**
	 * Timestep comprises all required information about a timestep.
	 * When the Timestep object describes timestep t, the contained correspondences are between
	 * timestep t-1 and t.
	 */
	struct Timestep {
		std::vector<count> partitionSizes;
		std::vector<Ownership> ownership;
		std::vector<Correspondence> correspondences;
	};

	class OwnershipAccessor {
	public:
		static count getNoOfParts(const DCGTimestepData& data) {
			return data.getNoOfParts();
		}

		static count owningPart(const Timestep& timestep, index part) {
			return timestep.ownership[part].largestPart;
		}

		static double ownershipMargin(const Timestep& timestep, index part) {
			return timestep.ownership[part].ownershipMargin;
		}
	};

	DCGTimestepData(const std::vector<std::vector<index>>& parts) : noOfParts(parts.size()) {}
	~DCGTimestepData() = default;

	inline void addTimestep(Timestep timestep) {
		this->timesteps.push_back(timestep);
	}

	inline count getNoOfParts() const {
		return this->noOfParts;
	}

	inline const std::vector<Timestep>& getTimesteps() const {
		return this->timesteps;
	}

	inline std::vector<Timestep>::iterator begin() {
		return this->timesteps.begin();
	}

	inline std::vector<Timestep>::iterator end() {
		return this->timesteps.end();
	}

	inline std::vector<Timestep>::const_iterator cbegin() const {
		return this->timesteps.cbegin();
	}

	inline std::vector<Timestep>::const_iterator cend() const {
		return this->timesteps.cend();
	}
protected:
	count noOfParts;
	std::vector<Timestep> timesteps;
};

class CorrespondencesExtractor {
public:
	CorrespondencesExtractor(Correspondences& correspondences);
	~CorrespondencesExtractor() = default;

	TimestepData::Correspondence extract(const std::vector<index>& parts) const;
protected:
	Correspondences& correspondences;
};

class LeafExpansion {
public:
	LeafExpansion() = default;
	virtual ~LeafExpansion() = default;

	virtual TimestepData::Timestep analyseFirst(const Partition& partition) const;
	virtual TimestepData::Timestep analyseStep(
		const Partition& partition1,
		const Partition& partition2
	) const;

protected:
	std::vector<count> calculateChildrenNo(const std::vector<index>& parents) const;
};

class DCGOwnershipExtractor {
public:
	DCGOwnershipExtractor(
		const DynamicCommunitiesGenerator& generator,
		const std::vector<std::vector<index>>& parts
	);
	~DCGOwnershipExtractor() = default;

	std::vector<DCGTimestepData::Ownership> extract(const Partition& partition) const;
	std::vector<DCGTimestepData::Ownership> extract(
		const Partition& partition,
		const std::vector<count>& subsetSizes
	) const;

protected:
	const GeneratorState generator;
	const std::vector<std::vector<index>>& parts;
	std::vector<index> subclusterToPart;

	std::vector<index> invertParts(const std::vector<std::vector<index>>& parts);
};

class DCGLeafExpansion {
public:
	DCGLeafExpansion(
		const DynamicCommunitiesGenerator& generator,
		const std::vector<std::vector<index>>& parts
	);
	virtual ~DCGLeafExpansion() = default;

	virtual DCGTimestepData::Timestep analyseFirst(const Partition& partition) const;
	virtual DCGTimestepData::Timestep analyseStep(
		const Partition& partition1,
		const Partition& partition2
	) const;
protected:
	DCGOwnershipExtractor ownershipExtractor;
	LeafExpansion leafExpansion;
};

/**
 * Namespace for HierarchicalTree data types.
 */
namespace HT {
	struct Result {
		std::vector<index> p;
		std::vector<index> pPrime;
		count cost;

		bool operator<(const HT::Result& rhs) const { return this->cost < rhs.cost; }
		bool operator>(const HT::Result& rhs) const { return this->cost > rhs.cost; }
	};

	template <typename QualityType>
	struct ResultSet {
		QualityType quality;
		std::vector<Result> results;

		bool operator<(const HT::ResultSet<QualityType>& rhs) const {
			return this->quality < rhs.quality;
		}
		bool operator>(const HT::ResultSet<QualityType>& rhs) const {
			return this->quality > rhs.quality;
		}
	};

	template <typename QualityType>
	struct Node {
		ResultSet<QualityType> resultSet;
		index parent;
		std::vector<index> children;
		count evaluatedChildren;
	};

	struct default_constructor {};
	struct correspondences_constructor {};

	template <class Objective>
	struct ObjectiveTraits {
			typedef typename Objective::QualityType QualityType;
			typedef typename Objective::Constructor Constructor;
	};

	template <class Objective>
	Objective instantiateObjective(Correspondences& c) {
		return instantiateObjective<Objective>(
			c, typename ObjectiveTraits<Objective>::Constructor()
		);
	}

	template <class Objective>
	Objective instantiateObjective(Correspondences& c, default_constructor) {
		return Objective();
	}

	template <class Objective>
	Objective instantiateObjective(Correspondences& c, correspondences_constructor) {
		return Objective(c);
	}
};

template<class Objective>
class TopDown {
public:
	typedef HT::Result Result;
	typedef HT::ResultSet<typename Objective::QualityType> ResultSet;
	typedef HT::Node<typename Objective::QualityType> Node;

	TopDown(Correspondences& c)
	: c(c),
		objective(HT::instantiateObjective<Objective>(c)),
		parents(c.gomoryHuParent),
		weights(c.cutWithGomoryHuParent),
		graph(GHGraph::build(c.gomoryHuParent, c.cutWithGomoryHuParent)),
		partSets(c.gomoryHuParent.size(), 0),
		resultNodes(2 * c.gomoryHuParent.size() - 1)
	{}
	~TopDown() = default;

	std::vector<Result> run() {
		#ifdef DEBUG_WRITE_HDT
		// DEBUG
		// Count timesteps
		static count timestepIndex = 1;
		++timestepIndex;

		std::ofstream file(std::to_string(timestepIndex).append("_hdt.dot").c_str());
		file << "digraph{";
		#endif /* DEBUG_WRITE_HDT */


		// Sort gomory hu tree indices by ascending weight
		std::vector<index> weightIndices = this->indexSortWeight();

		// Create root node of result set tree
		Node& rootNode = this->createRootNode();

		#ifdef DEBUG_WRITE_HDT
		// DEBUG
		file << "0[shape=box,label=\"c = 0\"];";
		#endif /* DEBUG_WRITE_HDT */

		for (auto it = weightIndices.cbegin(); it != weightIndices.cend() - 1; ++it) {
			// Iterate over weights in ascending order, ignore root's edge (non-existent)

			// Fetch current set the edge is in
			index parentSetIndex = this->partSets[*it];
			Node& parent = this->resultNodes[parentSetIndex];
			// Check if expansion stopped
			if (this->objective.expansionStopped(parent.resultSet))
				continue;

			// Assign parts a and b of edge (a -> b) new set indices
			this->partSets[*it] = ++this->maxSetIndex;
			this->partSets[this->parents[*it]] = ++this->maxSetIndex;

			Node& nodeA = this->propagateNode(*it, parentSetIndex);
			Node& nodeB = this->propagateNode(this->parents[*it], parentSetIndex);

			#ifdef DEBUG_WRITE_HDT
			// DEBUG
			file << parentSetIndex << "->" << this->partSets[*it] << ";";
			file << parentSetIndex << "->" << this->partSets[this->parents[*it]] << ";";
			if (nodeA.resultSet.results.front().p.size() == 1)
				file << this->partSets[*it] << "[label=\"c = " << nodeA.resultSet.results.front().cost << "\"]";
			else
				file << this->partSets[*it] << "[shape=box,label=\"c = " << nodeA.resultSet.results.front().cost << "\"]";
			if (nodeB.resultSet.results.front().p.size() == 1)
				file << this->partSets[this->parents[*it]] << "[label=\"c = " << nodeB.resultSet.results.front().cost << "\"]";
			else
				file << this->partSets[this->parents[*it]] << "[shape=box,label=\"c = " << nodeB.resultSet.results.front().cost << "\"]";
			#endif /* DEBUG_WRITE_HDT */
		}

		this->backtrack();

		#ifdef DEBUG_WRITE_HDT
		// DEBUG
		file << "}";
		#endif /* DEBUG_WRITE_HDT */

		return rootNode.resultSet.results;
	}

protected:
	Correspondences& c;
	Objective objective;

	const std::vector<index>& parents;
	const std::vector<count>& weights;

	GHGraph graph;

	std::vector<index> partSets;
	index maxSetIndex = 0;

	std::vector<Node> resultNodes;

	std::deque<index> backtrackQueue;

	std::vector<index> indexSortWeight() const {
		std::vector<index> weightIndices(this->graph.getSize(), 0);

		std::iota(weightIndices.begin(), weightIndices.end(), 0);

		std::sort(weightIndices.begin(), weightIndices.end(),
			[this](int i, int j) {
				return this->weights[i] < this->weights[j];
			}
		);

		return weightIndices;
	}

	Node& createRootNode() {
		// Create vector containing indices of all gomory hu tree parts
		std::vector<index> p(this->graph.getSize());
		std::iota(p.begin(), p.end(), 0);

		Node& rootNode = this->resultNodes[0];

		rootNode.resultSet.results.emplace_back();
		Result& rootResult = rootNode.resultSet.results.back();

		// Set parts of root node's result
		rootResult.p.swap(p);

		if (this->objective.calculatePPrime(rootResult))
			// If pPrime should be calculated for the result, calculate it
			this->calculatePPrime(rootResult);

		if (this->objective.isEligible(rootResult))
			// If the result is eligible to be a result set, calculate its quality
			rootNode.resultSet.quality = this->objective.calculateQuality(rootResult);

		if (rootResult.p.size() == 1 || this->objective.expansionStopped(rootNode.resultSet))
			// If only one node in result or expansion has stopped add the node to the
			//  backtrack queue
			this->backtrackQueue.push_back(0);

		// Set parent to the size of the resultNodes vector to mark root
		rootNode.parent = this->resultNodes.size();

		return rootNode;
	}

	Node& propagateNode(index l, index parentSetIndex) {
		index setIndex = this->partSets[l];

		// Create new node for set of l
		Node& node = this->resultNodes[setIndex];

		// Create new result for node of set of l
		node.resultSet.results.emplace_back();
		Result& result = node.resultSet.results.back();

		// Run bfsExpand on l
		result.cost = this->bfsExpand(l, parentSetIndex, result.p);

		if (this->objective.calculatePPrime(result))
			// If pPrime should be calculated for the result, calculate it
			this->calculatePPrime(result);

		if (this->objective.isEligible(result))
			// If the result is eligible to be a result set, calculate its quality
			node.resultSet.quality = this->objective.calculateQuality(result);

		if (result.p.size() == 1 || this->objective.expansionStopped(node.resultSet))
			// If only one node in result or expansion has stopped add the node to the
			//  backtrack queue
			this->backtrackQueue.push_back(setIndex);

		// Add to parent
		node.parent = parentSetIndex;
		this->resultNodes[parentSetIndex].children.push_back(setIndex);

		return node;
	}

	void backtrack() {
		#ifdef DEBUG_TOPDOWN_BACKTRACK
		std::cout << "backtrack()" << std::endl;
		#endif /* DEBUG_TOPDOWN_BACKTRACK */
		while (!this->backtrackQueue.empty()) {
			Node& node = this->resultNodes[this->backtrackQueue.front()];
			this->backtrackQueue.pop_front();

			if (node.children.size() == 0) {
				// Use own ResultSet
				#ifdef DEBUG_TOPDOWN_BACKTRACK
				std::cout << "    Leaf, using own" << std::endl;
				std::cout << "     " << node.resultSet.results.front().p.size() << " : " << node.resultSet.results.front().pPrime.size() << std::endl;
				#endif /* DEBUG_TOPDOWN_BACKTRACK */
			} else {
				ResultSet combined = {
					this->objective.combineQuality(
						this->resultNodes[node.children[0]].resultSet,
						this->resultNodes[node.children[1]].resultSet
					),
					std::vector<Result>(this->resultNodes[node.children[0]].resultSet.results)
				};
				combined.results.insert(
					combined.results.end(),
					// Breaking from c++11 to c++14:
					// combined.results.cend(),
					this->resultNodes[node.children[1]].resultSet.results.cbegin(),
					this->resultNodes[node.children[1]].resultSet.results.cend()
				);

				if (this->objective.isEligible(node.resultSet.results.front())
					&& this->objective.replaceCombined(combined, node.resultSet))
				{
					// Replace combined ResultSets of children by own ResultSet
					#ifdef DEBUG_TOPDOWN_BACKTRACK
					std::cout << "    Replaced by own" << std::endl;
					std::cout << "     " << node.resultSet.results.front().p.size() << " : " << node.resultSet.results.front().pPrime.size() << std::endl;
					#endif /* DEBUG_TOPDOWN_BACKTRACK */
				} else {
					#ifdef DEBUG_TOPDOWN_BACKTRACK
					std::cout << "    Combining children" << std::endl;
					std::cout << "     |.| = " << combined.results.size() << std::endl;
					#endif /* DEBUG_TOPDOWN_BACKTRACK */
					node.resultSet = std::move(combined);
				}
			}

			if (node.parent == this->resultNodes.size())
				// Skip for root node
				continue;

			Node& parent = this->resultNodes[node.parent];
			if (++parent.evaluatedChildren == parent.children.size())
				this->backtrackQueue.push_back(node.parent);
		}

		#ifdef DEBUG_TOPDOWN_BACKTRACK
		std::cout << "/" << std::endl;
		#endif /* DEBUG_TOPDOWN_BACKTRACK */
	}

	count bfsExpand(index from, index expandInto, std::vector<index>& nodes) {
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

	void calculatePPrime(Result& result) {
		for (index i = 0; i < this->c.cardPartition2; ++i) {
			count sum = 0;

			for (index part : result.p) {
				sum += this->c.distributions[part][i];
			}

			if (2 * sum > this->c.cardinalityOfCluster2[i])
				result.pPrime.push_back(i);
		}
	}

};

template <class Algorithm>
class HierarchicalTree {
public:
	HierarchicalTree() = default;
	virtual ~HierarchicalTree() = default;

	virtual TimestepData::Timestep analyseFirst(const Partition& partition) const {
		return {
			partitionSubsetSizes(partition),
			std::vector<TimestepData::Correspondence>()
		};
	}

	virtual TimestepData::Timestep analyseStep(
		const Partition& partition1,
		const Partition& partition2
	) const {
		TimestepData::Timestep timestep;
		timestep.partitionSizes = partitionSubsetSizes(partition2);

		Correspondences c;
		c.detect(2, partition1, partition2);

		#ifdef DEBUG_WRITE_GHT
		// DEBUG
		// Count timesteps
		static count timestepIndex = 1;
		++timestepIndex;

		// DEBUG
		// Write out gomory hu tree
		std::ofstream file(std::to_string(timestepIndex).append("_ght.dot").c_str());
		file << "digraph" << "{";
		for (auto it = c.gomoryHuParent.cbegin(); it != c.gomoryHuParent.cend(); ++it) {
			if (*it < c.gomoryHuParent.size()) {
				file << "p" << "_" << it - c.gomoryHuParent.cbegin();
				file << "->" << "p" << "_" << *it;
				file << "[" << "label=" << c.cutWithGomoryHuParent[it - c.gomoryHuParent.cbegin()] << "]";
				file << ";";
			} else {
				file << "p" << "_" << it - c.gomoryHuParent.cbegin() << ";";
			}
		}
		file << "}";
		file.close();
		#endif /* DEBUG_WRITE_GHT */


		Algorithm algorithm(c);
		for (auto result : algorithm.run()) {
			timestep.correspondences.push_back(
				this->extractCorrespondence(result.p, result.pPrime, c)
			);
		}

		return timestep;
	}

protected:
	TimestepData::Correspondence extractCorrespondence(
		std::vector<index>& p,
		std::vector<index>& pPrime,
		Correspondences& c
	) const {

		std::vector<index> pSizes(p.size(), 0);
		std::vector<index> pPrimeSizes(pPrime.size(), 0);

		std::vector<std::vector<count>> pToPPrimeSizes(p.size(),
			std::vector<count>(pPrime.size()));

		count intersection = 0;
		count sizePPrime = 0;
		count sizeP = 0;

		for (auto pPrimeIt = pPrime.cbegin(); pPrimeIt != pPrime.cend(); ++pPrimeIt) {
			for (auto pIt = p.cbegin(); pIt != p.cend(); ++pIt) {
				intersection += c.distributions[*pIt][*pPrimeIt];

				pSizes[pIt - p.cbegin()] += c.distributions[*pIt][*pPrimeIt];
				pPrimeSizes[pPrimeIt - pPrime.cbegin()] += c.distributions[*pIt][*pPrimeIt];

				pToPPrimeSizes[pIt - p.cbegin()][pPrimeIt - pPrime.cbegin()] =
					c.distributions[*pIt][*pPrimeIt];
			}

			sizePPrime += c.cardinalityOfCluster2[*pPrimeIt];
		}

		for (index pPart : p) {
			sizeP += c.cardinalityOfCluster1[pPart];
		}

		return {
			p,
			pPrime,
			pSizes,
			pPrimeSizes,
			pToPPrimeSizes,
			intersection,
			sizeP + sizePPrime - intersection
		};
	}
};

template <class Algorithm>
class DCGHierarchicalTree {
public:
	DCGHierarchicalTree(
		const DynamicCommunitiesGenerator& generator,
		const std::vector<std::vector<index>>& parts
	) : ownershipExtractor(generator, parts) {
	}
	virtual ~DCGHierarchicalTree() = default;

	virtual DCGTimestepData::Timestep analyseFirst(const Partition& partition) const {
		TimestepData::Timestep simpleTimestep = this->hierarchicalTree.analyseFirst(partition);

		return {
			simpleTimestep.partitionSizes,
			this->ownershipExtractor.extract(partition),
			simpleTimestep.correspondences
		};
	}

	virtual DCGTimestepData::Timestep analyseStep(
		const Partition& partition1,
		const Partition& partition2
	) const {
		TimestepData::Timestep simpleTimestep = this->hierarchicalTree.analyseStep(
			partition1, partition2);

		return {
			simpleTimestep.partitionSizes,
			this->ownershipExtractor.extract(partition2),
			simpleTimestep.correspondences
		};
	}
protected:
	DCGOwnershipExtractor ownershipExtractor;
	HierarchicalTree<Algorithm> hierarchicalTree;
};


/**
 * SystemOptimum minimises the overall cost (sum) of the correspondences while ensuring that
 * every chosen correspondence is a 1:x, 2:x, y:1 or y:2 correspondence.
 *
 * In case there are only (!) edges with zero weight in the gomory hu tree, two random parts will
 * end up together in a correspondence. The SystemOptimumZero class mitigates this case, putting
 * the two parts in separate correspondences.
 *
 * The reason for this effect is, that further parts are always added to the correspondence as long
 * as the overall cost does not increase. The cost can only not increase, if the weight of the edge
 * is zero.
 */
class SystemOptimum {
public:
	typedef HT::default_constructor Constructor;
	typedef count QualityType;

	QualityType calculateQuality(HT::Result result) const {
		return result.cost;
	}

	QualityType combineQuality(HT::ResultSet<QualityType> a, HT::ResultSet<QualityType> b) const {
		return a.quality + b.quality;
	}

	bool expansionStopped(HT::ResultSet<QualityType> r) const {
		return false;
	}

	bool isEligible(HT::Result r) const {
		return r.p.size() <= 2 || r.pPrime.size() <= 2;
	}

	bool calculatePPrime(HT::Result r) const {
		return true;
	}

	virtual bool replaceCombined(HT::ResultSet<QualityType> combined, HT::ResultSet<QualityType> parent) const {
		return parent.quality <= combined.quality;
	}
};

/**
 * SystemOptimumZero minimises the overall cost (sum) of the correspondences while mitigating an
 * issue with gomory hu trees containing only zero weight edges.
 *
 * See the SystemOptimum class for more details.
 */
class SystemOptimumZero : public SystemOptimum {
public:
	typedef HT::default_constructor Constructor;
	typedef count QualityType;

	bool replaceCombined(HT::ResultSet<QualityType> combined, HT::ResultSet<QualityType> parent) const override {
		return this->SystemOptimum::replaceCombined(combined, parent)
			&& combined.quality != 0;

		// Alternative implementation:
		// return parent.quality < combined.quality;
	}
};

class IndividualOptimumStrict {
public:
	typedef HT::default_constructor Constructor;
	typedef count QualityType;

	QualityType calculateQuality(HT::Result result) const {
		return result.cost;
	}

	QualityType combineQuality(HT::ResultSet<QualityType> a, HT::ResultSet<QualityType> b) const {
		return std::min(a.quality, b.quality);
	}

	bool expansionStopped(HT::ResultSet<QualityType> r) const {
		return false;
	}

	bool isEligible(HT::Result r) const {
		return r.p.size() <= 2 || r.pPrime.size() <= 2;
	}

	bool calculatePPrime(HT::Result r) const {
		return true;
	}

	bool replaceCombined(HT::ResultSet<QualityType> combined, HT::ResultSet<QualityType> parent) const {
		return parent.quality < combined.quality;
	}
};

class IndividualOptimumWeakened {
public:
	typedef HT::default_constructor Constructor;
	typedef count QualityType;

	QualityType calculateQuality(HT::Result result) const {
		return result.cost;
	}

	QualityType combineQuality(HT::ResultSet<QualityType> a, HT::ResultSet<QualityType> b) const {
		return std::min(a.quality, b.quality);
	}

	bool expansionStopped(HT::ResultSet<QualityType> r) const {
		return false;
	}

	bool isEligible(HT::Result r) const {
		return r.p.size() <= 2 || r.pPrime.size() <= 2;
	}

	bool calculatePPrime(HT::Result r) const {
		return true;
	}

	virtual bool replaceCombined(HT::ResultSet<QualityType> combined, HT::ResultSet<QualityType> parent) const {
		return parent.quality <= combined.quality;
	}
};

class IndividualOptimumWeakenedZero : public IndividualOptimumWeakened {
public:
	typedef HT::default_constructor Constructor;
	typedef count QualityType;

	bool replaceCombined(HT::ResultSet<QualityType> combined, HT::ResultSet<QualityType> parent) const override {
		return this->IndividualOptimumWeakened::replaceCombined(combined, parent)
			&& combined.quality != 0;
	}
};

/**
 * CheapestSetsGenerator returns all sets contained in the decomposition hierarchy of the gomory hu
 * tree in ascending cost order. The root result (the entire tree) is omitted.
 *
 * It is implemented as a forward iterator.
 */
class CheapestSetsGenerator {
public:
	CheapestSetsGenerator(Correspondences& c);
	CheapestSetsGenerator(Correspondences& c, std::vector<index> partSets, count rootSet);
	/**
	 * Constructs an ended CheapestSetsGenerator.
	 */
	CheapestSetsGenerator(Correspondences& c, bool);
	~CheapestSetsGenerator() = default;

	inline HT::Result& operator*() {
		if (!this->reachedEnd && this->resultQueue.empty())
			this->advance();

		return this->refResultSet(this->resultQueue.top());
	}

	inline HT::Result* operator->() {
		return &**this;
	}

	inline CheapestSetsGenerator& operator++() {
		if (this->ended())
			return *this;

		if (this->resultQueue.empty())
			this->advance();

		if (!this->resultQueue.empty())
			this->resultQueue.pop();

		if (this->resultQueue.empty())
			this->advance();

		return *this;
	}


	inline bool ended() const {
		return this->reachedEnd && this->resultQueue.empty();
	}

	inline CheapestSetsGenerator operator++(int) {
		CheapestSetsGenerator tmp(*this);
		++*this;
		return tmp;
	}

	friend inline bool operator==(const CheapestSetsGenerator& lhs, const CheapestSetsGenerator& rhs) {
		if (lhs.ended() && rhs.ended())
			return true;
		else
			return false;
	}

	friend inline bool operator!=(const CheapestSetsGenerator& lhs, const CheapestSetsGenerator& rhs) {
		return !(lhs == rhs);
	}
protected:
	bool reachedEnd = false;

	Correspondences& c;

	const std::vector<index>& parents;
	const std::vector<count>& weights;

	std::vector<index> weightIndices;
	std::vector<index>::const_iterator weightIt;
	std::vector<index>::const_iterator weightItEnd;

	std::vector<index> partSets;
	index minSetIndex;

	std::vector<HT::Result> sets;
	std::priority_queue<index> resultQueue;

	GHGraph graph;

	std::vector<index> indexSortWeight() const;

	HT::Result& createRootResult(index rootSet);

	/**
	 * Advance advances the generator, i.e. explores the tree until the next cheapest result is
	 * available or the tree is exhausted.
	 *
	 * That means when advance returns either reachedEnd is true or resultQueue is non-empty.
	 */
	void advance();

	HT::Result& propagateResult(index l, index parentSetIndex);

	count bfsExpand(index from, index expandInto, std::vector<index>& nodes);

	/**
	 * Returns a reference to the result from this->set with the passed in set index.
	 *
	 * Enlarges this->sets if necessary.
	 */
	HT::Result& refResultSet(int set);
};

/**
 * SmallestMutual is a correspondences exploration algorithm extracting the cheapest smallest
 * mutual correspondences from the gomory hu tree.
 */
class SmallestMutual {
public:
	SmallestMutual(Correspondences& c);
	~SmallestMutual() = default;

	std::vector<HT::Result> run();
protected:
	Correspondences& c;

	std::vector<index> partSets;
	index maxSetIndex = 0;

	HT::Result createRootResult() const;

	HT::ResultSet<count> exploreTree(index setIndex, const HT::Result& self, bool isSelfMutual);

	HT::Result buildInverseResult(
		const HT::Result& other,
		const std::vector<index>& superSet
	);

	/**
	 * Calculates pPrime from p.
	 */
	void calculatePPrime(HT::Result& result) const;
	/**
	 * Calculates pPrime from p; includes all "maybe" (zero-impact) optimal partners in pPrime.
	 *
	 * Useful if as many parts as possible should be used.
	 */
	void calculatePPrime(HT::Result& result, int) const;

	bool isMutual(const HT::Result& result) const;
};

/**
 * RecursiveMutual is a correspondences exploration algorithm extracting the smallest mutual
 * correspondences from the gomory hu tree.
 */
class RecursiveMutual {
public:
	RecursiveMutual(Correspondences& c);
	~RecursiveMutual() = default;

	std::vector<HT::Result> run();
protected:
	Correspondences& c;

	const std::vector<index>& parents;
	const std::vector<count>& weights;

	GHGraph graph;

	std::vector<index> partSets;
	index maxSetIndex = 0;

	std::vector<index> weightIndices;

	std::vector<index> indexSortWeight() const;

	HT::Result createRootResult() const;

	count bfsExpand(index from, index expandInto, std::vector<index>& nodes);

	HT::ResultSet<count> processTree(index setIndex, const HT::Result result);

	HT::Result buildResult(index l, index parentSetIndex);
	HT::Result buildInverseResult(
		index l,
		index parentSetIndex,
		const HT::Result& other,
		const std::vector<index>& superSet
	);

	void calculatePPrime(HT::Result& result) const;
	void calculateInversePPrime(
		HT::Result& result,
		const std::vector<index>& superSet,
		const std::vector<index>& other
	) const;
	bool isMutual(const HT::Result& result) const;
};


namespace Strings {
	const std::string space = " ";
	const std::string arrow = "->";
	const std::string equals = "=";
	const std::string colon = ":";
	const std::string semicolon = ";";
	const std::string comma = ",";
	const std::string underscore = "_";
	const std::string quote = "\"";

	const std::string opCurly = "{";
	const std::string clCurly = "}";
	const std::string opBracket = "[";
	const std::string clBracket = "]";

	const std::string digraph = "digraph";
	const std::string rankSame = "rank=same;";
	const std::string rankdirLR = "rankdir=LR;";
	const std::string shapePlaintext = "shape=plaintext";
	const std::string shapeSquare = "shape=square";
	const std::string fixedsizeTrue = "fixedsize=true";
	const std::string styleFilled = "style=filled";
	const std::string styleFilledSolid = "style=\"filled,solid\"";

	const std::string labelEquals = "label=";
	const std::string heightEquals = "height=";
	const std::string widthEquals = "width=";
	const std::string colorEquals = "color=";
	const std::string fillcolorEquals = "fillcolor=";
	const std::string penwidthEquals = "penwidth=";

	const std::string p = "p";
	const std::string c = "c";

	const std::string nodes = "nodes";
	const std::string links = "links";
	const std::string value = "value";
	const std::string rank = "rank";
	const std::string spec = "spec";
	const std::string source = "source";
	const std::string target = "target";
}

template <class Data>
class DotOutput {
public:
	virtual void aggregate(std::ostream& stream) {
		stream << Strings::digraph << Strings::opCurly << Strings::rankdirLR;

		this->writeTimesteps(stream);

		this->writeNodes(stream);

		this->writeEdges(stream);

		stream << Strings::clCurly;
	}

	virtual void aggregateToFile(const char* filename) {
		std::ofstream file(filename);

		this->aggregate(file);

		file.close();
	}

	virtual inline void setData(Data& data) {
		this->data = &data;
	}

protected:
	Data* data = NULL;

	virtual void writeTimesteps(std::ostream& stream) const {
		// Write timestep order
		stream << Strings::opCurly;

		for (count i = 0; i < this->data->getTimesteps().size(); ++i) {
			if (i != 0)
				stream << Strings::arrow;
			stream << i + 1;
		}
		stream << Strings::semicolon;

		stream << Strings::clCurly;
	}

	virtual void writeNodes(std::ostream& stream) const {
		// Write nodes
		for (auto ts = this->data->cbegin(); ts != this->data->cend(); ++ts) {
			// Enable subgraph clusters with the next line
			// stream << "subgraph cluster" << ts - this->data->cbegin() + 1 << Strings::opCurly;
			stream << Strings::opCurly;

			stream << Strings::rankSame;

			// Shapes are drawn with their size (:= square area) proportional to the size of the cluster
			//  A cluster of mean size is drawn with width=0.75 and height=0.5 (the default dot uses for
			//  nodes). The width / height proportion is kept constant for all nodes at 0.75/0.5 = 3/2.
			count sum = std::accumulate(ts->partitionSizes.cbegin(), ts->partitionSizes.cend(), 0);
			double mean = static_cast<double>(sum) / static_cast<double>(ts->partitionSizes.size());
			double coeff = 1.0 / (4.0 * mean);

			index timestep = ts - this->data->cbegin() + 1;
			stream << timestep << Strings::opBracket << Strings::shapePlaintext << Strings::clBracket << Strings::semicolon;

			for (auto it = ts->partitionSizes.cbegin(); it != ts->partitionSizes.cend(); ++it) {
				stream << Strings::p << Strings::underscore << timestep << Strings::underscore;
				stream << it - ts->partitionSizes.cbegin();

				double height = std::sqrt(coeff * (*it));

				stream << Strings::opBracket;

				this->writeNodeAttributes(stream, *ts, it - ts->partitionSizes.cbegin(), height);

				stream << Strings::clBracket;

				stream << Strings::semicolon;
			}

			stream << Strings::clCurly;
		}
	}

	virtual void writeNodeAttributes(
		std::ostream& stream,
		const typename Data::Timestep& timestep,
		index part,
		double height
	) const {
		stream << Strings::shapeSquare << Strings::comma;
		stream << Strings::labelEquals << part << Strings::comma;
		stream << Strings::widthEquals << 1.25 * height << Strings::comma;
		stream << Strings::heightEquals << height;
	}

	virtual void writeEdges(std::ostream& stream) const {
		// Write edges
		for (auto ts = this->data->cbegin() + 1; ts != this->data->cend(); ++ts) {
			index timestep = ts - this->data->cbegin() + 1;

			for (auto corr = ts->correspondences.cbegin(); corr != ts->correspondences.cend(); ++corr) {
				if (corr->p.size() == 0 || corr->pPrime.size() == 0)
					continue;

				stream << Strings::opCurly;

				for (auto pIt = corr->p.cbegin(); pIt != corr->p.cend(); ++pIt) {
					stream << Strings::p << Strings::underscore << timestep - 1 << Strings::underscore << *pIt << Strings::semicolon;
				}

				stream << Strings::clCurly << Strings::arrow << Strings::opCurly;

				for (auto pPrimeIt = corr->pPrime.cbegin(); pPrimeIt != corr->pPrime.cend(); ++pPrimeIt) {
					stream << Strings::p << Strings::underscore << timestep << Strings::underscore << *pPrimeIt << Strings::semicolon;
				}

				stream << Strings::clCurly;

				double jaccard = static_cast<double>(corr->intersectionSize) / static_cast<double>(corr->unionSize);
				stream << Strings::opBracket << Strings::penwidthEquals << jaccard * 10 << Strings::clBracket;
			}
		}
	}
};

template <class Data, class OwnershipAccessor = typename Data::OwnershipAccessor>
class OwnershipDotOutput : public DotOutput<Data> {
public:
	virtual void aggregate(std::ostream& stream) override {
		// Propagate the palette, that is a base hue value for each part
		this->palette = std::vector<int16_t>(OwnershipAccessor::getNoOfParts(*this->data));
		// Use 360 degree hue spectrum with wrap-around to find colours
		double increment = 360.0 / static_cast<double>(OwnershipAccessor::getNoOfParts(*this->data));

		double runningHue = 0.0;
		for (auto it = this->palette.begin(); it != this->palette.end(); ++it) {
			*it = static_cast<int16_t>(runningHue);
			runningHue += increment;
		}

		this->DotOutput<Data>::aggregate(stream);
	}
protected:
	std::vector<int16_t> palette;

	virtual void writeNodeAttributes(
		std::ostream& stream,
		const typename Data::Timestep& timestep,
		index part,
		double height
	) const override {
		this->DotOutput<Data>::writeNodeAttributes(stream, timestep, part, height);
		stream << Strings::comma;

		stream << Strings::styleFilled << Strings::comma;

		stream << Strings::fillcolorEquals << Strings::quote;
		stream << static_cast<double>(this->palette[OwnershipAccessor::owningPart(timestep, part)]) / 360.0;
		stream << Strings::comma;
		stream << 0.75 * 2 * OwnershipAccessor::ownershipMargin(timestep, part) + 0.1;
		stream << Strings::comma << 1.0;
		stream << Strings::quote;
	}
};

template <class Data>
class JSONOutput {
public:
	virtual void aggregate(std::ostream& stream) {
		stream << Strings::opCurly;

		stream << Strings::quote << Strings::nodes << Strings::quote << Strings::colon;
		this->writeNodes(stream);

		stream << Strings::comma;

		stream << Strings::quote << Strings::links << Strings::quote << Strings::colon;
		this->writeLinks(stream);

		stream << Strings::clCurly;
	}

	virtual void aggregateToFile(const char* filename) {
		std::ofstream file(filename);

		this->aggregate(file);

		file.close();
	}

	virtual void setData(Data& data) {
		this->data = &data;
	}

protected:
	Data* data = NULL;

	virtual void writeNodes(std::ostream& stream) const {
		stream << Strings::opBracket;

		bool first = true;

		// Write nodes
		for (auto ts = this->data->cbegin(); ts != this->data->cend(); ++ts) {
			index timestep = ts - this->data->cbegin() + 1;

			for (auto it = ts->partitionSizes.cbegin(); it != ts->partitionSizes.cend(); ++it) {
				if (first)
					first = false;
				else
					stream << Strings::comma;

				stream << Strings::opCurly;

				this->writePartNodeAttributes(stream, *ts, timestep, it -ts->partitionSizes.cbegin(), *it);

				stream << Strings::clCurly;
			}

			for (auto it = ts->correspondences.cbegin(); it != ts->correspondences.cend(); ++it) {
				if (first)
					first = false;
				else
					stream << Strings::comma;

				stream << Strings::opCurly;

				this->writeCorrespondenceNodeAttributes(stream, *ts, timestep,
					it - ts->correspondences.cbegin());

				stream << Strings::clCurly;
			}
		}

		stream << Strings::clBracket;
	}

	virtual void writePartNodeAttributes(
		std::ostream& stream,
		const typename Data::Timestep& timestep,
		index timestepIndex,
		index part,
		count size
	) const {
		stream << Strings::quote << Strings::spec << Strings::quote << Strings::colon;
		stream << Strings::quote << Strings::p << Strings::underscore << timestepIndex;
		stream << Strings::underscore << part << Strings::quote;

		stream << Strings::comma;

		stream << Strings::quote << Strings::rank << Strings::quote << Strings::colon;
		stream << 2 * (timestepIndex - 1);

		stream << Strings::comma;

		stream << Strings::quote << Strings::value << Strings::quote << Strings::colon << size;
	}

	virtual void writeCorrespondenceNodeAttributes(
		std::ostream& stream,
		const typename Data::Timestep& timestep,
		index timestepIndex,
		index corresIndex
	) const {
		stream << Strings::quote << Strings::spec << Strings::quote << Strings::colon;
		stream << Strings::quote << Strings::c << Strings::underscore << timestepIndex;
		stream << Strings::underscore << corresIndex << Strings::quote;

		stream << Strings::comma;

		stream << Strings::quote << Strings::rank << Strings::quote << Strings::colon;
		stream << 2 * (timestepIndex - 1) - 1;
	}

	virtual void writeLinks(std::ostream& stream) const {
		stream << Strings::opBracket;

		bool first = true;

		for (auto ts = this->data->cbegin(); ts != this->data->cend(); ++ts) {
			index timestepIndex = ts - this->data->cbegin() + 1;

			for (auto it = ts->correspondences.cbegin(); it != ts->correspondences.cend(); ++it) {
				index corresIndex = it - ts->correspondences.cbegin();

				for (auto pIt = it->p.cbegin(); pIt != it->p.cend(); ++pIt) {
					if (first)
						first = false;
					else
						stream << Strings::comma;

					stream << Strings::opCurly;

					this->writePLinkAttributes(stream, *ts, timestepIndex, *pIt, corresIndex,
						it->pSizes[pIt - it->p.cbegin()]);

					stream << Strings::clCurly;
				}

				for (auto pPrimeIt = it->pPrime.cbegin(); pPrimeIt != it->pPrime.cend(); ++pPrimeIt) {
					if (first)
						first = false;
					else
						stream << Strings::comma;

					stream << Strings::opCurly;

					this->writePPrimeLinkAttributes(stream, *ts, timestepIndex, corresIndex, *pPrimeIt,
						it->pPrimeSizes[pPrimeIt - it->pPrime.cbegin()]);

					stream << Strings::clCurly;
				}
			}
		}

		stream << Strings::clBracket;
	}

	virtual void writePLinkAttributes(
		std::ostream& stream,
		const typename Data::Timestep& timestep,
		index timestepIndex,
		index part,
		index corresIndex,
		count value
	) const {
		stream << Strings::quote << Strings::source << Strings::quote << Strings::colon;
		stream << Strings::quote << Strings::p << Strings::underscore << timestepIndex - 1;
		stream << Strings::underscore << part << Strings::quote;

		stream << Strings::comma;

		stream << Strings::quote << Strings::target << Strings::quote << Strings::colon;
		stream << Strings::quote << Strings::c << Strings::underscore << timestepIndex;
		stream << Strings::underscore << corresIndex << Strings::quote;

		stream << Strings::comma;

		stream << Strings::quote << Strings::value << Strings::quote << Strings::colon;
		stream << value;
	}

	virtual void writePPrimeLinkAttributes(
		std::ostream& stream,
		const typename Data::Timestep& timestep,
		index timestepIndex,
		index corresIndex,
		index part,
		count value
	) const {
		stream << Strings::quote << Strings::source << Strings::quote << Strings::colon;
		stream << Strings::quote << Strings::c << Strings::underscore << timestepIndex;
		stream << Strings::underscore << corresIndex << Strings::quote;

		stream << Strings::comma;

		stream << Strings::quote << Strings::target << Strings::quote << Strings::colon;
		stream << Strings::quote << Strings::p << Strings::underscore << timestepIndex;
		stream << Strings::underscore << part << Strings::quote;

		stream << Strings::comma;

		stream << Strings::quote << Strings::value << Strings::quote << Strings::colon;
		stream << value;
	}
};

template <class Data, class OwnershipAccessor = typename Data::OwnershipAccessor>
class OwnershipJSONOutput : public JSONOutput<Data> {
protected:
	virtual void writeNodes(std::ostream& stream) const override {
		this->JSONOutput<Data>::writeNodes(stream);
		stream << Strings::comma;

		stream << Strings::quote << "clusters" << Strings::quote << Strings::colon;
		stream << Strings::opCurly;

		stream << Strings::quote << "number" << Strings::quote << Strings::colon;
		stream << OwnershipAccessor::getNoOfParts(*this->data);

		stream << Strings::clCurly;
	}

	virtual void writePartNodeAttributes(
		std::ostream& stream,
		const typename Data::Timestep& timestep,
		index timestepIndex,
		index part,
		count size
	) const override {
		this->JSONOutput<Data>::writePartNodeAttributes(
			stream, timestep, timestepIndex, part, size);
		stream << Strings::comma;

		stream << Strings::quote << "cluster" << Strings::quote << Strings::colon;
		stream << OwnershipAccessor::owningPart(timestep, part);
	}
};

template <class Data, class OwnershipAccessor = typename Data::OwnershipAccessor>
class EdgeCutEvaluation {
public:
	EdgeCutEvaluation() = default;
	virtual ~EdgeCutEvaluation() = default;

	virtual count evaluate() {
		if (data == NULL)
			return 0;

		count noOfParts = OwnershipAccessor::getNoOfParts(*this->data);

		this->partToPartCutSum.assign(noOfParts, std::vector<count>(noOfParts, 0));

		this->cutSum = 0;
		this->edgeSum = 0;
		this->cuts = 0;
		this->edges = 0;

		for (auto ts = this->data->cbegin(); ts != this->data->cend(); ++ts) {
			++this->timesteps;

			for (auto it = ts->correspondences.cbegin(); it != ts->correspondences.cend(); ++it) {
				for (auto pIt = it->p.cbegin(); pIt != it->p.cend(); ++pIt) {
					for (auto pPrimeIt = it->pPrime.cbegin(); pPrimeIt != it->pPrime.cend(); ++pPrimeIt) {
						if (OwnershipAccessor::owningPart(*(ts - 1), *pIt) != OwnershipAccessor::owningPart(*ts, *pPrimeIt)) {
							this->cutSum += it->pToPPrimeSizes[pIt - it->p.cbegin()][pPrimeIt - it->pPrime.cbegin()];
							++this->cuts;

							this->partToPartCutSum[OwnershipAccessor::owningPart(*(ts - 1), *pIt)][OwnershipAccessor::owningPart(*ts, *pPrimeIt)]
							 	+= it->pToPPrimeSizes[pIt - it->p.cbegin()][pPrimeIt - it->pPrime.cbegin()];
						}

						this->edgeSum += it->pToPPrimeSizes[pIt - it->p.cbegin()][pPrimeIt - it->pPrime.cbegin()];
						++this->edges;
					}
				}
			}
		}

		this->evaluated = true;

		return this->cutSum;
	}

	virtual inline void setData(Data& data) {
		this->data = &data;
		this->evaluated = false;
	}

	virtual inline count getCutSum() {
		return this->cutSum;
	}
	virtual inline count getEdgeSum() {
		return this->edgeSum;
	}
	virtual inline count getCuts() {
		return this->cuts;
	}
	virtual inline count getEdges() {
		return this->edges;
	}
	virtual inline count getTimesteps() {
		return this->timesteps;
	}
	virtual inline std::vector<std::vector<count>>& getPartToPartCutSum() {
		return this->partToPartCutSum;
	}
	virtual inline count getPartToPartCutSum(index partA, index partB) {
		return this->partToPartCutSum[partA][partB];
	}
protected:
	Data* data = NULL;

	bool evaluated = false;

	count cutSum;
	count edgeSum;

	count cuts;
	count edges;

	count timesteps;

	std::vector<std::vector<count>> partToPartCutSum;
};


} /* namespace NetworKit */

#endif /* TRACKING_H_ */