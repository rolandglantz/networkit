/*
 * Tracking.h
 *
 *  Created on: July 8, 2017
 *      Author: Paul Skopnik
 */

#ifndef TRACKING_H_
#define TRACKING_H_

#include <iterator>
#include <ostream>
#include <vector>

#include "../Globals.h"
#include "../structures/Partition.h"
#include "../generators/DynamicCommunitiesGenerator.h"
#include "Correspondences.h"

namespace NetworKit {

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
		EdgeIterator(const GHGraph& graph, index node, bool b) : graph(&graph), node(node) {
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
								+ node.childrenIndex + node.totalChildren) {
						this->state = after;
					}

					break;
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

	static GHGraph build(const std::vector<index>& parents, const std::vector<count>& weights);
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
	Simple() = default;
	~Simple() = default;

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
	 */
	struct Correspondence {
		std::vector<index> p;
		std::vector<index> pPrime;
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
	 */
	typedef TimestepData::Correspondence Correspondence;
	// struct Correspondence {
	// 	std::vector<index> p;
	// 	std::vector<index> pPrime;
	// };

	/**
	 * Comprises ownership information about a cluster.
	 * largestPart is the index of the cluster holding the largest stake of the cluster,
	 * while stake expresses the relative size of that stake (fraction of cluster owned).
	 */
	struct Ownership {
		index largestPart;
		double stake;
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
	TimestepData::Correspondence extractCorrespondence(
		Correspondences& c,
		const std::vector<index>& parts
	) const;
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


class DotOutputer {
public:
	virtual void aggregate(std::ostream& stream);

	virtual void aggregateToFile(const char* filename);

	virtual void setData(TimestepData& data);

protected:
	TimestepData* data = NULL;

	virtual void writeTimesteps(std::ostream& stream) const;

	virtual void writeNodes(std::ostream& stream) const;

	virtual void writeEdges(std::ostream& stream) const;
};

class DCGDotOutputer : public DotOutputer {
public:
	virtual void aggregate(std::ostream& stream);

	virtual void setData(DCGTimestepData& data);

protected:
	DCGTimestepData* data = NULL;
	std::vector<int16_t> palette;

	virtual void writeTimesteps(std::ostream& stream) const;

	virtual void writeNodes(std::ostream& stream) const;

	virtual void writeEdges(std::ostream& stream) const;
};

template<class Data>
class DotOutput {
public:
	virtual void aggregateToFile(const char* filename);

	virtual void aggregate(std::ostream& stream);

	virtual inline void setData(Data& data) {
		this->data = &data;
	}
protected:
	Data* data = NULL;

	virtual void writeTimesteps(std::ostream& stream) const;

	virtual void writeNodes(std::ostream& stream) const;

	virtual void writeEdges(std::ostream& stream) const;
};

template <>
class DotOutput<DCGTimestepData> {
};

} /* namespace NetworKit */

#endif /* TRACKING_H_ */
