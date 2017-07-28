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

/**
 * Polyfill for initializer-list constructor for std::set introduced in C++14.
 */
template<typename T>
inline std::set<T> set(std::initializer_list<T> l) {
	std::set<T> s = l;
	return s;
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

/**
 * TODO
 * To be removed
 */
std::vector<count> calculateChildrenNo(const std::vector<index>& parents) {
	std::vector<count> childrenNo(parents.size(), 0);

	for (auto it = parents.cbegin(); it != parents.cend(); ++it) {
		if (*it < parents.size())
			++childrenNo[*it];
	}

	return childrenNo;
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

namespace Strings {
	const std::string space = " ";
	const std::string arrow = "->";
	const std::string equals = "=";
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
		partition.subsetSizes(),
		std::vector<TimestepData::Correspondence>()
	};
}

TimestepData::Timestep LeafExpansion::analyseStep(
	const Partition& partition1,
	const Partition& partition2
) const {
	TimestepData::Timestep timestep;
	timestep.partitionSizes = partition1.subsetSizes();

	std::vector<TimestepData::Correspondence>& correspondences = timestep.correspondences;

	Correspondences c;
	c.detect(2, partition1, partition2);

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
					correspondences.push_back(this->extractCorrespondence(c, parts));
				} else {
					// The "cut" is cheaper when excluding the parent
					//  -> only include the leaf in the correspondence
					correspondences.push_back(this->extractCorrespondence(c, std::vector<index>(1, it - subtreeSizes.cbegin())));
				}
			} else {
				// Parent node has more than one child -> ignore parent, only include leaf in the correspondence
				correspondences.push_back(this->extractCorrespondence(c, std::vector<index>(1, it - subtreeSizes.cbegin())));
			}
		} else if (*it == partition2.upperBound()) {
			// Root node
			std::vector<count> childrenNo = calculateChildrenNo(c.gomoryHuParent);
			if (childrenNo[it - subtreeSizes.cbegin()] == 0) {
				// The root node is the only node (only one partition)
				//  -> search for correspondence with this one node
				correspondences.push_back(this->extractCorrespondence(c, std::vector<index>(1, it - subtreeSizes.cbegin())));
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
						correspondences.push_back(this->extractCorrespondence(c, parts));
					} else {
						// The "cut" is cheaper excluding the child
						//  -> only include root node in correspondence
						correspondences.push_back(this->extractCorrespondence(c, std::vector<index>(1, it - subtreeSizes.cbegin())));
					}
				} else {
					// Child has no more than one child, there is no single "cut"
					//  -> only include root node in correspondence
					correspondences.push_back(this->extractCorrespondence(c, std::vector<index>(1, it - subtreeSizes.cbegin())));
				}
			// } else {
			// 	// Ignore as there is no single "cut" using the root
			}
		}
	}

	return timestep;
}

TimestepData::Correspondence LeafExpansion::extractCorrespondence(Correspondences& c, const std::vector<index>& parts) const {
	TimestepData::Correspondence correspondence = {parts, std::vector<index>(0)};

	for (index i = 0; i < c.cardPartition2; ++i) {
		count sum = 0;

		for (index part : parts) {
			sum += c.distributions.at(part).at(i);
			// sum += c.distributions[part][i];
		}

		if (2 * sum > c.cardinalityOfCluster2[i])
			correspondence.pPrime.push_back(i);
	}

	return correspondence;
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
	return this->extract(partition, partition.subsetSizes());
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

		auto partIt = std::max_element(partCounters.cbegin(), partCounters.cend());
		ownership[i] = {
			static_cast<index>(partIt - partCounters.cbegin()),
			static_cast<double>(*partIt) / subsetSizes[i]
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


void DotOutputer::aggregate(std::ostream& stream) {
	stream << Strings::digraph << Strings::opCurly << Strings::rankdirLR;

	this->writeTimesteps(stream);

	this->writeNodes(stream);

	this->writeEdges(stream);

	stream << Strings::clCurly;
}

void DotOutputer::setData(TimestepData& data) {
	this->data = &data;
}

void DotOutputer::aggregateToFile(const char* filename) {
	std::ofstream file(filename);

	this->aggregate(file);

	file.close();
}

void DotOutputer::writeTimesteps(std::ostream& stream) const {
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

void DotOutputer::writeNodes(std::ostream& stream) const {
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

		count timestep = ts - this->data->cbegin() + 1;
		stream << timestep << Strings::opBracket << Strings::shapePlaintext << Strings::clBracket << Strings::semicolon;

		for (auto it = ts->partitionSizes.cbegin(); it != ts->partitionSizes.cend(); ++it) {
			stream << "p" << Strings::underscore << timestep << Strings::underscore << it - ts->partitionSizes.cbegin();

			double height = std::sqrt(coeff * (*it));
			stream << Strings::opBracket << Strings::shapeSquare << Strings::comma;
			stream << Strings::labelEquals << it - ts->partitionSizes.cbegin() << Strings::comma;
			stream << Strings::widthEquals << 1.25 * height << Strings::comma;
			stream << Strings::heightEquals << height << Strings::clBracket;

			stream << Strings::semicolon;
		}

		stream << Strings::clCurly;
	}
}

void DotOutputer::writeEdges(std::ostream& stream) const {
	// Write edges
	for (auto ts = this->data->cbegin() + 1; ts != this->data->cend(); ++ts) {
		count timestep = ts - this->data->cbegin() + 1;

		for (auto corr = ts->correspondences.cbegin(); corr != ts->correspondences.cend(); ++corr) {
			stream << Strings::opCurly;

			for (auto pIt = corr->p.cbegin(); pIt != corr->p.cend(); ++pIt) {
				stream << "p" << Strings::underscore << timestep - 1 << Strings::underscore << *pIt << Strings::semicolon;
			}

			stream << Strings::clCurly << Strings::arrow << Strings::opCurly;

			for (auto pPrimeIt = corr->pPrime.cbegin(); pPrimeIt != corr->pPrime.cend(); ++pPrimeIt) {
				stream << "p" << Strings::underscore << timestep << Strings::underscore << *pPrimeIt << Strings::semicolon;
			}

			stream << Strings::clCurly;
		}
	}
}

void DCGDotOutputer::aggregate(std::ostream& stream) {
	// Propagate the palette, that is a base hue value for each part
	this->palette = std::vector<int16_t>(this->data->getNoOfParts());
	// Use 360 degree hue spectrum with wrap-around to find colours
	double increment = 360.0 / static_cast<double>(this->data->getNoOfParts());

	double runningHue = 0.0;
	for (auto it = this->palette.begin(); it != this->palette.end(); ++it) {
		*it = static_cast<int16_t>(runningHue);
		runningHue += increment;
	}

	this->DotOutputer::aggregate(stream);
}

void DCGDotOutputer::setData(DCGTimestepData& data) {
	this->data = &data;
}

void DCGDotOutputer::writeTimesteps(std::ostream& stream) const {
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

void DCGDotOutputer::writeNodes(std::ostream& stream) const {
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

		count timestep = ts - this->data->cbegin() + 1;
		stream << timestep << Strings::opBracket << Strings::shapePlaintext << Strings::clBracket << Strings::semicolon;

		for (auto it = ts->partitionSizes.cbegin(); it != ts->partitionSizes.cend(); ++it) {
			stream << "p" << Strings::underscore << timestep << Strings::underscore << it - ts->partitionSizes.cbegin();

			double height = std::sqrt(coeff * (*it));
			stream << Strings::opBracket << Strings::shapeSquare << Strings::comma;
			stream << Strings::labelEquals << it - ts->partitionSizes.cbegin() << Strings::comma;
			stream << Strings::widthEquals << 1.25 * height << Strings::comma;
			stream << Strings::heightEquals << height << Strings::comma;
			stream << Strings::fillcolorEquals << Strings::quote;
			stream << static_cast<double>(this->palette[ts->ownership[it - ts->partitionSizes.cbegin()].largestPart]) / 360.0 << Strings::comma;
			stream << 0.75 * 2 * (ts->ownership[it - ts->partitionSizes.cbegin()].stake - 0.5) + 0.1 << Strings::comma << 1.0;
			stream << Strings::quote << Strings::comma;
			stream << Strings::styleFilled << Strings::clBracket;

			stream << Strings::semicolon;
		}

		stream << Strings::clCurly;
	}
}

void DCGDotOutputer::writeEdges(std::ostream& stream) const {
	// Write edges
	for (auto ts = this->data->cbegin() + 1; ts != this->data->cend(); ++ts) {
		count timestep = ts - this->data->cbegin() + 1;

		for (auto corr = ts->correspondences.cbegin(); corr != ts->correspondences.cend(); ++corr) {
			stream << Strings::opCurly;

			for (auto pIt = corr->p.cbegin(); pIt != corr->p.cend(); ++pIt) {
				stream << "p" << Strings::underscore << timestep - 1 << Strings::underscore << *pIt << Strings::semicolon;
			}

			stream << Strings::clCurly << Strings::arrow << Strings::opCurly;

			for (auto pPrimeIt = corr->pPrime.cbegin(); pPrimeIt != corr->pPrime.cend(); ++pPrimeIt) {
				stream << "p" << Strings::underscore << timestep << Strings::underscore << *pPrimeIt << Strings::semicolon;
			}

			stream << Strings::clCurly;
		}
	}
}

} /* namespace NetworKit */
