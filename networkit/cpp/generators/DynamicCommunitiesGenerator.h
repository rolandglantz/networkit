/*
 * DynamicCommunitiesGenerator.h
 *
 *  Created on: May 26, 2017
 *      Author: Paul Skopnik
 */

#ifndef DYNAMICCOMMUNITIESGENERATOR_H_
#define DYNAMICCOMMUNITIESGENERATOR_H_

#define DYNAMICCOMMUNITIESGENERATOR_VALIDATOR
// #define DYNAMICCOMMUNITIESGENERATOR_VALIDATE

// Enable one of the following definitions to activate splits and merges

// #define DYNAMICCOMMUNITIESGENERATOR_MERGE_ORIG
// #define DYNAMICCOMMUNITIESGENERATOR_MERGE_1STEP // Not Implemented
#define DYNAMICCOMMUNITIESGENERATOR_MERGE_2STEP

// #define DYNAMICCOMMUNITIESGENERATOR_SPLIT_ORIG
// #define DYNAMICCOMMUNITIESGENERATOR_SPLIT_LARGEST
#define	DYNAMICCOMMUNITIESGENERATOR_SPLIT_3STEP
// #define	DYNAMICCOMMUNITIESGENERATOR_SPLIT_2STEP // Not Implemented
// #define	DYNAMICCOMMUNITIESGENERATOR_SPLIT_1STEP // Not Implemented

#define DOUBLE_MARGIN_OF_ERROR 1.0 / 32.0

#include <vector>
#include <stdexcept>

#include "../Globals.h"
#include "../structures/Partition.h"

namespace NetworKit {

template<typename T>
class Matrix {
public:
	typedef typename std::vector<T>::iterator iterator;
	typedef typename std::vector<T>::const_iterator const_iterator;
	struct Index {
		index i, j;
	};

	Matrix() = default;
	virtual ~Matrix() = default;

	T& operator[](const Index ind) {
		return this->at(ind);
	}

	virtual T& at(const Index ind) = 0;

	virtual T get(index i, index j) const = 0;
	virtual void set(index i, index j, T value) = 0;

	virtual count getN() const = 0;

	virtual iterator begin() = 0;
	virtual iterator end() = 0;

	virtual const_iterator cbegin() const = 0;
	virtual const_iterator cend() const = 0;

	virtual Index indexFromIterator(const_iterator it) const = 0;
	virtual Index indexFromPosition(index bufIndex) const = 0;
};

template<typename T>
class SymmetricMatrix : public Matrix<T> {
public:
	typedef typename Matrix<T>::iterator iterator;
	typedef typename Matrix<T>::const_iterator const_iterator;
	typedef typename Matrix<T>::Index Index;

	SymmetricMatrix() : n(0), buf(0) {};
	SymmetricMatrix(count n) : n(n), buf((n - 1)*n / 2) {}
	SymmetricMatrix(count n, const T& value) : n(n), buf((n - 1)*n / 2, value) {};
	~SymmetricMatrix() = default;

	virtual T& at(const Index ind) {
		index bufIndex = this->calcBufIndex(ind.i, ind.j);
		return this->buf[bufIndex];
	}

	virtual T get(index i, index j) const {
		index bufIndex = this->calcBufIndex(i, j);
		return this->buf[bufIndex];
	}

	virtual void set(index i, index j, T value) {
		index bufIndex = this->calcBufIndex(i, j);
		this->buf[bufIndex] = value;
	}

	virtual count getN() const {
		return this->n;
	}

	virtual iterator begin() {
		return this->buf.begin();
	}
	virtual iterator end() {
		return this->buf.end();
	}

	virtual const_iterator cbegin() const {
		return this->buf.cbegin();
	}
	virtual const_iterator cend() const {
		return this->buf.cend();
	}

	virtual Index indexFromIterator(const_iterator it) const {
		if (this->cbegin() > it || this->cend() < it) {
			throw(std::domain_error(
				"SymmetricMatrix::indexFromIterator: Iterator seems to be invalid"
			));
		}

		return this->indexFromPosition(it - this->cbegin());
	}

	virtual Index indexFromPosition(index bufIndex) const {
		SymmetricMatrix<T>::Index ind = {1,0};

		//     1 2 3 4
		//   + - - - -
		// 1 | x 0 1 2
		// 2 | / x 3 4
		// 3 | / / x 5
		// 4 | / / / x

		while (bufIndex >= this->n - ind.i) {
			// Not in the right row yet:
			// Remove one row, increment i.
			bufIndex -= this->n - ind.i;
			++ind.i;
		}
		// In the right row now:
		// Calculate j component of index.
		ind.j = bufIndex + ind.i + 1;

		return ind;
	}

protected:
	count n;
	std::vector<T> buf;

	virtual index calcBufIndex(index i, index j) const {
		if (i < 1 || j < 1)
			throw(std::domain_error(
				"SymmetricMatrix: Invalid index (<1)"
			));
		if (i == j)
			throw(std::domain_error(
				"SymmetricMatrix: Accessing invalid index: i == j"
			));
		if (i > n || j > n)
			throw(std::domain_error(
				"SymmetricMatrix: Invalid index (>n)"
			));

		if (j < i) {
			index tmp = j;
			j = i;
			i = tmp;
		}

		// Calculate index:
		//   (i - 1) * this->n : Start of row i index in regular matrix
		// + (j - 1)           : Index of column j in row i in regular matrix
		// - (i * (i + 1) / 2) : Correction due to the efficient nature of buf
		//                        (also removes the + 1 from the column)

		return (i - 1) * this->n + (j - 1) - (i * (i + 1) / 2);
	};
};

template<typename T>
class DiagonalSymmetricMatrix : public SymmetricMatrix<T> {
public:
	typedef typename SymmetricMatrix<T>::iterator iterator;
	typedef typename SymmetricMatrix<T>::const_iterator const_iterator;
	typedef typename SymmetricMatrix<T>::Index Index;

	DiagonalSymmetricMatrix() {};
	DiagonalSymmetricMatrix(count n) : SymmetricMatrix<T>(n + 1) {}
	DiagonalSymmetricMatrix(count n, const T& value) : SymmetricMatrix<T>(n + 1, value) {};
	~DiagonalSymmetricMatrix() = default;

	virtual count getN() const {
		return this->SymmetricMatrix<T>::getN() - 1;
	}

	virtual Index indexFromPosition(index bufIndex) const {
		Index ind = this->SymmetricMatrix<T>::indexFromPosition(bufIndex);
		ind.j -= 1;
		return ind;
	}

protected:
	virtual index calcBufIndex(index i, index j) const {
		if (j < i) {
			index tmp = j;
			j = i;
			i = tmp;
		}

		return this->SymmetricMatrix<T>::calcBufIndex(i, j + 1);
	}
};

class SubclusterEdges {
public:
	typedef SymmetricMatrix<int8_t>::Index Index;

	SubclusterEdges(SymmetricMatrix<double> affinities);
	~SubclusterEdges() = default;

	count getN() const;

	void selectEdge(index i, index j);
	void unselectEdge(index i, index j);

	double getAvailableEdgesWeight() const;
	double getSelectedEdgesWeight() const;

	count getAvailableEdgesNumber() const;
	count getSelectedEdgesNumber() const;

	Index availableEdgeAtDistribution(double affinity) const;
protected:
	SymmetricMatrix<double> affinities;
	SymmetricMatrix<int8_t> availableEdges;

	double availableEdgesWeight, selectedEdgesWeight;
	count availableEdgesNumber, selectedEdgesNumber;
};

class AffinitiesGenerator {
public:
	typedef SymmetricMatrix<double> Affinities;

	AffinitiesGenerator() = default;
	~AffinitiesGenerator() = default;

	Affinities halfHalf(count k, double w) const;
	Affinities halfHalf(count k, double w, std::vector<std::vector<index>>& parts) const;
};

class DynamicCommunitiesGenerator {
	friend class GeneratorState;

public:
	struct Parameters {
		SymmetricMatrix<double> affinities;
		count l;

		count n;
		double alpha_v;
	};

	struct Individual {
		index subcluster;
		index homeSubcluster;
	};

	struct Subcluster {
		index parent;
		index cluster;

		// TODO
		//  Switch to a size approach, significantly simplifying the implementation
		// count size;
		// index preOrder;
		index postOrder;
	};

	DynamicCommunitiesGenerator(const Parameters parameters);
	~DynamicCommunitiesGenerator() = default;

	Partition next();

protected:
	const Parameters parameters;

	index timestep;
	bool preGenerated;

	std::vector<Individual> individuals;

	std::vector<std::vector<index>> clusters;

	SubclusterEdges subclusterEdges;
	std::vector<Subcluster> subclusters;

	void preGenerate();

	void performSplit();
	void performMerge();

	void addSubclusterEdge(index i, index j);
	void removeSubclusterEdge(index i, index j);

	void makeClusterRoot(index i);
	std::vector<std::vector<index>> extractPartialTrees(index i);
	void rejoinPartialTrees(std::vector<std::vector<index>> partialTrees);

	void performIndividualMoves();
};

/**
 * GeneratorState provides safe access to a Generator's internal state.
 */
class GeneratorState {
public:
	GeneratorState(const DynamicCommunitiesGenerator& generator) : generator(generator) {};
	~GeneratorState() = default;

	inline DynamicCommunitiesGenerator::Parameters getParameters() const {
		return this->generator.parameters;
	}

	inline index getTimestep() const {
		return this->generator.timestep;
	}

	inline bool getPreGenerated() const {
		return this->generator.preGenerated;
	}

	inline const std::vector<DynamicCommunitiesGenerator::Subcluster>& getSubclusters() const {
		return this->generator.subclusters;
	}

	inline DynamicCommunitiesGenerator::Subcluster getSubcluster(index i) const {
		return this->generator.subclusters.at(i);
	}

	inline index getSubclusterParent(index i) const {
		return this->generator.subclusters.at(i).parent;
	}

	inline index getSubclusterCluster(index i) const {
		return this->generator.subclusters.at(i).cluster;
	}

	inline index getSubclusterPostOrder(index i) const {
		return this->generator.subclusters.at(i).postOrder;
	}

	inline const std::vector<std::vector<index>>& getClusters() const {
		return this->generator.clusters;
	}

	inline const std::vector<index>& getCluster(index i) const {
		return this->generator.clusters.at(i);
	}

	inline const std::vector<DynamicCommunitiesGenerator::Individual>& getIndividuals() const {
		return this->generator.individuals;
	}

	inline DynamicCommunitiesGenerator::Individual getIndividual(index v) const {
		return this->generator.individuals.at(v);
	}

	inline index getIndividualSubcluster(index v) const {
		return this->generator.individuals.at(v).subcluster;
	}

	inline index getIndividualHomeSubcluster(index v) const {
		return this->generator.individuals.at(v).homeSubcluster;
	}
protected:
	const DynamicCommunitiesGenerator& generator;
};

class Tracer {
public:
	virtual void performMerge(GeneratorState gS, index i, index j) = 0;
	virtual void performedMerge(GeneratorState gS, index i, index j) = 0;

	virtual void performSplit(GeneratorState gS, index i, index j) = 0;
	virtual void performedSplit(GeneratorState gS, index i, index j) = 0;

};

#ifdef DYNAMICCOMMUNITIESGENERATOR_VALIDATOR

class GeneratorValidator {
public:
	GeneratorValidator(DynamicCommunitiesGenerator& generator) : state(generator) {}
	~GeneratorValidator() = default;

	/**
	 * Checks the state of the generator.
	 * If there are inconsistencies an exception (std::runtime_error) is thrown.
	 */
	void validate() {
		this->validate0Cluster();

		this->validateSubclustersLength();
		this->validateSubclusterCluster();
		this->validateSubclusterParent();

		Forest forest = this->buildSubclusterForest(this->state.getClusters(), 1);
		this->validateClusters(forest);
		this->validateParents(forest);
		this->validateOrder(forest, this->state.getClusters());

		this->validateIndividuals();
	}

	/**
	 * Checks the passed in list of partial trees for validity.
	 * If there are inconsistencies an exception (std::runtime_error) is thrown.
	 */
	void validatePartialTrees(const std::vector<std::vector<index>>& partialTrees) {
		Forest forest = this->buildSubclusterForest(partialTrees, 0);

		this->validatePartialTreesParents(partialTrees, forest);
		this->validateOrder(forest, partialTrees);
	}

	void validateClusterNum();

protected:
	GeneratorState state;

	struct TreeNode {
		index ind;
		index parent;
		std::vector<index> children;
	};

	struct Forest {
		std::vector<index> roots;
		std::vector<TreeNode> nodes;
	};

	class DFS {
	public:
		DFS(const Forest& forest,
			const std::vector<DynamicCommunitiesGenerator::Subcluster>& subclusters,
			const std::vector<index>& preOrders)
			: postOrder(-1), preOrder(-1), forest(forest), subclusters(subclusters),
				preOrders(preOrders) {}

		void startAt(index node);

	protected:
		index postOrder, preOrder;
		const Forest& forest;
		const std::vector<DynamicCommunitiesGenerator::Subcluster>& subclusters;
		const std::vector<index>& preOrders;

		void process(index node);
	};

	void validate0Cluster();

	void validateSubclustersLength();
	void validateSubclusterCluster();
	void validateSubclusterParent();

	void validatePartialTreesParents(const std::vector<std::vector<index>>& partialTrees,
		Forest& forest);

	Forest buildSubclusterForest(const std::vector<std::vector<index>>& clusters, index offset);
	void validateClusters(Forest& forest);
	void validateParents(Forest& forest);
	void validateOrder(Forest& forest, const std::vector<std::vector<index>>& clusters);

	void validateIndividuals();
};

#endif /* DYNAMICCOMMUNITIESGENERATOR_VALIDATOR */

} /* namespace NetworKit */

#endif /* DYNAMICCOMMUNITIESGENERATOR_H_ */
