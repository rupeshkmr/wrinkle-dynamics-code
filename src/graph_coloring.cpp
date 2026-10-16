#include <argus/graph_coloring.hpp>
#include <array>
#include <iostream>
#include <random>
#include <chrono>
#include <queue>
#include <memory>
#include <set>
namespace argus
{
	class GraphColor
	{
	protected:
		/* Internal Members */
		std::vector<std::vector<int>> graph; // node x neighbors
		std::vector<int> graph_colors;
		std::vector<std::vector<int>> categories;

	public:
		/* Constructors */
		GraphColor();
		explicit GraphColor(const std::vector<std::vector<int>>& graph);

		/* Mutators */
		virtual std::vector<int>& color() = 0;
		void set_graph(const std::vector<std::vector<int>>& new_graph) { this->graph = new_graph; }
		void modify_graph(int node, const std::vector<int>& neighbors) { this->graph[node] = neighbors; }

		/* Accessors */
		unsigned size() { return this->graph.size(); }
		bool is_colored();
		std::vector<int>& get_coloring() { return this->graph_colors; }
		int get_color(int node) { return graph_colors[node]; }
		int get_num_colors();
		bool is_valid();

		/* Print functions */
		void convertToColoredCategories();

		/* Functions for coloring balancing */
		float findLargestSmallestCategories(int& biggestCategory, int& smallestCategory);
		// return the category id of the changable node, not the node id, -1 if not changable
		int findChangableNodeInCategory(int sourceColor, int destinationColor);
		void changeColor(int sourceColor, int categoryId, int destinationColor);
		bool changable(int node, int destinationColor);
		void balanceColoredCategories(float goalMaxMinRatio = 1.5);
	};

	/*
	 * MCS (Register Allocation via Coloring of Chordal Graphs - Magno et al.)
	 */
	class Mcs : public GraphColor
	{
	public:
		/* Constructors */
		// Mcs(const Graph& graph) : GraphColor(graph) {}
		explicit Mcs(const std::vector<std::vector<int>>& graph)
		    : GraphColor(graph)
		{
		}

		/* Mutators */
		std::vector<int>& color();
	};

	GraphColor::GraphColor(const std::vector<std::vector<int>>& inGraph)
	{
		graph = inGraph;
		for (size_t i = 0; i < graph.size(); i++)
		{
			graph_colors.push_back(-1);
		}
	}

	int GraphColor::get_num_colors()
	{
		int numColors = 0;
		for (size_t i = 0; i < graph.size(); i++)
		{
			if (graph_colors[i] + 1 > numColors)
			{
				numColors = graph_colors[i] + 1;
			}
		};
		return numColors;
	}

	bool GraphColor::is_valid()
	{
		if (this->graph_colors.size() == 0 || this->graph.size() != this->graph_colors.size())
		{
			return false;
		}
		for (size_t i = 0; i < graph.size(); i++)
		{
			if (graph_colors[i] == -1)
			{
				return false;
			}
			for (size_t j = 0; j < graph[i].size(); j++)
			{
				int neiNode = graph[i][j];
				if (graph_colors[i] == graph_colors[neiNode])
				{
					return false;
				}
			}
		}
		return true;
	}

	void GraphColor::convertToColoredCategories()
	{
		categories.clear();
		size_t numColors = get_num_colors();
		categories.resize(numColors);

		for (int iV = 0; iV < size(); iV++)
		{
			int color = get_color(iV);
			// in this
			categories[color].push_back(iV);
		}
	}

	float GraphColor::findLargestSmallestCategories(int& biggestCategory, int& smallestCategory)
	{
		if (categories.size() == 0)
		{
			biggestCategory = -1;
			smallestCategory = -1;

			return 1;
		}

		size_t maxSize = categories[0].size();
		biggestCategory = 0;
		size_t minSize = categories[0].size();
		smallestCategory = 0;

		for (size_t iColor = 0; iColor < categories.size(); iColor++)
		{
			if (maxSize < categories[iColor].size())
			{
				biggestCategory = iColor;
				maxSize = categories[iColor].size();
			}

			if (minSize > categories[iColor].size())
			{
				smallestCategory = iColor;
				minSize = categories[iColor].size();
			}
		}

		return float(categories[biggestCategory].size()) / float(categories[smallestCategory].size());
	}

	int GraphColor::findChangableNodeInCategory(int sourceColor, int destinationColor)
	{
		auto& sourceCategory = categories[sourceColor];
		for (size_t iNode = 0; iNode < sourceCategory.size(); iNode++)
		{
			if (changable(sourceCategory[iNode], destinationColor))
			{
				return iNode;
			}
		}
		return -1;
	}

	void GraphColor::changeColor(int sourceColor, int categoryId, int destinationColor)
	{
		int nodeId = categories[sourceColor][categoryId];
		graph_colors[nodeId] = destinationColor;

		if (categories.size())
		{
			categories[sourceColor].erase(categories[sourceColor].begin() + categoryId);
			categories[destinationColor].push_back(nodeId);
		}
	}

	bool GraphColor::changable(int node, int destinationColor)
	{
		// loop through node and see if it has destinationColor
		for (size_t i = 0; i < graph[node].size(); i++)
		{
			int neiId = graph[node][i];
			if (graph_colors[neiId] == destinationColor)
			{
				return false;
			}
		}
		return true;
	}

	void GraphColor::balanceColoredCategories(float goalMaxMinRatio)
	{
		float maxMinRatio = -1.f;

		do
		{
			int biggestCategory = -1, smallestCategory = -1;

			maxMinRatio = findLargestSmallestCategories(biggestCategory, smallestCategory);

			// find a availiable vertex from the biggest category to move to the smallest category
			int changableId = findChangableNodeInCategory(biggestCategory, smallestCategory);
			if (changableId == -1)
			{
				for (size_t iColor = 0; iColor < categories.size(); iColor++)
				{
					if (iColor == biggestCategory || iColor == smallestCategory)
					{
						continue;
					}

					changableId = findChangableNodeInCategory(iColor, smallestCategory);

					if (changableId != -1)
					{
						biggestCategory = iColor;

						break;
					}
				}
			}

			if (changableId == -1)
			{
				throw std::runtime_error("The graph is not opimizable anymore, terminated with a max/min ratio: " + std::to_string(maxMinRatio));
			}
			changeColor(biggestCategory, changableId, smallestCategory);

			// change the color of changable id

		} while (maxMinRatio > goalMaxMinRatio);

		// std::cout << "The graph optimization terminated with a max/min ratio: " << maxMinRatio << std::endl;
	}

	std::vector<int>& Mcs::color()
	{
		// std::list<int> temp_graph;
		std::vector<int> temp_graph;
		for (size_t i = 0; i < this->graph.size(); i++)
		{
			temp_graph.push_back(i);
		}

		std::vector<int> weight(temp_graph.size());
		std::queue<int> ordering;

		// std::cout << "Initializing.\n";

		// Initially set the weight of each node to 0
		for (int i = 0; i < weight.size(); ++i)
		{
			weight[i] = 0;
		}

		// std::cout << "Working through all the nodes in the graph to update maximum weight.\n";

		// Work through all the nodes in the graph, choosing the node
		// with maximum weight, then add that node to the queue. Increase
		// the weight of the queued nodes neighbors by 1. Continue until
		// every node in the graph has been added to the queue

		int percentage = 0;

		std::vector<int> coloringOrder(this->graph.size());
		for (int iNode = 0; iNode < this->graph.size(); iNode++)
		{
			coloringOrder[iNode] = iNode;
		}

		// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine e(0);
		std::shuffle(std::begin(coloringOrder), std::end(coloringOrder), e);

		for (int i = 0; i < this->graph.size(); i++)
		{
			int iNode = coloringOrder[i];

			int max_weight = -1;
			int max_vertex = -1;

			// Out of the remaining nodes, find the node with the highest weight
			;
			int maxWId;
			for (int j = 0; j < temp_graph.size(); j++)
			{
				int nodeId = temp_graph[j];
				if (nodeId < 0)
				{
					continue;
				}
				if (weight[nodeId] > max_weight)
				{
					max_weight = weight[nodeId];
					max_vertex = nodeId;
					maxWId = j;
				}
			}
			if (max_vertex == -1)
			{
				throw std::runtime_error("Error: Could not find a max weight node in the graph (reason unknown)");
			}

			// Add highest weight node to the queue and increment all of its
			// neighbors weights by 1
			ordering.push(max_vertex);
			for (unsigned j = 0; j < this->graph[max_vertex].size(); j++)
			{
				weight[this->graph[max_vertex][j]] += 1;
			}

			// Remove the maximum weight node from the graph so that it won't
			// be accidentally added again
			// temp_graph.erase(graphIterMaxW); // 73010 used 3:21
			//*graphIterMaxW = -1;             // 73010 used 2:31
			temp_graph[maxWId] = -1; // 73010 used 2:10

			// std::cout << i << "  th iteration.\n";
			if (100.0 * (double)iNode / this->graph.size() > percentage)
			{
				percentage += 1;

				// update temp_graph to remove negative ids
				std::vector<int> temp_graph_new;
				for (size_t j = 0; j < temp_graph.size(); j++)
				{
					if (temp_graph[j] > 0)
					{
						temp_graph_new.push_back(temp_graph[j]);
					}
				}
				temp_graph = std::move(temp_graph_new);
			}
		}

		int sizeOrdering = ordering.size();
		percentage = 0;
		// std::cout << "Work through the queue in order and color each node.\n";

		// Work through the queue in order and color each node
		while (!ordering.empty())
		{
			int color = 0;

			// Find the lowest possible graph_colors for this node between
			// its neighbors
			int min = ordering.front();

			// Thanks to Michael Kochte @ Universitaet Stuttgart for the below speedup snippit

			// Collect color numbers of neighbors
			std::set<int> colorset;
			for (unsigned i = 0; i < this->graph[min].size(); i++)
			{
				int col = this->graph_colors[this->graph[min][i]];
				colorset.insert(col);
			}

			// Sort and uniquify
			std::vector<int> colorvec;
			std::copy(colorset.begin(), colorset.end(), std::back_inserter(colorvec));
			std::sort(colorvec.begin(), colorvec.end());

			// Pick the lowest color not contained
			int newcolor = 0;
			for (unsigned i = 0; i < colorvec.size(); i++)
			{
				if (colorvec[i] == newcolor)
				{
					newcolor++;
				}
			}
			color = newcolor;

			this->graph_colors[min] = color;
			ordering.pop();

			if (100.0 * (double)(sizeOrdering - ordering.size()) / sizeOrdering > percentage)
			{
				percentage += 1;
			}
		}
		return this->graph_colors;
	}

	std::vector<std::vector<unsigned int>> graphColoringMcs(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& elements, const unsigned int numParticles)
	{
		unsigned int N = elements.cols();
		std::vector<std::set<int>> vertElements(numParticles);
		for (int i = 0; i < elements.rows(); i++)
		{
			for (int j = 0; j < elements.cols(); j++)
			{
				if (elements(i, j) != -1)
				{
					vertElements[elements(i, j)].insert(i);
				}
			}
		}
		std::vector<std::vector<int>> elementGraph(elements.rows());
		for (int i = 0; i < elements.rows(); i++)
		{
			std::set<int> neighbors;
			for (int j = 0; j < elements.cols(); j++)
			{
				if (elements(i, j) != -1)
				{
					for (auto& elem : vertElements[elements(i, j)])
					{
						if (elem != i)
						{
							neighbors.insert(elem);
						}
					}
				}
			}
			elementGraph[i] = std::vector<int>(neighbors.begin(), neighbors.end());
		}
		Mcs mcs(elementGraph);
		std::vector<int>& colors = mcs.color();
		std::vector<std::vector<unsigned int>> constraintGroups;
		constraintGroups.resize(mcs.get_num_colors());
		for (unsigned int i = 0; i < elements.rows(); i++)
		{
			constraintGroups[colors[i]].push_back(i);
		}
		return constraintGroups;
	}

	// template std::vector<std::vector<unsigned int>> graphColoringMcs<2u>(const Eigen::Matrix<Real, Eigen::Dynamic><2u>& elements, const unsigned int numParticles);
	// template std::vector<std::vector<unsigned int>> graphColoringMcs<3u>(const Eigen::Matrix<Real, Eigen::Dynamic><3u>& elements, const unsigned int numParticles);
	// template std::vector<std::vector<unsigned int>> graphColoringMcs<4u>(const Eigen::Matrix<Real, Eigen::Dynamic><4u>& elements, const unsigned int numParticles);
	// template std::vector<std::vector<unsigned int>> graphColoringMcs<32u>(const Eigen::Matrix<int, Eigen::Dynamic><32u>& elements, const unsigned int numParticles);
};
