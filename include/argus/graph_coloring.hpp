#ifndef GRAPH_COLORING_H
#define GRAPH_COLORING_H

#include <vector>
#include "common.hpp"

namespace argus
{
	// thanks to https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/blob/afa26c12594e18a0ed2c069fb992ee6ae225c482/Simulation/SimulationModel.cpp#L1033
	template <unsigned int N>
	std::vector<std::vector<unsigned int>> graphColoring(const Eigen::Matrix<int, Eigen::Dynamic, N>& elements, const unsigned int numParticles)
	{
		const unsigned int numConstraints = (unsigned int)elements.rows();
		const unsigned int numBodies = numParticles;
		std::vector<std::vector<unsigned int>> constraintGroups;

		// filter out -1 entries
		std::vector<std::vector<int>> elementsFiltered(numConstraints);
		for (unsigned int i = 0; i < numConstraints; i++)
		{
			for (unsigned int j = 0; j < N; j++)
			{
				if (elements(i, j) != -1)
				{
					elementsFiltered[i].push_back(elements(i, j));
				}
			}
		}

		// Maps in which group a particle is or 0 if not yet mapped
		std::vector<std::vector<unsigned int>> mapping;
		// std::vector<unsigned char*> mapping;

		for (unsigned int i = 0; i < numConstraints; i++)
		{
			// Vectori<N> constraint = elements.row(i);
			const std::vector<int>& constraint = elementsFiltered[i];

			bool addToNewGroup = true;
			for (unsigned int j = 0; j < constraintGroups.size(); j++)
			{
				bool addToThisGroup = true;

				for (unsigned int k = 0; k < constraint.size(); k++)
				{
					if (mapping[j][constraint[k]] != 0)
					{
						addToThisGroup = false;
						break;
					}
				}

				if (addToThisGroup)
				{
					constraintGroups[j].push_back(i);

					for (unsigned int k = 0; k < constraint.size(); k++)
						mapping[j][constraint[k]] = 1;

					addToNewGroup = false;
					break;
				}
			}
			if (addToNewGroup)
			{
				mapping.push_back(std::vector<unsigned int>(numBodies, 0));
				// memset(mapping[mapping.size() - 1], 0, sizeof(unsigned char)*numBodies);
				constraintGroups.resize(constraintGroups.size() + 1);
				constraintGroups[constraintGroups.size() - 1].push_back(i);
				for (unsigned int k = 0; k < constraint.size(); k++)
					mapping[constraintGroups.size() - 1][constraint[k]] = 1;
			}
		}

		// for (unsigned int i = 0; i < mapping.size(); i++)
		// {
		//      delete[] mapping[i];
		// }
		// mapping.clear();

		return constraintGroups;
	}

	std::vector<std::vector<unsigned int>> graphColoringMcs(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& elements, const unsigned int numParticles);
};

#endif // GRAPH_COLORING_H
