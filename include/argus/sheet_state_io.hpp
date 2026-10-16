#ifndef CHECKPOINTING_HPP
#define CHECKPOINTING_HPP
// Taken from https://github.com/gabrielliew/Save-and-Load-Eigen-Cpp-Matrices-Arrays-to-and-from-CSV-files/blob/header_only/save_load_eigen_csv.hpp

#include <Eigen/Dense>
#include <boost/filesystem.hpp>
#include <map>
#include "sheet.hpp"

namespace argus
{
	void initCheckpointing(Sheet& sheet);
	// VectorN convertSTDVecVecN(const std::vector<argus::Real>& ip);
	void saveSheetState(std::string saveDir, int frame, Sheet& sheet, std::map<std::string, int>* parametersToSave, nlohmann::json state);
	void saveSheetObj(std::string saveDir, int frame, Sheet& sheet, argus::Real& energy);
	// MatrixNN openData(std::string fileToOpen);
	bool loadSheetState(std::string loadDir, int frame, Sheet& sheet);
	void saveSheet(const std::string fileName, const Sheet& sheet);
};
#endif /* CHECKPOINTING */
