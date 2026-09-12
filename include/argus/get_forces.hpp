#ifndef GET_FORCES
#define GET_FORCES
#include "sheet_forces.hpp"
#include <nlohmann/json.hpp>

namespace argus
{
	std::optional<std::shared_ptr<SheetForces>> getSheetForce(nlohmann::json& config);
	std::optional<std::vector<std::shared_ptr<SheetForces>>> getSheetForcesVector(std::vector<nlohmann::json>& forces_data);
}
#endif
