#include <argus/get_forces.hpp>
#include <argus/gravity.hpp>
#include <argus/tfw_stretching.hpp>
#include <argus/tfw_bending.hpp>
#include <argus/tfw_gravity.hpp>
#include <argus/dphi_penalty.hpp>
#include <argus/thin_shell_stretching.hpp>
#include <argus/thin_shell_bending.hpp>
namespace argus
{

	std::optional<std::shared_ptr<argus::SheetForces>> getSheetForce(nlohmann::json& config)
	{
		if (config.contains("name") == false)
			return std::nullopt;

		std::string name = config["name"].get<std::string>();
		if (name == "tfwstretch")
			return std::make_shared<TFWStretch>();
		else if (name == "tfwbending")
			return std::make_shared<TFWBending>();
		else if (name == "thinshellbending")
			return std::make_shared<ThinShellBending>();
		else if (name == "thinshellstretching")
			return std::make_shared<ThinShellStretching>();
		else if (name == "tfwgravity")
		{
			if (config.contains("gravity") == false)
			{
				std::cerr << "Gravity value not provided with energy" << std::endl;
				return std::nullopt;
			}
			Vector3 gravity = config["gravity"].get<Vector3>();
			return std::make_shared<TFWGravity>(gravity);
		}
		else if (name == "gravity")
		{
			if (config.contains("gravity") == false)
			{
				std::cerr << "Gravity value not provided with energy" << std::endl;
				return std::nullopt;
			}
			Vector3 gravity = config["gravity"].get<Vector3>();
			return std::make_shared<Gravity>(gravity);
		}
		else if (name == "dphipenalty")
		{
			if (config.contains("threshold") == false)
			{
				std::cerr << "threshold value not provided for dphi penalty" << std::endl;
				return std::nullopt;
			}
			Real threshold = config["threshold"].get<Real>();
			return std::make_shared<DphiPenalty>(threshold);
		}
		else
		{
			std::cerr << "could not infer sheet force " << name << std::endl;
			return std::nullopt;
		}
	}

	std::optional<std::vector<std::shared_ptr<SheetForces>>> getSheetForcesVector(std::vector<nlohmann::json>& forces_data)
	{
		if (forces_data.size() == 0)
			return std::nullopt;
		std::vector<std::shared_ptr<SheetForces>> sheetForces;
		for (auto it = forces_data.begin(); it != forces_data.end(); it++)
		{
			auto retval = getSheetForce(*it);
			if (retval.has_value())
			{
				sheetForces.push_back(retval.value());
			}
		}
		if (sheetForces.size() == 0)
			return std::nullopt;
		else
			return sheetForces;
	}
}
