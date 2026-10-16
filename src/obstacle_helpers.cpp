#include <argus/obstacle_helpers.hpp>
#include <argus/obstacles.hpp>
namespace argus
{
	/* void solveCollisionVertex(
	    const Obstacle* obs,
	    VectorN& X,
	    VectorN& V,
	    const VectorN& Xold,
	    Real dt,
	    Real threshold,
	    const Real& frictionK,
	    std::vector<int>* contactVertices)
	{
	    int nverts = X.rows() / 3;

	    std::vector<uint8_t> contactFlags;
	    if (contactVertices)
	        contactFlags.assign(nverts, 0);

#ifdef USE_OMP
#pragma omp parallel for
#endif
	    for (int i = 0; i < nverts; i++)
	    {
	        Vector3 p = X.segment<3>(3 * i);

	        if (!obs->inBounds(p, threshold))
	        {
	            continue;
	        }

	        Vector3 gradient;

	        Real dist = obs->getSDF(p, &gradient);

	        if (dist < threshold)
	        {
	            if (gradient.squaredNorm() < 1e-12)
	                continue;
	            gradient.normalize();

	            Vector3 correction = (threshold - dist) * gradient;

	            Vector3 v_particle = (X.segment<3>(3 * i) - Xold.segment<3>(3 * i)) / dt;
	            X.segment<3>(3 * i) += correction;

	            Vector3 v_body = Vector3::Zero();
	            if (obs->hasExplicitVelocity())
	            {
	                v_body = obs->getVelocity(p);
	            }
	            // Clamp normal velocity — remove inward component, keep outward
	            Real v_n = v_particle.dot(gradient);
	            if (v_n < 0)
	                v_particle -= v_n * gradient; // zero out penetrating component

	            applyFriction(dt, v_particle, std::abs(frictionK * (threshold - dist)), gradient, v_body);

	            V.segment<3>(3 * i) = v_particle;

	            if (contactVertices)
	                contactFlags[i] = 1;
	        }
	    }

	    if (contactVertices)
	    {
	        for (int i = 0; i < nverts; i++)
	        {
	            if (contactFlags[i])
	                contactVertices->push_back(i);
	        }
	    }
	} */
	void solveCollisionVertex(
	    const Obstacle* obs,
	    VectorN& X,
	    VectorN& V,
	    const VectorN& Xold,
	    Real dt,
	    Real threshold,
	    const Real& frictionK, unsigned int substep,
	    std::vector<int>* contactVertices)
	{
		int nverts = X.rows() / 3;
		const Real thickness_zone = 2.0 * threshold;
		std::vector<int> clampverts = { 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 94, 98, 101, 122, 127, 132, 134, 137, 138, 141, 142, 145, 146, 147, 148, 149, 150, 158 };
		std::vector<uint8_t> contactFlags;
		if (contactVertices)
			contactFlags.assign(nverts, 0);
		Real time_offset = dt * substep;
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < nverts; i++)
		{
			Vector3 p;
			Vector3 pcurr = X.segment<3>(3 * i);
			Vector3 v_b = Vector3::Zero();
			if (obs->hasExplicitVelocity())
				v_b = obs->getVelocity(pcurr);

			p = pcurr - (v_b * time_offset);

			if (!obs->inBounds(p, threshold))
				continue;

			Vector3 gradient;
			Real dist = obs->getSDF(p, &gradient);
			// if (std::find(clampverts.begin(), clampverts.end(), i) != clampverts.end())
			// {
			// 	Vector3 correction = (threshold - dist) * gradient;
			// 	X.segment<3>(3 * i) += correction;
			// 	V.segment<3>(3 * i) *= 0;
			// 	continue;
			// }
			if (dist < threshold)
			{
				if (gradient.squaredNorm() < 1e-12)
					continue;
				gradient.normalize();

				Real depth = threshold - dist;
				Real alpha = std::min(1.0, depth / thickness_zone);

				Vector3 correction = alpha * depth * gradient;
				X.segment<3>(3 * i) += correction;

				Vector3 v_particle = V.segment<3>(3 * i);

				Vector3 v_body = v_b;
				Vector3 v_rel = v_particle - v_body;
				Vector3 v_rel_after_pos = (X.segment<3>(3 * i) - Xold.segment<3>(3 * i)) / dt - v_body;
				Real vn_pre = v_rel.dot(gradient);
				Real vn_post = v_rel_after_pos.dot(gradient);
				Real delta_vn = std::max(0.0, vn_post - vn_pre);
				if (delta_vn > 0)
				{
					Vector3 v_tan = v_rel_after_pos - vn_post * gradient;
					Real v_tan_norm = v_tan.norm();
					if (vn_post < 0)
					{
						vn_post = 0;
					}

					if (v_tan_norm > 1e-12)
					{
						Real friction_scale = std::max(1.0 - frictionK * delta_vn / v_tan_norm, 0.0);
						v_tan *= friction_scale;
					}

					Vector3 v_rel_final = vn_post * gradient + v_tan;
					Vector3 v_final = v_rel_final + v_body;

					V.segment<3>(3 * i) = v_final;
					X.segment<3>(3 * i) = Xold.segment<3>(3 * i) + v_final * dt;
				}
				else
				{
					V.segment<3>(3 * i) = (X.segment<3>(3 * i) - Xold.segment<3>(3 * i)) / dt;
				}

				if (contactVertices)
					contactFlags[i] = 1;
			}
			// project amplitudes
		}

		if (contactVertices)
		{
			for (int i = 0; i < nverts; i++)
			{
				if (contactFlags[i])
					contactVertices->push_back(i);
			}
		}
	}
	void solveCollisionQuadrature(
	    const Obstacle* obs,
	    Real dt,
	    VectorN& X,
	    VectorN& V,
	    const VectorN& Xold,
	    const VectorN& W,
	    const std::vector<trimesh::triangle_t>& faces,
	    Real threshold,
	    const Real& frictionK,
	    const probeWeights& weights,
	    std::vector<int>* contactVertices,
	    int probe_points)
	{
		int nfaces = (int)faces.size();
		int nverts = X.rows() / 3;
		int num_samples = (probe_points > 25) ? 25 : probe_points;

		std::vector<Matrix33> faceProjections(nfaces, Matrix33::Zero());
		std::vector<Vector3> faceWeights(nfaces, Vector3::Zero());

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
		for (int i = 0; i < nfaces; i++)
		{
			const auto& face = faces[i];
			int idx[3] = { (int)face.v[0], (int)face.v[1], (int)face.v[2] };

			if (W[idx[0]] == 0 && W[idx[1]] == 0 && W[idx[2]] == 0)
				continue;

			Vector3 x_curr[3] = { X.segment<3>(3 * idx[0]), X.segment<3>(3 * idx[1]), X.segment<3>(3 * idx[2]) };
			Vector3 x_old[3] = { Xold.segment<3>(3 * idx[0]), Xold.segment<3>(3 * idx[1]), Xold.segment<3>(3 * idx[2]) };

			for (int j = 0; j < num_samples; j++)
			{
				Real b[3] = { weights.weights[j][0], weights.weights[j][1], weights.weights[j][2] };

				Vector3 p_old = b[0] * x_old[0] + b[1] * x_old[1] + b[2] * x_old[2];
				Vector3 p_new = b[0] * x_curr[0] + b[1] * x_curr[1] + b[2] * x_curr[2];

				Vector3 p_min = p_old.cwiseMin(p_new);
				Vector3 p_max = p_old.cwiseMax(p_new);
				if (!obs->inBounds(0.5 * (p_min + p_max), (p_max - p_min).norm() + threshold))
					continue;

				Vector3 p_mid = 0.5 * (p_old + p_new);
				bool mid_safe = (obs->getSDF(p_mid) >= threshold);

				if (obs->getSDF(p_old) >= threshold && mid_safe && obs->getSDF(p_new) >= threshold)
					continue;

				Real t_min = 0.0;
				Real t_max = 1.0;

				if (!mid_safe)
				{
					t_max = 0.5;
				}

				for (int iter = 0; iter < 4; iter++)
				{
					Real t_test = 0.5 * (t_min + t_max);
					Vector3 p_test = (1 - t_test) * p_old + t_test * p_new;

					if (obs->getSDF(p_test) < threshold)
					{
						t_max = t_test;
					}
					else
					{
						t_min = t_test;
					}
				}

				Real t_hit = t_min;
				Vector3 p_verify = (1 - t_hit) * p_old + t_hit * p_new;
				Real final_res = obs->getSDF(p_verify);

				if (final_res < threshold)
				{
					t_hit = 0.0;
				}
				Vector3 p_hit = (1 - t_hit) * p_old + t_hit * p_new;

				Vector3 n_hit;
				obs->getSDF(p_hit, &n_hit);
				if (n_hit.squaredNorm() < 1e-12)
					continue;
				n_hit.normalize();

				Vector3 remaining_travel = p_new - p_hit;
				Vector3 sliding_travel = remaining_travel - remaining_travel.dot(n_hit) * n_hit;
				Vector3 p_target = p_hit + sliding_travel;

				Vector3 n_target;
				Real dist_target = obs->getSDF(p_target, &n_target);
				if (dist_target < threshold)
				{
					if (n_target.squaredNorm() > 1e-12)
					{
						n_target.normalize();
						p_target += n_target * (threshold - dist_target);
					}
				}

				Vector3 correction = p_target - p_new + n_hit * 1e-5;

				for (int k = 0; k < 3; k++)
				{
					if (W[idx[k]] > 0)
					{
						faceProjections[i].col(k) += b[k] * correction;
						faceWeights[i][k] += b[k];
					}
				}
			}
		}

		std::vector<Vector3> global_correction(nverts, Vector3::Zero());
		std::vector<Real> global_weight(nverts, 0.0);
		std::vector<uint8_t> contact_flags;
		if (contactVertices)
			contact_flags.assign(nverts, 0);

		for (int i = 0; i < nfaces; i++)
		{
			const auto& face = faces[i];
			int v0 = face.v[0];
			int v1 = face.v[1];
			int v2 = face.v[2];

			if (faceWeights[i][0] > 0)
			{
				global_correction[v0] += faceProjections[i].col(0);
				global_weight[v0] += faceWeights[i][0];
				if (contactVertices)
					contact_flags[v0] = 1;
			}
			if (faceWeights[i][1] > 0)
			{
				global_correction[v1] += faceProjections[i].col(1);
				global_weight[v1] += faceWeights[i][1];
				if (contactVertices)
					contact_flags[v1] = 1;
			}
			if (faceWeights[i][2] > 0)
			{
				global_correction[v2] += faceProjections[i].col(2);
				global_weight[v2] += faceWeights[i][2];
				if (contactVertices)
					contact_flags[v2] = 1;
			}
		}

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
		for (int i = 0; i < nverts; i++)
		{
			if (global_weight[i] > 1e-9)
			{
				Real softness = 0.1;
				Vector3 correction = global_correction[i] / (global_weight[i] + softness);
				if (correction.norm() > 1e-4)
				{
					Vector3 normal = correction.normalized();

					X.segment<3>(3 * i) += correction;

					Vector3 v_new = (X.segment<3>(3 * i) - Xold.segment<3>(3 * i)) / dt;

					Vector3 v_body = Vector3::Zero();
					if (obs->hasExplicitVelocity())
					{
						v_body = obs->getVelocity(X.segment<3>(3 * i));
					}
					applyFriction(dt, v_new, frictionK, normal, v_body);

					V.segment<3>(3 * i) = v_new;
				}
			}
		}

		if (contactVertices)
		{
			for (int i = 0; i < nverts; i++)
			{
				if (contact_flags[i])
				{
					contactVertices->push_back(i);
				}
			}
		}
	}
	void exportToVTK(Discregrid::CubicLagrangeDiscreteGrid* grid, const std::string& filename)
	{
		auto res = grid->resolution();
		auto domain = grid->domain();
		auto cell_size = grid->cellSize();

		std::cout << "Exporting grid to " << filename << "..." << std::endl;

		std::ofstream file(filename);
		if (!file.is_open())
		{
			std::cerr << "Error: Could not open file for writing!" << std::endl;
			return;
		}

		file << "# vtk DataFile Version 3.0\n";
		file << "Discregrid SDF Debug\n";
		file << "ASCII\n";
		file << "DATASET STRUCTURED_POINTS\n";
		file << "DIMENSIONS " << res[0] << " " << res[1] << " " << res[2] << "\n";
		file << "ORIGIN " << domain.min().x() << " " << domain.min().y() << " " << domain.min().z() << "\n";
		file << "SPACING " << cell_size[0] << " " << cell_size[1] << " " << cell_size[2] << "\n";
		file << "POINT_DATA " << res[0] * res[1] * res[2] << "\n";
		file << "SCALARS distance float 1\n";
		file << "LOOKUP_TABLE default\n";

		for (unsigned int k = 0; k < res[2]; ++k)
		{
			for (unsigned int j = 0; j < res[1]; ++j)
			{
				for (unsigned int i = 0; i < res[0]; ++i)
				{
					// Calculate world position for this node
					Eigen::Vector3d pos = domain.min() + Eigen::Vector3d(i * cell_size[0], j * cell_size[1], k * cell_size[2]);

					// Sample the SDF value (assuming field 0 is your SDF)
					double value = grid->interpolate(0, pos);

					file << value << "\n";
				}
			}
		}
		file.close();
		std::cout << "Done." << std::endl;
	}
}