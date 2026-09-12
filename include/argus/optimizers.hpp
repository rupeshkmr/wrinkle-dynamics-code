#ifndef OPTIMIZERS_HPP
#define OPTIMIZERS_HPP

#include "common.hpp"
#include <LBFGS.h>
#include <LBFGSB.h>

namespace argus
{
	class Optimizer
	{
	public:
		/** Must be called after all objects have been added to world
		 * else first frame computations will take more time
		 * assuming lbfgs solver, however the parameters will be similar for other static solvers
		 * the lower and upper bounds are for box constrained solvers
		 */

		virtual void init(Real tolerance, int max_iter) = 0;
		virtual Real solve(std::function<Real(VectorN&, VectorN&)>& func, VectorN& X, VectorN* lb = NULL, VectorN* ub = NULL) = 0;
		virtual ~Optimizer() = default;
	};
}

#endif /* OPTIMIZERS */
