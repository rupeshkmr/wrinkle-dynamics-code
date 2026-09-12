#ifndef LBFGS_WRAPPER_HPP
#define LBFGS_WRAPPER_HPP

#include "common.hpp"

#include "optimizers.hpp"

namespace argus
{
	class LBFGSWrapper : public Optimizer
	{
		Real m_epsilon;
		int m_max_iterations;

	public:
		/**
		 * Lower and upper bounds are for box constrained solver
		 */
		void init(Real tolerance, int max_iter) override;
		Real solve(std::function<Real(VectorN&, VectorN&)>& func, VectorN& X, VectorN* lb = NULL, VectorN* ub = NULL) override;
	};
}

#endif /* LBFGS_WRAPPER_HPP */
