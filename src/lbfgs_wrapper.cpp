#include <argus/lbfgs_wrapper.hpp>
namespace argus
{
	void LBFGSWrapper::init(Real tolerance, int max_iter)
	{
		m_epsilon = tolerance;
		m_max_iterations = max_iter;
	}
	// lower bounds and upper bounds provided for box solvers
	Real LBFGSWrapper::solve(std::function<Real(VectorN&, VectorN&)>& func, VectorN& X, VectorN* lb, VectorN* ub)
	{
		/**
		 * func should handle details about active parameters
		 */
		LBFGSpp::LBFGSParam param;
		param.max_iterations = m_max_iterations;
		param.epsilon = m_epsilon;
		LBFGSpp::LBFGSSolver solver(param);

		VectorN Xbak = X;
		VectorN td;
		Real tempE = func(X, td);
		// X will be updated inplace
		try
		{
			solver.minimize(func, X, tempE);
		}
		catch (const char* e)
		{
			std::cerr << "Error while optimization! " << e << std::endl;
			X = Xbak;
		}
		catch (const std::runtime_error& e)
		{ // specific
			std::cerr << "runtime_error while running LBFGS optimization: " << e.what() << '\n';
			X = Xbak;
		}
		catch (...)
		{
			std::cerr << "Some other error occurred while optimization using LBFGSpp " << std::endl;
			X = Xbak;
		}
		return tempE;
	}
}
