#include <argus/lbfgsb_wrapper.hpp>
namespace argus
{
	void LBFGSBWrapper::init(Real tolerance, int max_iter)
	{

		m_epsilon = tolerance;
		m_max_iterations = max_iter;
	}
	// lower bounds and upper bounds provided for box solvers
	Real LBFGSBWrapper::solve(std::function<Real(VectorN&, VectorN&)>& func, VectorN& X, VectorN* lb, VectorN* ub)
	{
		LBFGSpp::LBFGSBParam param;
		param.max_iterations = m_max_iterations;
		param.epsilon = m_epsilon;
		LBFGSpp::LBFGSBSolver solver(param);
		VectorN lowerB, upperB;
		if (lb == NULL || ub == NULL)
		{
#ifdef ARGUS_DEBUG
			// check if the bounds are defined
			std::cerr << "Warning! Upper and lower bounds not provided to box-constrained solver!" << std::endl;
#endif
			lowerB.resize(X.rows());
			lowerB.fill(-std::numeric_limits<Real>::infinity());
			upperB.resize(X.rows());
			upperB.fill(std::numeric_limits<Real>::infinity());
		}
		else
		{
			lowerB = *lb;
			upperB = *ub;
		}
		VectorN Xbak = X;
		VectorN td;
		Real tempE = func(X, td);
		// X will be updated inplace
		int iters = 0;
		try
		{
			iters = solver.minimize(func, X, tempE, lowerB, upperB);
		}
		catch (const std::runtime_error& e)
		{ // specific
			std::cerr << "runtime_error while running LBFGSB optimization: " << e.what() << '\n';
			X = Xbak;
		}
		catch (const char* e)
		{
			std::cerr << "Error while optimization! " << e << std::endl;
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
