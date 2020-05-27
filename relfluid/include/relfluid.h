#ifndef RELFLUID_H
#define RELFLUID_H

#include <ode.h>

class RelFluid : public ODE 
{
    private:
	// Variable labels
	static const unsigned int U_D = 0;
	static const unsigned int U_SX = 1;
	static const unsigned int U_TAU = 2;

	const double floor = 10e-6;
	const double k = 1.0;
	const double Gamma = 5.0/3.0;

    protected: 
	virtual void applyBoundaries(bool intermediate);
	virtual void rhs(const Grid &grid, double **u, double **dtu);
	virtual void doAfterBoundaries(bool intermediate);

    public:
	RelFluid(Domain &d, Solver &s);
	virtual ~RelFluid();
	virtual void initData();
};

#endif


