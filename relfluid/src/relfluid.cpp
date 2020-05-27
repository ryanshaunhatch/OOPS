#include <relfluid.h>
#include <operators.h>
#include <iostream>
#include <cmath>

// Constructor
RelFluid::RelFluid(Domain &d, Solver &s) : ODE(3,0) {
    if (d.getGhostPoints() < 2) {
        std::cerr<<"Warning: Domain has fewer ghost points than expected.\n";
    }

    domain = &d;
    solver = &s;

    // Set Parameters
    params = new Parameters();

    reallocateData();

}

// Destructor
RelFluid::~RelFluid() {
    delete params;
}

// Initial Data
RelFluid::initData() {
    double x0 = 0.0;
    double sigma = 0.5;

    for (auto it = data.begin(); it != data.end(); ++it) {
        const double *x = it->getGrid().getPoints();
	unsigned int nx = it->getGrid().getSize();
	double **u = it->getData();
	for (unsigned int i = 0; i < nx; i++) {
            u[V_RHO][i] = std::exp(-(x[i]-x0)*(x[i]-x0)/(sigma*sigma)) + floor;
            u[V_VX][i] = 0.0;
            u[V_P][i] = k * pow(u[V_RHO][i],Gamma);
	}
    }
}

/*-----------------------------------------------------------------------------------------------------
 * 	dt_D = - dx_(D*vx)
 * 	dt_Sx = - dx_(Sx*vx + P)
 * 	dt_Tau = - dx_(Sx - D*vx)
 *---------------------------------------------------------------------------------------------------*/

void RelFluid::rhs(const Grid &grid, double **u, double **dtu) {
        unsigned int nb = domain->getGhostPoints();
	double dx = grid.getSpacing();
	int shp = grid.getSize();
/*-----------------------------------------------------------------------------------------------------
 * 	We're going to need to solve these three equations in specific order.
 *----------------------------------------------------------------------------------------------------*/
	
}
void RelFluid::applyBoundaries(bool intermediate) {
	 unsigned int nb = domain->getGhostPoints();

	 auto left_it = data.begin();
	 auto right_it = --data.end(); 

	 double **left;
	 double **right; 

	 if (!intermediate) {
	 	left = left_it->getData();
		right =right_it->getData();
	 }
	 else {
		 left = left_it->getIntermediateData();
		 right = right_it->getIntermediateData();
	 }

	 unsigned int nr = right_it->getGrid().getSize();
	 double dx = right_it->getGrid().getSpacing();
/*-----------------------------------------------------------------------------------------------------
 * 	Then a for loop to apply the boundary conditions to all three quantities. 
 *----------------------------------------------------------------------------------------------------*/ 
}
