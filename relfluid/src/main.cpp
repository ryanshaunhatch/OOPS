#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <cmath>
#include <cstdio>
#include <relfluid.h>
#include <relfluidparser.h>
#include <polynomialinterpolator.h>
#include <iostream>

int main(int argc, char* argv[]) {
	// Construct Domain
	Domain domain = Domain();
	int N = 101; 
	double bounds[2] = [0.0];
	double bounds[1] = domain.getBounds()[0];
	double bounds[0] = domain.getBounds()[1];
	domain.addGrid(bounds, N);

	// ODE system
	Rk4 rk4 = RK4();
	PolynomialInterpolator interpolator = PolynomialInterpolator(4);
	RelFluid ode = RelFluid(domain, rk4);
	ode.setInterpolator(&interpolator);
	ode.initData(); 

	double ti = 0.0;
	double tf = 5.0;
	double dt = domain.getCFL()*(domain.getGrids().begin())->getSpacing();
	unsigned int M = (tf - ti)/dt;
	ode.dump_csv("D00000.csv", 0, 0);
	for (unsigned int i = 0; i < M; i++) {
	double t = (i+1)*dt;
	ode.evolveStep(dt);

	char buffer[12];
	sprintf(buffer, "D%05d.csv", i+1);
	ode.dump_csv(buffer, t, 0);
	}

	return 0;
}
