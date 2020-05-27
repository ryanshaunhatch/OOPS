#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <cmath>
#include <cstdio>
#include <wave.h>
#include <polynomialinterpolator.h>

int main(int argc, char* argv[]) {
      // Construct our domain and a grid to fit on it
      Domain domain = Domain();
      int N = 101;
      double bounds[2] = {0.0};
      bounds[0] = domain.getBounds()[0];
      bounds[1] = domain.getBounds()[1];
      domain.addGrid(bounds, N);

      // Set up our ODE system
      RK4 rk4 = RK4();
      PolynomialInterpolator interpolator = PolynomialInterpolator(4);

      Wave ode = Wave(domain, rk4);
      ode.setInterpolator(&interpolator);
      ode.initData();

      double ti = 0.0;
      double tf = 5.0;
      double dt = domain.getCFL()*(domain.getGrids().begin())->getSpacing();
      unsigned int M = (tf - ti)/dt;
      ode.dump_csv("phi00000.csv", 0, 0);
     // ode.output_frame(fnames[0], 0.0, it, 0);
     // ode.output_frame(fnames[1], 0.0, it, 1);
      for (unsigned int i = 0; i < M; i++) {
        double t = (i+1)*dt;
	ode.evolveStep(dt);

	char buffer[12];
	sprintf(buffer, "phi%05d.csv", i+1);
	ode.dump_csv(buffer, t, 0);
      }
      
      return 0;
}
