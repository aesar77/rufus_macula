#include "fluid_solver_2d.h"

int main() {
  int Nx = 258;
  int Ny = 258;
  double Lx = 1.0;
  double Ly = 1.0;
  double gamma = 5.0 / 3.0;
  double output_dt = 0.01;

  fluid_solver_2d solver(Lx, Ly, Nx, Ny, gamma, output_dt);

  // initial conditions
  std::vector<double> rho(Nx * Ny);
  std::vector<double> vx(Nx * Ny);
  std::vector<double> vy(Nx * Ny);
  std::vector<double> P(Nx * Ny);

  for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
      int idx = i + j * Nx;
      double x = (i - 1) * Lx / (Nx - 2);
      double y = (j - 1) * Ly / (Ny - 2);

      // parameters for a sound wave
      double rho0 = 1.0;
      double P0 = 2.5;
      double a = 0.45;
      double b = 0.55;
      double rho_d = 4.0;
      double v0 = 0.5;
      double k = 6*M_PI;
      double vsmall = 0.01;
      double sigma = 0.05;

      if ((a <= y) && (y <= b)){
        rho[idx] = rho_d;
        vx[idx] = v0;
        vy[idx] = vsmall*std::sin(k*x)*(std::exp(-((y-a)*(y-a))/(sigma*sigma)) + std::exp(-((y-b)*(y-b))/(sigma*sigma)));
        P[idx] = P0;

    } else {
        rho[idx] = rho0;
        vx[idx] = -v0;
        vy[idx] = vsmall*std::sin(k*x)*(std::exp(-((y-a)*(y-a))/(sigma*sigma)) + std::exp(-((y-b)*(y-b))/(sigma*sigma)));
        P[idx] = P0;       
    }
  }
}

  solver.init(rho, vx, vy, P);
  solver.solve(0.0, 4.0);

  return 0;

}
