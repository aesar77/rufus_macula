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

  double sigma = 0.05 / std::sqrt(2.0);
  for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
      int idx = i + j * Nx;
      double x = (i - 1) * Lx / (Nx - 2);
      double y = (j - 1) * Ly / (Ny - 2);

      // parameters for a sound wave
      double wave = std::sin(6.0 * M_PI * x);
      double amplitude = 0.001;
      double rho0 = 1.0;
      double P0 = 1.0;
      double c = std::sqrt(gamma * P0 / rho0);

      // initial values for primitive variables
      rho[idx] = 1.0 + amplitude * wave;
      vx[idx] = c * amplitude * wave / rho0;
      vy[idx] = 0.0;
      P[idx] = 1.0 + c * c * amplitude * wave;
    }
  }

  solver.init(rho, vx, vy, P);
  solver.solve(0.0, 1.0);

  return 0;

}
