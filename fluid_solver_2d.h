#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

double square(double x) { return x * x; }

class fluid_solver_2d {
public:
  fluid_solver_2d(double Lx0, double Ly0, int Nx0, int Ny0, double gamma0,
                  double output_dt0)
      : gamma(gamma0), Lx(Lx0), Ly(Ly0), Nx(Nx0), Ny(Ny0) {
    output_dt = output_dt0;
    dx = Lx / (Nx - 2);
    dy = Ly / (Ny - 2);

    rho.resize(Nx * Ny);
    vx.resize(Nx * Ny);
    vy.resize(Nx * Ny);
    P.resize(Nx * Ny);
    rho00.resize(Nx * Ny);

    mass.resize(Nx * Ny);
    mom_x.resize(Nx * Ny);
    mom_y.resize(Nx * Ny);
    energy.resize(Nx * Ny);

    rho_tmp.resize(Nx * Ny);
    vx_tmp.resize(Nx * Ny);
    vy_tmp.resize(Nx * Ny);
    P_tmp.resize(Nx * Ny);

    rho_Lx.resize(Nx * Ny);
    rho_Rx.resize(Nx * Ny);
    rho_Ly.resize(Nx * Ny);
    rho_Ry.resize(Nx * Ny);

    vx_Lx.resize(Nx * Ny);
    vx_Rx.resize(Nx * Ny);
    vx_Ly.resize(Nx * Ny);
    vx_Ry.resize(Nx * Ny);
    vy_Lx.resize(Nx * Ny);
    vy_Rx.resize(Nx * Ny);
    vy_Ly.resize(Nx * Ny);
    vy_Ry.resize(Nx * Ny);
    P_Lx.resize(Nx * Ny);
    P_Rx.resize(Nx * Ny);
    P_Ly.resize(Nx * Ny);
    P_Ry.resize(Nx * Ny);

    mass_flux_x.resize(Nx * Ny);
    mass_flux_y.resize(Nx * Ny);
    momx_flux_x.resize(Nx * Ny);
    momx_flux_y.resize(Nx * Ny);
    momy_flux_x.resize(Nx * Ny);
    momy_flux_y.resize(Nx * Ny);
    energy_flux_x.resize(Nx * Ny);
    energy_flux_y.resize(Nx * Ny);
  }

  ~fluid_solver_2d() {}

  void primitive_to_conserved() {
    // TODO: Compute conserved variables from primitive ones
    for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
        int n = i + j * Nx;
        mass[n] = rho[n] * dx * dy;
        mom_x[n] = rho[n] * vx[n] * dx * dy;
        mom_y[n] = rho[n] * vy[n] * dx * dy;
        energy[n] = rho[n] * (0.5 * (vx[n]*vx[n] + vy[n]*vy[n]) + P[n]/((gamma-1)*rho[n])) * dx * dy;
      }
    }
    periodic_boundary(mass);
    periodic_boundary(mom_x);
    periodic_boundary(mom_y);
    periodic_boundary(energy);
  }  

    
  void conserved_to_primitive() {
    // TODO: Compute primitive variables from conserved ones

     for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
        int n = i + j * Nx;
        rho[n] = mass[n] / (dx * dy);
        vx[n] = mom_x[n] / mass[n];
        vy[n] = mom_y[n] / mass[n];
        P[n] = rho[n] * (gamma-1) * (energy[n] / (rho[n] * dx * dy) - (0.5 * (vx[n]*vx[n] + vy[n]*vy[n])));


      }
    } 
    periodic_boundary(rho);
    periodic_boundary(vx);
    periodic_boundary(vy);
    periodic_boundary(P);
    // std::cout << "rho : " << rho[10000] <<  " vx : " << vx[10000] << " vy : " << vy[10000] << " P : " << P[10000] << std::endl;

  }

  void init(const std::vector<double> &rho0, const std::vector<double> &vx0,
            const std::vector<double> &vy0, const std::vector<double> &P0) {
    // TODO: Initialize the primitive variables using the given initial
    // condition

    for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
      int n = i + j * Nx;
      rho[n] = rho0[n];
      vx[n] = vx0[n];
      vy[n] = vy0[n];
      P[n] = P0[n];
      rho00[n] = rho0[n];
      //primitive_to_conserved();
      }
    }
      periodic_boundary(rho);
      periodic_boundary(vx);
      periodic_boundary(vy);
      periodic_boundary(P);

  }
  double find_dt() {
    // TODO: Find the optimal dt that satisfies the CFL condition, and return
    // its value
    double max;
    double max_old = 0;
    double CFL = 0.1;

    for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
        int n = i + j * Nx;
        max = (std::sqrt(gamma * (P[n]/rho[n])) + std::sqrt(vx[n]*vx[n] + vy[n]*vy[n]));

        if (max > max_old){
          max_old = max;
        }
      }
    }
  //std::cout << "max : " << max << std::endl;  
  double delta_t = CFL * std::min(dx, dy)/ max_old;

  return delta_t;
  }  

  void solve(double t0, double t_end) {
    // Solve the fluid equations, starting from t0 and stoping at t_end
    double t = t0;
    int n = 0; // n labels the output file
    while (t < t_end) {
      if (t >= output_dt * n) {
        output(n);
        n += 1;
      }
      double dt = find_dt();
      //std::cout << "deltat : " << dt << std::endl;

      step(dt);
      t += dt;
    }
  }

  void step(double dt) {
    // extrapolate a half step in time using primitive equations
    primitive_update(0.5 * dt);

    //std::cout << "primitive up " << " rho : " << rho[10000] <<  " vx : " << vx[10000] << " vy : " << vy[10000] << " P : " << P[10000] << std::endl;

    extrapolate_to_interface();

    // compute fluxes
    compute_fluxes();


    // update solultion from fluxes
    update_conserved(dt);


    // update primitive variables
    conserved_to_primitive();
    //std::cout << "cons to prim " << " rho : " << rho[10000] <<  " vx : " << vx[10000] << " vy : " << vy[10000] << " P : " << P[10000] << std::endl;

  }

  void periodic_boundary(std::vector<double> &f) {
    // TODO: apply periodic boundary conditions to an array f
    for (int j = 0; j < Ny; j++) {
      f[j * Nx] = f[j * Nx + Nx - 2];
      f[j * Nx + Nx - 1] = f[j * Nx + 1];
    }
    for (int i = 0; i < Nx; i++) {
      f[i] = f[(Ny - 2) * Nx + i];
      f[(Ny - 1) * Nx + i] = f[Nx + i];
    }
  } 
  void primitive_update(double dt) {
    // TODO: update the primitive variables using Euler equations in primitive
    // form using an FTCS scheme
    for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
        int n = i + j * Nx;
        rho_tmp[n] = rho[n] - dt * (vx[n] * ((rho[n+1]-rho[n-1])/(2*dx)) + vy[n] * ((rho[n+Nx]-rho[n-Nx])/(2*dy)) + rho[n] * ((vx[n+1]-vx[n-1])/(2*dx) + (vy[n+Nx]-vy[n-Nx])/(2*dy)));
        vx_tmp[n] = vx[n] - (dt * (vx[n] * ((vx[n+1]-vx[n-1])/(2*dx)) + vy[n] * ((vx[n+Nx]-vx[n-Nx])/(2*dy)) + ((P[n+1 ]-P[n-1 ])/(2*dx)))/rho[n]);
        vy_tmp[n] = vy[n] - (dt * (vx[n] * ((vy[n+1]-vy[n-1])/(2*dx)) + vy[n] * ((vy[n+Nx]-vy[n-Nx])/(2*dy)) + ((P[n+Nx]-P[n-Nx])/(2*dy)))/rho[n]);
        P_tmp[n] = P[n] - dt * (vx[n] * (((P[n+1]-P[n-1])/(2*dx)) + vy[n] * ((P[n+Nx]-P[n-Nx])/(2*dy)) + gamma * P[n] * ((vx[n+1]-vx[n-1])/(2*dx) + (vy[n+Nx]-vy[n-Nx])/(2*dy))));

      } 
    }
    for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
        int n = i + j * Nx;
        rho[n] = rho_tmp[n];
        vx[n] = vx_tmp[n];
        vy[n] = vy_tmp[n];
        P[n] = P_tmp[n];
      }
    }

    periodic_boundary(rho);
    periodic_boundary(vx);
    periodic_boundary(vy);
    periodic_boundary(P);

 }

  void extrapolate_to_interface() {
    // TODO: compute rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, and P_R here
    for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
        int n = i + j * Nx;

        rho_Lx[n] = rho[n] - 0.25 * (rho[n+1] - rho[n-1]);
        rho_Rx[n] = rho[n] + 0.25 * (rho[n+1] - rho[n-1]);

        rho_Ly[n] = rho[n] - 0.25 * (rho[n+Nx] - rho[n-Nx]);
        rho_Ry[n] = rho[n] + 0.25 * (rho[n+Nx] - rho[n-Nx]);

        vx_Lx[n] = vx[n] - 0.25 * (vx[n+1] - vx[n-1]);
        vx_Rx[n] = vx[n] + 0.25 * (vx[n+1] - vx[n-1]);

        vx_Ly[n] = vx[n] - 0.25 * (vx[n+Nx] - vx[n-Nx]);
        vx_Ry[n] = vx[n] + 0.25 * (vx[n+Nx] - vx[n-Nx]);

        vy_Lx[n] = vy[n] - 0.25 * (vy[n+1] - vy[n-1]);
        vy_Rx[n] = vy[n] + 0.25 * (vy[n+1] - vy[n-1]);

        vy_Ly[n] = vy[n] - 0.25 * (vy[n+Nx] - vy[n-Nx]);
        vy_Ry[n] = vy[n] + 0.25 * (vy[n+Nx] - vy[n-Nx]);

        P_Lx[n] = P[n] - 0.25 * (P[n+1] - P[n-1]);
        P_Rx[n] = P[n] + 0.25 * (P[n+1] - P[n-1]);

        P_Ly[n] = P[n] - 0.25 * (P[n+Nx] - P[n-Nx]);
        P_Ry[n] = P[n] + 0.25 * (P[n+Nx] - P[n-Nx]);

      }
    }
    periodic_boundary(rho_Lx);
    periodic_boundary(rho_Rx);
    periodic_boundary(rho_Ly);
    periodic_boundary(rho_Ry);
    periodic_boundary(vx_Lx);
    periodic_boundary(vx_Ly);
    periodic_boundary(vx_Rx);
    periodic_boundary(vx_Ry);
    periodic_boundary(vy_Lx);
    periodic_boundary(vy_Ly);
    periodic_boundary(vy_Rx);
    periodic_boundary(vy_Ry);
    periodic_boundary(P_Lx);
    periodic_boundary(P_Ly);
    periodic_boundary(P_Rx);
    periodic_boundary(P_Ry);



  }
  void compute_fluxes() {
    // TODO: compute the fluxes
    for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
        int n = i + j * Nx;
        double vmaxX, vmaxY;
        double vLx = (std::sqrt(gamma * (P_Rx[n]/rho_Rx[n])) + std::sqrt(vx_Rx[n]*vx_Rx[n] + vy_Rx[n]*vy_Rx[n]));
        double vRx = (std::sqrt(gamma * (P_Lx[n+1]/rho_Lx[n+1])) + std::sqrt(vx_Lx[n+1]*vx_Lx[n+1] + vy_Lx[n+1]*vy_Lx[n+1]));
        double vLy = (std::sqrt(gamma * (P_Ry[n]/rho_Ry[n])) + std::sqrt(vx_Ry[n]*vx_Ry[n] + vy_Ry[n]*vy_Ry[n]));
        double vRy = (std::sqrt(gamma * (P_Ly[n+Nx]/rho_Ly[n+Nx])) + std::sqrt(vx_Ly[n+Nx]*vx_Ly[n+Nx] + vy_Ly[n+Nx]*vy_Ly[n+Nx]));

        if (vLx < vRx){
          vmaxX = vRx;
        } else {
          vmaxX = vLx;
        }
        if (vLy < vRy){
          vmaxY = vRy;
        } else {
          vmaxY = vLy;
        }

        //mass fluxes and conserved values
        double mFy_L = rho_Ry[n]*vy_Ry[n];
        double mFy_R = rho_Ly[n+Nx]*vy_Ly[n+Nx];
        double mQy_R = rho_Ly[n+Nx];
        double mQy_L = rho_Ry[n];

        double mFx_L = rho_Rx[n]*vx_Rx[n];
        double mFx_R = rho_Lx[n+1]*vx_Lx[n+1];
        double mQx_R = rho_Lx[n+1];
        double mQx_L = rho_Rx[n];
        //momentum fluxes and conserved values
        //mom_x
        double momxFy_L = rho_Ry[n]*vx_Ry[n]*vy_Ry[n];
        double momxFy_R = rho_Ly[n+Nx]*vx_Ly[n+Nx]*vy_Ly[n+Nx];
        double momxQy_R = rho_Ly[n+Nx]*vx_Ly[n+Nx];
        double momxQy_L = rho_Ry[n]*vx_Ry[n];       

        double momxFx_L = rho_Rx[n]*vx_Rx[n]*vx_Rx[n] + P_Rx[n];
        double momxFx_R = rho_Lx[n+1]*vx_Lx[n+1]*vx_Lx[n+1] + P_Lx[n+1];
        double momxQx_R = rho_Lx[n+1]*vx_Lx[n+1];
        double momxQx_L = rho_Rx[n]*vx_Rx[n]; 

        //mom_y
        double momyFy_L = rho_Ry[n]*vy_Ry[n]*vy_Ry[n] + P_Ry[n];
        double momyFy_R = rho_Ly[n+Nx]*vy_Ly[n+Nx]*vy_Ly[n+Nx] + P_Ly[n+Nx];
        double momyQy_R = rho_Ly[n+Nx]*vy_Ly[n+Nx];
        double momyQy_L = rho_Ry[n]*vy_Ry[n];       

        double momyFx_L = rho_Rx[n]*vx_Rx[n]*vy_Rx[n];
        double momyFx_R = rho_Lx[n+1]*vx_Lx[n+1]*vy_Lx[n+1];
        double momyQx_R = rho_Lx[n+1]*vy_Lx[n+1];
        double momyQx_L = rho_Rx[n]*vy_Rx[n]; 

        //energy fluxes and conserved valuables
        double enFy_L_u =  P_Ry[n]/((gamma-1)*rho_Ry[n]);
        double enFy_L_U = rho_Ry[n]*(0.5*std::sqrt(vx_Ry[n]*vx_Ry[n] + vy_Ry[n]*vy_Ry[n]) + enFy_L_u);
        double enFy_L =   (enFy_L_U + P_Ry[n])*vy_Ry[n];

        double enFy_R_u =  P_Ly[n+Nx]/((gamma-1)*rho_Ly[n+Nx]);
        double enFy_R_U = rho_Ly[n+Nx]*(0.5*std::sqrt(vx_Ly[n+Nx]*vx_Ly[n+Nx] + vy_Ly[n+Nx]*vy_Ly[n+Nx]) + enFy_R_u);
        double enFy_R =   (enFy_R_U + P_Ly[n+Nx])*vy_Ly[n+Nx];

        double enQy_L = enFy_L_U;
        double enQy_R = enFy_R_U;

        double enFx_L_u =  P_Rx[n]/((gamma-1)*rho_Rx[n]);
        double enFx_L_U = rho_Rx[n]*(0.5 * std::sqrt(vx_Rx[n]*vx_Rx[n] + vy_Rx[n]*vy_Rx[n]) + enFx_L_u);
        double enFx_L =   (enFx_L_U + P_Rx[n])*vx_Rx[n];

        double enFx_R_u =  P_Lx[n+1]/((gamma-1)*rho_Lx[n+1]);
        double enFx_R_U = rho_Lx[n+1]*(0.5*std::sqrt(vx_Lx[n+1]*vx_Lx[n+1] + vy_Lx[n+1]*vy_Lx[n+1]) + enFx_R_u);
        double enFx_R =   (enFx_R_U + P_Lx[n+1])*vx_Lx[n+1];

        double enQx_L = enFx_L_U;
        double enQx_R = enFx_R_U;


        mass_flux_y[n] = 0.5 * ((mFy_L + mFy_R ) - vmaxY*(mQy_R - mQy_L));
        mass_flux_x[n] = 0.5 * ((mFx_L + mFx_R ) - vmaxX*(mQx_R - mQx_L));

        momy_flux_y[n] = 0.5 * ((momyFy_L + momyFy_R) - vmaxY*(momyQy_R - momyQy_L));
        momy_flux_x[n] = 0.5 * ((momyFx_L + momyFx_R) - vmaxX*(momyQx_R - momyQx_L));

        momx_flux_y[n] = 0.5 * ((momxFy_L + momxFy_R) - vmaxY*(momxQy_R - momxQy_L));
        momx_flux_x[n] = 0.5 * ((momxFx_L + momxFx_R) - vmaxX*(momxQx_R - momxQx_L));       

        energy_flux_y[n] = 0.5 * ((enFy_L + enFy_R) - vmaxY*(enQy_R - enQy_L));
        energy_flux_x[n] = 0.5 * ((enFx_L + enFx_R) - vmaxX*(enQx_R - enQx_L));

        //std::cout << "momxQx_R : " << momxQx_R <<  " momyQy_R : " << momyQx_R <<  " enQx_R : " << enQx_R << " enqy_R :" << enQy_R << std::endl;  


      }
    }
    periodic_boundary(mass_flux_x);
    periodic_boundary(mass_flux_y);
    periodic_boundary(momx_flux_x);
    periodic_boundary(momx_flux_y);
    periodic_boundary(momy_flux_x);
    periodic_boundary(momy_flux_y);
    periodic_boundary(energy_flux_x);
    periodic_boundary(energy_flux_y);
  }

  void update_conserved(double dt) {
    // TODO: update the conserved variables using the fluxes
    primitive_to_conserved();
    for (int j = 1; j < Ny-1; j++) {
      for (int i = 1; i < Nx-1; i++) {
        int n = i + j * Nx;

        mass[n] = ((rho[n] * dx * dy) - (mass_flux_x[n] - mass_flux_x[n-1])*dy*dt - (mass_flux_y[n] - mass_flux_y[n-Nx])*dx*dt);
        mom_x[n] = ((rho[n] * vx[n] * dx * dy) - (momx_flux_x[n] - momx_flux_x[n-1])*dy*dt - (momx_flux_y[n] - momx_flux_y[n-Nx])*dx*dt);
        mom_y[n] = ((rho[n] * vy[n] * dx * dy) - (momy_flux_x[n] - momy_flux_x[n-1])*dy*dt - (momy_flux_y[n] - momy_flux_y[n-Nx])*dx*dt);
        energy[n] = ((rho[n] * (0.5 * (vx[n]*vx[n] + vy[n]*vy[n]) + P[n]/((gamma-1)*rho[n])) * dx * dy) - (energy_flux_x[n] - energy_flux_x[n-1])*dy*dt - (energy_flux_y[n] - energy_flux_y[n-Nx])*dx*dt);
        
      }
    }
    periodic_boundary(mass);
    periodic_boundary(mom_x);
    periodic_boundary(mom_y);
    periodic_boundary(energy);
  }

  void output(int n) {
    std::ofstream outfile("output3d_rho_" + std::to_string(n) + ".csv");
    for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
        int idx = i + j * Nx;
        outfile << rho[idx];
        if (i != Nx - 2)
          outfile << ", ";
        else
          outfile << std::endl;
      }
    }
    outfile.close();
  }

  int Nx, Ny;
  double Lx, Ly;
  double dx, dy;
  double gamma, output_dt;
  std::vector<double> rho, vx, vy, P, rho00;                 // primitive variables
  std::vector<double> mass, mom_x, mom_y, energy;     // conserved variables
  // arrays to hold the results during primitive_update
  std::vector<double> rho_tmp, vx_tmp, vy_tmp, P_tmp;
  // arrays of fluxes for each conserved variable:
  std::vector<double> mass_flux_x, mass_flux_y;
  std::vector<double> momx_flux_x, momx_flux_y;
  std::vector<double> momy_flux_x, momy_flux_y;
  std::vector<double> energy_flux_x, energy_flux_y;
  // arrays for extrapolating to cell interfaces:
  std::vector<double> rho_Lx, rho_Ly, rho_Rx, rho_Ry;
  std::vector<double> vx_Lx, vx_Ly, vx_Rx, vx_Ry;
  std::vector<double> vy_Lx, vy_Ly, vy_Rx, vy_Ry;
  std::vector<double> P_Lx, P_Ly, P_Rx, P_Ry;
};
