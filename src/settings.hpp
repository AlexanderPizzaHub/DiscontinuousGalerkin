#pragma once
#include <cmath>

namespace settings
{
    namespace hyperbolic
    {
        double G(double x);
        double S(double x);
        double Truth(double x);
    }

    namespace shakhov
    {
        /*
        TEST CASE
        */
        // const double Kn0 = 0.645161;
        //  const double Kn0 = 20.0;
        // const double Kn0 = 0.0645161;
        // const double Kn0 = 0.00645161;
        const double Kn0 = 0.1;
        // const double Kn0=0.01;
        const double PI = 3.14159265358979323846;
        const double delta_0 = sqrt(PI) / (2.0 * Kn0);
        const double Pr = 2.0 / 3.0;
        // const double omega = 0.5;
        // const double omega = 0.74;
        const double omega = 0.74;

        // bc
        const double u_wall_l = 0.0;
        const double T_wall_l = 300.0 / 300.0;

        const double u_wall_r = 0.0;
        const double T_wall_r = 600.0 / 300.0;

        const double T_mach_l = 1.0;
        const double T_mach_r = 1.0;

        /*
        DISCRETIZE
        */
        const double xl = 0.0;
        const double xr = 1.0;

        const double vm = 10.0;

        const int Nx = 50;
        const int Nv = 50;

        /*
        SOLVER
        */
        const int dg_max_iter = 100;
        const double dg_tolerance = 1e-6;

        const int shakhov_max_iter = 10000;
        const double shakhov_tolerance = -1.0;

        const bool printinfo = true;    // print info during the simulation
        const bool write_result = true; // write result to file

        const int output_interval = 100;
    }
}