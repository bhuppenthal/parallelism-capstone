#ifndef HEAT_H
#define HEAT_H

// Heat Diffusion Equation
#define CALC_DTEMP(elem, left, right, up, down) (k_over_rho_c * DT * ((down - 2*elem + up)/(DY * DY) + (left - 2*elem + right)/(DX * DX)))

// Heat Diffusion Equation Constants
extern const float RHO;
extern const float C;
extern const float K;
extern float k_over_rho_c;                 // units of m^1/s NOTE: Cannot be a const (true for OpenMP too?)

extern const float DX;
extern const float DY;
extern const float DT;

#endif
