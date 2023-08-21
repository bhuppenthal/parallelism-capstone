#include "heat.h"

// Heat Diffusion Equation Constants
const float RHO = 8050;
const float C = 0.466;
const float K = 20;
float k_over_rho_c = K/(RHO*C);

const float DX = 1;
const float DY = 1;
const float DT = 1;
