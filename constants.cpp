#include "constants.h"

extern const int theDim[3] = { 30, 30, 30 };
extern const int theMillisecondsPerFrame = 10;
extern const float dt = 0.1;
extern const float theCellSize = 1;
extern const float fluidDensity = 2.0f;
extern const float theAirDensity = 1.0f;
extern const float theBoundConstant = (fluidDensity * theCellSize) / dt;
extern const float theBuoyancyAlpha = 0.08f; // Gravity portion of particles
extern const float theBuoyancyBeta = 0.87; // Buoyancy effect due to temperature diff
extern const float theBuoyancyAmbientTemperature = 0.0f;
extern const float theVorticityEpsilon = 0.1f;