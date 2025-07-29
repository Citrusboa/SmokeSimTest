#ifndef CONSTANTS_H
#define CONSTANTS_H

#define LERP(a, b, t) (1-t) * a + t * b


//#ifdef _DEBUG
//const int theDim[3] = { 4, 4, 1 };
//#else
//#endif

extern const int theDim[3];
extern const int theMillisecondsPerFrame;
extern const float dt;
extern const float theCellSize;
extern const float fluidDensity;
extern const float theAirDensity;
extern const float theBoundConstant;
extern const float theBuoyancyAlpha;
extern const float theBuoyancyBeta;
extern const float theBuoyancyAmbientTemperature;
extern const float theVorticityEpsilon;

#endif