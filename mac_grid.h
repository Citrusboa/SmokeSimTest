#pragma once
#ifndef MACGrid_H_
#define MACGrid_H_

#include "core.h"
#include "constants.h"

#include "linalg.h"
using namespace linalg::aliases;

#include "grid_data.h"
#include "grid_data_matrix.h"

#include <iostream>

class MACGrid {

public:
	MACGrid();
	~MACGrid();
	MACGrid(const MACGrid& orig);
	MACGrid& operator=(const MACGrid& orig);

	// member variables
	std::vector<float3> particles;
	std::vector<float3> particles_vel;

	void resetSim();
	void updateSources(int RES);
	void advectVelocity(float dt);
	void addExternalForces(float dt);
	void project(float dt);
	void project_simple(float dt);
	void advectTemperature(float dt);
	void advectDensity(float dt);
	void advectParticles(float dt);

	std::vector<float3>& getParticles();
	std::vector<float3>& getParticlesVelocity();

	void saveSmoke(const char* fileName);
	void saveParticles(std::string filename);
	void saveDensity(std::string filename);

	void copyDensityToBuffer(float* bufferOut, int resolution);
	void SimulateStep(float dt, int RES);

protected:
	// member variables
	GridDataX mU;		// X component of velocity, on X faces
	GridDataY mV;
	GridDataZ mW;
	GridData mP;
	GridData mD;
	GridData mT;

	GridData centeral_vel_x;
	GridData centeral_vel_y;
	GridData centeral_vel_z;

	GridData omegaX;
	GridData omegaY;
	GridData omegaZ;
	GridData omegaN;

	GridData omegaGX;
	GridData omegaGY;
	GridData omegaGZ;

	GridData vorConfFX;
	GridData vorConfFY;
	GridData vorConfFZ;

	GridData divergence;
	GridDataMatrix AMatrix;
	GridData precon;
	float3 sphereCenter;

	void initialize();

	// Simulation
	void computeBuoyancy(float dt);
	void computeVorticityConfinement(float dt);
	void computeCentralVel();
	void computeOmega();
	void computeOmegaGradient();
	void computeDivergence();
	void computeBound();

	// Rendering

	struct Cube { float3 pos; float4 color; float dist; };

	// GridData accessors
	enum Direction { X, Y, Z };
	float3 getVelocity(const float3& pt);
	float getVelocityX(const float3& pt);
	float getVelocityY(const float3& pt);
	float getVelocityZ(const float3& pt);
	float getTemperature(const float3& pt);
	float getDensity(const float3& pt);
	float3 getCenter(int i, int j, int k);

	void checkPressure(int& i, int& j, int& k, const GridData& p, float3& pLow, float3& pHigh);
	float3 getRewoundPosition(const float3& currPosition, const float dt);
	float3 clipToGrid(const float3& outsidePoint, const float3& insidePoint);
	float getSize(int dimension);
	int getCellIndex(int i, int j, int k);
	void getCellIndexReverse(int idx, int& i, int& j, int& k);
	int getNumberOfCells();
	bool isValidCell(int i, int j, int k);
	bool isValidFace(int dimension, int i, int j, int k);
	float3 getFacePosition(int dimension, int i, int j, int k);
	void calculateAMatrix();
	bool preconditionedConjugateGradient(const GridDataMatrix& A, GridData& p, const GridData& d, int maxIter, float tolerance);
	bool conjugateGradient(const GridDataMatrix& A, GridData& p, const GridData& d, int maxIter, float tolerance);
	void calculatePreconditioner(const GridDataMatrix& A);
	void applyPreconditioner(const GridData& r, const GridDataMatrix& A, GridData& a);

	float dotProduct(const GridData& v1, const GridData& v2);
	void add(const GridData& v1, const GridData& v2, GridData& result);
	void subtract(const GridData& v1, const GridData& v2, GridData& result);
	void multiply(const float scalar, const GridData& v, GridData& result);
	float maxMagnitude(const GridData& v);
	void apply(const GridDataMatrix& matrix, const GridData& v, GridData& result);

};
#endif