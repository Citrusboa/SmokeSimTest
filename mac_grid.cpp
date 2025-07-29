#include "mac_grid.h"

//#include <Eigen/Eigen>
//#include <Eigen/Sparse>
//#include <Eigen/Dense>
//#include <math.h>
//#include <map>
//#include <fstream>

// linear algebra 
#include "linalg.h"
using namespace linalg::aliases;

using namespace Eigen;
using namespace std;


// Global variables
MACGrid target;

// Indexing: x -> cols, z -> rows, y-> stacks

#define FOR_EACH_CELL \
	for(int k = 0; k < theDim[MACGrid::Z]; ++k) \
		for (int j = 0; j < theDim[MACGrid::Y]; ++j) \
			for (int i = 0; i < theDim[MACGrid::X]; ++i)

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; --k)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; --j) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; --i) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]+1; ++i) 

#define FOR_EACH_FACE_X \
   for(int k = 0; k < theDim[MACGrid::Z]; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]+1; ++i)

#define FOR_EACH_FACE_Y \
   for(int k = 0; k < theDim[MACGrid::Z]; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]; ++i)

#define FOR_EACH_FACE_Z \
   for(int k = 0; k < theDim[MACGrid::Z]+1; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]; ++i)

MACGrid::MACGrid() {
    initialize();
}

MACGrid::MACGrid(const MACGrid& orig) {
    mU = orig.mU;
    mV = orig.mV;
    mW = orig.mW;
    mP = orig.mP;
    mD = orig.mD;
    mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig) {
    if (&orig == this) return*this;

    mU = orig.mU;
    mV = orig.mV;
    mW = orig.mW;
    mP = orig.mP;
    mD = orig.mD;
    mT = orig.mT;

    return *this;
}

MACGrid::~MACGrid() {

}

void MACGrid::resetSim() {
    cout << "reset MAC grid" << endl;

    mU.initialize();
    mV.initialize();
    mW.initialize();
    mP.initialize();
    mD.initialize();
    mT.initialize();

    centeral_vel_x.initialize();
    centeral_vel_y.initialize();
    centeral_vel_z.initialize();

    omegaX.initialize();
    omegaY.initialize();
    omegaZ.initialize();
    omegaN.initialize();

    omegaGX.initialize();
    omegaGY.initialize();
    omegaGZ.initialize();

    vorConfFX.initialize();
    vorConfFY.initialize();
    vorConfFZ.initialize();

    divergence.initialize();

    calculateAMatrix();
    calculatePreconditioner(AMatrix);
    int l = 0;
}

void MACGrid::SimulateStep(float dt, int RES) {

    updateSources(RES);
    advectVelocity(dt);
    addExternalForces(dt);
    project(dt);
    //project_simple(dt);
    advectTemperature(dt);
    advectDensity(dt);
}

void MACGrid::copyDensityToBuffer(float* bufferOut, int resolution) {
    // check resolution is same as theDim
    assert(resolution == theDim[0]);

    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            for (int k = 0; k < resolution; ++k) {
                bufferOut[i * resolution * resolution + j * resolution + k] = mD(i, j, k);
            }
        }
    }
}

void MACGrid::initialize() {
    resetSim();

    // Set initial amount of particles in system
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                float3 cell_center(theCellSize * (i + 0.5), theCellSize * (j + 0.5), theCellSize * (k + 0.5));
                for (int p = 0; p < 10; p++) {
                    double a = ((float)rand() / RAND_MAX - 0.5) * theCellSize;
                    double b = ((float)rand() / RAND_MAX - 0.5) * theCellSize;
                    double c = ((float)rand() / RAND_MAX - 0.5) * theCellSize;
                    float3 shift(a, b, c);
                    float3 xp = cell_center + shift;
                    particles.push_back(xp);
                }
            }
        }
    }
}

void MACGrid::updateSources(int RES) {

    const int centerY = RES / 4;
    const int centerZ = RES / 2;
    float dens = (rand() % 1000) / 1000.0f;
    for (int i = 1; i <= (RES - 2); i++) {
        for (int j = 1; j <= (RES - 2); j++) {
            if (hypot(i - centerY, j - centerZ) < RES / 10) {
                this->mD(5, i, j) = dens;
                this->mU(5, i, j) = 2.0f;
            }
        }
    }
}

void MACGrid::advectVelocity(float dt) {
    target.mU = mU;
    target.mV = mV;
    target.mW = mW;

    // RK2 Integration
    FOR_EACH_FACE_X{
        float3 curPos = float3(i * theCellSize, (j + 0.5f) * theCellSize, (k + 0.5f) * theCellSize);
        float3 midPos = curPos - 0.5f * getVelocity(curPos) * dt;
        target.mU(i, j, k) = getVelocityX(curPos - getVelocity(midPos) * dt);
    };

    FOR_EACH_FACE_Y{
        float3 curPos = float3((i + 0.5f) * theCellSize, j * theCellSize, (k + 0.5f) * theCellSize);
        float3 midPos = curPos - 0.5f * getVelocity(curPos) * dt;
        target.mV(i, j, k) = getVelocityY(curPos - getVelocity(midPos) * dt);

    };

    FOR_EACH_FACE_Z{
        float3 curPos = float3((i + 0.5f) * theCellSize, (j + 0.5f) * theCellSize, k * theCellSize);
        float3 midPos = curPos - 0.5f * getVelocity(curPos) * dt;
        target.mV(i, j, k) = getVelocityZ(curPos - getVelocity(midPos) * dt);

    };

    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(float dt) {
    target.mT = mT;

    FOR_EACH_CELL{
        float3 currentptc = float3((i + 0.5) * theCellSize, (j + 0.5) * theCellSize, (k + 0.5) * theCellSize);
        float3 curr_vel = getVelocity(currentptc);
        float3 mid_pt = currentptc - 0.5f * curr_vel * dt;
        float3 old_pt = currentptc - getVelocity(mid_pt);
        float new_temp = getTemperature(old_pt);
        target.mT(i, j, k) = new_temp;
    };
    mT = target.mT;
}

void MACGrid::advectParticles(float dt) {
    particles_vel.resize(particles.size());

    for (size_t p = 0; p < particles.size(); ++p) {
        float3 currentPosition = particles[p];
        float3 currentVelocity = getVelocity(currentPosition);
        float3 nextPosition = currentPosition + currentVelocity * dt;
        float3 clippedNextPosition = clipToGrid(nextPosition, currentPosition);

        // Next
        float3 nexv = getVelocity(nextPosition);
        float3 nextVelocity = getVelocity(clippedNextPosition);
        float3 averageVelocity = (currentVelocity + nextVelocity) / 2.0f;
        float3 betterNextPosition = currentPosition + averageVelocity * dt;
        float3 clippedBetterNextPosition = clipToGrid(betterNextPosition, currentPosition);

        particles[p] = clippedBetterNextPosition;
        particles_vel[p] = averageVelocity;
        //particles[p] = clippedNextPosition;
        //particles_vel[p] = currentVelocity;
    }
}

void MACGrid::advectDensity(float dt) {
    target.mD = mD;
    FOR_EACH_CELL
    {
        float3 currentptc = float3((i + 0.5) * theCellSize,(j + 0.5) * theCellSize,(k + 0.5) * theCellSize);
        float3 curr_vel = getVelocity(currentptc);
        float3 mid_pt = currentptc - curr_vel * dt / 2;
        float3 old_pt = currentptc - getVelocity(mid_pt);
        float new_density = getDensity(old_pt);
        target.mD(i,j,k) = new_density;
    };

    mD = target.mD;
}

std::vector<float3>& MACGrid::getParticles() {
    return particles;
}

std::vector<float3>& MACGrid::getParticlesVelocity() {
    return particles_vel;
}

void MACGrid::computeBuoyancy(float dt) {
    target.mV = mV;

    FOR_EACH_CELL
    {
        float s = (mD(i,j,k) + mD(i,j - 1,k)) / 2;
        float T = (mT(i,j,k) + mT(i,j - 1,k)) / 2;
        float fbuoy = -theBuoyancyAlpha * s + theBuoyancyBeta * (T - theBuoyancyAmbientTemperature);
        //cout << fbuoy << endl;
        target.mV(i,j,k) += dt * fbuoy;
    };

    mV = target.mV;
}

void MACGrid::computeCentralVel() {
    FOR_EACH_CELL
    {
        centeral_vel_x(i,j,k) = (mU(i,j,k) + mU(i + 1,j,k)) / 2;
        centeral_vel_y(i,j,k) = (mV(i,j,k) + mV(i,j + 1,k)) / 2;
        centeral_vel_z(i,j,k) = (mW(i,j,k) + mW(i,j,k + 1)) / 2;
    };
}

void MACGrid::computeOmega() {
    //std::cout << "Before compute omega" << std::endl;
    FOR_EACH_CELL
    {
        omegaX(i,j,k) = (centeral_vel_z(i,j + 1,k) - centeral_vel_z(i,j - 1,k) - centeral_vel_y(i,j,k + 1) + centeral_vel_y(i,j,k - 1)) / (2 * theCellSize);
        omegaY(i,j,k) = (centeral_vel_x(i,j,k + 1) - centeral_vel_x(i,j,k - 1) - centeral_vel_z(i + 1,j,k) + centeral_vel_z(i - 1,j,k)) / (2 * theCellSize);
        omegaZ(i,j,k) = (centeral_vel_y(i + 1,j,k) - centeral_vel_y(i - 1,j,k) - centeral_vel_x(i,j + 1,k) + centeral_vel_x(i,j - 1,k)) / (2 * theCellSize);
        float3 normalO = float3(omegaX(i,j,k),omegaY(i,j,k),omegaZ(i,j,k));
        omegaN(i,j,k) = length(normalO);

    };
}

void MACGrid::computeOmegaGradient() {
    FOR_EACH_CELL
    {
        omegaGX(i,j,k) = (omegaN(i + 1,j,k) - omegaN(i - 1,j,k)) / (2 * theCellSize);
        omegaGY(i,j,k) = (omegaN(i,j + 1,k) - omegaN(i,j - 1,k)) / (2 * theCellSize);
        omegaGZ(i,j,k) = (omegaN(i,j,k + 1) - omegaN(i,j,k - 1)) / (2 * theCellSize);

    };
}

void MACGrid::computeVorticityConfinement(float dt) {
    target.mU = mU;
    target.mV = mV;
    target.mW = mW;

    computeCentralVel();
    computeOmega();
    computeOmegaGradient();

    FOR_EACH_CELL{
        float3 N = float3(omegaGX(i, j, k), omegaGY(i, j, k), omegaGZ(i, j, k));
    //N = N / (N.Length() + 1e-20);
    N = normalize(N);
    float3 Omega = float3(omegaX(i, j, k), omegaY(i, j, k), omegaZ(i, j, k));
    // TODO: Check what ^ does
    float3 confF = cross(N, Omega) * theCellSize * theVorticityEpsilon;

    vorConfFX(i, j, k) = confF[0];
    vorConfFY(i, j, k) = confF[1];
    vorConfFZ(i, j, k) = confF[2];

    FOR_EACH_FACE_X
    {
        target.mU(i, j, k) += dt * (vorConfFX(i - 1, j, k) + vorConfFX(i, j, k)) / fluidDensity * 2.0f;
    };

    FOR_EACH_FACE_Y
    {

        target.mV(i, j, k) += dt * (vorConfFY(i, j - 1, k) + vorConfFY(i, j, k)) / fluidDensity * 2.0f;

    };

    FOR_EACH_FACE_Z
    {
        target.mW(i, j, k) += dt * (vorConfFZ(i, j, k - 1) + vorConfFZ(i, j, k)) / fluidDensity * 2.0f;
    };

    // Then save the result to our object
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
    };
}

void MACGrid::addExternalForces(float dt) {
    std::cout << "Adding external forces" << std::endl;
    computeBuoyancy(dt);
    //computeVorticityConfinement(dt);
}

void MACGrid::computeDivergence() {
    FOR_EACH_CELL{

        float velLowX = mU(i, j, k);
        float velHighX = mU(i + 1, j, k);
        float velLowY = mV(i, j, k);
        float velHighY = mV(i, j + 1, k);
        float velLowZ = mW(i, j, k);
        float velHighZ = mW(i, j, k + 1);
        if (i == 0) {
            velLowX = 0.0;
        }
        if (j == 0)
        {
            velLowY = 0.0;
        }
        if (k == 0)
        {
            velLowZ = 0.0;
        }
        if (i == theDim[MACGrid::X] - 1)
        {
            velHighX = 0.0;
        }
        if (j == theDim[MACGrid::Y] - 1)
        {
            velHighY = 0.0;
        }
        if (k == theDim[MACGrid::Z] - 1)
        {
            velHighZ = 0.0;
        }
        divergence(i,j,k) = -((velHighX - velLowX) + (velHighY - velLowY) + (velHighZ - velLowZ)) / theCellSize;

    }
}

void MACGrid::checkPressure(int& i, int& j, int& k, const GridData& p, float3& minPressurebound, float3& maxPressurebound) {
    if (isValidFace(MACGrid::X, i, j, k)) {
        if (i - 1 >= 0)
        {
            minPressurebound[0] = p(i - 1, j, k);
        }
        if (i < theDim[MACGrid::X])
        {
            maxPressurebound[0] = p(i, j, k);
        }
        if (i - 1 < 0)
        {
            minPressurebound[0] = maxPressurebound[0] - theBoundConstant * (mU(i, j, k) - 0);
        }
        if (i >= theDim[MACGrid::X])
        {
            maxPressurebound[0] = minPressurebound[0] + theBoundConstant * (mU(i, j, k) - 0);
        }
    }
    if (isValidFace(MACGrid::Y, i, j, k)) {
        if (j - 1 >= 0)
        {
            minPressurebound[1] = p(i, j - 1, k);
        }
        if (j < theDim[MACGrid::Y])
        {
            maxPressurebound[1] = p(i, j, k);
        }
        if (j - 1 < 0)
        {
            minPressurebound[1] = maxPressurebound[1] - theBoundConstant * (mV(i, j, k) - 0);
        }
        if (j >= theDim[MACGrid::Y])
        {
            maxPressurebound[1] = minPressurebound[1] + theBoundConstant * (mV(i, j, k) - 0);
        }
    }
    if (isValidFace(MACGrid::Z, i, j, k)) {
        if (k - 1 >= 0)
        {
            minPressurebound[2] = p(i, j, k - 1);
        }

        if (k < theDim[MACGrid::Z])
        {
            maxPressurebound[2] = p(i, j, k);
        }

        if (k - 1 < 0)
        {
            minPressurebound[2] = maxPressurebound[2] - theBoundConstant * (mW(i, j, k) - 0);
        }

        if (k >= theDim[MACGrid::Z])
        {
            maxPressurebound[2] = minPressurebound[2] + theBoundConstant * (mW(i, j, k) - 0);
        }
    }
}

void MACGrid::computeBound() {
    FOR_EACH_FACE_X{
        if (i == 0 || i == theDim[MACGrid::X]) {
            mU(i, j, k) = 0.0f;
        }
    };
    FOR_EACH_FACE_Y{
        if (j == 0 || j == theDim[MACGrid::Y]) {
            mV(i, j, k) = 0.0f;
        }
    };
    FOR_EACH_FACE_Z{
        if (k == 0 || k == theDim[MACGrid::Z]) {
            mW(i, j, k) = 0.0f;
        }
    };

    FOR_EACH_CELL{
        if (i == 0 || i == theDim[MACGrid::X] - 1 ||
            j == 0 || j == theDim[MACGrid::Y] - 1 ||
            k == 0 || k == theDim[MACGrid::Z] - 1) {
            mP(i, j, k) = 0.0f;
        }
    }
}

void MACGrid::project_simple(float dt) {
    GridData p, p_new;
    p.initialize();
    p_new.initialize();

    const float scale = dt / (theAirDensity * theCellSize);
    // Compute divergence
    computeDivergence();

    // Solve Poisson Equation using Gauss Siedel
    for (int iter = 0; iter < 100; ++iter) {
        FOR_EACH_CELL{
            float sum = p(i - 1, j, k) + p(i + 1, j, k) +
                        p(i, j - 1, k) + p(i, j + 1, k) +
                        p(i, j, k - 1) + p(i, j, k + 1);
        p(i, j, k) = (sum - theCellSize * divergence(i, j, k)) / 6.0f;
        }
        computeBound();
    }

    // Apply pressure gradient to update velocities
    FOR_EACH_FACE_X{
        if (i > 0 && i < theDim[MACGrid::X]) {
            mU(i, j, k) -= scale * (p(i, j, k) - p(i - 1, j, k));
        }
    }
        FOR_EACH_FACE_Y{
            if (j > 0 && j < theDim[MACGrid::Y]) {
                mV(i, j, k) -= scale * (p(i, j, k) - p(i, j - 1, k));
            }
    }
        FOR_EACH_FACE_Z{
            if (k > 0 && k < theDim[MACGrid::Z]) {
                mW(i, j, k) -= scale * (p(i, j, k) - p(i, j, k - 1));
            }
    }
    computeBound();
}

void MACGrid::project(float dt) {
    target.mU = mU;
    target.mV = mV;
    target.mW = mW;
    target.mP = mP;
    int size = theDim[MACGrid::X] * theDim[MACGrid::Y] * theDim[MACGrid::Z];
    SparseMatrix<float> A(size, size);

    A.setZero();
    VectorXd b(size);
    float pho = 4.0f;
    const float constant = theAirDensity * (theCellSize * theCellSize) / dt;

    GridData p;
    p.initialize();

    divergence.initialize();
    computeDivergence();
    preconditionedConjugateGradient(AMatrix, p, divergence, 1000, 0.01f);
    FOR_EACH_CELL
    {
        p(i,j,k) *= constant;
        target.mP(i,j,k) = p(i,j,k);
    }

    mP = target.mP;

    FOR_EACH_FACE_X
    {
        float3 minPressurebound = float3(0.0, 0.0, 0.0);
        float3 maxPressurebound = float3(0.0, 0.0, 0.0);
        checkPressure(i, j, k, p, minPressurebound, maxPressurebound);
        target.mU(i,j,k) = mU(i,j,k) - (dt / theAirDensity) * (maxPressurebound[0] - minPressurebound[0]) / theCellSize;
    };
    FOR_EACH_FACE_Y
    {
        float3 minPressurebound = float3(0.0, 0.0, 0.0);
        float3 maxPressurebound = float3(0.0, 0.0, 0.0);
        checkPressure(i, j, k, p, minPressurebound, maxPressurebound);
        target.mV(i,j,k) = mV(i,j,k) - (dt / theAirDensity) * (maxPressurebound[1] - minPressurebound[1]) / theCellSize;
    };
    FOR_EACH_FACE_Z
    {
        float3 minPressurebound = float3(0.0, 0.0, 0.0);
        float3 maxPressurebound = float3(0.0, 0.0, 0.0);
        checkPressure(i, j, k, p, minPressurebound, maxPressurebound);
        target.mW(i,j,k) = mW(i,j,k) - (dt / theAirDensity) * (maxPressurebound[2] - minPressurebound[2]) / theCellSize;

    };

    computeBound();
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;

    // Check velocities along wall if debug
    //FOR_EACH_CELL{
    //    // Construct the vector of divergences d:
    //     float velLowX = mU(i,j,k);
    //     float velHighX = mU(i + 1,j,k);
    //     float velLowY = mV(i,j,k);
    //     float velHighY = mV(i,j + 1,k);
    //     float velLowZ = mW(i,j,k);
    //     float velHighZ = mW(i,j,k + 1);
    //     float divergence = ((velHighX - velLowX) + (velHighY - velLowY) + (velHighZ - velLowZ)) / theCellSize;
    //     if (abs(divergence) > 0.02) {
    //         std::cout << "WARNING: Divergent! ";
    //         std::cout << "Divergence: " << divergence;
    //         std::cout << "Cell: " << i << ", " << j << ", " << k;
    //     }
    //}
}

float3 MACGrid::getVelocity(const float3& pt) {
    float3 vel = { getVelocityX(pt), getVelocityY(pt), getVelocityZ(pt) };
    return vel;
}

float MACGrid::getVelocityX(const float3& pt)
{
    return mU.trilinear(pt);

    //return mU.interpolate(pt);
}

float MACGrid::getVelocityY(const float3& pt)
{
    return mV.trilinear(pt);
    //return mV.interpolate(pt);
}

float MACGrid::getVelocityZ(const float3& pt)
{
    return mW.trilinear(pt);
    //return mW.interpolate(pt);
}

float MACGrid::getTemperature(const float3& pt)
{
    return mT.trilinear(pt);
    //return mT.interpolate(pt);
}

float MACGrid::getDensity(const float3& pt)
{
    return mD.trilinear(pt);
    //return mD.interpolate(pt);
}

float3 MACGrid::getCenter(int i, int j, int k) {
    float xstart = theCellSize / 2.0;
    float ystart = theCellSize / 2.0;
    float zstart = theCellSize / 2.0;

    float x = xstart + i * theCellSize;
    float y = ystart + j * theCellSize;
    float z = zstart + k * theCellSize;
    return float3(x, y, z);
}

float3 MACGrid::getRewoundPosition(const float3& currentPosition, const float dt) {
    /*
    // EULER (RK1):
    float3 currentVelocity = getVelocity(currentPosition);
    float3 rewoundPosition = currentPosition - currentVelocity * dt;
    float3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
    return clippedRewoundPosition;
    */

    // HEUN / MODIFIED EULER (RK2):
    float3 currentVelocity = getVelocity(currentPosition);
    float3 rewoundPosition = currentPosition - currentVelocity * dt;
    float3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
    // Keep going...
    float3 rewoundVelocity = getVelocity(clippedRewoundPosition);
    float3 averageVelocity = (currentVelocity + rewoundVelocity) / 2.0f;
    float3 betterRewoundPosition = currentPosition - averageVelocity * dt;
    float3 clippedBetterRewoundPosition = clipToGrid(betterRewoundPosition, currentPosition);
    return clippedBetterRewoundPosition;
}

float3 MACGrid::clipToGrid(const float3& outsidePoint, const float3& insidePoint) {
    float3 clippedPoint = outsidePoint;

    for (int i = 0; i < 3; i++) {
        if (clippedPoint[i] < 0) {
            float3 distance = clippedPoint - insidePoint;
            float newDistanceI = 0 - insidePoint[i];
            float ratio = newDistanceI / distance[i];
            clippedPoint = insidePoint + distance * ratio;
        }
        if (clippedPoint[i] > getSize(i)) {
            float3 distance = clippedPoint - insidePoint;
            float newDistanceI = getSize(i) - insidePoint[i];
            float ratio = newDistanceI / distance[i];
            clippedPoint = insidePoint + distance * ratio;
        }
    }
    // Ensure point lies within grid;
    if (clippedPoint[0] < 0 || clippedPoint[1] < 0 || clippedPoint[2] < 0 || clippedPoint[0] > getSize(0) || clippedPoint[1] > getSize(1) || clippedPoint[2] > getSize(2)) {
        std::cout << "WARNING: Clipped point is outside grid!";
    }

    return clippedPoint;
}

float MACGrid::getSize(int dimension) {
    return theDim[dimension] * theCellSize;
}

int MACGrid::getCellIndex(int i, int j, int k) {
    return i * j * theDim[MACGrid::X] + k * theDim[MACGrid::Y] * theDim[MACGrid::X];
}

void MACGrid::getCellIndexReverse(int idx, int& i, int& j, int& k) {
    i = idx / theDim[MACGrid::X];
    j = idx / (theDim[MACGrid::X] * theDim[MACGrid::Y]);
    k = idx / (theDim[MACGrid::X] * theDim[MACGrid::Y] * theDim[MACGrid::Z]);
}


int MACGrid::getNumberOfCells()
{
    return theDim[MACGrid::X] * theDim[MACGrid::Y] * theDim[MACGrid::Z];
}


bool MACGrid::isValidCell(int i, int j, int k)
{
    if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
        return false;
    }
    if (i < 0 || j < 0 || k < 0) {
        return false;
    }
    return true;
}

bool MACGrid::isValidFace(int dimension, int i, int j, int k)
{
    if (dimension == 0) {
        if (i > theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
            return false;
        }
    }
    else if (dimension == 1) {
        if (i >= theDim[MACGrid::X] || j > theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
            return false;
        }
    }
    else if (dimension == 2) {
        if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k > theDim[MACGrid::Z]) {
            return false;
        }
    }
    if (i < 0 || j < 0 || k < 0) {
        return false;
    }
    return true;
}

float3 MACGrid::getFacePosition(int dimension, int i, int j, int k) {
    if (dimension == 0) {
        return float3(i * theCellSize, (j + 0.5) * theCellSize, (k + 0.5) * theCellSize);
    }
    else if (dimension == 1) {
        return float3((i + 0.5) * theCellSize, j * theCellSize, (k + 0.5) * theCellSize);
    }
    else if (dimension == 2) {
        return float3((i + 0.5) * theCellSize, (j + 0.5) * theCellSize, k * theCellSize);
    }

    return float3(0, 0, 0); //???
}

void MACGrid::calculateAMatrix() {

    FOR_EACH_CELL{

        int numFluidNeighbors = 0;
        if (i - 1 >= 0) {
            AMatrix.plusI(i - 1,j,k) = -1;
            numFluidNeighbors++;
        }
        if (i + 1 < theDim[MACGrid::X]) {
            AMatrix.plusI(i,j,k) = -1;
            numFluidNeighbors++;
        }
        if (j - 1 >= 0) {
            AMatrix.plusJ(i,j - 1,k) = -1;
            numFluidNeighbors++;
        }
        if (j + 1 < theDim[MACGrid::Y]) {
            AMatrix.plusJ(i,j,k) = -1;
            numFluidNeighbors++;
        }
        if (k - 1 >= 0) {
            AMatrix.plusK(i,j,k - 1) = -1;
            numFluidNeighbors++;
        }
        if (k + 1 < theDim[MACGrid::Z]) {
            AMatrix.plusK(i,j,k) = -1;
            numFluidNeighbors++;
        }
        // Set the diagonal:
        AMatrix.diag(i,j,k) = numFluidNeighbors;
    }
}

bool MACGrid::preconditionedConjugateGradient(const GridDataMatrix& A, GridData& p, const GridData& d, int maxIterations, float tolerance) {
    FOR_EACH_CELL
    {
        p(i,j,k) = 0.0; // Initial guess p = 0.
    }

    GridData r = d; // Residual vector.
    GridData z;
    z.initialize();
    applyPreconditioner(r, A, z); // aux vector

    GridData s = z;
    // TODO: check other dot product functions are called correctly
    float sigma = dotProduct(z, r);
    bool converged = false;

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        float rho = sigma;
        apply(A, s, z);
        float alpha = rho / dotProduct(z, s);

        GridData alphaTimesS;
        alphaTimesS.initialize();
        multiply(alpha, s, alphaTimesS);
        add(p, alphaTimesS, p);
        // aka: p += alpha * s

        GridData alphaTimesZ;
        alphaTimesZ.initialize();
        multiply(alpha, z, alphaTimesZ);
        subtract(r, alphaTimesZ, r);
        // aka: r -= alpha * z

        if (maxMagnitude(r) <= tolerance) {
            converged = true;
        }

        applyPreconditioner(r, A, z);

        float sigmaNew = dotProduct(z, r);
        float beta = sigmaNew / rho;

        GridData betaTimesS;
        betaTimesS.initialize();
        multiply(beta, s, betaTimesS);
        add(z, betaTimesS, s);
        // aka: s = z + beta * s

        sigma = sigmaNew;
    }
    if (!converged) {
        std::cout << "PCG did not converge";
    }
    return converged;
}

void MACGrid::calculatePreconditioner(const GridDataMatrix& A) {
    precon.initialize();
    const float tau = 0.97;
    FOR_EACH_CELL
    {
        {
            float Aii = A.plusI(i - 1,j,k) * precon(i - 1,j,k);
            float Ajj = A.plusJ(i,j - 1,k) * precon(i,j - 1,k);
            float Akk = A.plusK(i,j,k - 1) * precon(i,j,k - 1);
            float Aijk = Aii * Aii + Ajj * Ajj + Akk * Akk;

            float Aiii = Aii * (A.plusJ(i - 1,j,k) + A.plusK(i - 1,j,k)) * precon(i - 1,j,k);
            float Ajjj = Ajj * (A.plusI(i,j - 1,k) + A.plusK(i,j - 1,k)) * precon(i,j - 1,k);
            float Akkk = Akk * (A.plusI(i,j,k - 1) + A.plusJ(i,j,k - 1)) * precon(i,j,k - 1);
            float temp = Aijk + tau * (Aiii + Ajjj + Akkk);

            float e = A.diag(i,j,k) - temp;
            precon(i,j,k) = 1.0 / sqrt(e + 1e-30);
        }
    }
}

void MACGrid::applyPreconditioner(const GridData& r, const GridDataMatrix& A, GridData& z) {

    if (1) {

        // APPLY THE PRECONDITIONER:
        // Solve Lq = r for q:
        GridData q;
        q.initialize();
        FOR_EACH_CELL{
            //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
            float t = r(i, j, k) - A.plusI(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k)
                       - A.plusJ(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k)
                       - A.plusK(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
            q(i, j, k) = t * precon(i, j, k);
            //}
        }
            // Solve L^Tz = q for z:
            FOR_EACH_CELL_REVERSE{
            //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
            float t = q(i, j, k) - A.plusI(i, j, k) * precon(i, j, k) * z(i + 1, j, k)
                       - A.plusJ(i, j, k) * precon(i, j, k) * z(i, j + 1, k)
                       - A.plusK(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
            z(i, j, k) = t * precon(i, j, k);
            //}
        }
    }
    else {
        // Unpreconditioned CG: Bypass preconditioner:
        z = r;
        return;
    }
}

float MACGrid::dotProduct(const GridData& vector1, const GridData& vector2) {

    float result = 0.0f;

    FOR_EACH_CELL{
        result += vector1(i,j,k) * vector2(i,j,k);
    }

    return result;
}

void MACGrid::add(const GridData& vector1, const GridData& vector2, GridData& result) {


    FOR_EACH_CELL{
        if (vector1(i,j,k) != NAN && vector2(i,j,k) != NAN)
            result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
        else
            cout << i << "," << j << "" << k << endl;
    }

}

void MACGrid::subtract(const GridData& vector1, const GridData& vector2, GridData& result) {

    FOR_EACH_CELL{
        result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
    }

}

void MACGrid::multiply(const float scalar, const GridData& vector, GridData& result) {

    FOR_EACH_CELL{
        result(i,j,k) = scalar * vector(i,j,k);
    }

}

float MACGrid::maxMagnitude(const GridData& vector) {

    float result = 0.0;

    FOR_EACH_CELL{
        if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
    }

    return result;
}

void MACGrid::apply(const GridDataMatrix& matrix, const GridData& vector, GridData& result) {

    FOR_EACH_CELL{ // For each row of the matrix.

        float diag = 0;
        float plusI = 0;
        float plusJ = 0;
        float plusK = 0;
        float minusI = 0;
        float minusJ = 0;
        float minusK = 0;

        diag = matrix.diag(i,j,k) * vector(i,j,k);
        if (isValidCell(i + 1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i + 1,j,k);
        if (isValidCell(i,j + 1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j + 1,k);
        if (isValidCell(i,j,k + 1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k + 1);
        if (isValidCell(i - 1,j,k)) minusI = matrix.plusI(i - 1,j,k) * vector(i - 1,j,k);
        if (isValidCell(i,j - 1,k)) minusJ = matrix.plusJ(i,j - 1,k) * vector(i,j - 1,k);
        if (isValidCell(i,j,k - 1)) minusK = matrix.plusK(i,j,k - 1) * vector(i,j,k - 1);

        result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
    }

}

void MACGrid::saveSmoke(const char* fileName) {
    std::ofstream fileOut(fileName);
    if (fileOut.is_open()) {
        FOR_EACH_CELL{
            fileOut << mD(i,j,k) << std::endl;
        }
        fileOut.close();
    }
}

void MACGrid::saveParticles(std::string filename) {

}

void MACGrid::saveDensity(std::string filename) {

}
