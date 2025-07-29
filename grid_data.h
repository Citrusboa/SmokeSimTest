#pragma once
#ifndef GridData_H_
#define GridData_H_

#pragma warning(disable: 4244 4267 4996)

#include <vector>
#include "constants.h"

#include "linalg.h"
using namespace linalg::aliases;

#include <iostream>
#include <cfloat>
#include <chrono>

// GridData is container for storing data in a grid
// Columns indexed with i and increasing with x
// Rows indexed with j and increasing z
// Stacks are indexed with k and increasing y

// GridData world space goes from 90, 0, 0) to mMax, which is
// (theCellSize * theDim[0], theCellSize*theDim[1], theCellSize*theDim[2])

class GridData {
public:
	GridData();
	GridData(const GridData& orig);
	virtual ~GridData();
	virtual GridData& operator=(const GridData& orig);

	virtual void initialize(float defaultVal = 0.0f);

	virtual float& operator()(int i, int j, int k);
	virtual const float operator()(int i, int j, int k) const;
	virtual int3 getDim() const;
	virtual float trilinear(const float3& pt);
	virtual float interpolate(const float3& pt);

	float CINT(float qi_minus1, float qi, float qi_plus1, float qi_plus2, float x) const;

	std::vector<float>& data();

	virtual void getCell(const float3& pt, int& i, int& j, int& k);

protected:
	virtual float3 worldToSelf(const float3& pt) const;
	float defaultVal;
	float3 mMax;
	std::vector<float> mData;
};

class GridDataX : public GridData {
public:
	GridDataX();
	virtual ~GridDataX();
	virtual void initialize(float defaultVal = 0.0f);
	virtual float& operator()(int i, int j, int k);
	virtual const float operator()(int i, int j, int k) const;
	virtual float3 worldToSelf(const float3& pt) const;
};

class GridDataY : public GridData {
public:
	GridDataY();
	virtual ~GridDataY();
	virtual void initialize(float defaultVal = 0.0f);
	virtual float& operator()(int i, int j, int k);
	virtual const float operator()(int i, int j, int k) const;
	virtual float3 worldToSelf(const float3& pt) const;
};

class GridDataZ : public GridData {
public:
	GridDataZ();
	virtual ~GridDataZ();
	virtual void initialize(float defaultVal = 0.0f);
	virtual float& operator()(int i, int j, int k);
	virtual const float operator()(int i, int j, int k) const;
	virtual float3 worldToSelf(const float3& pt) const;
};


#endif