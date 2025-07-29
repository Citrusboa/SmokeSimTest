#include "grid_data.h"

GridData::GridData() : defaultVal(0.0f), mMax(0.0f, 0.0f, 0.0f) {

}

GridData::GridData(const GridData& orig) : defaultVal(orig.defaultVal) {
	mData = orig.mData;
	mMax = orig.mMax;
}

GridData::~GridData() {

}

std::vector<float>& GridData::data() {
	return mData;
}

GridData& GridData::operator=(const GridData& orig) {
	if (this == &orig) {
		return *this;
	}
	defaultVal = orig.defaultVal;
	mData = orig.mData;
	mMax = orig.mMax;
	return *this;
}

void GridData::initialize(float defaultVal) {
	this->defaultVal = defaultVal;
	mMax[0] = theCellSize * theDim[0];
	mMax[1] = theCellSize * theDim[1];
	mMax[2] = theCellSize * theDim[2];
	mData.resize(theDim[0] * theDim[1] * theDim[2], false);
	std::fill(mData.begin(), mData.end(), this->defaultVal);
}

int3 GridData::getDim() const {
	return (int3)(mMax / theCellSize);
}

float& GridData::operator()(int i, int j, int k) {
	static float default1 = 0.0f;
	default1 = defaultVal;

	if (i < 0 || j < 0 || k < 0 ||
		i > theDim[0] - 1 ||
		j > theDim[1] - 1 ||
		k > theDim[2] - 1) {
		return default1;
	}
	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * theDim[2];
	return mData[col + row + stack];
}

const float GridData::operator()(int i, int j, int k) const {
	static float default1 = 0.0f;
	default1 = defaultVal;

	if (i < 0 || j < 0 || k < 0 ||
		i > theDim[0] - 1 ||
		j > theDim[1] - 1 ||
		k > theDim[2] - 1) {
		return default1;
	}
	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * theDim[2];
	return mData[col + row + stack];
}

void GridData::getCell(const float3& pt, int& i, int& j, int& k) {
	float3 pos = worldToSelf(pt);
	i = (int)(pos[0] / theCellSize);
	j = (int)(pos[1] / theCellSize);
	k = (int)(pos[2] / theCellSize);
}

float GridData::trilinear(const float3& pt) {
	int i, j, k;
	float3 pos = worldToSelf(pt);
	getCell(pt, i, j, k);

	float scale = 1.0f / theCellSize;
	float fx = scale * (pos[0] - i * theCellSize);
	float fy = scale * (pos[1] - j * theCellSize);
	float fz = scale * (pos[2] - k * theCellSize);

	int3 dims = getDim();
	i = std::min(std::max(0, i), dims[0] - 2);
	j = std::min(std::max(0, j), dims[1] - 2);
	k = std::min(std::max(0, k), dims[2] - 2);

	float val[2][2][2];
	for (int dx = 0; dx <= 1; ++dx) {
		for (int dy = 0; dy <= 1; ++dy) {
			for (int dz = 0; dz <= 1; ++dz) {
				val[dx][dy][dz] = (*this)(i + dx, j + dy, k + dz);
			}
		}
	}
	auto lerp = [](float a, float b, float t) {
		return a * (1 - t) + b * t;
		};
	float c00 = lerp(val[0][0][0], val[0][0][1], fz);
	float c01 = lerp(val[0][1][0], val[0][1][1], fz);
	float c10 = lerp(val[1][0][0], val[1][0][1], fz);
	float c11 = lerp(val[1][1][0], val[1][1][1], fz);

	float c0 = lerp(c00, c01, fy);
	float c1 = lerp(c10, c11, fy);

	float c = lerp(c0, c1, fx);
	return c;
}

float GridData::interpolate(const float3& pt) {
	// Sharp cubic interpolation

	int i, j, k;
	float3 pos = worldToSelf(pt);
	getCell(pt, i, j, k);

	float scale = 1.0f / theCellSize;
	float fractx = scale * (pos[0] - i * theCellSize);
	float fracty = scale * (pos[1] - j * theCellSize);
	float fractz = scale * (pos[2] - k * theCellSize);

	float t[4][4];
	float u[4];
	float f;

#define ONE 1
#define EVAL(a, b, c) (*this)(a, b, c)
	for (int x = -1; x <= 2; ++x) {
		for (int y = -1; y <= 2; ++y) {
			t[x + ONE][y + ONE] = CINT(EVAL(i + x, j + y, k - 1), EVAL(i + x, j + y, k + 0), EVAL(i + x, j + y, k + 1), EVAL(i + x, j + y, k + 2), fractz);
		}
	}
#undef EVAL
	for (int x = -1; x <= 2; ++x) {
		u[x + ONE] = CINT(t[x + ONE][-1 + ONE], t[x + ONE][0 + ONE], t[x + ONE][1 + ONE], t[x + ONE][2 + ONE], fracty);
	}
	f = CINT(u[-1 + ONE], u[0 + ONE], u[1 + ONE], u[2 + ONE], fractx);
#undef ONE
	return f;
}

float GridData::CINT(float qi_minus1, float qi, float qi_plus1, float qi_plus2, float x) const {
	// slopes
	float di = (qi_plus1 - qi_minus1) / 2.0f;
	float di_plus1 = (qi_plus2 - qi) / 2.0f;

	float delta_q = qi_plus1 - qi;

	// constrain slope
	if (delta_q > 0) {
		if (di < 0) di = 0.0f;
		if (di_plus1 < 0) di_plus1 = 0.0f;
	}
	else if (delta_q < 0) {
		if (di > 0) di = 0.0f;
		if (di_plus1 > 0) di_plus1 = 0.0f;
	}

	// Hermite cubic;
	float q_x = qi + di * x + (3.0f * delta_q - 2.0f * di - di_plus1) * (x * x) + (-2.0f * delta_q + di + di_plus1) * (x * x * x);

	return q_x;
}

float3 GridData::worldToSelf(const float3& pt) const {
	float3 out;
	out[0] = std::min(std::max(0.0f, pt[0] - theCellSize * 0.5f), mMax[0]);
	out[1] = std::min(std::max(0.0f, pt[1] - theCellSize * 0.5f), mMax[1]);
	out[2] = std::min(std::max(0.0f, pt[2] - theCellSize * 0.5f), mMax[2]);
	return out;
}


/// <summary>
/// Grid Data X Class
/// </summary>
GridDataX::GridDataX() : GridData() {

}

GridDataX::~GridDataX() {

}

void GridDataX::initialize(float defaultVal) {
	GridData::initialize(defaultVal);
	mMax[0] = theCellSize * (theDim[0] + 1);
	mMax[1] = theCellSize * theDim[1];
	mMax[2] = theCellSize * theDim[2];
	mData.resize((theDim[0] + 1) * theDim[1] * theDim[2], false);
	std::fill(mData.begin(), mData.end(), defaultVal);
}

float& GridDataX::operator()(int i, int j, int k) {
	static float default1 = 0.0f;
	default1 = this->defaultVal;

	if (i < 0 || i > theDim[0]) return default1;
	if (j < 0) j = 0.0f;
	if (j > theDim[1] - 1) j = theDim[1] - 1;
	if (k < 0) k = 0.0f;
	if (k > theDim[2] - 1) k = theDim[2] - 1;

	int col = i;
	int row = k * (theDim[0] + 1);
	int stack = j * (theDim[0] + 1) * theDim[2];
	return mData[stack + row + col];
}

const float GridDataX::operator()(int i, int j, int k) const {
	static float default1 = 0.0f;
	default1 = this->defaultVal;

	if (i < 0 || i > theDim[0]) return default1;
	if (j < 0) j = 0.0f;
	if (j > theDim[1] - 1) j = theDim[1] - 1;
	if (k < 0) k = 0.0f;
	if (k > theDim[2] - 1) k = theDim[2] - 1;

	int col = i;
	int row = k * (theDim[0] + 1);
	int stack = j * (theDim[0] + 1) * theDim[2];
	return mData[stack + row + col];
}

float3 GridDataX::worldToSelf(const float3& pt) const {
	float3 out;
	out[0] = std::min(std::max(0.0f, pt[0]), mMax[0]);
	out[1] = std::min(std::max(0.0f, pt[1] - theCellSize * 0.5f), mMax[1]);
	out[2] = std::min(std::max(0.0f, pt[2] - theCellSize * 0.5f), mMax[2]);
	return out;
}

/// <summary>
/// Grid Data Y Class
/// </summary>

GridDataY::GridDataY() : GridData() {

}

GridDataY::~GridDataY() {

}


void GridDataY::initialize(float defaultVal) {
	GridData::initialize(defaultVal);
	mMax[0] = theCellSize * theDim[0];
	mMax[1] = theCellSize * (theDim[1] + 1);
	mMax[2] = theCellSize * theDim[2];
	mData.resize(theDim[0] * (theDim[1] + 1) * theDim[2], false);
	std::fill(mData.begin(), mData.end(), this->defaultVal);
}

float& GridDataY::operator()(int i, int j, int k) {
	static float default1 = 0.0f;
	default1 = defaultVal;

	if (j < 0 || j > theDim[1]) return default1;

	if (i < 0) i = 0.0f;
	if (i > theDim[0] - 1) i = theDim[0] - 1;
	if (k < 0) k = 0.0f;
	if (k > theDim[2] - 1) k = theDim[2] - 1;

	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * theDim[2];
	return mData[stack + row + col];
}

const float GridDataY::operator()(int i, int j, int k) const {
	static float default1 = 0.0f;
	default1 = defaultVal;

	if (j < 0 || j > theDim[1]) return default1;

	if (i < 0) i = 0.0f;
	if (i > theDim[0] - 1) i = theDim[0] - 1;
	if (k < 0) k = 0.0f;
	if (k > theDim[2] - 1) k = theDim[2] - 1;

	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * theDim[2];
	return mData[stack + row + col];
}

float3 GridDataY::worldToSelf(const float3& pt) const {
	float3 out;
	out[0] = std::min(std::max(0.0f, pt[0] - theCellSize * 0.5f), mMax[0]);
	out[1] = std::min(std::max(0.0f, pt[1]), mMax[1]);
	out[2] = std::min(std::max(0.0f, pt[2] - theCellSize * 0.5f), mMax[2]);
	return out;
}

/// <summary>
/// Grid Data Z Class
/// </summary>
GridDataZ::GridDataZ() : GridData() {

}

GridDataZ::~GridDataZ() {

}

void GridDataZ::initialize(float defaultVal) {
	GridData::initialize(defaultVal);
	mMax[0] = theCellSize * theDim[0];
	mMax[1] = theCellSize * theDim[1];
	mMax[2] = theCellSize * (theDim[2] + 1);
	mData.resize(theDim[0] * theDim[1] * (theDim[2] + 1), false);
	std::fill(mData.begin(), mData.end(), this->defaultVal);
}

float& GridDataZ::operator()(int i, int j, int k) {
	static float default1 = 0.0f;
	default1 = this->defaultVal;

	if (k < 0 || k > theDim[2]) return default1;

	if (i < 0) i = 0;
	if (i > theDim[0] - 1) i = theDim[0] - 1;
	if (j < 0) j = 0;
	if (j > theDim[1] - 1) j = theDim[1] - 1;

	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * (theDim[2] + 1);

	return mData[stack + row + col];
}

const float GridDataZ::operator()(int i, int j, int k) const {
	static float default1 = 0.0f;
	default1 = this->defaultVal;

	if (k < 0 || k > theDim[2]) return default1;

	if (i < 0) i = 0;
	if (i > theDim[0] - 1) i = theDim[0] - 1;
	if (j < 0) j = 0;
	if (j > theDim[1] - 1) j = theDim[1] - 1;

	int col = i;
	int row = k * theDim[0];
	int stack = j * theDim[0] * (theDim[2] + 1);

	return mData[stack + row + col];
}

float3 GridDataZ::worldToSelf(const float3& pt) const {
	float3 out;
	out[0] = std::min(std::max(0.0f, pt[0] - theCellSize * 0.5f), mMax[0]);
	out[1] = std::min(std::max(0.0f, pt[1] - theCellSize * 0.5f), mMax[1]);
	out[2] = std::min(std::max(0.0f, pt[2]), mMax[2]);
	return out;
}
