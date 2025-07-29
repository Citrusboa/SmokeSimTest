#pragma once
#ifndef SMOKE_SIM_H_
#define SMOKE_SIM_H_

#include "mac_grid.h"

#include <fstream>

class SmokeSim {
public:
	SmokeSim();
	virtual ~SmokeSim();

	virtual void reset();
	virtual void step();
	virtual void setRecording(bool on, int width, int height);
	virtual bool isRecording();

	int getTotalFrames();

	std::vector<float3>& getParticles();
	std::vector<float3>& getParticlesVelocity();

protected:
	//virtual void drawAxes();
	//virtual void grabScreen();

	MACGrid mGrid;
	bool mRecordEnabled;
	int mFrameNum;
	int mTotalFrameNum;

	int recordWidth;
	int recordHeight;
};

SmokeSim::SmokeSim() : mFrameNum(0), mTotalFrameNum(0), mRecordEnabled(true) {
	reset();
}

SmokeSim::~SmokeSim() {

}

void SmokeSim::reset() {
	mGrid.resetSim();
	mTotalFrameNum = 0;
}

void SmokeSim::step() {
	std::cout << "=============================================================" << mTotalFrameNum << std::endl;

	mGrid.updateSources();
	// Advect
	mGrid.advectVelocity(dt);
	mGrid.addExternalForces(dt);
	//mGrid.project(dt);
	mGrid.project_simple(dt);

	mGrid.advectTemperature(dt);
	mGrid.advectDensity(dt);
	//mGrid.advectParticles(dt);

	mTotalFrameNum++;
}

void SmokeSim::setRecording(bool on, int width, int height) {
	if (on && !mRecordEnabled) {
		mFrameNum = 0;
	}

	mRecordEnabled = on;

	recordWidth = width;
	recordHeight = height;
}

bool SmokeSim::isRecording() {
	return mRecordEnabled;
}

int SmokeSim::getTotalFrames() {
	return mTotalFrameNum;
}

std::vector<float3>& SmokeSim::getParticles() {
	return mGrid.getParticles();
}

std::vector<float3>& SmokeSim::getParticlesVelocity() {
	return mGrid.getParticlesVelocity();
}

#endif