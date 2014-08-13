#include "core.h"

Flip::Flip()
{
	Init();
}

void Flip::Reset()
{
	Object::Reset();
	/*FIXME: free memory and reinitialize*/
	//Init();
}

void Flip::Show()
{
	
	//draw cube
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	GLfloat vertices[][3] = { { 0.0f, 0.0f, 0.0f },
		{ 1.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 0.0f },
		{ 0.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f },
		{ 1.0f, 0.0f, 1.0f }, { 1.0f, 1.0f, 1.0f },
		{ 0.0f, 1.0f, 1.0f } };              // 頂点座標値 
	int faces[][4] = { { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 4, 5, 6, 7 },
		{ 0, 4, 7, 3 }, { 0, 1, 5, 4 }, { 0, 3, 2, 1 } };
	// 各面の頂点番号列
	GLfloat colors[][3] = { { 0.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f },
		{ 1.0f, 1.0f, 0.0f }, { 0.0f, 0.5f, 0.5f },
		{ 0.5f, 0.0f, 0.5f }, { 0.5f, 0.5f, 0.0f } };
	// 各面の描画色
	glBegin(GL_QUADS);                  // 四角形描画開始
	for (int i = 0; i < 6; i++) {
		//	glBegin(GL_QUAD_STRIP);                  // 四角形描画開始
		//	glColor3fv(colors[i]);
		for (int j = 0; j < 4; j++)
			glVertex3fv(vertices[faces[i][j]]);
		//	glEnd();
	}                                         // 四角形頂点列の指定
	glEnd();                               // 四角形描画終了
	
	glBegin(GL_POINTS);
	//glVertex3d(0.0, 0.0, 0.0);
	
	for (int i = 0 ; i < particles_.size() ; i++) {
		GLdouble x = particles_[i]->position_.x_;
		GLdouble y = particles_[i]->position_.y_;
		GLdouble z = particles_[i]->position_.z_;
		glVertex3d(x, y, z);
	}
	
	glEnd();
	
}

void Flip::_InsertToCell(Particle* particle)
{
	int i = (int)(particle->position_.x_ / grid_.cellSize_);
	int j = (int)(particle->position_.y_ / grid_.cellSize_);
	int k = (int)(particle->position_.z_ / grid_.cellSize_);

#if 0
	std::cout << i << " " << j << " " << k  << std::endl;
#endif
	grid_.cells_[i][j][k].particles_.push_back(particle);
}

void Flip::_GenerateFluidBlock(GLdouble xL, GLdouble xU, GLdouble yL, GLdouble yU, 
			GLdouble zL, GLdouble zU)
{
	GLdouble particleIntervalX = simualateAreaX_ / ((PARTICLE_PER_CELL)*grid_.nX_);
	GLdouble particleIntervalY = simualateAreaY_ / ((PARTICLE_PER_CELL)*grid_.nY_);
	GLdouble particleIntervalZ = simualateAreaZ_ / ((PARTICLE_PER_CELL)*grid_.nZ_);
	GLdouble startX = 1.25 * particleIntervalX;
	GLdouble startY = 1.25 * particleIntervalY;
	GLdouble startZ = 1.25 * particleIntervalZ;

	for (GLdouble i = startX ; i <= simualateAreaX_-startX;
			i += particleIntervalX) {
		for (GLdouble j = startY ; j <= simualateAreaY_-startY; 
				j += particleIntervalY) {
			for (GLdouble k = startZ ; k <= simualateAreaZ_-startZ;
					k += particleIntervalZ) {
				if (i > xL && i < xU && j > yL && j < yU && k > zL && k < zU) {
					Particle* particle = new Particle;
					particle->position_.x_ = i;
					particle->position_.y_ = j;
					particle->position_.z_ = k;
					particles_.push_back(particle);
#if 0
					std::cout << i << " " << j << " " << k << std::endl;
#endif
					_InsertToCell(particle);
				}
			}
		}
	}

	_MarkFluid();
}


void Flip::_MarkFluid()
{
	LOOP_FOR_CELLS(grid_.nX_, grid_.nY_, grid_.nZ_) {
		if(grid_.cells_[i][j][k].type_ != SOLID) {
			grid_.cells_[i][j][k].type_ = AIR;
		}
		if(!grid_.cells_[i][j][k].particles_.empty() &&
				grid_.cells_[i][j][k].type_ != SOLID) {
#if 0
			std::cout << "is fluid" << std::endl;
#endif
			grid_.cells_[i][j][k].type_ = FLUID;
		}
	}END_LOOP
}

void Flip::_InitGrid()
{
	grid_.nX_ = grid_.nY_ = grid_.nZ_ = N;

	Allocator<Cell> cellAlloc;
	Allocator<GLdouble> velAlloc;

	grid_.cells_ = cellAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_);
	grid_.velX_ = velAlloc.Alloc3D(grid_.nX_+1, grid_.nY_, grid_.nZ_);
	grid_.velY_ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_+1, grid_.nZ_);
	grid_.velZ_ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_+1);
	grid_.cellSize_ = simualateAreaX_ / grid_.nX_;

	savedGrid_.cells_ = cellAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_);
	savedGrid_.velX_ = velAlloc.Alloc3D(grid_.nX_+1, grid_.nY_, grid_.nZ_);
	savedGrid_.velY_ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_+1, grid_.nZ_);
	savedGrid_.velZ_ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_+1);
	savedGrid_.cellSize_ = simualateAreaX_ / grid_.nX_;
}

void Flip::_InitBoundary()
{
	//set all cell as air
	LOOP_FOR_CELLS(grid_.nX_,grid_.nY_,grid_.nZ_) {
		grid_.cells_[i][j][k].type_ = AIR;
	} END_LOOP

	GLdouble rate = 0.01;
	//init wall and wall normals
	boundary = new GridObject;
	Allocator<Normal> normAlloc;
	boundary->normals = normAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_);
	for (int i = 0 ; i < grid_.nX_ ; i++) {
		for (int j = 0 ; j < grid_.nY_ ; j++) {
			grid_.cells_[i][j][0].type_ = SOLID;
			grid_.cells_[i][j][grid_.nZ_-1].type_ = SOLID;

			boundary->normals[i][j][0].x_ += 0;
			boundary->normals[i][j][0].y_ += 0;
			boundary->normals[i][j][0].z_ += rate * grid_.cellSize_;

			boundary->normals[i][j][grid_.nZ_-1].x_ += 0;
			boundary->normals[i][j][grid_.nZ_-1].y_ += 0;
			boundary->normals[i][j][grid_.nZ_-1].z_ += -rate * grid_.cellSize_;
#if 0
			std::cout << i << " " << j <<" " << std::endl;
#endif
		}
	}

	for (int j = 0 ; j < grid_.nX_ ; j++) {
		for (int k = 0 ; k < grid_.nY_ ; k++) {
			grid_.cells_[0][j][k].type_ = SOLID;
			grid_.cells_[grid_.nX_-1][j][k].type_ = SOLID;

			boundary->normals[0][j][k].x_ += rate * grid_.cellSize_;
			boundary->normals[0][j][k].y_ += 0;
			boundary->normals[0][j][k].z_ += 0;

			boundary->normals[grid_.nX_-1][j][k].x_ += -rate * grid_.cellSize_;
			boundary->normals[grid_.nX_-1][j][k].y_ += 0;
			boundary->normals[grid_.nX_-1][j][k].z_ += 0;
		}
	}

	for (int i = 0 ; i < grid_.nX_ ; i++) {
		for (int k = 0 ; k < grid_.nY_ ; k++) {
			grid_.cells_[i][0][k].type_ = SOLID;
			grid_.cells_[i][grid_.nY_-1][k].type_ = SOLID;

			boundary->normals[i][0][k].x_ += 0;
			boundary->normals[i][0][k].y_ += rate * grid_.cellSize_;
			boundary->normals[i][0][k].z_ += 0;

			boundary->normals[i][grid_.nY_-1][k].x_ += 0;
			boundary->normals[i][grid_.nY_-1][k].y_ += -rate * grid_.cellSize_;
			boundary->normals[i][grid_.nY_-1][k].z_ += 0;
		}
	}
#if 0
	std::cout << "end boundary" << std::endl;
#endif

}

void Flip::Init()
{
	step_ = 0;
	simualateAreaX_ = simualateAreaY_ = simualateAreaZ_ = 1.0;

	//allocate memory
	_InitGrid();

	//init boundary and mark solid
	_InitBoundary();


	_GenerateFluidBlock(0.2, 0.4, 0.10, 0.6, 0.2, 0.8);
	_GenerateFluidBlock(WALL_THICKNESS, 1.0-WALL_THICKNESS, WALL_THICKNESS, 
			0.10, WALL_THICKNESS, 1.0-WALL_THICKNESS);
	//_GenerateFluidBlock(WALL_THICKNESS, 1.0-WALL_THICKNESS, WALL_THICKNESS, 
	//		0.6, WALL_THICKNESS, 1.0-WALL_THICKNESS);
#if 0
	LOOP_FOR_CELLS(grid_.nX_,grid_.nY_,grid_.nZ_) {
		if(grid_.cells_[i][j][k].type_ == FLUID) {
			std::cout << "is fluid" << std::endl;
		}
	} END_LOOP
#endif

#if 0
	LOOP_FOR_CELLS(grid_.nX_, grid_.nY_, grid_.nZ_) {
		if(grid_.cells_[i][j][k].type_ != SOLID)
			std::cout << i << " " << j << " " << k << " is not solid"<< std::endl;
	} END_LOOP
#endif

#if 0
	//init air
	LOOP_FOR_CELLS(grid_.nX_, grid_.nY_, grid_.nZ_) {
		if (grid_.cells_[i][j][k].type_ != SOLID && grid_.cells_[i][j][k].type_ != FLUID)
			grid_.cells_[i][j][k].type_ = AIR;
	} END_LOOP
#endif

}

void Flip::EnforceBoundary()
{
	LOOP_FOR_CELLS(grid_.nX_,grid_.nY_,grid_.nZ_) {
		if(grid_.cells_[i][j][k].type_ == SOLID) {
			grid_.velX_[i][j][k] = usolid;
			grid_.velX_[i+1][j][k] = usolid;
			grid_.velY_[i][j][k] = usolid;
			grid_.velY_[i][j+1][k] = usolid;
			grid_.velZ_[i][j][k] = usolid;
			grid_.velZ_[i][j][k+1] = usolid;
		}
	} END_LOOP
}

void Flip::Project()
{
	GLdouble dx = grid_.cellSize_;
	GLdouble div = 0;


	//compute right-hand side divergence
	int index;	//linear index
	int rhsLength = MAT2LINEAR(grid_.nX_-1,grid_.nY_-1,grid_.nZ_-1) + 1;
	Eigen::VectorXd rightHand(rhsLength);
	rightHand.setZero();
	LOOP_FOR_CELLS(grid_.nX_, grid_.nY_, grid_.nZ_) {
		if (grid_.cells_[i][j][k].type_ == FLUID) {
			index = MAT2LINEAR(i,j,k);
			div = - ( (grid_.velX_[i+1][j][k] - grid_.velX_[i][j][k])
				+ (grid_.velY_[i][j+1][k] - grid_.velY_[i][j][k]) 
				+ (grid_.velZ_[i][j][k+1] - grid_.velZ_[i][j][k]) ) / dx;
			rightHand(index) = div;
		}
	} END_LOOP
	
	/*FIXME: modify divergence to account for solid velocities*/
	LOOP_FOR_CELLS(grid_.nX_, grid_.nY_, grid_.nZ_) {
		index = MAT2LINEAR(i,j,k);
		if (grid_.cells_[i][j][k].type_ == FLUID) {
			if(grid_.cells_[i-1][j][k].type_ == SOLID) {
				rightHand(index) -= (grid_.velX_[i][j][k] - usolid) / dx;
			}
			if(grid_.cells_[i+1][j][k].type_ == SOLID) {
				rightHand(index) += (grid_.velX_[i+1][j][k] - usolid) / dx;
			}

			if(grid_.cells_[i][j-1][k].type_ == SOLID) {
				rightHand(index) -= (grid_.velY_[i][j][k] - usolid) / dx;
			}
			if(grid_.cells_[i][j+1][k].type_ == SOLID) {
				rightHand(index) += (grid_.velY_[i][j+1][k] - usolid) / dx;
			}

			if(grid_.cells_[i][j][k-1].type_ == SOLID) {
				rightHand(index) -= (grid_.velZ_[i][j][k] - usolid) / dx;
			}
			if(grid_.cells_[i][j][k+1].type_ == SOLID) {
				rightHand(index) += (grid_.velY_[i][j][k+1] - usolid) / dx;
			}
		}
	} END_LOOP
#if 1
	std::cout << std::endl << "rightHand: "<< std::endl << rightHand.transpose() << std::endl;
#endif


	//set coefficient matrix using sparse matrix
	std::vector<DoubleTriplet> coefficients;
	DoubleSpaMat coeffMat(rhsLength, rhsLength);
	coeffMat.setZero();
	int index0, indexIPlus1, indexJPlus1, indexKPlus1,
		indexIMinus1, indexJMinus1, indexKMinus1;	//linear index
	GLdouble diagIJK, iPlus1, jPlus1, kPlus1, 
			 iMinus1, jMinus1, kMinus1;	//value
	GLdouble scale = DT / (RHO*dx*dx);

#if 1
	std::cout << std::endl << "coefficients:" << std::endl;
#endif
	LOOP_FOR_CELLS(grid_.nX_, grid_.nY_, grid_.nZ_) {

		index0 = MAT2LINEAR(i,j,k);
#if 0
		std::cout << index0 << std::endl;
#endif
		indexIPlus1 = MAT2LINEAR(i+1,j,k);
		indexJPlus1 = MAT2LINEAR(i,j+1,k);
		indexKPlus1 = MAT2LINEAR(i,j,k+1);
		indexIMinus1 = MAT2LINEAR(i-1,j,k);
		indexJMinus1 = MAT2LINEAR(i,j-1,k);
		indexKMinus1 = MAT2LINEAR(i,j,k-1);
		diagIJK = iPlus1 = jPlus1 = kPlus1 
			= iMinus1 = jMinus1 = kMinus1 = 0.0;

		if(grid_.cells_[i][j][k].type_ != FLUID)
			continue;

		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i+1][j][k].type_ == FLUID) {
			diagIJK += scale;
			iPlus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i+1][j][k].type_ == AIR) {
			diagIJK += scale;
		}


		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i][j+1][k].type_ == FLUID) {
			diagIJK += scale;
			jPlus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i][j+1][k].type_ == AIR) {
			diagIJK += scale;
		}


		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i][j][k+1].type_ == FLUID) {
			diagIJK += scale;
			kPlus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i][j][k+1].type_ == AIR) {
			diagIJK += scale;
		}


		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i-1][j][k].type_ == FLUID) {
			diagIJK += scale;
			iMinus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i-1][j][k].type_ == AIR) {
			diagIJK += scale;
		}

		
		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i][j-1][k].type_ == FLUID) {
			diagIJK += scale;
			jMinus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i][j-1][k].type_ == AIR) {
			diagIJK += scale;
		}


		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i][j][k-1].type_ == FLUID) {
			diagIJK += scale;
			kMinus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i][j][k-1].type_ == AIR) {
			diagIJK += scale;
		}

		coefficients.push_back(DoubleTriplet(index0, index0, diagIJK) );
		coefficients.push_back(DoubleTriplet(index0, indexIPlus1, iPlus1) );
		coefficients.push_back(DoubleTriplet(index0, indexJPlus1, jPlus1) );
		coefficients.push_back(DoubleTriplet(index0, indexKPlus1, kPlus1) );
		coefficients.push_back(DoubleTriplet(index0, indexIMinus1, iMinus1) );
		coefficients.push_back(DoubleTriplet(index0, indexJMinus1, jMinus1) );
		coefficients.push_back(DoubleTriplet(index0, indexKMinus1, kMinus1) );
#if 1
		std::cout << index0 << " " << indexIPlus1 << " " << indexJPlus1 << 
			" " << indexKPlus1 << " " << indexIMinus1  << " " 
			<< " " << indexJMinus1 << " " << indexKMinus1 << std::endl;
		std::cout << diagIJK << " " << iPlus1 << " " << jPlus1 << " " << kPlus1 << 
			" " << iMinus1 << " " << jMinus1 << " " << kMinus1 << std::endl;
#endif
	} END_LOOP
	//set coefficient matrix
	coeffMat.setFromTriplets(coefficients.begin(), coefficients.end() );
	
	Eigen::VectorXd pressure(rhsLength);
	Eigen::ConjugateGradient<DoubleSpaMat> solver;
	solver.compute(coeffMat);
	if(solver.info() != Eigen::Success) {
		std::cerr << "ConjugateGradient solver error!" << std::endl;
		exit(-1);
	}
	pressure = solver.solve(rightHand);
#if 1
	std::cout << std::endl << "pressure:"<< std::endl << pressure.transpose() << std::endl;
#endif

	for(int l = 0 ; l < pressure.size() ; l++) {
		//change linear index to 3D index
		int tmp = l / grid_.nY_;
		int i = tmp / grid_.nX_;
		int j = tmp % grid_.nX_;
		int k = l - (i*grid_.nX_*grid_.nY_+j*grid_.nY_);
#if 0
		std::cout << i << " " << j << " " << k << std::endl;
#endif
		grid_.cells_[i][j][k].press_ = pressure[l];
	}
	
	
	//update velocities
	//FIXME:set boundary vel to zero
	scale = DT/(RHO*dx);

	//update x vel;
	LOOP_FOR_XVELS {
		if(grid_.cells_[i][j][k].type_ == FLUID) {
			grid_.velX_[i][j][k] -= scale*grid_.cells_[i][j][k].press_;
			grid_.velX_[i+1][j][k] += scale*grid_.cells_[i][j][k].press_;
		}
	} END_LOOP

	//update y vel;
	LOOP_FOR_YVELS {
		if(grid_.cells_[i][j][k].type_ == FLUID) {
			grid_.velY_[i][j][k] -= scale*grid_.cells_[i][j][k].press_;
			grid_.velY_[i][j+1][k] += scale*grid_.cells_[i][j][k].press_;
		}
	} END_LOOP

	//update z vel;
	LOOP_FOR_ZVELS {
		if(grid_.cells_[i][j][k].type_ == FLUID) {
			grid_.velZ_[i][j][k] -= scale*grid_.cells_[i][j][k].press_;
			grid_.velZ_[i][j][k+1] += scale*grid_.cells_[i][j][k].press_;
		}
	} END_LOOP

#if 0
		std::cout << grid_.vel_[i][j][k].x_ << " " <<
			grid_.vel_[i][j][k].y_ << " " <<
			grid_.vel_[i][j][k].z_ << std::endl;
#endif
}

void Flip::AddExtForce()
{
	LOOP_FOR_PARTICLES(particles_) {
		p->velocity_.y_ += -DT*GRAVITY;
	} END_LOOP
}

void Flip::SaveGridVel()
{
	LOOP_FOR_XVELS {
		savedGrid_.velX_[i][j][k] = grid_.velX_[i][j][k];
	} END_LOOP

	LOOP_FOR_YVELS {
		savedGrid_.velY_[i][j][k] = grid_.velY_[i][j][k];
	} END_LOOP

	LOOP_FOR_ZVELS {
		savedGrid_.velZ_[i][j][k] = grid_.velZ_[i][j][k];
	} END_LOOP
}

void Flip::ExtrapolateVelocity()
{
}

void Flip::SolvePICFLIP()
{
	Allocator<GLdouble> velAlloc;
	GLdouble*** deltaVelX = velAlloc.Alloc3D(grid_.nX_+1, grid_.nY_, grid_.nZ_);
	GLdouble*** deltaVelY = velAlloc.Alloc3D(grid_.nX_, grid_.nY_+1, grid_.nZ_);
	GLdouble*** deltaVelZ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_+1);

#if 1
	std::cout << std::endl << "deltaVelY: " << std::endl;
#endif
	LOOP_FOR_XVELS {
		deltaVelX[i][j][k] = grid_.velX_[i][j][k] - savedGrid_.velX_[i][j][k];
	} END_LOOP

	LOOP_FOR_YVELS {
		deltaVelY[i][j][k] = grid_.velY_[i][j][k] - savedGrid_.velY_[i][j][k];
#if 1
		std::cout << deltaVelY[i][j][k] << " ";
#endif
	} END_LOOP
#if 0
	std::cout << std::endl;
#endif

	LOOP_FOR_ZVELS {
		deltaVelZ[i][j][k] = grid_.velZ_[i][j][k] - savedGrid_.velZ_[i][j][k];
	} END_LOOP

	//trilinear interpolation
#if 0
	std::cout << std::endl << "i j k" << std::endl;
#endif
	GLdouble tx, ty, tz;
	LOOP_FOR_PARTICLES(particles_) {
		int i = (int)(p->position_.x_ / grid_.cellSize_);
		int j = (int)(p->position_.y_ / grid_.cellSize_);
		int k = (int)(p->position_.z_ / grid_.cellSize_);

		tx = (p->position_.x_ - i*grid_.cellSize_) / grid_.cellSize_;
		ty = (p->position_.y_ - j*grid_.cellSize_) / grid_.cellSize_;
		tz = (p->position_.z_ - k*grid_.cellSize_) / grid_.cellSize_;

		p->velocity_.x_ += tx*deltaVelX[i+1][j][k] + (1-tx)*deltaVelX[i][j][k];
		p->velocity_.y_ += ty*deltaVelY[i][j+1][k] + (1-ty)*deltaVelY[i][j][k];
		p->velocity_.z_ += tz*deltaVelZ[i][j][k+1] + (1-tz)*deltaVelZ[i][j][k];
#if 0
		//if(std::abs(deltaVelX[i][j+1][k]-0.0588) > 1e-3)
		std::cout << i << " " << j << " " << k << std::endl;
#endif
#if 0
		std::cout << tx << " " << ty << " " << tz << std::endl;
#endif
#if 0
		std::cout << p->velocity_.x_ << " " << p->velocity_.y_ << " " <<
			p->velocity_.z_ << std::endl;
#endif
	} END_LOOP
}

void Flip::ComputeDensity()
{
}

/*FIXME: Loop for particles is better?*/
void Flip::TransferParticleVelToCell()
{
	int i, j, k;
	GLdouble weightedVel;
	GLdouble weight, totalWeight;
	ParticleList pl;
	GLdouble pX, pY, pZ;
	GLdouble cellCenterX, cellCenterY, cellCenterZ;
	
	//interpolate x velocity to grid
#if 0 
	std::cout <<"x vel " << std::endl;
#endif
	LOOP_FOR_XVELS {
		pl = _GetNeighborParticles(i, j, k, 1, 2, 2);
		weight = totalWeight = 0;
		weightedVel = 0;
		LOOP_FOR_PARTICLES(pl) {
			cellCenterX = (i)*grid_.cellSize_;
			cellCenterY = (j+0.5)*grid_.cellSize_;
			cellCenterZ = (k+0.5)*grid_.cellSize_;
			pX = p->position_.x_;
			pY = p->position_.y_;
			pZ = p->position_.z_;
			Eigen::Vector3d subVec(pX-cellCenterX, 
					pY-cellCenterY,	pZ-cellCenterZ);
			weight = _ReciprocalKernal(subVec.norm());
			weightedVel += p->velocity_.x_ * weight;
			totalWeight += weight;
		} END_LOOP
		if(weight != 0) {
			grid_.velX_[i][j][k] = weightedVel / totalWeight;
#if 0
			std::cout << weightedVel << " " << weight << " "
				<< grid_.velX_[i][j][k] << std::endl;
#endif
		}
	}END_LOOP

#if 0 
	std::cout << std::endl <<"transfering y vel " << std::endl;
#endif
	//interpolate y velocity to grid
	LOOP_FOR_YVELS {
		pl = _GetNeighborParticles(i, j, k, 2, 1, 2);
		weight = totalWeight = 0;
		weightedVel = 0;
		LOOP_FOR_PARTICLES(pl) {
			cellCenterX = (i+0.5)*grid_.cellSize_;
			cellCenterY = (j)*grid_.cellSize_;
			cellCenterZ = (k+0.5)*grid_.cellSize_;
			pX = p->position_.x_;
			pY = p->position_.y_;
			pZ = p->position_.z_;
			Eigen::Vector3d subVec(pX-cellCenterX, 
					pY-cellCenterY,	pZ-cellCenterZ);
			weight = _ReciprocalKernal(subVec.norm());
			weightedVel += p->velocity_.y_ * weight;
			totalWeight += weight;
		} END_LOOP
		if(weight != 0) {
			grid_.velY_[i][j][k] = weightedVel / totalWeight;
#if 0
			std::cout << grid_.velY_[i][j][k] << std::endl;
#endif
		}
	}END_LOOP

		//interpolate z velocity to grid
	LOOP_FOR_ZVELS {
		pl = _GetNeighborParticles(i, j, k, 2, 2, 1);
		weight = totalWeight = 0;
		weightedVel = 0;
		LOOP_FOR_PARTICLES(pl) {
			cellCenterX = (i+0.5)*grid_.cellSize_;
			cellCenterY = (j+0.5)*grid_.cellSize_;
			cellCenterZ = (k)*grid_.cellSize_;
			pX = p->position_.x_;
			pY = p->position_.y_;
			pZ = p->position_.z_;
			Eigen::Vector3d subVec(pX-cellCenterX, 
					pY-cellCenterY,	pZ-cellCenterZ);
			weight = _ReciprocalKernal(subVec.norm());
			weightedVel += p->velocity_.z_ * weight;
			totalWeight += weight;
		} END_LOOP
		if(weight != 0)
			grid_.velZ_[i][j][k] = weightedVel / totalWeight;

	} END_LOOP

}

ParticleList Flip::_GetNeighborParticles(int i, int j, int k, int w, int h, int d)
{
	ParticleList res;
	for( int si=i-w; si<=i+w-1; si++ )
		for( int sj=j-h; sj<=j+h-1; sj++ ) 
			for( int sk=k-d; sk<=k+d-1; sk++ ) {
		if(si < 0 || si > grid_.nX_-1 || 
				sj < 0 || sj > grid_.nY_-1 || sk < 0 || sk > grid_.nZ_-1) 
			continue;
		ParticleList pl = grid_.cells_[si][sj][sk].particles_;
		for( int l=0; l<pl.size(); l++ ) { 
			Particle *p = pl[l];
			res.push_back(p);
		}
	}
	return res;
}

GLdouble Flip::_ReciprocalKernal(GLdouble d)
{
	return KERNEL_H*KERNEL_H/std::max(d,1e-5) - 1.0;
}

void Flip::MatchParticlesToCell()
{
	LOOP_FOR_CELLS(grid_.nX_, grid_.nY_, grid_.nZ_) {
		grid_.cells_[i][j][k].particles_.clear();
	} END_LOOP

	LOOP_FOR_PARTICLES(particles_) {
		_InsertToCell(p);
	} END_LOOP

	_MarkFluid();

}

void Flip::MixPICFLIP()
{

}

void Flip::AdvectParticles()
{
	LOOP_FOR_PARTICLES(particles_) {
#if 0
		std::cout << p->velocity_.x_ << " " << p->velocity_.y_
			<< " " << p->velocity_.z_ << std::endl;
#endif
		p->position_.x_ += DT*p->velocity_.x_;
		p->position_.y_ += DT*p->velocity_.y_;
		p->position_.z_ += DT*p->velocity_.z_;
	} END_LOOP
}

void Flip::ParticleCollisionDetection()
{
	/*FIXME: using more flexible collision detection*/
	LOOP_FOR_PARTICLES(particles_) {
		int i, j, k;
		do {
			p->position_.x_ = std::max(p->position_.x_,0.0);
			p->position_.y_ = std::max(p->position_.y_,0.0);
			p->position_.z_ = std::max(p->position_.z_,0.0);
		//	p->position_.x_ = isnan(p->position_.x_) ? 0.0 : p->position_.x_;
		//	p->position_.y_ = isnan(p->position_.y_) ? 0.0 : p->position_.y_;
		//	p->position_.z_ = isnan(p->position_.z_) ? 0.0 : p->position_.z_;
			i = (int)(p->position_.x_ / grid_.cellSize_);
			j = (int)(p->position_.y_ / grid_.cellSize_);
			k = (int)(p->position_.z_ / grid_.cellSize_);
			i = (i >= grid_.nX_ ? grid_.nX_-1 : i);
			j = (j >= grid_.nY_ ? grid_.nY_-1 : j);
			k = (k >= grid_.nZ_ ? grid_.nZ_-1 : k);
#if 0
			std::cout << i << " " << j << " " << k << std::endl;
#endif
			if(grid_.cells_[i][j][k].type_ == SOLID) {
				p->position_.x_ += boundary->normals[i][j][k].x_;
				p->position_.y_ += boundary->normals[i][j][k].y_;
				p->position_.z_ += boundary->normals[i][j][k].z_;
			//	p->velocity_.x_ += 10* boundary->normals[i][j][k].x_;
			//	p->velocity_.y_ += 10* boundary->normals[i][j][k].y_;
			//	p->velocity_.z_ += 10* boundary->normals[i][j][k].z_;
			//	p->position_.x_ += p->velocity_.x_ * DT;
			//	p->position_.y_ += p->velocity_.y_ * DT;
			//	p->position_.z_ += p->velocity_.z_ * DT;
			//	if (boundary->normals[i][j][k].x_ != 0)
			//		p->velocity_.x_ = 0;
			//	if (boundary->normals[i][j][k].y_ != 0)
			//		p->velocity_.y_ = 0;
			//	if (boundary->normals[i][j][k].z_ != 0)
			//		p->velocity_.z_ = 0;
			}
		}while(i < 1 || j < 1 || k < 1 || 
				i > grid_.nX_-2 || j > grid_.nY_-2 || k > grid_.nZ_-2);
	} END_LOOP
}

void Flip::GenerateSurface()
{

}

void Flip::SimulateStep()
{
    
	//timestep++;
	step_++;

	//add externaal force, e.g. gravity
	AddExtForce();

	//match particles to corresponding cells
	MatchParticlesToCell();
	
	ComputeDensity();
	
	//transfer particle velocities to cell
	TransferParticleVelToCell();
#if 1
	std::cout << std::endl << "after TransferParticleVelToCell, grid_.velY_: " << std::endl;
	LOOP_FOR_YVELS {
		std::cout << grid_.velY_[i][j][k] <<  "   ";
	} END_LOOP
	std::cout << std::endl;
#endif

	SaveGridVel();

	//solve pressure, enforce boundaries, extrapolate velocities
	EnforceBoundary();
	Project();
	EnforceBoundary();
#if 1
	std::cout << std::endl << "after Project: grid_.velY_: " << std::endl;
	LOOP_FOR_YVELS {
		std::cout << grid_.velY_[i][j][k] << "    ";
	} END_LOOP
	std::cout << std::endl;
#endif

#if 1
	std::cout << std::endl << "before SolvePICFLIP: particle.velY: " << std::endl;
	LOOP_FOR_PARTICLES(particles_) {
		std::cout << p->velocity_.y_ << "    ";
	} END_LOOP
	std::cout << std::endl;
#endif

    SolvePICFLIP();
#if 1
	std::cout << std::endl << "after SolvePICFLIP: particle.velY: " << std::endl;
	LOOP_FOR_PARTICLES(particles_) {
		std::cout << p->velocity_.y_ << "    ";
	} END_LOOP
	std::cout << std::endl;
#endif

	//mix pic and flip velocities
	MixPICFLIP();
	
	AdvectParticles();

	ParticleCollisionDetection();
	
	//genenrate surface mesh
	GenerateSurface();
}
