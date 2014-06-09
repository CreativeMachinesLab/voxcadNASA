/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VX_Voxel.h"
#include "VX_Material.h"
#include "VX_Link.h"
#include <assert.h>

CVX_Voxel::CVX_Voxel(CVX_Material* material) 
{
	for (int i=0; i<6; i++) links[i]=NULL;
	mat = material;
	pos = linMom = angMom = Vec3D<>(0,0,0);
	dofFixed = boolStates = 0;
	//...etc

}

CVX_Voxel::~CVX_Voxel(void)
{
}

void CVX_Voxel::addLinkInfo(linkDirection direction, CVX_Link* link)
{
	links[direction] = link;
}

void CVX_Voxel::removeLinkInfo(linkDirection direction)
{
	links[direction]=NULL;
}

void CVX_Voxel::replaceMaterial(CVX_Material* newMaterial)
{
	if (newMaterial != NULL){
		linMom *= newMaterial->_mass/mat->_mass; //adjust momuntums to keep velocity constant accorss material change
		angMom *= newMaterial->_momentInertia/mat->_momentInertia;

		mat = newMaterial;

		for (int i=0; i<6; i++)	if (links[i]) links[i]->updateProperties();
		setFloorStaticFriction(false);
	}
}

bool CVX_Voxel::isYielded() const
{
	for (int i=0; i<6; i++){
		if (links[i] && links[i]->isYielded(isNegative((linkDirection)i))) return true;
	}
	return false;
}

bool CVX_Voxel::isFailed() const
{
	for (int i=0; i<6; i++){
		if (links[i] && links[i]->isFailed(isNegative((linkDirection)i))) return true;
	}
	return false;
}

void CVX_Voxel::setFixed(dofComponent dof, float displacement)
{
	dofSet(dofFixed, dof, true);
	if (dof && X_TRANSLATE) pos.x+=displacement;
	if (dof && Y_TRANSLATE) pos.y+=displacement;
	if (dof && Z_TRANSLATE) pos.z+=displacement;
	if (dof && X_ROTATE) orient = CQuat<>(displacement, Vec3D<>(1,0,0))*orient; //test me!
	if (dof && Y_ROTATE) orient = CQuat<>(displacement, Vec3D<>(0,1,0))*orient; //test me!
	if (dof && Z_ROTATE) orient = CQuat<>(displacement, Vec3D<>(0,0,1))*orient; //test me!

}

void CVX_Voxel::setAllFixed(const Vec3D<>& translation, const Vec3D<>& rotation)
{
	dofSetAll(dofFixed, true);
	pos += translation;
	CQuat<> rot(rotation);
//	rot.FromRotationVector(rotation);
	orient = rot*orient; //test me!
}

Vec3D<> CVX_Voxel::cornerNegative() const
{
	return -mat->nominalSize()*Vec3D<>(	1 + (links[1]?links[1]->axialStrain(true):0),
										1 + (links[3]?links[3]->axialStrain(true):0),
										1 + (links[5]?links[5]->axialStrain(true):0));
}

Vec3D<> CVX_Voxel::cornerPositive() const
{
	return mat->nominalSize()*Vec3D<>(	1 + (links[0]?links[0]->axialStrain(false):0),
										1 + (links[2]?links[2]->axialStrain(false):0),
										1 + (links[4]?links[4]->axialStrain(false):0));
}

//http://klas-physics.googlecode.com/svn/trunk/src/general/Integrator.cpp (reference)
void CVX_Voxel::timeStep(float dt)
{
	previousDt = dt;
	if (dt == 0.0f) return;

	if (dofIsAllSet(dofFixed)){ //if this voxel is fixed, no need to change its location
		haltMotion();
		return;
	}

	//Translation
	Vec3D<> curForce = force();
	linMom += curForce*dt;

	if (dofIsSet(dofFixed, X_TRANSLATE)){linMom.x=0;}
	if (dofIsSet(dofFixed, Y_TRANSLATE)){linMom.y=0;}
	if (dofIsSet(dofFixed, Z_TRANSLATE)){linMom.z=0;}

	Vec3D<double> translate(linMom*(dt*mat->_massInverse)); //movement of the voxel this timestep

	if (isFloorEnabled() && floorPenetration() > 0){ //we must catch a slowing voxel here since it all boils down to needing access to the dt of this timestep.
		vfloat work = curForce.x*translate.x + curForce.y*translate.y; //F dot disp
		if(kineticEnergy() + work < 0) setFloorStaticFriction(true); //this checks for a change of direction according to the work-energy principle

		if (isFloorStaticFriction()){ //if we're in a state of static friction, zero out all horizontal motion
			linMom.x = linMom.y = 0;
			translate.x = translate.y = 0;
		}
	}

	pos += translate;

	//Rotation
	Vec3D<> curMoment = moment();
	angMom += curMoment*dt;

	if (dofIsSet(dofFixed, X_ROTATE)){angMom.x=0;}
	if (dofIsSet(dofFixed, Y_ROTATE)){angMom.y=0;}
	if (dofIsSet(dofFixed, Z_ROTATE)){angMom.z=0;}

	orient = CQuat<>(angularVelocity()*dt)*orient; //update the orientation




//	if(pSim->StatToCalc & CALCSTAT_KINE) KineticEnergy = 0.5*Mass*Vel.Length2() + 0.5*Inertia*AngVel.Length2(); //1/2 m v^2
//	if(pSim->StatToCalc & CALCSTAT_PRESSURE) {
////		vfloat VolumetricStrain = GetVoxelStrain(AXIS_X) + GetVoxelStrain(AXIS_Y) + GetVoxelStrain(AXIS_Z);
//		vfloat VolumetricStrain = strain(false).x+strain(false).y+strain(false).z;
//		Pressure = - _pMat->GetElasticMod()*VolumetricStrain/(3*(1-2*_pMat->GetPoissonsRatio())); //http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf
//	}
//
//	updatePoissonsstrain();
//
//

}

Vec3D<> CVX_Voxel::force()
{
	//forces from internal bonds
	Vec3D<> totalForce(0,0,0);
	for (int i=0; i<6; i++){ 
		if (links[i]) totalForce += links[i]->force(isNegative((linkDirection)i)); //total force in LCS
	}
	totalForce = orient.RotateVec3D(totalForce);

	//other forces
	totalForce += extForce; //external forces
	totalForce -= velocity()*mat->globalDampingTranslateC(); //global damping
	if (isGravityEnabled()) totalForce.z -= mat->_mass*9.80665; //gravity, according to f=mg
	if (isFloorEnabled()) addFloorForces(&totalForce); //adds forces from interacting with the floor
	//Collisions: todo

	return totalForce;
}

Vec3D<> CVX_Voxel::moment()
{
	//moments from internal bonds
	Vec3D<> totalMoment(0,0,0);
	for (int i=0; i<6; i++){ 
		if (links[i]) totalMoment += links[i]->moment(isNegative((linkDirection)i)); //total force in LCS
	}
	totalMoment = orient.RotateVec3D(totalMoment);
	
	//other moments
	totalMoment += extMoment; //external moments
	totalMoment -= angularVelocity()*mat->globalDampingRotateC(); //global damping
	return totalMoment;
}


void CVX_Voxel::addFloorForces(Vec3D<>* pTotalForce)
{
	float CurPenetration = floorPenetration(); //for now use the average.

	if (CurPenetration>0){ 
		Vec3D<> vel = velocity();
		Vec3D<> horizontalVel(vel.x, vel.y, 0);
		
		float normalForce = mat->penetrationStiffness()*CurPenetration;
		pTotalForce->z += normalForce - mat->internalDampingTranslateC()*vel.z; //in the z direction: k*x-C*v - spring and damping

		if (isFloorStaticFriction()){ //If this voxel is currently in static friction mode (no lateral motion) 
			assert(horizontalVel.Length2() == 0);
			float surfaceForce = sqrt(pTotalForce->x*pTotalForce->x + pTotalForce->y*pTotalForce->y);
			if (surfaceForce > mat->muStatic*normalForce) setFloorStaticFriction(false); //if we're breaking static friction, leave the forces as they currently have been calculated to initiate motion this time step
			//otherwise we'll freeze horizontal motion in the timestep function
		}
		else { //else if this voxel is in kinetic friction mode (non-zero lateral motion)
			assert(horizontalVel.Length2() != 0);
			float surfaceVelAngle = atan2(vel.y, vel.x); //angle of sliding along floor...
			*pTotalForce -=  mat->muKinetic*normalForce*horizontalVel.Normalized(); //add a friction force opposing velocity according to the normal force and the kinetic coefficient of friction
		}
	}
}

//Vec3D<> CVX_Voxel::calcStressTensor() //sums all forces acting on this voxel
//{
	
//}

Vec3D<> CVX_Voxel::strain(bool tensionStrain) const
{
	//if no connections in the positive and negative directions of a particular axis, strain is zero
	//if one connection in positive or negative direction of a particular axis, strain is that strain - ?? and force or constraint?
	//if connections in both the positive and negative directions of a particular axis, strain is the average. 
	Vec3D<> intStrRet(0,0,0); //intermediate strain return value. axes according to linkAxis enum
	int numBondAxis[3] = {0}; //number of bonds in this axis (0,1,2). axes according to linkAxis enum
	for (int i=0; i<6; i++){ //cycle through link directions
		if (links[i]){
			int axis = toAxis((linkDirection)i);
			intStrRet[axis] += links[i]->axialStrain(isNegative((linkDirection)i));
			numBondAxis[axis]++;
		}
	}
	for (int i=0; i<3; i++){ //cycle through axes
		if (numBondAxis[i]==2) intStrRet[i] *= 0.5f; //average
		if (tensionStrain && numBondAxis[i]==1){ //if just one bond
			if (!dofIsSet(dofFixed, (dofComponent)i) && extForce[i] == 0) intStrRet[i]=0; //if no other external means of providing tension, zero out strain.
		}
	}

	return intStrRet;
}

Vec3D<> CVX_Voxel::poissonsStrain() const
{
	vfloat mu = mat->poissonsRatio();
	Vec3D<> tensileStrain = strain(true);
	return Vec3D<>(	tensileStrain.x + pow(1+tensileStrain.y + tensileStrain.z, -mu)-1, 
					tensileStrain.y + pow(1+tensileStrain.x + tensileStrain.z, -mu)-1, 
					tensileStrain.z + pow(1+tensileStrain.x + tensileStrain.y, -mu)-1);
}

float CVX_Voxel::axialStrain(linkAxis axis) const
{
	switch (axis){
	case X_AXIS: return poissonsStrain().x;
	case Y_AXIS: return poissonsStrain().y;
	case Z_AXIS: return poissonsStrain().z;
	default: return 0.0f;
	}
}


float CVX_Voxel::transverseStrainSum(linkAxis axis) const
{
	if (mat->poissonsRatio() == 0) return 0;
	else switch (axis){
	case X_AXIS: return poissonsStrain().y+poissonsStrain().z;
	case Y_AXIS: return poissonsStrain().x+poissonsStrain().z;
	case Z_AXIS: return poissonsStrain().x+poissonsStrain().y;
	default: return 0.0f;
	}
}

float CVX_Voxel::transverseArea(linkAxis axis) const
{
	float size = mat->nominalSize();
	if (mat->poissonsRatio() == 0) return size*size;
	else switch (axis){
	case X_AXIS: return size*size*(1+poissonsStrain().y)*(1+poissonsStrain().z);
	case Y_AXIS: return size*size*(1+poissonsStrain().x)*(1+poissonsStrain().z);
	case Z_AXIS: return size*size*(1+poissonsStrain().x)*(1+poissonsStrain().y);
	default: return size*size;
	}
}