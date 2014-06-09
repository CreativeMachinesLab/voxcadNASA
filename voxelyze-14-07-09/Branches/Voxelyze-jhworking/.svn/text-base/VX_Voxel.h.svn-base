/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#ifndef VX_VOXEL_H
#define VX_VOXEL_H

class CVX_Link;
#include "Vec3D.h"
#include "VX_Enums.h"
#include "VX_Material.h" //needed for inline of some "get" functions

//A voxel has a local coordinate system (LCS) that always stays centered on the center of the voxel and oriented with the axes of the cube.
//the global coordinate system (GCS) is just that - the global coordinate system.
class CVX_Voxel
{
public:
	CVX_Voxel(CVX_Material* material); //a voxel cannot exist without a material
	~CVX_Voxel();

	void replaceMaterial(CVX_Material* newMaterial); 
	CVX_Material* material() {return mat;}

	void timeStep(float dt);

	//physical location
	Vec3D<> position() const {return pos;} //center of voxel (GCS)
	Vec3D<> size() const {return cornerPositive()-cornerNegative();} // in local voxel coordinates system (not guaranteed to be centered on center
	Vec3D<> cornerNegative() const; // in local voxel coordinates system
	Vec3D<> cornerPositive() const; // in local voxel coordinates system
	bool isInterior() const {return links[0] && links[1] && links[2] && links[3] && links[4] && links[5];} //returns true if this voxel has no exposed faces
	bool isSurface() const {return !isInterior();} //returns true if this voxel has at least one exposed face

	Vec3D<> nominalSize() const {return mat->size()*(1+tempRelative*mat->alphaCTE);} //accounts for temperature and external actuation
	float nominalSize(linkAxis axis) const {return mat->size()[axis]*(1+tempRelative*mat->alphaCTE);} //accounts for temperature and external actuation
	float nominalSizeAverage() const {Vec3D<> nomSize=nominalSize(); return (nomSize.x+nomSize.y+nomSize.z)/3.0f;}

	CQuat<> orientation() const {return orient;} //global CS
	float orientationAngle() const {return orient.Angle();} //global CS
	Vec3D<> orientationAxis() const {return orient.AxisNormalized();} //global CS

	Vec3D<> velocity() const {return linMom*mat->_massInverse;} //global CS
	float velocityMagnitude() const {return linMom.Length()*mat->_massInverse;}
	Vec3D<> angularVelocity() const {return angMom*mat->_momentInertiaInverse;} //global CS
	float angularVelocityMagnitude() const {return angMom.Length()*mat->_momentInertiaInverse;} //global CS
	float kineticEnergy() const {return 0.5*(mat->_massInverse*linMom.Length2() + mat->_momentInertiaInverse*angMom.Length2());}
	float volumetricStrain() const {return strain(false).x+strain(false).y+strain(false).z;}
	float pressure() const {return -mat->youngsModulus()*volumetricStrain()/(3*(1-2*mat->poissonsRatio()));} //http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf

	//material state
	bool isYielded() const; //returns true if any of the links are yielded
	bool isFailed() const; //returns true if any of the links are failed

	//@ voxel level for heat diffusion experiments later
	void temperatureRelative(); //returns the temperature relative to zero-expansion temperature
	void setTemperatureRelative(); //determines the amount of scaling from the temperature


	Vec3D<> externalForce() const {return extForce;}
	void setExternalForce(const Vec3D<>& force) {extForce = force;} //in N
	Vec3D<> externalMoment() const {return extForce;}
	void setExternalMoment(const Vec3D<>& moment) {extMoment = moment;} //in N-m

	void haltMotion(){linMom = angMom = Vec3D<>(0,0,0);} //stops all motion

	void setFixed(dofComponent dof, float displacement=0.0f); //applies the specified displacement then fixes the specified dof.
	void setUnfixed(dofComponent dof) {dofSet(dofFixed, dof, false);}
	bool isFixed(dofComponent dof) {dofIsSet(dofFixed, dof);}
	
	//convenience
	void setAllFixed(const Vec3D<>& translation = Vec3D<>(0,0,0), const Vec3D<>& rotation = Vec3D<>(0,0,0));
	void setAllUnfixed() {dofSetAll(dofFixed, false);}
	bool isAllFixed() {dofIsAllSet(dofFixed);}

	bool isGravityEnabled() const {return boolStates & GRAVITY_ENABLED;}
	bool isFloorEnabled() const {return boolStates & FLOOR_ENABLED;}
	bool isFloorStaticFriction() const {return boolStates & FLOOR_STATIC_FRICTION;}
	float floorPenetration() const {return nominalSizeAverage()-pos.z;} //returns positive for interference

	Vec3D<> force(); //sum of the current forces on this voxel
	Vec3D<> moment(); //sum of the current moments on this voxel

private:
	typedef int voxState;
	enum voxFlags { //default of each should be zero for easy clearing
		FLOOR_ENABLED = 1<<2, //interact with a floor at z=0?
		FLOOR_STATIC_FRICTION = 1<<3, //is the voxel in a state of static friction with the floor?
		GRAVITY_ENABLED = 1<<4 //should this voxel respond to gravity?
	};

	CVX_Material* mat;

	void addLinkInfo(linkDirection direction, CVX_Link* link); //adds the information about a link connected to this voxel in the specified direction
	void removeLinkInfo(linkDirection direction); //removes the information about a link connected to this voxel in the specified direction
	CVX_Link* links[6]; //links in the 6 cardinal directions according to linkDirection enumeration



	//voxel state
	Vec3D<vfloat> pos;					//current center position (meters) (GCS)
	Vec3D<vfloat> linMom;				//current linear momentum (kg*m/s) (GCS)
	CQuat<vfloat> orient;				//current orientation (GCS)
	Vec3D<vfloat> angMom;				//current angular momentum (kg*m^2/s) (GCS)

	voxState boolStates;				//single int to store many boolean state values as bit flags according to 
	void setFloorStaticFriction(bool active) {active? boolStates |= FLOOR_STATIC_FRICTION : boolStates &= !FLOOR_STATIC_FRICTION;}

	dofObject dofFixed;
	
	//potential inputs affecting this local voxel
	Vec3D<float> extForce, extMoment; //External force, moment applied to this voxel (N, N-m) if relevant DOF are unfixed

	float tempRelative; //0 is no expansion

	void addFloorForces(Vec3D<>* pTotalForce); //modifies pTotalForce to include the object's interaction with a floor. This should be calculated as the last step of sumForce so that pTotalForce is complete.

	float dampingMultiplier() {return 2*mat->_sqrtMass*mat->zetaInternal/previousDt;} //an arbitrary multiplier used in internal bond damping

	Vec3D<> strain(bool tensionStrain) const; //LCS returns voxel strain. if tensionStrain true and no actual tension in that
	Vec3D<> poissonsStrain() const;
	float axialStrain(linkAxis axis) const;
	float transverseStrainSum(linkAxis axis) const;
	float transverseArea(linkAxis axis) const;

	void eulerStep(float dt); //execute an euler time step at the specified dt
	float previousDt; //remember the duration of the last timestep of this voxel
	friend class CVoxelyze; //give access to private members directly
	friend class CVX_Link;
};

//References:
//http://gafferongames.com/game-physics/physics-in-3d/

//Poissons ratio (volumetric expansion)
//http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf


#endif //VX_VOXEL_H