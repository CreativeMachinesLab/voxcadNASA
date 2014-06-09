/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VX_Link.h"
#include "VX_Voxel.h"
#include <assert.h>

//we can do better...
#ifdef PREC_LOW
	static const vfloat SA_BOND_BEND_RAD = 0.1; //Amount for small angle bond calculations
#elif defined PREC_HIGH
	static const vfloat SA_BOND_BEND_RAD = 0.02; //Amount for small angle bond calculations
#elif defined PREC_MAX 
	static const vfloat SA_BOND_BEND_RAD = 0.002; //Amount for small angle bond calculations
#else //defined PREC_MED 
	static const vfloat SA_BOND_BEND_RAD = 0.05; //Amount for small angle bond calculations
#endif

static const vfloat SA_BOND_EXT_PERC = 1.30; //Amount for small angle bond calculations
//end we can do better


CVX_Link::CVX_Link(CVX_Voxel* voxel1, CVX_Voxel* voxel2, linkDirection direction)
{
	assert(voxel1 != NULL);
	assert(voxel2 != NULL);

	if (isPositive(direction)){
		pVNeg=voxel1;
		pVPos=voxel2;
	}
	else {
		pVNeg=voxel2;
		pVPos=voxel1;
	}
	axis = toAxis(direction);

	boolStates=0;
	strainNeg=0.0f;
	strainPos=0.0f;
	smallAngle=true;

	updateProperties();
}

CVX_Link::~CVX_Link(void)
{

}

float CVX_Link::avgTransverseStrainSum()
{
	return 0.5*(pVNeg->transverseStrainSum(axis)+pVPos->transverseStrainSum(axis));
}

float CVX_Link::avgTransverseArea()
{
	return 0.5*(pVNeg->transverseArea(axis)+pVPos->transverseArea(axis));
}

CQuat<> CVX_Link::orientLink(float restLength) //updates pos2, angle1, angle2, and smallAngle
{
	pos2 = toAxisX(Vec3D<double>(pVPos->position() - pVNeg->position())); //digit truncation happens here...
	angle1 = toAxisX(pVNeg->orientation());
	angle2 = toAxisX(pVPos->orientation());

	CQuat<> totalRot = angle1.Conjugate(); //keep track of the total rotation of this bond (after toAxisX())
	pos2 = totalRot.RotateVec3D(pos2);
	angle2 = totalRot*angle2;
	angle1 = CQuat<>(); //zero for now...

	//small angle approximation?
	float SmallTurn = (abs(pos2.z)+abs(pos2.y))/pos2.x;
	float ExtendPerc = pos2.x/restLength;
	if (!smallAngle && angle2.IsSmallAngle() && SmallTurn < SA_BOND_BEND_RAD && ExtendPerc < SA_BOND_EXT_PERC){
		smallAngle = true;
		setBoolState(LOCAL_VELOCITY_VALID, false);
	}
	else if (smallAngle && (!angle2.IsSmallishAngle() || SmallTurn > VEC3D_HYSTERESIS_FACTOR*SA_BOND_BEND_RAD || ExtendPerc > VEC3D_HYSTERESIS_FACTOR*SA_BOND_EXT_PERC)){
		smallAngle = false;
		setBoolState(LOCAL_VELOCITY_VALID, false);
	}

	if (smallAngle)	{ //Align so Angle1 is all zeros
		pos2.x -= restLength; //only valid for small angles
	}
	else { //Large angle. Align so that Pos2.y, Pos2.z are zero.
		angle1.FromAngleToPosX(pos2); //get the angle to align Pos2 with the X axis
		totalRot = angle1 * totalRot; //update our total rotation to reflect this
		angle2 = angle1 * angle2; //rotate angle2
		pos2 = Vec3D<>(pos2.Length() - restLength, 0, 0); 
	}

	angle1v = angle1.ToRotationVector();
	angle2v = angle2.ToRotationVector();
	return totalRot;
}

bool CVX_Link::isYielded() const
{
	return (pVNeg->mat->isYielded(maxStrainNeg()) || pVPos->mat->isYielded(maxStrainPos()));
}

bool CVX_Link::isYielded(bool positiveEnd) const
{
	return positiveEnd ? pVPos->mat->isYielded(maxStrainPos()) : pVNeg->mat->isYielded(maxStrainNeg());
}

bool CVX_Link::isFailed() const
{
	return (pVNeg->mat->isFailed(maxStrainNeg()) || pVPos->mat->isFailed(maxStrainPos()));
}

bool CVX_Link::isFailed(bool positiveEnd) const
{
	return positiveEnd ? pVPos->mat->isFailed(maxStrainPos()) : pVNeg->mat->isFailed(maxStrainNeg());
}

void CVX_Link::updateForces()
{
	Vec3D<double> oldPos2 = pos2, oldAngle1v = angle1v, oldAngle2v = angle2v; //remember the positions/angles from last timestep to calculate velocity

	float restLength = 0.5*(pVNeg->nominalSize(axis) + pVPos->nominalSize(axis));
	CQuat<> rotation = orientLink(restLength); //sets pos2, angle1, angle2

	Vec3D<double> dPos2 = pos2-oldPos2; //deltas for local damping
	Vec3D<double> dAngle1 = angle1v-oldAngle1v;
	Vec3D<double> dAngle2 = angle2v-oldAngle2v;

	float curStress = updateStrain(pos2.x/restLength);
	currentA1 = isLinear() && isAreaConstant() || pos2.x==0 ? a1 : curStress*avgTransverseArea()/pos2.x; // calculate the linear stiffness

	//Beam equations. All relevant terms are here, even though some are zero for small angle and others are zero for large angle (profiled as negligible performance penalty)
	forceNeg = Vec3D<> (	currentA1*pos2.x,
							b1*pos2.y - b2*(angle1v.z + angle2v.z),
							b1*pos2.z + b2*(angle1v.y + angle2v.y)); //Use Curstress instead of -a1*Pos2.x to account for non-linear deformation 
	forcePos = -forceNeg;

	momentNeg = Vec3D<> (	a2*(angle1v.x - angle2v.x),
							b2*pos2.z + b3*(2*angle1v.y + angle2v.y),
							-b2*pos2.y + b3*(2*angle1v.z + angle2v.z));
	momentPos = Vec3D<> (	a2*(angle2v.x - angle1v.x),
							b2*pos2.z + b3*(angle1v.y + 2*angle2v.y),
							-b2*pos2.y + b3*(angle1v.z + 2*angle2v.z));

	//local damping:
	if (isLocalVelocityValid()){ //if we don't have the basis for a good damping calculation, don't do any damping.
		assert(pVNeg->previousDt == pVPos->previousDt); //in the future this may change, but for now assume synchronous timestep.
		Vec3D<> posCalc(	(a1==currentA1?sqA1:sqrt(currentA1))*dPos2.x,
							sqB1*dPos2.y - sqB2xFMp*(dAngle1.z+dAngle2.z), //units: sqrt(N*m)/s
							sqB1*dPos2.z + sqB2xFMp*(dAngle1.y+dAngle2.y));
		Vec3D<> angCalc(	sqA2xIp*(dAngle2.x - dAngle1.x), //units: sqrt(N-m)*
							sqB2xFMp*dPos2.z + sqB3xIp*(dAngle1.y + dAngle2.y),
							-sqB2xFMp*dPos2.y + sqB3xIp*(dAngle1.z + dAngle2.z));
		
		forceNeg += pVNeg->dampingMultiplier() * posCalc;
		forcePos -= pVPos->dampingMultiplier() * posCalc;

		momentNeg += 0.5*pVNeg->dampingMultiplier()*Vec3D<>(-angCalc.x, angCalc.y+sqB3xIp*dAngle1.y, angCalc.z+sqB3xIp*dAngle1.z);
		momentPos += 0.5*pVPos->dampingMultiplier()*Vec3D<>(angCalc.x, angCalc.y+sqB3xIp*dAngle2.y, angCalc.z+sqB3xIp*dAngle2.z);
	}
	else setBoolState(LOCAL_VELOCITY_VALID, true); //we're good for next go-around unless something changes

	//	transform forces and moments to local voxel coordinates
	if (!smallAngle){
		forceNeg = angle1.RotateVec3DInv(forceNeg);
		momentNeg = angle1.RotateVec3DInv(momentNeg);
	}
	forcePos = angle2.RotateVec3DInv(forcePos);
	momentPos = angle2.RotateVec3DInv(momentPos);

	toAxisOriginal(&forceNeg);
	toAxisOriginal(&forcePos);
	toAxisOriginal(&momentNeg);
	toAxisOriginal(&momentPos);


}


void CVX_Link::scaleStrains(float axialStrain)
{
	double strainFactor;
	if (strainNeg+strainPos == 0.0) strainFactor = 0.0;
	else strainFactor = 2.0*axialStrain/(strainNeg+strainPos); //new axialStain / previos axialStrain ->double because this could be a source of float roundoff
	strainNeg *= strainFactor;
	strainPos *= strainFactor;
}

float CVX_Link::updateStrain(float axialStrain)
{
	float returnStress;
	if (isLinear()){
		if (isHomogenous()) strainNeg = strainPos = axialStrain;
		else { //non homogenous (use springs in series equations with some shortcuts to use elastic moduli. since all we need is the correct ratio before scaleStrains(), further shortcuts taken.)
			strainNeg = E/pVNeg->mat->E; // = Epos/(Eneg + Epos), because E = (Eneg*Epos/(Eneg+Epos))*2
			strainPos = E/pVPos->mat->E;
			scaleStrains(axialStrain);
		}

		if (isAreaConstant()) returnStress = E*axialStrain;
		else returnStress = equalizeStress(axialStrain);

	}
	else { //at least one material is nonlinear
		if (axialStrain > maxStrain){ //if new territory on the stress/strain curve
			if (isHomogenous() && isAreaConstant()){
				strainNeg = strainPos = axialStrain;
				returnStress = pVNeg->mat->stress(strainNeg);
			}
			else returnStress = equalizeStress(axialStrain); //nonHomogenous or non-constant areas
		}
		else { //backed off a non-linear material, therefore in linear region.
			float relativeStrain = axialStrain-strainOffset; // treat the material as linear with a strain offset according to the maximum plastic deformation
			if (isAreaConstant()) returnStress = E*relativeStrain;
			else returnStress = equalizeStress(relativeStrain); //correctly handle changes due to transverse strains
		}
	}

	//keep track of maximum (even for linear materials, in case they have a yield/fail point specified)
	if (axialStrain > maxStrain){ //if new territory on the stress/strain curve
		maxStrain = axialStrain; //remember this maximum for easy reference
		maxStrainRatio = (1.0+strainNeg)/(1.0+axialStrain); //dividing point = L*0.5*maxStrainRatio
		strainOffset = isLinear() ? 0.0f : maxStrain-returnStress/E; //precalculate strain offset for when we back off
	}

	return returnStress;
}

float CVX_Link::equalizeStress(float axialStrain, bool forceLinear) //the long (but thorough) way: updates strainNeg and strainPos according to the provided axial strain. returns current stress as well (MPa)
{
	CVX_Material *matNeg = pVNeg->mat, *matPos = pVPos->mat;
	float tssNeg = pVNeg->transverseStrainSum(axis), tssPos = pVPos->transverseStrainSum(axis);

	scaleStrains(axialStrain);
	float stressNeg = matNeg->stress(strainNeg, tssNeg, forceLinear);
	float stressPos = matPos->stress(strainPos, tssPos, forceLinear);
	if (fabs((stressNeg - stressPos)/stressPos) > 0.001f){ //if division of strain isn't quite accurate and stresses are off by more than 0.1%, update the relative srains a bit
		strainNeg *= 2*stressPos/(stressNeg+stressPos);
		strainPos *= 2*stressNeg/(stressNeg+stressPos);
		scaleStrains(axialStrain);
		stressNeg = matNeg->stress(strainNeg, tssNeg, forceLinear);
		stressPos = matPos->stress(strainPos, tssPos, forceLinear);
	}
	return 0.5*(stressNeg+stressPos);
}

void CVX_Link::updateProperties() //updates a1, b1, etc. based on the materials of the linked voxels
{
	CVX_Material *matNeg = pVNeg->mat, *matPos = pVPos->mat;
	setBoolState(HOMOGENOUS, matNeg == matPos); //note if both materials are the same
	setBoolState(LINEAR, matNeg->isModelLinear() && matPos->isModelLinear()); //note if material is linear
	setBoolState(AREACONSTANT, matNeg->poissonsRatio()==0.0f && matPos->poissonsRatio()==0.0f); //note if our cross-sectional area will change with strain
	setBoolState(LOCAL_VELOCITY_VALID, false);

	float L = matNeg->nominalSize();
	float E1= matNeg->youngsModulus(), E2 = matPos->youngsModulus();
	float v1= matNeg->poissonsRatio(), v2 = matPos->poissonsRatio();

	E = (E1*E2/(E1+E2))*2; //x2 derived from case of equal stiffness: E1*E1/(E1+E1) = 0.5*E1 
	float v = (v1+v2)/2; //poissons ratio

	//stiffnesses terms
	a1 = E*L; //EA/L : Units of N/m
	currentA1 = a1; //instantaneous axial stiffness : Units of N/m
	a2 = E * L*L*L / (12.0f*(1+v)); //GJ/L : Units of N-m
	b1 = E*L; //12EI/L^3 : Units of N/m
	b2 = E*L*L/2.0f; //6EI/L^2 : Units of N (or N-m/m: torque related to linear distance)
	b3 = E*L*L*L/6.0f; //2EI/L : Units of N-m
	
	//damping sqrt(mk) terms (with sqrt(m) factored out)
	sqA1=sqrt(a1); //Units of sqrt(N/m)
	sqA2xIp=sqrt(a2*L*L/6.0f); //m*Ip = I : Units of sqrt(N-m^3) = m*sqrt(N-m)
	sqB1=sqrt(b1); //units of sqrt(N/m)
	sqB2xFMp=sqrt(b2*L/2); //m*FMp = FM : Units of sqrt(N-m)
	sqB3xIp=sqrt(b3*L*L/6.0f); //m*Ip = I : Units of sqrt(N-m^3) = m*sqrt(N-m)

	//reset the max strains (i.e. erase any plastic deformation in case materials are switched)
	maxStrain=0.0f;
	maxStrainRatio = 0.5f;
	strainOffset=0.0f;

}

float CVX_Link::strainEnergy() const
{
	return	forceNeg.x*forceNeg.x/(2.0f*a1) + //Tensile strain
			momentNeg.x*momentNeg.x/(2.0*a2) + //Torsion strain
			(momentNeg.z*momentNeg.z - momentNeg.z*momentPos.z +momentPos.z*momentPos.z)/(3.0*b3) + //Bending Z
			(momentNeg.y*momentNeg.y - momentNeg.y*momentPos.y +momentPos.y*momentPos.y)/(3.0*b3); //Bending Y
}
