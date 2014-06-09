/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VX_Material.h"
#include <assert.h>

//if any of these change, it has implications on save files.
#define DEFAULT_YOUNGSMODULUS 1000000.0f
//#define DEFAULT_YIELDSTRESS -1.0f
//#define DEFAULT_FAILURESTRESS -1.0f
#define DEFAULT_POISSONSRATIO 0.0f
#define DEFAULT_DENSITY 1000.0f
//#define DEFAULT_CTE 0.0f;
//#define DEFAULT_STATICFRICTION 0.0f;
//#define DEFAULT_KINETICFRICTION 0.0f;

CVX_Material::CVX_Material(void)
{
	r = -1;
	g = -1;
	b = -1;
	a = -1;
	linear = true;
	E = DEFAULT_YOUNGSMODULUS;
	sigmaYield = -1.0;
	sigmaFail = -1.0;
	epsilonYield = -1.0;
	epsilonFail = -1.0;
	nu = DEFAULT_POISSONSRATIO;
	rho = DEFAULT_DENSITY;
	alphaCTE = 0.0;
	muStatic = 0.0;
	muKinetic = 0.0;
	zetaInternal = 0.1f;
	zetaGlobal = 0.1f;
	nomSize=1.0;
	extScale=Vec3D<>(1.0, 1.0, 1.0);
	updateDerived();
}

CVX_Material::~CVX_Material(void)
{

}

CVX_Material& CVX_Material::operator=(const CVX_Material& vIn)
{
	error = vIn.error;
	myName = vIn.myName;
	r = vIn.r;
	g = vIn.g;
	b = vIn.b;
	a = vIn.a;
	linear = vIn.linear;
	E = vIn.E;
	sigmaYield = vIn.sigmaYield;
	sigmaFail = vIn.sigmaFail;
	epsilonYield = vIn.epsilonYield;
	epsilonFail = vIn.epsilonFail;
	strainData = vIn.strainData;
	stressData = vIn.stressData;
	nu = vIn.nu;
	rho = vIn.rho;
	alphaCTE = vIn.alphaCTE;
	muStatic = vIn.muStatic;
	muKinetic = vIn.muKinetic;
	zetaInternal = vIn.zetaInternal;
	zetaGlobal = vIn.zetaGlobal;
	nomSize=vIn.nomSize;
	_mass=vIn._mass;
	_massInverse=vIn._massInverse;
	_sqrtMass=-vIn._sqrtMass;
	_firstMoment=vIn._firstMoment;
	_momentInertia=vIn._momentInertia;
	_momentInertiaInverse=vIn._momentInertiaInverse;
	_2xSqMxExS=vIn._2xSqMxExS;
	_2xSqIxExSxSxS=vIn._2xSqIxExSxSxS;

	return *this;
}

//float CVX_Material::stress(float strain)
//{
//	if (isFailed(strain)) return 0.0f; //if a failure point is set and exceeded, we've broken!
//	if (strain <= strainData[1] || linear) return E*strain; //for compression/first segment and linear materials, simple calculation
//
//	int DataCount = modelDataPoints();
//	for (int i=2; i<DataCount; i++){ //go through each segment in the material model (skipping the first segment because it has already been handled.
//		if (strain <= strainData[i] || i==DataCount-1){ //if in the segment ending with this point (or if this is the last point extrapolate out
//			float Perc = (strain-strainData[i-1])/(strainData[i]-strainData[i-1]);
//			return stressData[i-1] + Perc*(stressData[i]-stressData[i-1]);
//		}
//	}
//
//	assert(false); //should never reach this point
//	return 0.0f;
//}

float CVX_Material::stress(float strain, float transverseStrainSum, bool forceLinear)
{
	//reference: http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf page 10
	if (isFailed(strain)) return 0.0f; //if a failure point is set and exceeded, we've broken!
	
	if (strain <= strainData[1] || linear || forceLinear){ //for compression/first segment and linear materials (forced or otherwise), simple calculation
		if (nu==0.0f) return E*strain;
		else return eHat()*((1-nu)*strain + nu*transverseStrainSum); 
	}

	//the non-linear feature with non-zero poissons ratio is currently experimental
	int DataCount = modelDataPoints();
	for (int i=2; i<DataCount; i++){ //go through each segment in the material model (skipping the first segment because it has already been handled.
		if (strain <= strainData[i] || i==DataCount-1){ //if in the segment ending with this point (or if this is the last point extrapolate out
			float Perc = (strain-strainData[i-1])/(strainData[i]-strainData[i-1]);
			float basicStress = stressData[i-1] + Perc*(stressData[i]-stressData[i-1]);
			if (nu==0.0f) return basicStress;
			else { //accounting for volumetric effects
				float modulus = (stressData[i]-stressData[i-1])/(strainData[i]-strainData[i-1]);
				float modulusHat = modulus/((1-2*nu)*(1+nu));
				float effectiveStrain = basicStress/modulus; //this is the strain at which a simple linear stress strain line would hit this point at the definied modulus
				float effectiveTransverseStrainSum = transverseStrainSum*(effectiveStrain/strain);
				return modulusHat*((1-nu)*effectiveStrain + nu*effectiveTransverseStrainSum);
			}
		}
	}

	assert(false); //should never reach this point
	return 0.0f;
}


float CVX_Material::modulus(float strain)
{
	if (isFailed(strain)) return 0.0f; //if a failure point is set and exceeded, we've broken!
	if (strain <= strainData[1] || linear) return E; //for compression/first segment and linear materials, simple calculation

	int DataCount = modelDataPoints();
	for (int i=2; i<DataCount; i++){ //go through each segment in the material model (skipping the first segment because it has already been handled.
		if (strain <= strainData[i] || i==DataCount-1) return (stressData[i]-stressData[i-1])/(strainData[i]-strainData[i-1]); //if in the segment ending with this point
	}
	assert(false); //we should never reach this point in the function
	return 0.0f;

}

void CVX_Material::setColor(int red, int green, int blue, int alpha)
{
	setRed(red);
	setGreen(green);
	setBlue(blue);
	setAlpha(alpha);
}

void CVX_Material::setRed(int red)
{
	if (red>255) red=255;
	if (red<0) red=0;
	r = red;
}

void CVX_Material::setGreen(int green)
{
	if (green>255) green=255;
	if (green<0) green=0;
	g = green;
}

void CVX_Material::setBlue(int blue)
{
	if (blue>255) blue=255;
	if (blue<0) blue=0;
	b = blue;
}

void CVX_Material::setAlpha(int alpha)
{
	if (alpha>255) alpha=255;
	if (alpha<0) alpha=0;
	a = alpha;
}

/*! The arrays are assumed to be of equal length.
The first data point is assumed to be [0,0] and need not be provided.
At least 1 non-zero data point must be provided.
The inital segment from [0,0] to the first strain and stress value is interpreted as young's modulus.
The slope of the stress/strain curve should never exceed this value in subsequent segments.
The last data point is assumed to represent failure of the material. The 0.2% offset method is used to calculate the yield point.

Restrictions on pStrainValues:
	- The values must be positive and increasing in order.
	- Strains are defined in absolute numbers according to delta l / L.

Restrictions on pStressValues:
	- The values must be positive and increasing in order.

Special cases:
	- 1 data point (linear): Yield and failure are assumed to occur simultaneously at the single data point.
	- 2 data points (bilinear): Yield is taken as the first data point, failure at the second.

*/
bool CVX_Material::setModel(int dataPointCount, float* pStrainValues, float* pStressValues)
{
	//Pre-checks
	if (*pStrainValues==0 && *pStressValues==0) { //if first data point is 0,0, ignore it
		pStrainValues++; //advance the pointers...
		pStressValues++;
		dataPointCount--; //decrement the count
	}
	if (dataPointCount<=0){
		error = "Not enough data points";
		return false;
	}
	if (*pStrainValues<=0 || *pStressValues<=0){
		error = "First stress and strain data points negative or zero";
		return false; 
	}

	//Copy the data into something more usable (and check for monotonically increasing)
	std::vector<float> tmpStrainData;
	std::vector<float> tmpStressData;
	tmpStrainData.push_back(0); //add in the zero data point (required always)
	tmpStressData.push_back(0);
	float sweepStrain = 0.0f, sweepStress = 0.0f;
	for (int i=0; i<dataPointCount; i++){
		float thisStrain = *(pStrainValues+i); //grab the values
		float thisStress = *(pStressValues+i);

		if (thisStrain <= sweepStrain){
			error = "Out of order strain data";
			return false;
		}

		if (thisStress <= sweepStress){
			error = "Stress data is not monotonically increasing";
			return false;
		}

		if (i>0 && (thisStress-sweepStress)/(thisStrain-sweepStrain) > tmpStressData[0]/tmpStrainData[0]){
			error = "Slope of stress/strain curve should never exceed that of the first line segment (youngs modulus)";
			return false;
		}

		sweepStrain = thisStrain;
		sweepStress = thisStress;

		tmpStrainData.push_back(thisStrain); //add to the temporary vector
		tmpStressData.push_back(thisStress);
	}

	assert(tmpStrainData.size() == dataPointCount+1 && tmpStressData.size() == dataPointCount+1); //sizes match up? (+1 to include the zero point)

	//at this point, we know we have valid data and will return true
	strainData = tmpStrainData;
	stressData = tmpStressData;
	E=stressData[1]/strainData[1]; //youngs modulus is the inital slope
	sigmaFail = stressData[stressData.size()-1]; //failure stress is the highest stress data point
	epsilonFail = strainData[strainData.size()-1]; //failure strain is the highest strain data point
	linear = (dataPointCount==1);

	if (dataPointCount == 1 || dataPointCount == 2){ //linear or blinear
		sigmaYield = stressData[1];
		epsilonYield = strainData[1];
	}
	else { //.2% (0.002) strain offset to find a good yield point...
		setYieldFromData();
	}

	return updateDerived();
}

/*! Specified Young's modulus and failure stress must both be positive.
Yield stress is interpreted as identical to failure stress. If failure stress is not specified an arbitrary data point consistent with the specified Young's modulus is added to the model.
*/
bool CVX_Material::setModelLinear(float youngsModulus, float failureStress)
{
	if (youngsModulus<=0){
		error = "Young's modulus must be positive";
		return false;
	}
	if (failureStress != -1.0f && failureStress<=0){
		error = "Failure stress must be positive";
		return false;
	}

	float tmpfailureStress = failureStress; //create a dummy failure stress if none was provided
	if (tmpfailureStress == -1) tmpfailureStress = 1000000;
	float tmpfailStrain = tmpfailureStress/youngsModulus;

	strainData.clear();
	stressData.clear();
	strainData.push_back(0); //add in the zero data point (required always)
	stressData.push_back(0);
	strainData.push_back(tmpfailStrain);
	stressData.push_back(tmpfailureStress);

	linear=true;
	E=youngsModulus;
	sigmaYield = failureStress; //yield and failure are one in the same here.
	sigmaFail = failureStress;
	epsilonYield = (failureStress == -1) ? -1 : tmpfailStrain;
	epsilonFail = (failureStress == -1) ? -1 : tmpfailStrain;
	return updateDerived();
}

/*! Specified Young's modulus, plastic modulus, yield stress, and failure stress must all be positive.
Plastic modulus must be less than Young's modulus and failure stress must be greater than the yield stress.
*/
bool CVX_Material::setModelBilinear(float youngsModulus, float plasticModulus, float yieldStress, float failureStress)
{
	if (youngsModulus<=0){
		error = "Young's modulus must be positive";
		return false;
	}
	if (plasticModulus<=0 || plasticModulus>=youngsModulus){
		error = "Plastic modulus must be positive but less than Young's modulus";
		return false;
	}
	if (yieldStress<=0){
		error = "Yield stress must be positive";
		return false;
	}
	if (failureStress != -1.0f && failureStress <= yieldStress){
		error = "Failure stress must be positive and greater than the yield stress";
		return false;
	}

	float yieldStrain = yieldStress / youngsModulus;
	float tmpfailureStress = failureStress; //create a dummy failure stress if none was provided
	if (tmpfailureStress == -1) tmpfailureStress = 3*yieldStress;

	float tM = plasticModulus;
	float tB = yieldStress-tM*yieldStrain; //y-mx=b
	float tmpfailStrain = (tmpfailureStress-tB)/tM; // (y-b)/m = x

	strainData.clear();
	strainData.push_back(0); //add in the zero data point (required always)
	strainData.push_back(yieldStrain);
	strainData.push_back(tmpfailStrain);

	stressData.clear();
	stressData.push_back(0);
	stressData.push_back(yieldStress);
	stressData.push_back(tmpfailureStress);

	linear = false;
	E=youngsModulus;
	sigmaYield = yieldStress;
	sigmaFail = failureStress;
	epsilonYield = yieldStrain;
	epsilonFail = failureStress==-1 ? -1 : tmpfailStrain;
	return updateDerived();

}


bool CVX_Material::setYieldFromData(float percentStrainOffset)
{
	sigmaYield = -1.0f; //assume we fail until we succeed.
	epsilonYield = -1.0f; //assume we fail until we succeed.

	float oM = E; //the offset line slope (y=Mx+B)
	float oB = (-percentStrainOffset/100*oM); //offset line intercept (100 factor turns percent into absolute

	float tmpYieldStr = sigmaFail; //assume yield stress is the failure stress value unless we find a better one...
	
	assert(strainData.size() == stressData.size());
	assert(strainData.size() > 2); // more than 2 data points (more than bilinear)
	int dataPoints = strainData.size()-1;
	for (int i=1; i<dataPoints-1; i++){
		float x1=strainData[i];
		float x2=strainData[i+1];
		float y1=stressData[i];
		float y2=stressData[i+1];

		float tM = (y2-y1)/(x2-x1); //temporary slope
		float tB = y1-tM*x1; //temporary intercept

		if (oM!=tM){ //if not parallel lines...
			float xIntersect = (tB-oB)/(oM-tM);
			if (xIntersect>x1 && xIntersect<x2){ //if intersects at this segment...
				float percentBetweenPoints = (xIntersect-x1)/(x2-x1);
				sigmaYield = y1+percentBetweenPoints*(y2-y1);
				epsilonYield = xIntersect;
				return true;
			}
		}
	}
	return false;
}

void CVX_Material::setPoissonsRatio(float poissonsRatio)
{
	if (poissonsRatio < 0) poissonsRatio = 0;
	if (poissonsRatio >= 0.5 ) poissonsRatio = 0.5-FLT_EPSILON*2; //exactly 0.5 will still cause problems, but it can get very close.
	nu = poissonsRatio;
}

void CVX_Material::setDensity(float density)
{
	if (density <= 0) density = FLT_MIN; //density of exactly 0 will cause problems, but can get as close as desired.
	rho = density;
	updateDerived();
}

void CVX_Material::setStaticFriction(float staticFrictionCoefficient)
{
	if (staticFrictionCoefficient <= 0) staticFrictionCoefficient = 0;
	muStatic = staticFrictionCoefficient;
}

void CVX_Material::setKineticFriction(float kineticFrictionCoefficient)
{
	if (kineticFrictionCoefficient <= 0) kineticFrictionCoefficient = 0;
	muKinetic = kineticFrictionCoefficient;
}

void CVX_Material::setInternalDamping(float zeta)
{
	if (zeta <= 0) zeta = 0;
	zetaInternal = zeta;
}

void CVX_Material::setGlobalDamping(float zeta)
{
	if (zeta <= 0) zeta = 0;
	zetaGlobal = zeta;
}

bool CVX_Material::setNominalSize(float size)
{
	if (size <= 0) size = FLT_MIN;
	nomSize=size;
	return updateDerived(); //update derived quantities
}

void CVX_Material::setExternalScaleFactor(Vec3D<> factor)
{
	if (factor.x <= 0) factor.x = FLT_MIN;
	if (factor.y <= 0) factor.y = FLT_MIN;
	if (factor.z <= 0) factor.z = FLT_MIN;
	extScale = factor;
}

bool CVX_Material::updateDerived() 
{
	float volume = nomSize*nomSize*nomSize;
	_mass = volume*rho; 
	_momentInertia = _mass*nomSize*nomSize / 6.0f; //simple 1D approx
	_firstMoment = _mass*nomSize / 2.0f;

	if (volume==0 || _mass==0 || _momentInertia==0){
		_massInverse = _sqrtMass = _momentInertiaInverse = _2xSqMxExS = _2xSqIxExSxSxS = 0.0f; //zero everything out
		return false;
	}

	_massInverse = 1.0f / _mass;
	_sqrtMass = sqrt(_mass);
	_momentInertiaInverse = 1.0f / _momentInertia;
	_2xSqMxExS = 2.0f*sqrt(_mass*E*nomSize);
	_2xSqIxExSxSxS = 2.0f*sqrt(_momentInertia*E*nomSize*nomSize*nomSize);
	return true;
}
