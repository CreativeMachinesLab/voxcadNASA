/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "Voxelyze.h"
#include "VX_Material.h"
#include "VX_Voxel.h"
#include "VX_Link.h"
#include <assert.h>

#define DEFAULT_LATTICESIZE 0.001 //1mm

CVoxelyze::CVoxelyze(void)
{
	latSize=DEFAULT_LATTICESIZE;
}

CVoxelyze::~CVoxelyze(void)
{
	clear();
}

bool CVoxelyze::doTimeStep(float dt)
{
	if (dt==0) return true;
	else if (dt<0) dt = recommendedTimeStep();

	//Euler integration:
	bool Diverged = false;
//#pragma omp parallel for
	for (std::list<CVX_Link*>::const_iterator it=linksList.begin(); it != linksList.end(); it++){ //for each link
		(*it)->updateForces();
		if ((*it)->axialStrain() > 100) Diverged = true; //catch divergent condition! (if any thread sets true we will fail, so don't need mutex...
	}
	if (Diverged) return false;

//#pragma omp parallel for
//	for (int i=0; i<iT; i++) BondArrayCollision[i].UpdateBond();

//#pragma omp parallel for
	for (std::list<CVX_Voxel*>::iterator it=voxelsList.begin(); it != voxelsList.end(); it++){
		(*it)->timeStep(dt);
	}

	currentTime += dt;
	return true;
}

float CVoxelyze::recommendedTimeStep() const
{
	//find the largest natural frequency (sqrt(k/m)) that anything in the simulation will experience, then multiply by 2*pi and invert to get the optimally largest timestep that should retain stability
	float MaxFreq2 = 0.0f; //maximum frequency in the simulation

	for (std::list<CVX_Link*>::const_iterator it=linksList.begin(); it != linksList.end(); it++){ //for each link
		CVX_Link* pL = (*it);
		float m1 = pL->pVNeg->mat->_mass,  m2 = pL->pVPos->mat->_mass;
		float thisMaxFreq2 = pL->axialStiffness()/(m1<m2?m1:m2);
		if (thisMaxFreq2 > MaxFreq2) MaxFreq2 = thisMaxFreq2;
	}
	for (std::vector<CVX_Material*>::const_iterator it=matList.begin(); it != matList.end(); it++){
		float thisMaxFreq2 = (*it)->penetrationStiffness()/(*it)->_mass;
		if (thisMaxFreq2 > MaxFreq2) MaxFreq2 = thisMaxFreq2;
	}

	if (MaxFreq2 <= 0.0f) return 0.0f;
	else return 1.0/(sqrt(MaxFreq2)*2.0*3.1415926f); //convert to time... (seconds)
}

void CVoxelyze::resetTime()
{
	currentTime=0.0f;

	for (int iz=indexMinZ(); iz<=indexMaxZ(); iz++){
		for (int iy=indexMinY(); iy<=indexMaxY(); iy++){
			for (int ix=indexMinX(); ix<=indexMaxX(); ix++){
				CVX_Voxel* pV = voxels(ix, iy, iz);
				if (pV){
					pV->pos = Vec3D<vfloat>(ix*latSize, iy*latSize, iz*latSize); //set position to nominal
					pV->orient = CQuat<>(); //no angular displacement
					pV->haltMotion(); //zeros linMom and angMom
					pV->setFloorStaticFriction(false);
					pV->tempRelative=0.0f;
					pV->previousDt=0.0f;
				}
			}
		}
	}

	for (std::list<CVX_Link*>::iterator it=linksList.begin(); it != linksList.end(); it++){ //for each link
		CVX_Link* pL = (*it);
		pL->updateProperties();
	}
}

void CVoxelyze::clear() //deallocates and returns everything to defaults
{
	//delete and remove links
	for (std::list<CVX_Link*>::iterator it = linksList.begin(); it!=linksList.end(); it++) delete *it;
	for (int i=0; i<3; i++)	links[i].clear();
	linksList.clear();

	//delete and remove voxels
	for (std::list<CVX_Voxel*>::iterator it = voxelsList.begin(); it!=voxelsList.end(); it++) delete *it;
	voxelsList.clear();
	voxels.clear();

	//delete and remove materials
	for (int i=0; i<matList.size(); i++) delete matList[i];
	matList.clear();

	latSize = DEFAULT_LATTICESIZE;
	currentTime = 0.0f;
}

CVX_Material* CVoxelyze::addMaterial()
{
	CVX_Material* pMat;
	try {
		pMat = new CVX_Material();
		matList.push_back(pMat);
	}
	catch (std::bad_alloc&){
		return NULL;
	}
	return pMat; 
}

bool CVoxelyze::removeMaterial(CVX_Material* toRemove)
{
	int matIndex = exists(toRemove);
	if (matIndex == -1) return false;

	//remove all voxels that use this material
	for (int k=indexMinZ(); k<=indexMaxZ(); k++){
		for (int j=indexMinY(); j<=indexMaxY(); j++){
			for (int i=indexMinX(); i<=indexMaxX(); i++){
				if (voxels(i, j, k)->material() == toRemove) removeVoxel(i, j, k);
			}
		}
	}

	//remove the material
	delete toRemove;
	matList.erase(matList.begin()+matIndex); //remove from list
	assert(exists(toRemove) == -1); //the material should no longer exist.

	return true;
}

bool CVoxelyze::replaceMaterial(CVX_Material* replaceMe, CVX_Material* replaceWith)
{
	if (!exists(replaceMe) || !exists(replaceWith)) return false;
	
	//switch all voxel references
	for (std::list<CVX_Voxel*>::iterator it = voxelsList.begin(); it!=voxelsList.end(); it++){ //remove from the list
		if ((*it)->material()==replaceMe) (*it)->replaceMaterial(replaceWith);
	}
	return true;
}


CVX_Voxel* CVoxelyze::addVoxel(CVX_Material* newVoxelMaterial, int xIndex, int yIndex, int zIndex) //creates a new voxel if there isn't one here. Otherwise
{
	//if "adding" a null material, we're actually deleting.
	if (newVoxelMaterial == NULL){
		removeVoxel(xIndex, yIndex, zIndex);
		return NULL;
	}

	//if there's already a voxel here, just replace it with the new material
	CVX_Voxel* pV = voxels(xIndex, yIndex, zIndex);
	if (pV != NULL)	pV->replaceMaterial(newVoxelMaterial);
	else { //make a new voxel!
		try {
			pV = new CVX_Voxel(newVoxelMaterial);
			voxels.addValue(xIndex, yIndex, zIndex, pV); //add to the array
			voxelsList.push_back(pV);
			pV->pos = Vec3D<vfloat>(xIndex*latSize, yIndex*latSize, zIndex*latSize); //set initial voxel location (extrapolate?)
		}
		catch (std::bad_alloc&){
			return NULL;
		}

		//add any possible links utilizing this voxel
		for (int i=0; i<6; i++){ //from X_POS to Z_NEG (0-5 enums)
			addLink(xIndex, yIndex, zIndex, (linkDirection)i); 
		}


	}
	return pV;
}

void CVoxelyze::removeVoxel(int xIndex, int yIndex, int zIndex)
{
	const CVX_Voxel* pV = voxel(xIndex, yIndex, zIndex);
	if (pV==NULL) return; //no voxel exists here.
	delete pV;
	voxels.removeValue(xIndex, yIndex, zIndex); //remove from the array
	for (std::list<CVX_Voxel*>::iterator it = voxelsList.begin(); it!=voxelsList.end(); it++){ //remove from the list
		if (*it == pV) voxelsList.erase(it);
	}

	//make sure no references are left in the list This should be compiled away in release
	for (std::list<CVX_Voxel*>::iterator it = voxelsList.begin(); it!=voxelsList.end(); it++) assert(*it != pV); 

	//remove any links to this voxel
	for (int i=0; i<6; i++){ //from X_POS to Z_NEG (0-5 enums)
		removeLink(xIndex, yIndex, zIndex, (linkDirection)i); 
	}

	//remove any collisions involving this voxel
}

CVX_Link* CVoxelyze::link(int xIndex, int yIndex, int zIndex, linkDirection direction) const
{
	return links[toAxis(direction)](
		xIndex+xIndexLinkOffset(direction),
		yIndex+yIndexLinkOffset(direction),
		zIndex+zIndexLinkOffset(direction));
}

CVX_Link* CVoxelyze::addLink(int xIndex, int yIndex, int zIndex, linkDirection direction)
{
	CVX_Link* pL = link(xIndex, yIndex, zIndex, direction);
	if (pL){return pL;} //if a link already exists... well, then it should be up to date.

	//ensure that there are voxels at both ends of the link
	CVX_Voxel* voxel1 = voxels(xIndex, yIndex, zIndex);
	CVX_Voxel* voxel2 = voxels(
		xIndex+xIndexVoxelOffset(direction),
		yIndex+yIndexVoxelOffset(direction),
		zIndex+zIndexVoxelOffset(direction));
	if (voxel1 == NULL || voxel2 == NULL) return NULL; //if no voxel at either position, don't make a link

	//make the link and add it to the array+list
	try {
		pL = new CVX_Link(voxel1, voxel2, direction);	//make the new link
		linksList.push_back(pL);							//add to the list
		links[toAxis(direction)].addValue(
			xIndex + xIndexLinkOffset(direction),
			yIndex + yIndexLinkOffset(direction),
			zIndex + zIndexLinkOffset(direction), pL);
	}
	catch (std::bad_alloc&){
		return NULL;
	}
	//Add reference to this link to the relevant voxels
	voxel1->addLinkInfo(direction, pL);
	voxel2->addLinkInfo(toOpposite(direction), pL);
}

void CVoxelyze::removeLink(int xIndex, int yIndex, int zIndex, linkDirection direction)
{
	const CVX_Link* pL = link(xIndex, yIndex, zIndex, direction);
	if (pL==NULL) return; //no link here to see!

	//remove the reference in the appropriate link 3d array
	links[linkAxis(direction)].removeValue( 
		xIndex + xIndexLinkOffset(direction),
		yIndex + yIndexLinkOffset(direction),
		zIndex + zIndexLinkOffset(direction)); 

	//remove the reference in the list
	for (std::list<CVX_Link*>::iterator it = linksList.begin(); it!=linksList.end(); it++){ //remove from the list
		if (*it == pL) linksList.erase(it);
	}
	
	//make sure no references are left in the list. This should be compiled away in release
	for (std::list<CVX_Link*>::iterator it = linksList.begin(); it!=linksList.end(); it++) assert(*it != pL); 

	//remove the reference to this link from one voxel (if it exists)
	CVX_Voxel* voxel1 = voxels(xIndex, yIndex, zIndex);
	if (voxel1) voxel1->removeLinkInfo(direction);

	//remove the reference to this link from the other voxel (if it exists)
	CVX_Voxel* voxel2 = voxels(
		xIndex+xIndexVoxelOffset(direction),
		yIndex+yIndexVoxelOffset(direction),
		zIndex+zIndexVoxelOffset(direction));
	if (voxel2) voxel2->removeLinkInfo(toOpposite(direction));

	delete pL;
}


int CVoxelyze::exists(const CVX_Material* toCheck) //returns the index in materialList if a material exists at this location, -1 otherwise.
{
	int voxMatCount = matList.size();
	for (int i=0; i<voxMatCount; i++){
		if (matList[i] == toCheck) return i;
	}
	return -1;
}

void CVoxelyze::setLatticeSize(float latticeSize) //sets the voxel size.
{
	float scaleFactor = latticeSize/latSize; //scaling factor
	latSize = latticeSize;
	
	//update materials
	for (std::vector<CVX_Material*>::iterator it=matList.begin(); it != matList.end(); it++){
		(*it)->setNominalSize(latticeSize);
	}

	//update voxels
	for (std::list<CVX_Voxel*>::iterator it=voxelsList.begin(); it != voxelsList.end(); it++){
		CVX_Voxel* pV = (*it);
		pV->pos *= scaleFactor;
		pV->haltMotion(); //stop motion to avoid weird huge kinetic energy disparities
		pV->setFloorStaticFriction(false);
	}

	//update links
	for (std::list<CVX_Link*>::iterator it=linksList.begin(); it != linksList.end(); it++){
		CVX_Link* pL = (*it);
		pL->updateProperties();
	}
}