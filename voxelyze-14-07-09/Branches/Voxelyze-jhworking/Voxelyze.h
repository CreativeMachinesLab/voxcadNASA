//*******************************************************************************
//Copyright (c) 2014, Jonathan Hiller (Cornell University)
//If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"
//
//This file is part of Voxelyze.
//Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
//*******************************************************************************/

#ifndef VOXELYZE_H
#define VOXELYZE_H

#include "VX_Enums.h"
#include "Array3D.h"
#include <vector> //delete if PIMPL'd
#include <list> //delete if PIMPL'd


class CVX_Material;
class CVX_Voxel;
class CVX_Link;
class CVX_Collision;

class CVoxelyze {
public:
	CVoxelyze(void);
	~CVoxelyze(void);
	CVoxelyze(const CVoxelyze& VIn) {*this = VIn;} //copy constructor
	CVoxelyze& operator=(const CVoxelyze& VIn); //equal operator

	void clear(); //deallocates and returns everything to defaults

	bool doTimeStep(float dt = -1.0f); //timestep. -1.0f = calculate the recommended timestep
	float recommendedTimeStep() const; //returns recommended time step
	void resetTime(); //resets simulation

	CVX_Material* addMaterial();
	bool removeMaterial(CVX_Material* toRemove); //if replaceWith is not null, replace all uses of toRemove material with replaceWith.
	bool replaceMaterial(CVX_Material* replaceMe, CVX_Material* replaceWith);

//	CVX_Material* material(int index);
//	int materialCount() const; //returns the total number of materials
//	CVX_Material* materialList() const; //returns pointer to the beginning of a null-terminated list of material pointers

	CVX_Voxel* voxel(int xIndex, int yIndex, int zIndex) const {return voxels(xIndex, yIndex, zIndex);}
	const std::list<CVX_Voxel*>* voxelList() const {return &voxelsList;}

	CVX_Voxel* addVoxel(CVX_Material* newVoxelMaterial, int xIndex, int yIndex, int zIndex); //creates a new voxel if there isn't one here. Otherwise
	void removeVoxel(int xIndex, int yIndex, int zIndex);

	int indexMinX() const {return voxels.minIndices().x;}
	int indexMaxX() const {return voxels.maxIndices().x;}
	int indexMinY() const {return voxels.minIndices().y;}
	int indexMaxY() const {return voxels.maxIndices().y;}
	int indexMinZ() const {return voxels.minIndices().z;}
	int indexMaxZ() const {return voxels.maxIndices().z;}

//	CVX_Voxel* voxelsList() const; //null-terminated list
//	int voxelCount() const; //total voxels
//	int voxelCount(CVX_Material* pMaterialToCount) const; //count of voxels of this material

	//hide?
	CVX_Link* link(int xIndex, int yIndex, int zIndex, linkDirection direction) const;
	const std::list<CVX_Link*>* linkList() const {return &linksList;}
//	CVX_Link* LinkList() const;
//	int LinkCount() const;

	CVX_Collision* CollisionList() const;
	int CollisionCount() const;

	void setLatticeSize(float latticeSize); //sets the voxel size.
	float latticeSize() const {return latSize;}

private:
	float latSize; //lattice size
	float currentTime; //current time of the simulation in seconds

	std::vector<CVX_Material*> matList; //up-to-date list of all voxel materials existing in this simulation
	int exists(const CVX_Material* toCheck); //returns the index in materialList if a material exists at this location, -1 otherwise.

//	CVX_Voxel* voxel(int xIndex, int yIndex, int zIndex);
//	CVX_Link* link(int xIndex, int yIndex, int zIndex, linkDirection direction) const;

	CArray3D<CVX_Voxel*> voxels; //main voxel array 3D lookup
	std::list<CVX_Voxel*> voxelsList; //main list of existing voxels (no particular order) (always kept syncd with voxels)

	CArray3D<CVX_Link*> links[3]; //main link arrays in the X[0], Y[1] and Z[2] directions. (0,0,0) is the bond pointting in the positive direction from voxel (0,0,0)
	std::list<CVX_Link*> linksList; //main list of all existing links (no particular order) (always kept syncd with voxels)

	CVX_Link* addLink(int xIndex, int yIndex, int zIndex, linkDirection direction); //adds a link (if one isn't already present) and updates parameters
	void removeLink(int xIndex, int yIndex, int zIndex, linkDirection direction); //removes just the link and all references to it in connected voxels

	int xIndexLinkOffset(linkDirection direction) const {return (direction == X_NEG) ? -1 : 0;} //the link X index offset from a voxel index and link direction
	int yIndexLinkOffset(linkDirection direction) const {return (direction == Y_NEG) ? -1 : 0;} //the link Y index offset from a voxel index and link direction
	int zIndexLinkOffset(linkDirection direction) const {return (direction == Z_NEG) ? -1 : 0;} //the link Z index offset from a voxel index and link direction
	int xIndexVoxelOffset(linkDirection direction) const {return (direction == X_NEG) ? -1 : ((direction == X_POS) ? 1 : 0);} //the voxel X index offset of a voxel across a link in the specified direction
	int yIndexVoxelOffset(linkDirection direction) const {return (direction == Y_NEG) ? -1 : ((direction == Y_POS) ? 1 : 0);} //the voxel Y index offset of a voxel across a link in the specified direction
	int zIndexVoxelOffset(linkDirection direction) const {return (direction == Z_NEG) ? -1 : ((direction == Z_POS) ? 1 : 0);} //the voxel Z index offset of a voxel across a link in the specified direction
	//updateCollisionList()
};


#endif //VOXELYZE_H