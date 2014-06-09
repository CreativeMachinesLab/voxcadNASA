/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#ifndef VX_VOXMATERIAL_H
#define VX_VOXMATERIAL_H

#include "VX_Material.h"

//!< VoxType is a VX_Material + voxel size.
class CVX_VoxMaterial : public CVX_Material {
	private:
	CVX_VoxMaterial(); //!< Default Constructor
	virtual ~CVX_VoxMaterial(void); //!< Destructor. Virtual so we can just keep track of CVX_Material pointers.
	CVX_VoxMaterial(const CVX_VoxMaterial& vIn) {*this = vIn;} //!< Copy constructor
	CVX_VoxMaterial& operator=(const CVX_VoxMaterial& vIn); //!< Equals operator

	bool setNominalSize(float nominalSize);

//	float nominalSize();
//	float mass();
//	float inertia();
//	float massInverse();
//	float inertiaInverse();

//	private:
	float nomSize; //nominal size (i.e. lattice dimension) (m)

	//derived quantities to cache
	bool updateDerived(); //updates all the derived quantities cache
	float _mass; //The mass of this voxel (kg)
	float _massInverse; //1/Mass (1/kg)
	float _firstMoment; //1st moment "inertia" (needed for certain calculations) (kg*m)
	float _momentInertia; //mass moment of inertia (i.e. rotational "mass") (kg*m^2)
	float _momentInertiaInverse; //1/Inertia (1/(kg*m^2))
	float _2xSqMxExS; //needed for quick damping calculations (Kg*m/s)
	float _2xSqIxExSxSxS; //needed for quick rotational damping calculations (Kg*m/s^2)


	//catch the virtual functions from CVX_Material indicated relevant parameters changed
	void changedE(){updateDerived();}
	void changedRho(){updateDerived();}

	friend class CVoxelyze; //give the main simulation class full access
	friend class CVX_Voxel; //give our voxel class direct access to all the members for quick access
};



#endif //VX_VOXMATERIAL_H