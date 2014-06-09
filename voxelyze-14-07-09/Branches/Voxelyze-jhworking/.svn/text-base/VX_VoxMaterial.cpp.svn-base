/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VX_VoxMaterial.h"


CVX_VoxMaterial::CVX_VoxMaterial()
{
	//CVX_Material constructor called first

	nomSize=1.0;
	updateDerived();
}

CVX_VoxMaterial::~CVX_VoxMaterial(void)
{

}

CVX_VoxMaterial& CVX_VoxMaterial::operator=(const CVX_VoxMaterial& vIn)
{
	CVX_Material::operator=(vIn); //set base CVX_Material class variables equal

	nomSize=vIn.nomSize;
	_mass=vIn._mass;
	_massInverse=vIn._massInverse;
	_firstMoment=vIn._firstMoment;
	_momentInertia=vIn._momentInertia;
	_momentInertiaInverse=vIn._momentInertiaInverse;
	_2xSqMxExS=vIn._2xSqMxExS;
	_2xSqIxExSxSxS=vIn._2xSqIxExSxSxS;

	return *this;
}

bool CVX_VoxMaterial::setNominalSize(float nominalSize)
{
	nomSize=nominalSize;
	return updateDerived(); //update derived quantities
}

bool CVX_VoxMaterial::updateDerived()
{
	float volume = nomSize*nomSize*nomSize;
	_mass = volume*rho; 
	_momentInertia = _mass*nomSize*nomSize / 6.0; //simple 1D approx
	_firstMoment = _mass*nomSize / 2.0;

	if (volume==0 || _mass==0 || _momentInertia==0){
		_massInverse = _momentInertiaInverse = _2xSqMxExS = _2xSqIxExSxSxS = 0; //zero everything out
		return false;
	}

	_massInverse = 1.0 / _mass;
	_momentInertiaInverse = 1 / _momentInertia;
	_2xSqMxExS = 2*sqrt(_mass*E*nomSize);
	_2xSqIxExSxSxS = 2*sqrt(_momentInertia*E*nomSize*nomSize*nomSize);
}
