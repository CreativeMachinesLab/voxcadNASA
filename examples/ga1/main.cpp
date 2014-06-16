/*===================================================================
Nick Cheney
2014-06-13
Naive GA Implemelntation: Template for GA design with voxelyze.
(direct encoding hillclimber)
===================================================================*/
#include <iostream>
#include <fstream>
#include <stdlib.h>    
#include <stdio.h>
#include <sstream>
#include "Individual.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Parameters: (best moved to seperate config file.  More params in "Individual.cpp")
int num_x_voxel = 5;
int num_z_voxel = 5; // max size of voxel creature
int num_y_voxel = 5;

int entriesMutated = 4; // average number of genotype entrieds mutated each generation

int keepVxaEvery = 0; // 0:never keep, -1: keep only if new best, 1: keep all VXAs, [n>0]: keep every [n]th generation

int numGenerations = 5;

int randomSeed = 0;
//////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{	
	srand(randomSeed);
	
	/////////////////////////////////////////////////////////////////////
	// Initialize random individual, and assign as current best
	Individual parent(num_x_voxel, num_z_voxel, num_z_voxel);

	parent.writeVxaFile();

	cout << "starting voxelyze... " << endl;
	popen("./voxelize -f genome.vxa","r");

	parent.collectFitness();

	//////////////////////////////////////////////////////////////////////
	// Keep first VXA (if keeping any at all)
	ostringstream mvFileCmd;
	if (keepVxaEvery != 0)
	{
		mvFileCmd << "mv genome.vxa VoxelyzeGA_Gen=" << 0 << "_Fit=" << parent.getFitness() << ".vxa";
	}
	else
	{
		mvFileCmd << "rm genome.vxa";
	}

	FILE* mvFile = popen(mvFileCmd.str().c_str(),"r");
	pclose(mvFile);

	///////////////////////////////////////////////////////////////////////
	// Initialize file to log the fitness over time
	ofstream logging;
	logging.open("champFile.txt");
	logging << "Soft robot evolution with Voxelyze.  Direct-encoding single hillclimber.  Robot size: (" << num_x_voxel << "x" << num_y_voxel << "x" << num_z_voxel << ").  On average, " << entriesMutated << " genes mutationed per generation. Random seed: " << randomSeed  << endl << endl;
	logging << "Gen" << "\t" << "Best Fitness" << "\t" << "Current Gen Fitness" << endl;
	logging << 0 << "\t" << parent.getFitness() << "\t" << parent.getFitness() << endl;
	logging.close();
	/////////////////////////////////////////////////////////////////////////

	// ITERATIVE LOOP
	for (int genNum=1; genNum < numGenerations; genNum++)
	{
		cout << endl;
		/////////////////////////////////////////////////////////////////////
		// Create mutated copy of current best
		Individual child = parent;

		child.mutate(entriesMutated);
		/////////////////////////////////////////////////////////////////////
		// Evaluate child
		child.writeVxaFile();

		cout << "starting voxelyze... " << endl;;
		popen("./voxelize -f genome.vxa","r");

		child.collectFitness();

		//////////////////////////////////////////////////////////////////////
		// Replace parent with child, if child is more fit
		bool keepVxa = false;
		if (child.getFitness() >= parent.getFitness()) // also allow horrizontal mutations
		{
			cout << "NEW BEST! (" << child.getFitness() << " better than " << parent.getFitness() << ")" << endl;
			parent = child;
			if (keepVxaEvery < 0) { keepVxa = true; }
		}
		else
		{
			cout << "no good... (" << child.getFitness() << " worse than " << parent.getFitness() << ")" << endl;
		}
		//////////////////////////////////////////////////////////////////////
		// Keep or remove Vxa Files
		if (keepVxaEvery > 0 and genNum % keepVxaEvery == 0) { keepVxa = true;}

		ostringstream mvFileCmd;
		if (keepVxa)
		{
			mvFileCmd << "mv genome.vxa VoxelyzeGA_Gen=" << genNum << "_Fit=" << child.getFitness() << ".vxa";
		}
		else
		{
			mvFileCmd << "rm genome.vxa";
		}

		FILE* mvFile = popen(mvFileCmd.str().c_str(),"r");
		pclose(mvFile);

		///////////////////////////////////////////////////////////////////////
		// Record fitness of current generation
		ofstream logging;
		logging.open("champFile.txt", ios_base::app);
		logging << genNum << "\t" << parent.getFitness() << "\t" << child.getFitness() << endl;
		logging.close();
		/////////////////////////////////////////////////////////////////////////
	}
}