/*===================================================================
Nick Cheney
2014-06-13
Naive GA Implemelntation: Template for GA design with voxelyze.
(direct encoding hillclimber)

Individual.cpp
===================================================================*/
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <math.h>
#include "Individual.h"

////////////////////////////////////////////////
// Parameters: (best moved to seperate config file.  More params in "main.cpp")
float simulationTime = 0.5; // "Temp Period" is length of actuation cycle.  With 0.025 Temp Period (below), 0.5 sec = 20 actuation cycles.
bool usingStiffMaterial = false;  // runs 3-4 times faster with no stiff material (soft materials that vary too much in stiffness cause ill-conditioned matrices in Voxelyze)
float evaluationTimeOut = 60.0;  // timeout for voxelyze evaluation (can hang in certain conditions), the best length to wait will depend heavily on the size, material, and evaluation sim time of the robots evolved.
///////////////////////////////////////////////

Individual::Individual(int _num_x_voxels, int _num_y_voxels, int _num_z_voxels)
{
	fitness = -999;
	num_x_voxels = _num_x_voxels;
	num_y_voxels = _num_y_voxels;
	num_z_voxels = _num_z_voxels;
	genotypeLength = num_x_voxels * num_y_voxels * num_z_voxels;
	for (int i=0; i<genotypeLength; i++) 
	{ 
		genotypeShape.push_back(rand()%2);
		genotypeMuscleOrSupport.push_back(rand()%2);
		genotypeMuscleType.push_back(rand()%2);
		genotypeSupportType.push_back(rand()%2); 
	}
}

void Individual::setFitness(float _fitness)
{
	fitness = _fitness;
}

float Individual::getFitness()
{
	return fitness;
}

// note: this genome is more efficiently implemented as a matrix, rather than 4 vectors.  But seeing how this is a tutorial, the naming is more straightforward with this implementation
void Individual::setGenotypeShapeAt(int _index, int _value){ genotypeShape[_index] = _value; }
void Individual::setGenotypeMuscleOrSupportAt(int _index, int _value){ genotypeMuscleOrSupport[_index] = _value; }
void Individual::setGenotypeMuscleTypeAt(int _index, int _value){ genotypeMuscleType[_index] = _value; }
void Individual::setGenotypeSupportTypeAt(int _index, int _value){ genotypeSupportType[_index] = _value; }

vector<int> Individual::getGenotypeShape() { return genotypeShape; }
vector<int> Individual::getGenotypeMuscleOrSupport() { return genotypeMuscleOrSupport; }
vector<int> Individual::getGenotypeMuscleType() { return genotypeMuscleType; }
vector<int> Individual::getGenotypeSupportType() { return genotypeSupportType; }

// note, to best find values/meaning of unknown parameters below, grep in "Voxleyze" folder, or manually adjust in VoxCad GUI.  Default values here correspond to those used in GECCO 2013 paper (besides varaible stiffness of hard support)
// If running more than one thread in same filesystem, change output file name to a unique individual ID with "FitnessFileName"
void Individual::writeVxaFile()
{
  	ofstream myfile;
	myfile.open ("genome.vxa");
	
	myfile << "\
<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n\
<VXA Version=\"1.0\">\n\
<Simulator>\n\
<Integration>\n\
<Integrator>0</Integrator>\n\
<DtFrac>0.9</DtFrac>\n\
</Integration>\n\
<Damping>\n\
<BondDampingZ>1</BondDampingZ>\n\
<ColDampingZ>0.8</ColDampingZ>\n\
<SlowDampingZ>0.01</SlowDampingZ>\n\
</Damping>\n\
<Collisions>\n\
<SelfColEnabled>1</SelfColEnabled>\n\
<ColSystem>3</ColSystem>\n\
<CollisionHorizon>2</CollisionHorizon>\n\
</Collisions>\n\
<Features>\n\
<FluidDampEnabled>0</FluidDampEnabled>\n\
<PoissonKickBackEnabled>0</PoissonKickBackEnabled>\n\
<EnforceLatticeEnabled>0</EnforceLatticeEnabled>\n\
</Features>\n\
<SurfMesh>\n\
<CMesh>\n\
<DrawSmooth>1</DrawSmooth>\n\
<Vertices/>\n\
<Facets/>\n\
<Lines/>\n\
</CMesh>\n\
</SurfMesh>\n\
<StopCondition>\n\
<StopConditionType>2</StopConditionType>\n\
<StopConditionValue>" << simulationTime << "</StopConditionValue>\n\
</StopCondition>\n\
<GA>\n\
<WriteFitnessFile>1</WriteFitnessFile>\n\
<FitnessFileName> fitness.xml </FitnessFileName>\n\
</GA>\n\
</Simulator>\n\
<Environment>\n\
<Fixed_Regions>\n\
<NumFixed>0</NumFixed>\n\
</Fixed_Regions>\n\
<Forced_Regions>\n\
<NumForced>0</NumForced>\n\
</Forced_Regions>\n\
<Gravity>\n\
<GravEnabled>1</GravEnabled>\n\
<GravAcc>-27.468</GravAcc>\n\
<FloorEnabled>1</FloorEnabled>\n\
</Gravity>\n\
<Thermal>\n\
<TempEnabled>1</TempEnabled>\n\
<TempAmp>39</TempAmp>\n\
<TempBase>25</TempBase>\n\
<VaryTempEnabled>1</VaryTempEnabled>\n\
<TempPeriod>0.025</TempPeriod>\n\
</Thermal>\n\
</Environment>\n\
<VXC Version=\"0.93\">\n\
<Lattice>\n\
<Lattice_Dim> 0.001</Lattice_Dim>\n\
<X_Dim_Adj>1</X_Dim_Adj>\n\
<Y_Dim_Adj>1</Y_Dim_Adj>\n\
<Z_Dim_Adj>1</Z_Dim_Adj>\n\
<X_Line_Offset>0</X_Line_Offset>\n\
<Y_Line_Offset>0</Y_Line_Offset>\n\
<X_Layer_Offset>0</X_Layer_Offset>\n\
<Y_Layer_Offset>0</Y_Layer_Offset>\n\
</Lattice>\n\
<Voxel>\n\
<Vox_Name>BOX</Vox_Name>\n\
<X_Squeeze>1</X_Squeeze>\n\
<Y_Squeeze>1</Y_Squeeze>\n\
<Z_Squeeze>1</Z_Squeeze>\n\
</Voxel>\n\
<Palette>\n";
	
	myfile << "\
<Material ID=\"1\">\n\
<MatType>0</MatType>\n\
<Name>Soft</Name>\n\
<Display>\n\
<Red>0</Red>\n\
<Green>1</Green>\n\
<Blue>1</Blue>\n\
<Alpha>1</Alpha>\n\
</Display>\n\
<Mechanical>\n\
<IsConductive>1</IsConductive>\n\
<MatModel>0</MatModel>\n\
<Elastic_Mod>1e+007</Elastic_Mod>\n\
<Plastic_Mod>0</Plastic_Mod>\n\
<Yield_Stress>0</Yield_Stress>\n\
<FailModel>0</FailModel>\n\
<Fail_Stress>0</Fail_Stress>\n\
<Fail_Strain>0</Fail_Strain>\n\
<Density>2e+006</Density>\n\
<Poissons_Ratio>0.35</Poissons_Ratio>\n\
<CTE>0.01</CTE>\n\
<uStatic>1</uStatic>\n\
<uDynamic>0.5</uDynamic>\n\
</Mechanical>\n\
</Material>\n";
	
	myfile << "\
<Material ID=\"2\">\n\
<MatType>0</MatType>\n\
<Name>Hard</Name>\n\
<Display>\n\
<Red>0</Red>\n\
<Green>0</Green>\n\
<Blue>1</Blue>\n\
<Alpha>1</Alpha>\n\
</Display>\n\
<Mechanical>\n\
<IsConductive>1</IsConductive>\n\
<MatModel>0</MatModel>\n\
<Elastic_Mod> 1e+008 </Elastic_Mod>\n\
<Plastic_Mod>0</Plastic_Mod>\n\
<Yield_Stress>0</Yield_Stress>\n\
<FailModel>0</FailModel>\n\
<Fail_Stress>0</Fail_Stress>\n\
<Fail_Strain>0</Fail_Strain>\n\
<Density>2e+006</Density>\n\
<Poissons_Ratio>0.35</Poissons_Ratio>\n\
<CTE>0.01</CTE>\n\
<uStatic>1</uStatic>\n\
<uDynamic>0.5</uDynamic>\n\
</Mechanical>\n\
</Material>\n";
	
	myfile << "\
<Material ID=\"3\">\n\
<MatType>0</MatType>\n\
<Name>Active_-</Name>\n\
<Display>\n\
<Red>0</Red>\n\
<Green>1</Green>\n\
<Blue>0</Blue>\n\
<Alpha>1</Alpha>\n\
</Display>\n\
<Mechanical>\n\
<IsConductive>1</IsConductive>\n\
<MatModel>0</MatModel>\n\
<Elastic_Mod>1e+007</Elastic_Mod>\n\
<Plastic_Mod>0</Plastic_Mod>\n\
<Yield_Stress>0</Yield_Stress>\n\
<FailModel>0</FailModel>\n\
<Fail_Stress>0</Fail_Stress>\n\
<Fail_Strain>0</Fail_Strain>\n\
<Density>2e+006</Density>\n\
<Poissons_Ratio>0.35</Poissons_Ratio>\n\
<CTE>-0.01</CTE>\n\
<uStatic>1</uStatic>\n\
<uDynamic>0.5</uDynamic>\n\
</Mechanical>\n\
</Material>\n";
	
	myfile << "\
<Material ID=\"4\">\n\
<MatType>0</MatType>\n\
<Name>Active_+</Name>\n\
<Display>\n\
<Red>1</Red>\n\
<Green>0</Green>\n\
<Blue>0</Blue>\n\
<Alpha>1</Alpha>\n\
</Display>\n\
<Mechanical>\n\
<IsConductive>1</IsConductive>\n\
<MatModel>0</MatModel>\n\
<Elastic_Mod>1e+007</Elastic_Mod>\n\
<Plastic_Mod>0</Plastic_Mod>\n\
<Yield_Stress>0</Yield_Stress>\n\
<FailModel>0</FailModel>\n\
<Fail_Stress>0</Fail_Stress>\n\
<Fail_Strain>0</Fail_Strain>\n\
<Density>2e+006</Density>\n\
<Poissons_Ratio>0.35</Poissons_Ratio>\n\
<CTE>0.01</CTE>\n\
<uStatic>1</uStatic>\n\
<uDynamic>0.5</uDynamic>\n\
</Mechanical>\n\
</Material>\n";

	myfile << "\
</Palette>\n\
<Structure Compression=\"ASCII_READABLE\">\n\
<X_Voxels>" << num_x_voxels << "</X_Voxels>\n\
<Y_Voxels>" << num_y_voxels << "</Y_Voxels>\n\
<Z_Voxels>" << num_z_voxels << "</Z_Voxels>\n\
<Data>\n\
<Layer><![CDATA[";

	int v=0;
	for (int z=0; z<num_z_voxels; z++)
	{
		for (int y=0; y<num_y_voxels; y++)
		{
			for (int x=0; x<num_x_voxels; x++)
			{
				if (genotypeShape[v] > 0)
				{
					if (genotypeMuscleOrSupport[v] > 0)
					{
						if (genotypeMuscleType[v] > 0)
						{
							myfile << "3";
						}
						else
						{
							myfile << "4";
						}
					}
					else
					{
						if (genotypeSupportType[v] > 0 or not usingStiffMaterial)
						{
							myfile << "1";
						}
						else
						{
							myfile << "2";
						}
					}
				}
				else
				{
					myfile << "0";
				}
				v++;
			}
		}
		if (z<num_z_voxels-1) { myfile << "]]></Layer>\n<Layer><![CDATA["; }
	}

	myfile << "]]></Layer>\n\
</Data>\n\
</Structure>\n\
</VXC>\n\
</VXA>";
  		
  	myfile.close();

}

void Individual::collectFitness()
{
	clock_t startEvaluationTime = clock();
	while (true) // wait for voxelyze to finish
	{
		ifstream infile("fitness.xml"); // find and record fitness value using line by line file parsing.  Could also use xml parsing.
		if (infile.is_open())
		{
			sleep(1);  // if file is not done being written before it is acessed, it will not read properly.  To increase efficiency: write to tmp.xml then have voxelyze call "mv tmp.xml fitness.xml" after file write is completed.
			string line;
			while (std::getline(infile, line))
			{
				std::size_t foundx = line.find("<FinalCOM_Dist>"); // name of returned fitness (change/add within voxelyze "SimGA" file.  In other branch, this is editied to not include the z displacement, which I believe it does here)
			    if (foundx!=std::string::npos)
			    {
			    	fitness = fabs(atof(line.substr(foundx+strlen("<FinalCOM_Dist>"),line.find("</")-(foundx+strlen("<FinalCOM_Dist>"))).c_str()));
			    	fitness = fitness;
			    }
			}

			cout << "voxelyze finished in " << float(clock() - startEvaluationTime)/CLOCKS_PER_SEC <<  " seconds." << endl;
			infile.close();
			FILE* rmCmd = popen("rm fitness.xml","r");
			pclose(rmCmd);
			break;
		}

		if (float(clock() - startEvaluationTime)/CLOCKS_PER_SEC  > evaluationTimeOut) // just in case voxelyze crashes, catch it and continue 
		{
			cout << "Voxelyze timed out" << endl;
			fitness = 0;
			break;

		}
	}
}

void Individual::mutate(int _entriesMutated)
{
	float mutationRate = (float) _entriesMutated / ((float)genotypeLength * 4); // sets mutation rate to target a number of mutated entries
	bool mutated = false;
	while (not mutated) // ensures child is not exact copy of parent
	{
		for (int i=0; i<genotypeLength; i++)
		{
			if ((double) rand() / RAND_MAX < mutationRate)
			{
				setGenotypeShapeAt( i, (getGenotypeShape()[i] + 1 ) % 2 ); // flip bit
				mutated = true;
			}

			if ((double) rand() / RAND_MAX < mutationRate)
			{
				setGenotypeMuscleOrSupportAt( i, (getGenotypeMuscleOrSupport()[i] + 1 ) % 2 ); // flip bit
				mutated = true;
			}

			if ((double) rand() / RAND_MAX < mutationRate)
			{
				setGenotypeMuscleTypeAt( i, (getGenotypeMuscleType()[i] + 1 ) % 2 ); // flip bit
				mutated = true;
			}

			if (usingStiffMaterial) // don't count changes here if not using the stiff material
			{
				if ((double) rand() / RAND_MAX < mutationRate)
				{
					setGenotypeSupportTypeAt( i, (getGenotypeSupportType()[i] + 1 ) % 2 ); // flip bit
					mutated = true;
				}
			}
		}	
	}
}
