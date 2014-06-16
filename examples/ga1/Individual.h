/*===================================================================
Nick Cheney
2014-06-13
Naive GA Implemelntation: Template for GA design with voxelyze.
(direct encoding hillclimber)

Individual.h
===================================================================*/
#include <vector>

using namespace std;

class Individual
{
	protected:
		float fitness;
		int num_x_voxels;
		int num_y_voxels;
		int num_z_voxels;
		int genotypeLength;
		vector<int> genotypeShape;
		vector<int> genotypeMuscleOrSupport;
		vector<int> genotypeMuscleType;
		vector<int> genotypeSupportType;

	public:
		Individual(int _num_x_voxels, int _num_y_voxels, int _num_z_voxels);

		void setFitness(float _fitness);
		float getFitness();

		// note: this genome is more efficiently implemented as a matrix, rather than 4 vectors.  But seeing how this is a tutorial, the naming is more straightforward with this implementation
		void setGenotypeShapeAt(int _index, int _value);
		void setGenotypeMuscleOrSupportAt(int _index, int _value);
		void setGenotypeMuscleTypeAt(int _index, int _value);
		void setGenotypeSupportTypeAt(int _index, int _value);

		vector<int> getGenotypeShape();
		vector<int> getGenotypeMuscleOrSupport();
		vector<int> getGenotypeMuscleType();
		vector<int> getGenotypeSupportType();

		void writeVxaFile();

		void collectFitness();

		void mutate(int entriesMutated);
};