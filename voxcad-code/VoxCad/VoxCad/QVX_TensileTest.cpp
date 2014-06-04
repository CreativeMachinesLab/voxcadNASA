#include "QVX_TensileTest.h"
#include "QVX_Interfaces.h"
#include <QFileDialog>
#include <QFile>
#include <QDataStream>
#include <QMessageBox>

QVX_TensileTest::QVX_TensileTest()
{
	NumStep = 0;
	ConvThresh = 0.0;
	OutFilePath = "";
	TestRunning = false;

	IsBasicTensile = false;
	CSArea = 0.0;
	IniLength = 0.0;
	TensileAxis = XAXIS;

	CurTick = 0;
	TotalTick = 1;
	ProgressMessage = "Initializing...";
	CancelFlag = false;
	TensileThread.LinkProgress(&CurTick, &TotalTick, &ProgressMessage, &CancelFlag);

	QObject::connect(&TensileThread, SIGNAL(CallFunc(QString*)), this, SLOT(RunTensileTest(QString*)), Qt::DirectConnection);
}

QVX_TensileTest::~QVX_TensileTest()
{

}

void QVX_TensileTest::BeginTensileTest(QVX_Environment* pEnvIn, int NumStepIn, double ConvThreshIn, double MixRadiusIn, MatBlendModel ModelIn, double PolyExpIn)
{
	//set up local copy of the environment and object...

//	if (MixRadius != 0)	RenderMixedObject(pEnvIn->pObj, &LocObj, MixRadius); //if there's any mixing, we need to make it...
	//else 
	LocObj = *(pEnvIn->pObj); //otherwise just simulate the original object.
	LocEnv = *pEnvIn; //set local environment params to those of the input
	LocEnv.pObj = &LocObj; //make sure local environment points to local object

	//set up blending if desired
	BlendingEnabled = false;
	if (MixRadiusIn != 0){ //if some blending...
		//check for 2 or fewer materials!
		if (LocObj.GetNumLeafMatInUse() > 2){
			QMessageBox::warning(NULL, "Warning", "Currently blending only supported with 2 or fewer materials. Aborting.");
			return;
		}
		BlendingEnabled = true;
		MixRadius = MixRadiusIn;
		BlendModel = ModelIn;
		PolyExp = PolyExpIn;
	}
	
	NumStep = NumStepIn;
	ConvThresh = ConvThreshIn;
	CurTick = 0;
	TotalTick = NumStep;
	CancelFlag = false;

	std::string Message = "";

	pEnvIn->EnableFloor(false);
	pEnvIn->EnableGravity(false);
	pEnvIn->EnableTemp(false);

	Import(&LocEnv, NULL, &Message);
	EnableSelfCollision(false);

	EnableEquilibriumMode(true);

	SetStopConditionType(SC_MIN_MAXMOVE);
	SetStopConditionValue(ConvThreshIn);
//	SetBondDampZ(1.0);
//	SetSlowDampZ(0.03);

	EnablePlasticity(false);
	EnableFailure(false);

	if (!DoBCChecks()) return;

	CurViewMode = RVM_VOXELS;
	CurViewCol = RVC_STRAIN_EN;
	CurViewVox = RVV_DEFORMED;

	//QString DispMesg;

	OutFilePath = QFileDialog::getSaveFileName(NULL, "Save Tensile Test Results", "", "TXT Files (*.txt)");

	emit StartExternalGLUpdate(33);

	TestRunning = true;
	if (OutFilePath != "") TensileThread.Execute(false);
//	RunTensileTest(&DispMesg);
	if (ProgressMessage != "") QMessageBox::warning(NULL, "warning", QString::fromStdString(ProgressMessage));

}


void QVX_TensileTest::RunTensileTest(QString* pDispMessage)
{
	QFile File(OutFilePath);
		
	if (!File.open(QIODevice::WriteOnly | QIODevice::Text)) {
		ProgressMessage = "Could not open file. Aborting.";
		return;
	}
	QTextStream out(&File); 

//	double CurMaxDisp = 9e9;
	std::string IntMsg;

	vfloat StepPercAmt = 1.0/NumStep;
	int VoxCount = NumVox();

	int count = 0;
	int MinimumPerStep = 5;

	int NumBCs = pEnv->GetNumBCs();
	for (int j=0; j<NumBCs; j++){
		if (pEnv->GetBC(j)->Displace.Length2() == 0 ) continue; //if zero displacement, continue
		out << "Disp (m)" << "\t" << "Force (N)" << "\t";
		if (IsBasicTensile){out << "Strain (%)" << "\t" << "Stress (MPa)" << "\t" << "Modulus (MPa)" << "\t";}
	}
	
	out << "\n";

	double LastStress; //so we can calculate modulus at each point...
	double LastStrain;

	for (int i=0; i<NumStep; i++){
		ProgressMessage = "Performing tensile test...";
		for (int j=0; j<VoxCount; j++) VoxArray[j].ScaleExternalInputs((i+1)*StepPercAmt);
		//wiat to settle between timesteps...
		int LastBroken = -1;
		ClearHistories(); //needed so it doesn't immediately "converge" on the next time step

		while (NumBroken() != LastBroken){ //if one breaks, repeat the settling undtil we're done breaking...
			LastBroken = NumBroken();
			EnableFailure(false);
	
//			SetSlowDampZ(CurDamp);

			count = 0;
//			bool LastUnder = false; //were we under the threshhold last time?
//			while (!LastUnder || CurMaxDisp > ConvThresh){
			while (!StopConditionMet()){
	//			if (CurMaxDisp < ConvThresh) LastUnder = true;

				for (int i=0; i<MinimumPerStep; i++){
					if (CancelFlag) break;

					if (!TimeStep(&IntMsg)){ //always do at least 5 steps...
						ProgressMessage = "Tensile test failed. \n \n" + IntMsg;
						CancelFlag = true; //get outta here!
					}
//					CurMaxDisp = SS.NormObjDisp;
				}
				if (CancelFlag) break;
				
				count+=MinimumPerStep;
				if (count > 20000){
					QString Num = QString::number(MaxMoveHistory[0]*1000);
					ProgressMessage = "Simulation not converging at " + Num.toStdString() + ". \n Retry with larger threshold.";
				}
			}
			if (CancelFlag) break;

			EnableFailure(true); //do one step to apply breaking and re-settle as needed...
			if (!TimeStep(&IntMsg)){
				ProgressMessage = "Tensile test failed. \n \n" + IntMsg;
				CancelFlag = true; //get outta here!
			}
		}

		for (int j=0; j<NumBCs; j++){
			CVX_FRegion* pThisBC = pEnv->GetBC(j);
			if (pThisBC->Displace.Length2() != 0 ){ //if non-zero displacement
				double CurDisp = pThisBC->Displace.Length()*(i+1.0)/((double)NumStep);
				double tmp2 = -GetSumForceDir(pThisBC);
				out << CurDisp << "\t" << tmp2 << "\t";
				if (IsBasicTensile){ //only two materials, only one with displacement, so we should only ever enter here once!!
					double ThisStress = tmp2/CSArea;
					double ThisStrain = CurDisp/IniLength;
					out << ThisStrain*100 << "\t" << ThisStress/1e6 << "\t";
					if (i!=0) out << (ThisStress-LastStress)/(ThisStrain-LastStrain)/1e6 << "\t";
					else out << "" << "\t";
					LastStress = ThisStress;
					LastStrain = ThisStrain;
				}
			}
		}
		
		out << "\n";

//		for (int k=0; k<VoxArray.size(); k++){
//			VoxArray[k].ExternalDisp *= (i+2.0)/(i+1.0);
//		}

		CurTick = i+1;
	}
	int stop = 1;

	
	File.close();
	ProgressMessage = ""; //flag to not display message boc on return...
	TestRunning = false;
	emit StopExternalGLUpdate();

}

bool QVX_TensileTest::DoBCChecks(void) //returns true if passes all checks
{
	int NumBCs = pEnv->GetNumBCs();
	std::vector<bool> FullFixed(NumBCs, false);
	int NumFullFixed = 0;
	for (int i=0; i<NumBCs; i++){if (IS_ALL_FIXED(pEnv->GetBC(i)->DofFixed)){ FullFixed[i] = true; NumFullFixed++;}}


	//Are there at least two fixed regions?
	if (NumFullFixed < 2){
		QMessageBox::critical(NULL, "Tensile test error", "At least two fixed boundary conditions are necessary. Aborting.");
		return false;

	}

	//(Calculate maximum presribed deflection)
	double MaxPrescDisp = 0.0;
	for (int i=0; i<NumBCs; i++){
		if (pEnv->GetBC(i)->Displace.Length() > MaxPrescDisp) MaxPrescDisp = pEnv->GetBC(i)->Displace.Length();
	}

	//Does at least one fixed region have a displacement attached to it?
	if (MaxPrescDisp == 0.0){
		QMessageBox::critical(NULL, "Tensile test error", "At least one fixed boundary conditions must have non-zero displacement. Aborting.");
		return false;
	}

	//Is the displacement reasonable? (warn if large displacement...?)
	if (MaxPrescDisp/NumStep > pEnv->pObj->GetLatticeDim()){
		if (QMessageBox::warning(NULL, "Tensile test warning", "Displacement steps are larger than the voxel size. This may cause instability. Do you want to continue?", QMessageBox::Yes | QMessageBox::No)!=QMessageBox::Yes) return false;
	}

	//Are there any forces? (warn that they will be constant throughout simulation)
	//if (pEnv->Forced.size() != 0) QMessageBox::warning(NULL, "Tensile test warning", "Forcing boundary conditions present. Applied forces will be constant throughout tensile test.");

	//checks for a basic, full tensile test:
	
	IsBasicTensile = true; //assume it is until proven false
	int DisplacedBCInd = 0; //which BC is fixed? (can only be 0 or 1)
	//Are all voxel present?
	if (LocObj.GetNumVox() != LocObj.GetStArraySize()) IsBasicTensile = false;
	//are there exactly two BC's, both of them fully fixed?
	if (NumFullFixed != 2 || NumBCs != 2) IsBasicTensile = false;
	//is at least one not moving?
	if (!(pEnv->GetBC(0)->Displace.Length() == 0 || pEnv->GetBC(1)->Displace.Length()==0)) IsBasicTensile = false;
	//is at least one moving?
	if (!(pEnv->GetBC(0)->Displace.Length() != 0 || pEnv->GetBC(1)->Displace.Length()!=0)) IsBasicTensile = false;
	//record which one is fixed
	if (pEnv->GetBC(0)->Displace.Length() == 0) DisplacedBCInd = 1;
	//else DisplacedBCInd already 0...
	//is box shaped regions?
	if (!pEnv->GetBC(0)->IsBox()) IsBasicTensile = false;
	if (!pEnv->GetBC(1)->IsBox()) IsBasicTensile = false;


	//is the BCdisplaced in only one dimension?

	//check if in X Direction...
	if (pEnv->GetBC(DisplacedBCInd)->Displace.x != 0){
		if (pEnv->GetBC(DisplacedBCInd)->Displace.y != 0) IsBasicTensile = false;
		if (pEnv->GetBC(DisplacedBCInd)->Displace.z != 0) IsBasicTensile = false;
		for (int i=0; i<2; i++){
			if (!(pEnv->GetBC(i)->GetRegion()->X == 0.0 || pEnv->GetBC(i)->GetRegion()->X > 1-1.0/LocObj.GetVXDim())) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->Y != 0.0) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->Z != 0.0) IsBasicTensile = false;

			if (pEnv->GetBC(i)->GetRegion()->dX >= 1.0/LocObj.GetVXDim()) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->dY != 1.0) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->dZ != 1.0) IsBasicTensile = false;
		}

		if (IsBasicTensile) {
			TensileAxis = XAXIS;
			Vec3D<> Size = LocObj.GetWorkSpace();
			CSArea = Size.y*Size.z;
			IniLength = Size.x*(((double)LocObj.GetVXDim()-1.0)/LocObj.GetVXDim());
		}
	}

	//check if in Y Direction...
	if (pEnv->GetBC(DisplacedBCInd)->Displace.y != 0){
		if (pEnv->GetBC(DisplacedBCInd)->Displace.x != 0) IsBasicTensile = false;
		if (pEnv->GetBC(DisplacedBCInd)->Displace.z != 0) IsBasicTensile = false;
		for (int i=0; i<2; i++){
			if (!(pEnv->GetBC(i)->GetRegion()->Y == 0.0 || pEnv->GetBC(i)->GetRegion()->Y > 1-1.0/LocObj.GetVYDim())) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->X != 0.0) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->Z != 0.0) IsBasicTensile = false;

			if (pEnv->GetBC(i)->GetRegion()->dY >= 1.0/LocObj.GetVYDim()) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->dX != 1.0) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->dZ != 1.0) IsBasicTensile = false;
		}

		if (IsBasicTensile) {
			TensileAxis = YAXIS;
			Vec3D<> Size = LocObj.GetWorkSpace();
			CSArea = Size.x*Size.z;
			IniLength = Size.y*((LocObj.GetVYDim()-1)/LocObj.GetVYDim());
		}
	}

	//check if in Z Direction...
	if (pEnv->GetBC(DisplacedBCInd)->Displace.z != 0){
		if (pEnv->GetBC(DisplacedBCInd)->Displace.y != 0) IsBasicTensile = false;
		if (pEnv->GetBC(DisplacedBCInd)->Displace.z != 0) IsBasicTensile = false;
		for (int i=0; i<2; i++){
			if (!(pEnv->GetBC(i)->GetRegion()->Z == 0.0 || pEnv->GetBC(i)->GetRegion()->Z > 1-1.0/LocObj.GetVZDim())) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->Y != 0.0) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->X != 0.0) IsBasicTensile = false;

			if (pEnv->GetBC(i)->GetRegion()->dZ >= 1.0/LocObj.GetVZDim()) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->dY != 1.0) IsBasicTensile = false;
			if (pEnv->GetBC(i)->GetRegion()->dX != 1.0) IsBasicTensile = false;
		}

		if (IsBasicTensile) {
			TensileAxis = ZAXIS;
			Vec3D<> Size = LocObj.GetWorkSpace();
			CSArea = Size.x*Size.y;
			IniLength = Size.z*((LocObj.GetVZDim()-1)/LocObj.GetVZDim());
		}
	}


	return true;
}

void QVX_TensileTest::RenderMixedObject(CVX_Object* pSrcObj, CVX_Object* pDestObj, double MixRadius)
{
//	MixRadius=MixRadiusIn;
	//set up palette
	int NumContMat = 100;
	double AvgPoisson = 0;
	int NumAvgd = 0;
	double MinEMod = 9e99, MaxEMod = 0;

	//get minimum and maximum stiffnesses of materials in simulation (use average poissons ratio)
	for (int i=1; i<pSrcObj->Palette.size(); i++){
		//store min/max elastic modulus
		double ThisEMod = pSrcObj->Palette[i].GetElasticMod();
		if (ThisEMod > MaxEMod) MaxEMod = ThisEMod;
		if (ThisEMod < MinEMod) MinEMod = ThisEMod;

		//average poissons
		AvgPoisson += pSrcObj->Palette[i].GetPoissonsRatio();
		NumAvgd++;
	}
	AvgPoisson /= NumAvgd;

	*pDestObj = *pSrcObj; //we want everything except the palette and the structure pattern...
	pDestObj->ClearPalette();

	//add back in a gradient palette:
	for (int i=0; i<NumContMat; i++){
		QString Name = "GM" + QString::number(i);
		pDestObj->AddMat(Name.toStdString(), MinEMod + (MaxEMod-MinEMod)*(((double)i)/(NumContMat-1)), AvgPoisson); //NumContMat, inclusive of endpoints
	}

	//set the material index in the new object
	int x, y, z, ThisInd;
	Vec3D<> EnvMult, BasePos, ThisPos;
	double ThisDist, ThisWeight;
	for (int i=0; i<pSrcObj->GetStArraySize(); i++){
		if (pSrcObj->GetMat(i) == 0) continue; //skip this one if it is empty...

		pSrcObj->GetXYZNom(&x, &y, &z, i);
		BasePos = pSrcObj->GetXYZ(i);
		EnvMult = pSrcObj->Lattice.GetDimAdj();
		int xLook = 3*MixRadius/EnvMult.x+1;
		int yLook = 3*MixRadius/EnvMult.y+1;
		int zLook = 3*MixRadius/EnvMult.z+1;

		//MixRadius is the standard deviation of the gaussian mixing weighting . (multiply by GetLatticeDim to get the real-valued standard dev)
		double AvgEMod = 0.0;
		double TotalWeight = 0.0;

		//up to three sigma...
		for (int ix = x-xLook; ix<=x+xLook; ix++){
			for (int jy = y-yLook; jy<=y+yLook; jy++){
				for (int kz = z-zLook; kz<=z+zLook; kz++){
					ThisInd = pSrcObj->GetIndex(ix, jy, kz);
					if (ThisInd == -1) continue; //if invalid location then skip
					if (pSrcObj->GetMat(ThisInd) == 0) continue; //skip this one if it is empty (no material...
					
					ThisPos = pSrcObj->GetXYZ(ThisInd);
					ThisDist = ((ThisPos-BasePos).Length())/pSrcObj->GetLatticeDim();

					ThisWeight = exp(-0.5*ThisDist*ThisDist/(MixRadius*MixRadius));
					double ThisLocEMod = pSrcObj->GetLeafMat(ThisInd)->GetElasticMod();

					AvgEMod += ThisWeight*ThisLocEMod;
					TotalWeight += ThisWeight;
				}
			}
		}

		AvgEMod /= TotalWeight;

		//convert this stiffness to a material index..
		int MatIndex = 1 + (int)((AvgEMod-MinEMod)/(MaxEMod-MinEMod)*(NumContMat-1)+0.5);
		pDestObj->SetMat(i, MatIndex);
	}
}
