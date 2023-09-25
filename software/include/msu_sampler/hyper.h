#ifndef __HYPER_H__
#define __HYPER_H__
#include <cstdio>
#include "classdefs.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/randy.h"
#include "msu_eos/eos.h"
#include "msu_commonutils/misc.h"
using namespace std;

//class Csampler;
//class CbulkQuantities;

	// ------------------------
	// Information about hyper-surface element
	//
// ------------------------
class Chyper{
public:
	//Chyper();
	//~Chyper();
	double T0;
	int firstcall;
	double sigma,rhoB,rhoII,rhoS,rhoQ;
	double muB,muII,muS,nhadrons,Rshear,Rbulk,RTbulk,epsilon,P,dedT;
	FourVector dOmega; // hyper volume
	double udotdOmega; // dot product of u and dOmega
	FourVector u;  // collective velocity
	FourVector r; // position
	FourVector qmu; // not used, but read in from Chun's hydro
	double tau,eta; // Bjorken tau and spatial rapidity
	FourTensor pitilde; // shear tensor
	double biggestpitilde; // largest eigenvalue
	double PItilde; // bulk tensor correction
	bool Rvisc_calculated;
	bool epsilon_calculated;
	Eigen::Matrix<double,4,4> chi4,chi4inv,chi4BQS,chi4BQSinv;
			
	void CalcBiggestpitilde();
	double GetEntropyDensity();
	void SetSampler(Csampler *samplerptr);
	//double pixx,pixy,pixz,piyy,piyz,pizz;
	void FillOutShearTensor(double &pixx,double &pixy,double &pixz,double &piyy,double &piyz,double &pizz);
	void Print2D();  // prints out info for 2D element (Bjorken symm)
	void Print();    // prints out info for 3D element
	void Copy(Chyper *oldhyper);
	void Init();
	Csampler *sampler;
	static char *message;
};

#endif
