#ifndef __PART_H__
#define __PART_H__
#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/misc.h"
#include "msu_sampler/resonances.h"
#include "msu_sampler/classdefs.h"
#include "msu_sampler/eos.h"
using namespace std;

// ----------------------------
// Each particle has following information in Cpart
// CpartList is a vector of such particles, 'nparts' is the number of particles
// Rather than adding memory one particle at a time, particles are added in increments of 'nparts_blocksize'
// Increasing particles done through AddPart -- which checks to see if vector size needs to be increased
// Only particles 0 through nparts-1 are to be considered
// ----------------------------

class Cpart{
public:
	Cpart();
	~Cpart();
	int pid;
	double msquared;
	FourVector p,r;
	void Print();
	double GetMass();
	void AddPart(int pid,FourVector &p,FourVector &r);
	void Setp0();
	void SetMsquared();
	void Boost(FourVector &u);
	void BoostP(FourVector &u);
	void BoostR(FourVector &u);
};


class CpartList{
public:
	double SE[4][4];
	CpartList(CparameterMap *parmap,CresList *reslist);
	~CpartList();
		int nparts,nparts_blocksize;  // increases array by nparts_blocksize when needed
		vector<Cpart> partvec;
		Cpart *GetPart();
		void Clear(); // set nparts=0 and frees memory
		void Reset();     // just sets nparts=0
		void WriteParts(string filename);
		void CountResonances();
		long long int CountResonances(int pid);
		void IncrementSpectra(int pid,double dp,vector<double> &spectra);
		void IncrementMassDist(int pid,double dm,vector<double> &massdist);
		double SumEnergy();
		double SumEnergy(int pid);
		void SumSETensor();
		void AddPart(int pidset,FourVector &pset,FourVector &rset);
		static CresList *reslist;
	};
#endif
