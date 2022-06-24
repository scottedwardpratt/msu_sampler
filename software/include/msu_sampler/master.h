#ifndef __MASTER_SAMPLER_H__
#define __MASTER_SAMPLER_H__
#include <list>
#include "classdefs.h"
#include "meanfield.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/misc.h"
#include "msu_eos/resonances.h"
#include "msu_sampler/hyper.h"
#include "msu_sampler/sampler.h"
#include "msu_eos/eos.h"
using namespace std;

// --------------------------------
// This object contains array of sampler objects, and information used across all samplers, e.g. resonance info lists
// sampler objects have unique temperature and sigma field
// It also has list of hyper-volume elements
//
// Typical Usage:
// CmasterSampler ms(parametermap);
// ms.partlist=new CpartList(&parmap);
// ms.ReadHyper();
// nparts+=ms.MakeEvent();
//
// ---------------------------------

class CmasterSampler{
public:
	CmasterSampler(CparameterMap *parmapin);
	CmasterSampler(string parsfilename); // Constructor
	~CmasterSampler();
	char *message;
	CparameterMap *parmap;
	CresList *reslist;
	Crandy *randy;
	double TFmax;
	string MEANFIELD;
	double SIGMAFmin,SIGMAFmax;
	int NTF,NSIGMAF;
	double DELTF,DELSIGMAF;
	bool SETMU0,CALCMU;
	double RESWIDTH_ALPHA;
	double YMAX; // only for 2D sampler
	long long int NEVENTS,NEVENTS_TOT;
	int nelements;
	bool FINDT;  // true if you need to find T(epsilon) in hyper elements

	list<Chyper *> hyperlist;
	vector<vector<Csampler *>> sampler;  // array of samplers indexed by T and sigma
	CpartList *partlist;

	void ReadHyper_BEST_Binary3D();
	void ReadHyper_OSU_2D();
	int MakeEvent(); // returns number of parts
	Csampler* ChooseSampler(Chyper *hyper);
	void ChooseSampler(double Tf,double sigmaf,int &itf,int &isigmaf);
	void MakeDummyHyper(int nhyper);
	void GetPitilde(FourTensor &pivisc,FourTensor &pitilde,FourVector &u);
	static CmeanField *meanfield;
	void TransformPiTotz(FourTensor &piMline, const double cosh_eta,
		const double sinh_eta);
	void ClearHyperList();
};

#endif
