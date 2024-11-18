#ifndef __MASTER_SAMPLER_H__
#define __MASTER_SAMPLER_H__
#include <list>
#include "msu_sampler/classdefs.h"
#include "msu_sampler/meanfield.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/misc.h"
#include "msu_eos/resonances.h"
#include "msu_sampler/hyper.h"
#include "msu_sampler/part.h"
#include "msu_sampler/sampler.h"
#include "msu_eos/eos.h"
using namespace std;
namespace NMSUPratt{

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
		void DeleteHyperElements();
		char *message;
		CparameterMap *parmap;
		CresList *reslist;
		Crandy *randy;
		double TFmax,TFmin;
		string MEANFIELD;
		double SIGMAFmin,SIGMAFmax;
		int NTF,NSIGMAF;
		double DELTF,DELSIGMAF;
		double FUGACITY_TAU0,FUGACITY_GAMMA_0,FUGACITY_TAU_EQ;
		bool SETMU0,CALCMU,VARY_FUGACITY;
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
		void ReadHyper_Duke_2D(int run_number,string qual);
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

}

#endif
