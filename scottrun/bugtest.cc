#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"
//#include "msu_erhosampler/erhosampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(int argc,char *argv[]){
	int ranseed=atoi(argv[1]);
	
	Crandy *randy=new Crandy(-ranseed);
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/parameters.txt");
	int NEVENTS_TOT=parmap.getI("NEVENTS_TOT",10);
	
	Chyper *hyper1=new Chyper();
	
	CresList *reslist=new CresList(&parmap);
	CresInfo::randy=randy;
	
	Chyper *hyper2=new Chyper();
	
	
	printf("check a\n");
	Csampler *sampler=new Csampler(0.15,0.093,&parmap,reslist,randy);
		printf("check b\n");
		int a=2,b=2;
		printf("%d+%d=%d\n",a,b,a+b);
		
		
	Chyper *hyper3=new Chyper();
	
	printf("SUCCESS!!!!!!!!!!!!!!!\n");
	
	/*
	printf("check b\n");
	//Chyper *hyper=new Chyper();
	//CcorrVsEtaScott corrvseta;
	printf("check c\n");
	CcorrVsEtaOlehAlt corrvseta;
	printf("check d\n");
	CcorrVsY corrvsy(&parmap);
	printf("check e\n");
	
	//Csampler *sampler=new Csampler(T,0.093,&parmap,reslist,randy);
	
	CpartList *partlista=new CpartList(&parmap,reslist);
	CpartList *partlistb=new CpartList(&parmap,reslist);
	
	NMSU_ERhoSampler::FillOutHyperBjorken(hyper,T,tau,A,deleta,rhoB,rhoQ);
	hyper->sampler=sampler;
	
	printf("check d\n");
	sampler->CalcDensitiesMu0();
	sampler->GetNHMu0();
	sampler->GetMuNH(hyper);
	sampler->CalcChi4BQS(hyper);
	
	for(ievent=0;ievent<NEVENTS_TOT;ievent++){
		sampler->partlist=partlista;
		sampler->MakeParts(hyper);
		
		partlista->SetEQWeightVec(hyper);
		NMSU_ERhoSampler::DecayParts(randy,partlista);
		
		sampler->partlist=partlistb;
		sampler->MakeParts(hyper);
		
		partlistb->SetEQWeightVec(hyper);
		NMSU_ERhoSampler::DecayParts(randy,partlistb);
		
		corrvsy.Increment(partlista,partlistb,&corrvseta);
		partlista->Clear();
		partlistb->Clear();
		
		if(10*(ievent+1)%NEVENTS_TOT==0){
			CLog::Info("finished "+to_string(((ievent+1)*100)/NEVENTS_TOT)+" percent\n");
		}
	}
	
	printf("howdy a\n");
	corrvsy.WriteResults(to_string(ranseed));
	
	printf("howdy b\n");
	
	//delete partlista;
	
	printf("howdy c\n");
	//delete partlistb;
	printf("howdy d\n");
	*/
	
	return 0;
}
