#include "msu_sampler/master.h"
#include "msu_commonutils/constants.h"

//#define __TEST_PITILDE_TRACE__

using namespace std;
CmeanField *CmasterSampler::meanfield=NULL;

CmasterSampler::CmasterSampler(CparameterMap *parmapin){
	parmap=parmapin;
	randy=new Crandy(-1234);
	reslist=new CresList(parmap);
	partlist=new CpartList(parmap,reslist);
	NEVENTS=0;
	RESWIDTH_ALPHA=parmap->getD("SAMPLER_RESWIDTH_ALPHA",0.5);
	TFmax=0.250;
	DELTF=parmap->getD("SAMPLER_DELTF",0.001);
	NTF=lrint(TFmax/DELTF);
	NSIGMAF=parmap->getD("SAMPLER_NSIGMAF",1);
	if(NSIGMAF<1)
		NSIGMAF=1;
	SIGMAFmin=parmap->getD("SAMPLER_SIGMAFMIN",93.0);
	SIGMAFmax=parmap->getD("SAMPLER_SIGMAFMAX",93.0);
	SETMU0=parmap->getB("SAMPLER_SETMU0",false);
	CALCMU=parmap->getB("SAMPLER_CALCMU",false);
	FINDT=parmap->getB("SAMPLER_FINDT", false);
	NEVENTS_TOT=parmap->getLongI("SAMPLER_NEVENTS_TOT",1);
	NEVENTS=0; // running count of events
	DELSIGMAF=(SIGMAFmax-SIGMAFmin)/double(NSIGMAF);
	if(NSIGMAF==0)
		DELSIGMAF=0.0;
	MEANFIELD=parmap->getS("SAMPLER_MEANFIELD","simple");
	if(MEANFIELD=="simple")
		meanfield=new CmeanField_Simple(parmap);
	else{
		printf("Don't recognize SAMPLER_MEANFIELD from parameter map=%s\n",MEANFIELD.c_str());
		exit(1);
	}
	Csampler::randy=randy;
	CresInfo::randy=randy;
	Csampler::mastersampler=this;
	Csampler::reslist=reslist;
	Csampler::parmap=parmap;
	Csampler::CALCMU=CALCMU;
	Csampler::bose_corr=parmap->getB("SAMPLER_BOSE_CORR",false);
	Csampler::n_bose_corr=parmap->getI("SAMPLER_N_BOSE_CORR",1);
	Csampler::BJORKEN_2D=parmap->getB("SAMPLER_BJORKEN_2D",false);
	Csampler::BJORKEN_YMAX=parmap->getD("SAMPLER_BJORKEN_YMAX",1.0);
	Csampler::USE_POLE_MASS=parmap->getB("SAMPLER_USE_POLE_MASS",false);
	Csampler::INCLUDE_BULK_VISCOSITY=parmap->getB("SAMPLER_INCLUDE_BULK_VISCOSITY",false);
	Csampler::INCLUDE_SHEAR_VISCOSITY=parmap->getB("SAMPLER_INCLUDE_SHEAR_VISCOSITY",false);
	int it,isigma;
	hyperlist.clear();
	sampler.resize(NTF+1);
	for(it=0;it<=NTF;it++){
		sampler[it].resize(NSIGMAF);
		for(isigma=0;isigma<NSIGMAF;isigma++){
			sampler[it][isigma]=nullptr;
		}
	}
}

CmasterSampler::~CmasterSampler(){
	for(auto p : hyperlist)
		delete p;
	hyperlist.clear();

	Csampler *sampleri;
	int iT,isigma;
	for(iT=0;iT<NTF;iT++){
		for(isigma=0;isigma<NSIGMAF;isigma++){
			sampleri=sampler[iT][isigma];
			delete sampleri;
		}
		sampler[iT].clear();
	}
	sampler.clear();
	delete reslist;
	delete randy;
}

int CmasterSampler::MakeEvent(){
	int np,nparts=0;
	Chyper *hyper;
	Csampler *samplerptr=nullptr,*sampler_findT=nullptr;
	//double Omega0Sum=0.0;
	partlist->nparts=0;
	list<Chyper *>::iterator it;

	for(it=hyperlist.begin();it!=hyperlist.end();it++){
		hyper=*it;
		if(hyper->firstcall){
			if(FINDT){
				if(sampler_findT==nullptr)
					sampler_findT=new Csampler(0.140,hyper->sigma);
				sampler_findT->GetTfMuNH(hyper);
				CALCMU=false;
			}
			samplerptr=ChooseSampler(hyper);
			hyper->SetSampler(samplerptr);
			if(samplerptr->FIRSTCALL){
				samplerptr->GetNHMu0();
				samplerptr->CalcDensitiesMu0();
				samplerptr->FIRSTCALL=false;
			}
			if(CALCMU)
				samplerptr->GetMuNH(hyper);
			hyper->firstcall=false;
			samplerptr->CalcNHadronsEpsilonP(hyper);
		}
		np=hyper->sampler->MakeParts(hyper);
		nparts+=np;
		//printf("np=%d,nparts=%d\n",np,nparts);
		//Omega0Sum+=hyper->dOmega[0];
	}
	if(sampler_findT!=nullptr)
		delete sampler_findT;
	NEVENTS+=1;
	return nparts;
}

Csampler* CmasterSampler::ChooseSampler(Chyper *hyper){
	double T,sigma,del;
	int it,isigma;
	T=hyper->T0;
	it=floorl(T/DELTF);
	del=T-it*DELTF;
	if(randy->ran()<del/DELTF)
		it+=1;
	if(it<0)
		it=0;
	if(it>=NTF){
		it=NTF-1;
		printf("WARNING in CmasterSampler::ChooseSampler\n");
		printf("your temperature, T=%g, is higher than TFmax=%g\n",T,TFmax);
	}
	if(NSIGMAF==0){
		isigma=0;
	}
	else{
		sigma=hyper->sigma;
		isigma=floorl((sigma-SIGMAFmin)/DELSIGMAF);
		del=sigma-(SIGMAFmin+isigma*DELSIGMAF);
		if(randy->ran()<del/DELSIGMAF)
			isigma+=1;
		if(isigma<0)
			isigma=0;
		if(isigma>=NSIGMAF)
			isigma=NSIGMAF-1;
	}
	if(sampler[it][isigma]==nullptr){
		printf("making Csampler object, DELTF=%g, it=%d, Tf-T0=%g\n",DELTF,it,it*DELTF-hyper->T0);
		sampler[it][isigma]=new Csampler(it*DELTF,SIGMAFmin+(isigma+0.5)*DELSIGMAF);
	}
	//eprintf("hyper->T0=%g, it=%d, T=%g, sampler->Tf=%g\n",hyper->T0,it,T,sampler[it][isigma]->Tf);
	return sampler[it][isigma];
}

void CmasterSampler::MakeDummyHyper(int nhyper){
	Chyper *hyper;
	for(int i=0;i<nhyper;i++){
		hyper=new Chyper();
		hyperlist.push_back(hyper);
	}
}

void CmasterSampler::ClearHyperList(){
	list<Chyper *>::iterator it;
	for(it=hyperlist.begin();it!=hyperlist.end();it++)
		delete *it;
	hyperlist.clear();
}

void CmasterSampler::GetPitilde(FourTensor &pivisc,FourTensor &pitilde,FourVector &u){
	int alpha,beta;
	double picontract=0.0;
	double pivec[4]={0.0};
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			pitilde[alpha][beta]=0.0;
			picontract+=pivisc[alpha][beta]*u[alpha]*u[beta];
			pivec[alpha]+=pivisc[alpha][beta]*u[beta];
		}
	}
	// Test Tr pivisc
#ifdef __TEST_PITILDE_TRACE__
	double trace=picontract/(u[0]*u[0]);
	for(alpha=1;alpha<4;alpha++){
		trace-=pivisc[alpha][alpha];
	}
	printf("-----------------------------\n");
	printf("before, trace=%g, pi_00=%g=?%g, u=(%g,%g,%g,%g)\n",trace,picontract/(u[0]*u[0]),pivisc[0][0],u[0],u[1],u[2],u[3]);
	printf("pivisc=\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			printf("%10.7f ",pivisc[alpha][beta]);
		}
		printf("\n");
	}
#endif
	//
	double x=u[0]*(1.0+u[0]);
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			pitilde[alpha][beta]=pivisc[alpha][beta]-(u[alpha]*pivec[beta]/x)-(pivec[alpha]*u[beta]/x)+picontract*u[alpha]*u[beta]/(x*x);
		}
	}
#ifdef __TEST_PITILDE_TRACE__
	// Test Tr pitilde
	trace=0.0;
	for(alpha=1;alpha<4;alpha++){
		trace-=pitilde[alpha][alpha];
	}
	printf("after, trace=%g\n",trace);
	if(fabs(picontract)>0.1)
		Misc::Pause();
#endif
}

void CmasterSampler::TransformPiTotz(FourTensor &piMline, const double cosh_eta,const double sinh_eta) {
    double piCart[4][4];
    piCart[0][0] = (piMline[0][0]*cosh_eta*cosh_eta
                    + 2.*piMline[0][3]*cosh_eta*sinh_eta
                    + piMline[3][3]*sinh_eta*sinh_eta);
    piCart[0][1] = piMline[0][1]*cosh_eta + piMline[1][3]*sinh_eta;
    piCart[0][2] = piMline[0][2]*cosh_eta + piMline[2][3]*sinh_eta;
    piCart[0][3] = (piMline[0][0]*cosh_eta*sinh_eta
                    + piMline[0][3]*(cosh_eta*cosh_eta + sinh_eta*sinh_eta)
                    + piMline[3][3]*sinh_eta*cosh_eta);
    piCart[1][0] = piCart[0][1];
    piCart[1][1] = piMline[1][1];
    piCart[1][2] = piMline[1][2];
    piCart[1][3] = piMline[0][1]*sinh_eta + piMline[1][3]*cosh_eta;
    piCart[2][0] = piCart[0][2];
    piCart[2][1] = piCart[1][2];
    piCart[2][2] = piMline[2][2];
    piCart[2][3] = piMline[0][2]*sinh_eta + piMline[2][3]*cosh_eta;
    piCart[3][0] = piCart[0][3];
    piCart[3][1] = piCart[1][3];
    piCart[3][2] = piCart[2][3];
    piCart[3][3] = (piMline[0][0]*sinh_eta*sinh_eta
                    + 2.*piMline[0][3]*cosh_eta*sinh_eta
                    + piMline[3][3]*cosh_eta*cosh_eta);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            piMline[i][j] = piCart[i][j];
}
