#include "msu_sampler/master.h"
#include "msu_commonutils/constants.h"

//#define __TEST_PITILDE_TRACE__

using namespace std;
CmeanField *CmasterSampler::meanfield=NULL;

CmasterSampler::CmasterSampler(CparameterMap *parmapin){
	message=new char[CLog::CHARLENGTH];
	parmap=parmapin;
	randy=new Crandy(-1234);
	reslist=new CresList(parmap);
	partlist=new CpartList(parmap,reslist);
	NEVENTS=0;
	RESWIDTH_ALPHA=parmap->getD("MSU_SAMPLER_RESWIDTH_ALPHA",0.5);
	TFmax=0.250;
	DELTF=parmap->getD("MSU_SAMPLER_DELTF",0.001);
	NTF=lrint(TFmax/DELTF);
	NSIGMAF=parmap->getD("MSU_SAMPLER_NSIGMAF",1);
	if(NSIGMAF<1)
		NSIGMAF=1;
	SIGMAFmin=parmap->getD("MSU_SAMPLER_SIGMAFMIN",93.0);
	SIGMAFmax=parmap->getD("MSU_SAMPLER_SIGMAFMAX",93.0);
	SETMU0=parmap->getB("MSU_SAMPLER_SETMU0",false);
	CALCMU=parmap->getB("MSU_SAMPLER_CALCMU",false);
	FINDT=parmap->getB("MSU_SAMPLER_FINDT", false);
	NEVENTS_TOT=parmap->getLongI("MSU_SAMPLER_NEVENTS_TOT",1);
	NEVENTS=0; // running count of events
	DELSIGMAF=(SIGMAFmax-SIGMAFmin)/double(NSIGMAF);
	if(NSIGMAF==0)
		DELSIGMAF=0.0;
	MEANFIELD=parmap->getS("MSU_SAMPLER_MEANFIELD","simple");
	if(MEANFIELD=="simple")
		meanfield=new CmeanField_Simple(parmap);
	else{
		snprintf(message,CLog::CHARLENGTH,"Don't recognize MSU_SAMPLER_MEANFIELD from parameter map=%s\n",MEANFIELD.c_str());
		CLog::Info(message);
	}
	Csampler::randy=randy;
	CresInfo::randy=randy;
	Csampler::mastersampler=this;
	Csampler::reslist=reslist;
	Csampler::parmap=parmap;
	Csampler::CALCMU=CALCMU;
	Csampler::bose_corr=parmap->getB("MSU_SAMPLER_BOSE_CORR",false);
	Csampler::n_bose_corr=parmap->getI("MSU_SAMPLER_N_BOSE_CORR",1);
	Csampler::BJORKEN_2D=parmap->getB("MSU_SAMPLER_BJORKEN_2D",false);
	Csampler::BJORKEN_YMAX=parmap->getD("MSU_SAMPLER_BJORKEN_YMAX",1.0);
	Csampler::USE_POLE_MASS=parmap->getB("MSU_SAMPLER_USE_POLE_MASS",false);
	Csampler::INCLUDE_BULK_VISCOSITY=parmap->getB("MSU_SAMPLER_INCLUDE_BULK_VISCOSITY",false);
	Csampler::INCLUDE_SHEAR_VISCOSITY=parmap->getB("MSU_SAMPLER_INCLUDE_SHEAR_VISCOSITY",false);
	Csampler::NSAMPLE=parmap->getI("MSU_SAMPLER_NSAMPLE",1);
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
	Csampler *samplerptr=nullptr,*MSU_SAMPLER_findT=nullptr;
	//double Omega0Sum=0.0;
	partlist->nparts=0;
	list<Chyper *>::iterator it;


	for(it=hyperlist.begin();it!=hyperlist.end();it++){
		hyper=*it;
		if(hyper->firstcall){
			if(FINDT){
				if(MSU_SAMPLER_findT==nullptr)
					MSU_SAMPLER_findT=new Csampler(0.140,hyper->sigma);
				MSU_SAMPLER_findT->GetTfMuNH(hyper);
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
	}
	if(MSU_SAMPLER_findT!=nullptr)
		delete MSU_SAMPLER_findT;
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
		snprintf(message,CLog::CHARLENGTH,"WARNING in CmasterSampler::ChooseSampler\n");
		CLog::Info(message);
		snprintf(message,CLog::CHARLENGTH,"your temperature, T=%g, is higher than TFmax=%g\n",T,TFmax);
		CLog::Info(message);
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
		snprintf(message,CLog::CHARLENGTH,"making Csampler object, DELTF=%g, it=%d, Tf-T0=%g\n",DELTF,it,it*DELTF-hyper->T0);
		CLog::Info(message);
		sampler[it][isigma]=new Csampler(it*DELTF,SIGMAFmin+(isigma+0.5)*DELSIGMAF);
	}
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
	snprintf(message,CLog::CHARLENGTH,"-----------------------------\n");
	CLog::Info(message);
	snprintf(message,"before, trace=%g, pi_00=%g=?%g, u=(%g,%g,%g,%g)\n",trace,picontract/(u[0]*u[0]),pivisc[0][0],u[0],u[1],u[2],u[3]);
	CLog::Info(message);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"pivisc=\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			snprintf(message,CLog::CHARLENGTH,"%10.7f ",pivisc[alpha][beta]);
			CLog::Info(message);
		}
		snprintf(message,CLog::CHARLENGTH,"\n");

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
	snprintf(message,CLog::CHARLENGTH,"after, trace=%g\n",trace);
	CLog::Info(message);
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
