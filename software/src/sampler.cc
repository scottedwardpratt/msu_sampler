#include "msu_sampler/sampler.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/constants.h"

using namespace std;

Crandy* Csampler::randy=NULL;
CresList *Csampler::reslist=NULL;
CparameterMap *Csampler::parmap=NULL;
CmasterSampler *Csampler::mastersampler=NULL;
bool Csampler::CALCMU=true;
bool Csampler::SETMU0=false;
int Csampler::n_bose_corr=1;
bool Csampler::bose_corr=false;
bool Csampler::BJORKEN_2D=false;
double Csampler::BJORKEN_YMAX=1.0;
bool Csampler::USE_POLE_MASS=false;
bool Csampler::INCLUDE_BULK_VISCOSITY=false;
bool Csampler::INCLUDE_SHEAR_VISCOSITY=false;
int Csampler::NSAMPLE=1;
char *Csampler::message=new char[300];

// Constructor
Csampler::Csampler(double Tfset,double sigmafset){
	FIRSTCALL=true;
	SETMU0=false;
	Tf=Tfset;
	sigmaf=sigmafset;
	if(!bose_corr)
		n_bose_corr=1;
	int nres=reslist->resmap.size();
	if(reslist->GetResInfoPtr(22)->pid==22)
		nres-=1;
	density0i.resize(nres);
	epsilon0i.resize(nres);
	P0i.resize(nres);
	sfdens0imap.resize(nres);
	if(bose_corr){
		pibose_P0.resize(n_bose_corr+1);
		pibose_epsilon0.resize(n_bose_corr+1);
		pibose_dedt0.resize(n_bose_corr+1);
		pibose_dens0.resize(n_bose_corr+1);
	}
	forMU0_calculated=false;
}

// Destructor
Csampler::~Csampler(){
	density0i.clear();
	epsilon0i.clear();
	P0i.clear();
	sfdens0imap.clear();
	pibose_dens0.clear();
	pibose_P0.clear();
	pibose_epsilon0.clear();
	pibose_dedt0.clear();
}

// Calculates array of densities for each resonance with mu=0, to be used in sampling
void Csampler::CalcDensitiesMu0(){
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	double Pi,epsiloni,densi,dedti,p4overE3i,Ji;
	int ires,a,b;
	nhadrons0=P0=epsilon0=0.0;
	chi0.setZero();
	sigma0.setZero();

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			// note this does not change mass due to sigma!=93 MeV
			MSU_EOS::GetEpsilonPDens_OneSpecies(Tf,resinfo,epsiloni,Pi,densi,dedti,p4overE3i,Ji,USE_POLE_MASS);
			ires=resinfo->ires;
			density0i[ires]=densi;
			epsilon0i[ires]=epsiloni;
			P0i[ires]=Pi;
			nhadrons0+=density0i[ires];
			epsilon0+=epsiloni;
			P0+=Pi;
			p4overE30+=p4overE3i;
			if(resinfo->decay && !USE_POLE_MASS)
				CalcSFDensMap(resinfo,Tf,sfdens0imap[ires]);
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					chi0(a,b)+=densi*resinfo->q[a]*resinfo->q[b];
					sigma0(a,b)+=Ji*resinfo->q[a]*resinfo->q[b];
				}
			}
		}
	}
	if(bose_corr){
		for(int nbose=2;nbose<=n_bose_corr;nbose++){
			MSU_EOS::freegascalc_onespecies(Tf/nbose,0.138,pibose_epsilon0[nbose],pibose_P0[nbose],pibose_dens0[nbose],
			pibose_dedt0[nbose]);
			p4overE3i=MSU_EOS::Getp4overE3(Tf/nbose,0.138,pibose_dens0[nbose]);
			nhadrons0+=3.0*pibose_dens0[nbose];
			epsilon0+=3.0*pibose_epsilon0[nbose];
			P0+=3.0*pibose_P0[nbose];
			p4overE30+=3.0*p4overE3i;
			Ji=3*MSU_EOS::GetJi(Tf/double(nbose),0.138,pibose_dens0[nbose]);
		}
	}
}

// Calculates factors (depend only on T) used for Newton's method to get muB, muII, muS from rhoB, rhoII, rhoS
void Csampler::GetNHMu0(){
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	double Pi,epsiloni,densi,dedti,p4overE3i,Ji,m2;
	int B,S,II,Q;

	nh0_b0i0s0=nh0_b0i2s0=nh0_b0i1s1=0.0;
	nh0_b1i0s1=nh0_b1i0s3=nh0_b1i1s0=nh0_b1i1s2=nh0_b1i2s1=nh0_b1i3s0=0.0;
	nh0_b2i0s0=0.0;

	eh0_b0i0s0=eh0_b0i2s0=eh0_b0i1s1=0.0;
	eh0_b1i0s1=eh0_b1i0s3=eh0_b1i1s0=eh0_b1i1s2=eh0_b1i2s1=eh0_b1i3s0=0.0;
	eh0_b2i0s0=0.0;

	dedth0_b0i0s0=dedth0_b0i2s0=dedth0_b0i1s1=0.0;
	dedth0_b1i0s1=dedth0_b1i0s3=dedth0_b1i1s0=dedth0_b1i1s2=dedth0_b1i2s1=dedth0_b1i3s0=0.0;
	dedth0_b2i0s0=0.0;
	
	p4overE3h0_b0i0s0=p4overE3h0_b0i2s0=p4overE3h0_b0i1s1=0.0;
	p4overE3h0_b1i0s1=p4overE3h0_b1i0s3=p4overE3h0_b1i1s0=p4overE3h0_b1i1s2=p4overE3h0_b1i2s1=p4overE3h0_b1i3s0=0.0;
	p4overE3h0_b2i0s0=0.0;
	
	eEbarh0_b0i0s0=eEbarh0_b0i2s0=eEbarh0_b0i1s1=0.0;
	eEbarh0_b1i0s1=eEbarh0_b1i0s3=eEbarh0_b1i1s0=eEbarh0_b1i1s2=eEbarh0_b1i2s1=eEbarh0_b1i3s0=0.0;
	eEbarh0_b2i0s0=0.0;
	
	m2densh0_b0i0s0=m2densh0_b0i2s0=m2densh0_b0i1s1=0.0;
	m2densh0_b1i0s1=m2densh0_b1i0s3=m2densh0_b1i1s0=m2densh0_b1i1s2=m2densh0_b1i2s1=m2densh0_b1i3s0=0.0;
	m2densh0_b2i0s0=0.0;

	if(bose_corr){
		densi=p4overE3i=epsiloni=dedti=0.0;
		for(int nbose=2;nbose<=n_bose_corr;nbose++){
			MSU_EOS::freegascalc_onespecies(Tf/nbose,0.138,pibose_epsilon0[nbose],pibose_P0[nbose],pibose_dens0[nbose],
			pibose_dedt0[nbose]);
			epsiloni+=pibose_epsilon0[nbose];
			densi+=pibose_dens0[nbose];
			dedti+=pibose_dedt0[nbose];
			p4overE3i+=MSU_EOS::Getp4overE3(Tf/nbose,0.138,pibose_dens0[nbose]);
		}
		m2=0.138*0.138;
		nh0_b0i0s0+=densi;
		eh0_b0i0s0+=epsiloni;
		dedth0_b0i0s0+=dedti;
		p4overE3h0_b0i0s0+=p4overE3i;
		eEbarh0_b0i0s0+=epsiloni*epsiloni/densi;
		m2densh0_b0i0s0+=m2*densi;
		nh0_b0i2s0+=densi;
		eh0_b0i2s0+=epsiloni;
		dedth0_b0i2s0+=dedti;
		p4overE3h0_b0i2s0+=p4overE3i;
		eEbarh0_b0i2s0+=epsiloni*epsiloni/densi;
		m2densh0_b0i2s0+=m2*densi;
	}

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			MSU_EOS::GetEpsilonPDens_OneSpecies(Tf,resinfo,epsiloni,Pi,densi,dedti,p4overE3i,Ji,USE_POLE_MASS);
			m2=resinfo->mass*resinfo->mass;

			B=resinfo->baryon;
			S=resinfo->strange;
			Q=resinfo->charge;
			II=2*Q-B-S; // II is 2*I3
			B=abs(B);
			S=abs(S);
			II=abs(II);

			if(B==0 && II==0 && S==0){
				nh0_b0i0s0+=densi;
				eh0_b0i0s0+=epsiloni;
				dedth0_b0i0s0+=dedti;
				p4overE3h0_b0i0s0+=p4overE3i;
				eEbarh0_b0i0s0+=epsiloni*epsiloni/densi;
				m2densh0_b0i0s0+=m2*densi;
			}
			else if(B==0 && II==1 && S==1){
				nh0_b0i1s1+=densi;
				eh0_b0i1s1+=epsiloni;
				dedth0_b0i1s1+=dedti;
				p4overE3h0_b0i1s1+=p4overE3i;
				eEbarh0_b0i1s1+=epsiloni*epsiloni/densi;
				m2densh0_b0i1s1+=m2*densi;
			}
			else if(B==0 && II==2 && S==0){
				nh0_b0i2s0+=densi;
				eh0_b0i2s0+=epsiloni;
				dedth0_b0i2s0+=dedti;
				p4overE3h0_b0i2s0+=p4overE3i;
				eEbarh0_b0i2s0+=epsiloni*epsiloni/densi;
				m2densh0_b0i2s0+=m2*densi;
			}
			else if(B==1 && II==0 && S==1){
				nh0_b1i0s1+=densi;
				eh0_b1i0s1+=epsiloni;
				dedth0_b1i0s1+=dedti;
				p4overE3h0_b1i0s1+=p4overE3i;
				eEbarh0_b1i0s1+=epsiloni*epsiloni/densi;
				m2densh0_b1i0s1+=m2*densi;
			}
			else if(B==1 && II==0 && S==3){
				nh0_b1i0s3+=densi;
				eh0_b1i0s3+=epsiloni;
				dedth0_b1i0s3+=dedti;
				p4overE3h0_b1i0s3+=p4overE3i;
				eEbarh0_b1i0s3+=epsiloni*epsiloni/densi;
				m2densh0_b1i0s3+=m2*densi;
			}
			else if(B==1 && II==1 && S==0){
				nh0_b1i1s0+=densi;
				eh0_b1i1s0+=epsiloni;
				dedth0_b1i1s0+=dedti;
				p4overE3h0_b1i1s0+=p4overE3i;
				eEbarh0_b1i1s0+=epsiloni*epsiloni/densi;
				m2densh0_b1i1s0+=m2*densi;
			}
			else if(B==1 && II==1 && S==2){
				nh0_b1i1s2+=densi;
				eh0_b1i1s2+=epsiloni;
				dedth0_b1i1s2+=dedti;
				p4overE3h0_b1i1s2+=p4overE3i;
				eEbarh0_b1i1s2+=epsiloni*epsiloni/densi;
				m2densh0_b1i1s2+=m2*densi;
			}
			else if(B==1 && II==2 && S==1){
				nh0_b1i2s1+=densi;
				eh0_b1i2s1+=epsiloni;
				dedth0_b1i2s1+=dedti;
				p4overE3h0_b1i2s1+=p4overE3i;
				eEbarh0_b1i2s1+=epsiloni*epsiloni/densi;
				m2densh0_b1i2s1+=m2*densi;
			}
			else if(B==1 && II==3 && S==0){
				nh0_b1i3s0+=densi;
				eh0_b1i3s0+=epsiloni;
				dedth0_b1i3s0+=dedti;
				p4overE3h0_b1i3s0+=p4overE3i;
				eEbarh0_b1i3s0+=epsiloni*epsiloni/densi;
				m2densh0_b1i3s0+=m2*densi;
			}
			else if(B==2 && II==0 && S==0){ //deuteron, not currently used in calculations
				nh0_b2i0s0+=densi;
				eh0_b2i0s0+=epsiloni;
				dedth0_b2i0s0+=dedti;
				p4overE3h0_b2i0s0+=p4overE3i;
				eEbarh0_b2i0s0+=epsiloni*epsiloni/densi;
				m2densh0_b2i0s0+=m2*densi;
			}
			else{
				snprintf(message,CLog::CHARLENGTH,"NO BIS Match!!!! B=%d, II=%d, S=%d\n",B,II,S);
				CLog::Fatal(message);
			}
		}
	}
	forMU0_calculated=true;
}

void Csampler::CalcNHadrons(Chyper *hyper){
	int ires,II3;
	CresInfo *resinfo;
	double mutot=0;
	CresMassMap::iterator iter;
	
	double nhadrons=0.0;
	for(iter=reslist->massmap.begin();iter!=reslist->massmap.end();++iter){
		resinfo=iter->second;
		if(resinfo->pid!=22){
			ires=resinfo->ires;
			II3=2.0*resinfo->charge-resinfo->baryon-resinfo->strange;
			mutot=hyper->muB*resinfo->baryon+hyper->muII*II3+hyper->muS*resinfo->strange;
			nhadrons+=density0i[ires]*exp(mutot);
		}
	}
	hyper->nhadrons=nhadrons;
}

/*
double Csampler::CalcLambdaF(double muB,double muII,double muS,double Pf){
	int ires=0,nbose;
	double Ipp=0.0,xx,lambdaf;
	double lambdafact,mutot,I3;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	if(mastersampler->SETMU0){
		return lambda0;
	}
	else{
		for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
			resinfo=rpos->second;
			if(resinfo->pid!=22){
				I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
				mutot=muB*resinfo->baryon+muII*I3+muS*resinfo->strange;
				Ipp+=lambda0i[ires]*exp(mutot);
				ires+=1;
			}
		}
		if(bose_corr){
			for(nbose=2;nbose<=n_bose_corr;nbose++){
				xx=exp(muII);
				Ipp+=pibose_lambda0[nbose]*(pow(xx,nbose)+pow(xx,-nbose)+1.0);
			}
		}
		if(mastersampler->SETMU0)
			lambdafact=2.0*P0-4.0*Ipp;
		else
			lambdafact=2.0*Pf-4.0*Ipp;
		lambdaf=lambdafact;
		return lambdaf;
	}
}
*/

// Calculates muB, muII, muS from rhoB, rhoII, rhoS, also calculates nhadronsf (factors above must be calculated first)
void Csampler::GetMuNH(Chyper *hyper){
	GetMuNH(hyper->rhoB,hyper->rhoII,hyper->rhoS,hyper->muB,hyper->muII,hyper->muS,hyper->nhadrons);
}

void Csampler::GetMuNH(double rhoBtarget,double rhoIItarget,double rhoStarget,double &muB,double &muII,double &muS,double &nhadrons){
	Eigen::MatrixXd A(3,3);
	Eigen::VectorXd mu(3),dmu(3),drho(3);

	//3D Newton's Method
	// Here rhoII refers to rho_u-rho_d = 2*I3 and mu[1]=muII/2
	double rhoB,rhoS,rhoII=0.0;

	double drhoB_dmuB,drhoB_dmuS,drhoB_dmuII;
	double drhoII_dmuS,drhoII_dmuII;
	double drhoS_dmuS;

	double xB,xI,xS,xxB,xxI,xxS,dmumag;

	mu[0]=muB;
	mu[1]=muII;
	mu[2]=muS;
	do{
		xB=exp(mu[0]);
		xI=exp(mu[1]);
		xS=exp(mu[2]);
		xxB=1.0/xB;
		xxI=1.0/xI;
		xxS=1.0/xS;
		
		nhadrons=nh0_b0i0s0+0.5*nh0_b0i2s0*(xI*xI+xxI*xxI)
			+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
					+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
						+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
							+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
								+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
									+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
										+0.5*nh0_b2i0s0*(xB*xB+xxB*xxB);
		
		rhoB=0.5*nh0_b1i0s1*(xB*xxS-xxB*xS)
			+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS-xxB*xS*xS*xS)
				+0.25*nh0_b1i1s0*(xB-xxB)*(xI+xxI)
					+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(xI*xI+xxI*xxI)
							+0.25*nh0_b1i3s0*(xB-xxB)*(xI*xI*xI+xxI*xxI*xxI);
		drhoB_dmuB=0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
			+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
				+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
					+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
							+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);
		drhoB_dmuII=0.25*nh0_b1i1s0*(xB-xxB)*(xI-xxI)
			+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI-xxI)
				+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(2*xI*xI-2*xxI*xxI)
					+0.25*nh0_b1i3s0*(xB-xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);
		drhoB_dmuS=0.5*nh0_b1i0s1*(-xB*xxS-xxB*xS)
			+0.5*nh0_b1i0s3*(-3*xB*xxS*xxS*xxS-3*xxB*xS*xS*xS)
				+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS-2*xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(-xB*xxS-xxB*xS)*(xI*xI+xxI*xxI);

		rhoII=0.5*nh0_b0i2s0*(2*xI*xI-2*xxI*xxI)
			+0.25*nh0_b0i1s1*(xI-xxI)*(xS+xxS)
				+0.25*nh0_b1i1s0*(xB+xxB)*(xI-xxI)
					+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI-xxI)
						+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI)
							+0.25*nh0_b1i3s0*(xB+xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);
		
		drhoII_dmuII=0.5*nh0_b0i2s0*(4*xI*xI+4*xxI*xxI)
			+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
					+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(4*xI*xI+4*xxI*xxI)
							+0.25*nh0_b1i3s0*(xB+xxB)*(9*xI*xI*xI+9*xxI*xxI*xxI);
		if (bose_corr) {
			for(int nbose=2;nbose<=n_bose_corr;nbose++){
				rhoII+=pibose_dens0[nbose]*2.0*sinh(2.0*nbose*mu[1]);
				drhoII_dmuII+=pibose_dens0[nbose]*nbose*4.0*cosh(2.0*nbose*mu[1]);
			}
		}
		drhoII_dmuS=0.25*nh0_b0i1s1*(xI-xxI)*(xS-xxS)
			+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI-xxI)
				+0.25*nh0_b1i2s1*(-xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI);

		rhoS=0.25*nh0_b0i1s1*(xI+xxI)*(xS-xxS)
			+0.5*nh0_b1i0s1*(-xB*xxS+xxB*xS)
				+0.5*nh0_b1i0s3*(-3*xB*xxS*xxS*xxS+3*xxB*xS*xS*xS)
					+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(-xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);
		drhoS_dmuS=0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*nh0_b1i0s3*(9*xB*xxS*xxS*xxS+9*xxB*xS*xS*xS)
					+0.25*nh0_b1i1s2*(4*xB*xxS*xxS+4*xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

		drho[0]=rhoBtarget-rhoB;
		drho[1]=rhoIItarget-rhoII;
		drho[2]=rhoStarget-rhoS;

		A(0,0)=drhoB_dmuB;
		A(0,1)=A(1,0)=drhoB_dmuII;
		A(0,2)=A(2,0)=drhoB_dmuS;
		A(1,1)=drhoII_dmuII;
		A(1,2)=A(2,1)=drhoII_dmuS;
		A(2,2)=drhoS_dmuS;
		dmu=A.colPivHouseholderQr().solve(drho);
		dmumag=sqrt(dmu[0]*dmu[0]+dmu[1]*dmu[1]+dmu[2]*dmu[2]);
		if(dmumag>2.0){
			dmu[0]=2.0*dmu[0]/dmumag;
			dmu[1]=2.0*dmu[1]/dmumag;
			dmu[2]=2.0*dmu[2]/dmumag;
		}
		mu[0]+=dmu[0]; mu[1]+=dmu[1]; mu[2]+=dmu[2];
		
	}while(fabs(drho[0])>1.0E-10 || fabs(drho[1])>1.0E-10 || fabs(drho[2])>1.0E-10);
	muB=mu[0];
	muII=mu[1];
	muS=mu[2];
}

// Same as above, but also calculates T, also uses epsilon (uses GetEpsilonRhoDerivatives to get factors )
void Csampler::GetTfMuNH(Chyper *hyper){
	GetTfMuNH(hyper->epsilon,hyper->rhoB,hyper->rhoII,hyper->rhoS,hyper->muB,hyper->muII,hyper->muS);
	sigmaf=hyper->sigma;
	hyper->T0=Tf;
}

void Csampler::GetTfMuNH(double epsilontarget,double rhoBtarget,double rhoIItarget,double rhoStarget,double &muB,double &muII,double &muS){
	//4D Newton's Method
	// Here rhoII refers to rho_u-rho_d = 2*I3 and mu[1]=muII/2
	double epsilon,rhoB,rhoS,rhoII=0;
	GetNHMu0();
	Eigen::MatrixXd A(4,4);
	Eigen::VectorXd dmu(4),drho(4);
	double cmb,smb;
	int ntries=0;
	do{
		ntries+=1;
		if(ntries>30000){
			snprintf(message,CLog::CHARLENGTH,"FAILURE, ntries=%d\n",ntries);
			CLog::Fatal(message);
		}
		smb=sinh(muB);
		cmb=cosh(muB);

		GetEpsilonRhoDerivatives(muB,muII,muS,epsilon,rhoB,rhoII,rhoS,A);
		for(int i=0;i<4;i++){
			A(i,1)=A(i,1)/cmb;
		}
		drho[0]=epsilontarget-epsilon;
		drho[1]=rhoBtarget-rhoB;
		drho[2]=rhoIItarget-rhoII;
		drho[3]=rhoStarget-rhoS;
		dmu=A.colPivHouseholderQr().solve(drho);
		Tf+=dmu[0];
		GetNHMu0();
		smb+=dmu[1];
		muB=asinh(smb);
		muII+=dmu[2];
		muS+=dmu[3];
	}while(fabs(drho[0])>1.0E-4 || fabs(drho[1])>1.0E-6 || fabs(drho[2])>1.0E-6 || fabs(drho[3])>1.0E-6);
}

void Csampler::GetEpsilonRhoDerivatives(double muB,double muII,double muS,double &epsilon,double &rhoB,double &rhoII,double &rhoS,Eigen::MatrixXd &A){
	double xB,xI,xS,xxB,xxI,xxS;

	double drhoB_dT,drhoB_dmuB,drhoB_dmuS,drhoB_dmuII;
	double drhoS_dT,drhoS_dmuB,drhoS_dmuS,drhoS_dmuII;
	double drhoII_dT,drhoII_dmuB,drhoII_dmuS,drhoII_dmuII;
	double de_dT,de_dmuB,de_dmuII,de_dmuS;

	if(!forMU0_calculated)
		GetNHMu0();

	xB=exp(muB);
	xI=exp(muII);
	xS=exp(muS);
	xxB=1.0/xB;
	xxI=1.0/xI;
	xxS=1.0/xS;

	epsilon=eh0_b0i0s0+0.5*eh0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*eh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*eh0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*eh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*eh0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*eh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*eh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*eh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
									+0.5*eh0_b2i0s0*(xB*xB+xxB*xxB);

	de_dT=dedth0_b0i0s0+0.5*dedth0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*dedth0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*dedth0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*dedth0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*dedth0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*dedth0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*dedth0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*dedth0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
									+0.5*dedth0_b2i0s0*(xB*xB+xxB*xxB);


	de_dmuB=0.5*eh0_b1i0s1*(xB*xxS-xxB*xS)
		+0.5*eh0_b1i0s3*(xB*xxS*xxS*xxS-xxB*xS*xS*xS)
			+0.25*eh0_b1i1s0*(xB-xxB)*(xI+xxI)
				+0.25*eh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI+xxI)
					+0.25*eh0_b1i2s1*(xB*xxS-xxB*xS)*(xI*xI+xxI*xxI)
						+0.25*eh0_b1i3s0*(xB-xxB)*(xI*xI*xI+xxI*xxI*xxI)
							+eh0_b2i0s0*(xB*xB-xxB*xxB);

	de_dmuII=0.5*eh0_b0i2s0*(2*xI*xI-2*xxI*xxI)
		+0.25*eh0_b0i1s1*(xI-xxI)*(xS+xxS)
			+0.25*eh0_b1i1s0*(xB+xxB)*(xI-xxI)
				+0.25*eh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI-xxI)
					+0.25*eh0_b1i2s1*(xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI)
						+0.25*eh0_b1i3s0*(xB+xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);

	de_dmuS=0.25*eh0_b0i1s1*(xI+xxI)*(xS-xxS)
		+0.5*eh0_b1i0s1*(-xB*xxS+xxB*xS)
			+0.5*eh0_b1i0s3*(-3*xB*xxS*xxS*xxS+3*xxB*xS*xS*xS)
				+0.25*eh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI+xxI)
					+0.25*eh0_b1i2s1*(-xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

	rhoB=0.5*nh0_b1i0s1*(xB*xxS-xxB*xS)
		+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS-xxB*xS*xS*xS)
			+0.25*nh0_b1i1s0*(xB-xxB)*(xI+xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(xI*xI+xxI*xxI)
						+0.25*nh0_b1i3s0*(xB-xxB)*(xI*xI*xI+xxI*xxI*xxI)
							+nh0_b2i0s0*(xB*xB-xxB*xxB);

	drhoB_dT=de_dmuB/(Tf*Tf);

	drhoB_dmuB=0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
		+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
			+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
						+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
							+2.0*nh0_b2i0s0*(xB*xB-xxB*xxB);

	drhoB_dmuII=0.25*nh0_b1i1s0*(xB-xxB)*(xI-xxI)
		+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI-xxI)
			+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(2*xI*xI-2*xxI*xxI)
				+0.25*nh0_b1i3s0*(xB-xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);

	drhoB_dmuS=0.5*nh0_b1i0s1*(-xB*xxS-xxB*xS)
		+0.5*nh0_b1i0s3*(-3*xB*xxS*xxS*xxS-3*xxB*xS*xS*xS)
			+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS-2*xxB*xS*xS)*(xI+xxI)
				+0.25*nh0_b1i2s1*(-xB*xxS-xxB*xS)*(xI*xI+xxI*xxI);

	rhoII=0.5*nh0_b0i2s0*(2*xI*xI-2*xxI*xxI)
		+0.25*nh0_b0i1s1*(xI-xxI)*(xS+xxS)
			+0.25*nh0_b1i1s0*(xB+xxB)*(xI-xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI-xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI)
						+0.25*nh0_b1i3s0*(xB+xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);
	
	drhoII_dT=de_dmuII/(Tf*Tf);

	drhoII_dmuB=drhoB_dmuII;

	drhoII_dmuII=0.5*nh0_b0i2s0*(4*xI*xI+4*xxI*xxI)
		+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(4*xI*xI+4*xxI*xxI)
						+0.25*nh0_b1i3s0*(xB+xxB)*(9*xI*xI*xI+9*xxI*xxI*xxI);

	drhoII_dmuS=0.25*nh0_b0i1s1*(xI-xxI)*(xS-xxS)
		+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI-xxI)
			+0.25*nh0_b1i2s1*(-xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI);

	rhoS=0.25*nh0_b0i1s1*(xI+xxI)*(xS-xxS)
		+0.5*nh0_b1i0s1*(-xB*xxS+xxB*xS)
			+0.5*nh0_b1i0s3*(-3*xB*xxS*xxS*xxS+3*xxB*xS*xS*xS)
				+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(-xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

	drhoS_dT=de_dmuS/(Tf*Tf);
	drhoS_dmuB=drhoB_dmuS;
	drhoS_dmuII=drhoII_dmuS;
	drhoS_dmuS=0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
		+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
			+0.5*nh0_b1i0s3*(9*xB*xxS*xxS*xxS+9*xxB*xS*xS*xS)
				+0.25*nh0_b1i1s2*(4*xB*xxS*xxS+4*xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

	if (bose_corr) {
		double ch,sh;
		for(int nbose=2;nbose<=n_bose_corr;nbose++){
			ch=cosh(nbose*muII);
			sh=sinh(nbose*muII);
			epsilon+=pibose_epsilon0[nbose]*(2.0*ch+1.0);
			de_dT+=pibose_dedt0[nbose]*(2.0*ch+1.0);
			de_dmuII+=pibose_epsilon0[nbose]*4.0*nbose*sh;
			rhoII+=pibose_dens0[nbose]*2.0*sh;
			drhoII_dT+=pibose_epsilon0[nbose]*4.0*nbose*sh/(Tf*Tf);
			drhoII_dmuII+=pibose_dens0[nbose]*nbose*4.0*ch;
		}
	}

	A(0,0)=de_dT;
	A(0,1)=de_dmuB;
	A(0,2)=de_dmuII;
	A(0,3)=de_dmuS;

	A(1,0)=drhoB_dT;
	A(1,1)=drhoB_dmuB;
	A(1,2)=drhoB_dmuII;
	A(1,3)=drhoB_dmuS;

	A(2,0)=drhoII_dT;
	A(2,1)=drhoII_dmuB;
	A(2,2)=drhoII_dmuII;
	A(2,3)=drhoII_dmuS;

	A(3,0)=drhoS_dT;
	A(3,1)=drhoS_dmuB;
	A(3,2)=drhoS_dmuII;
	A(3,3)=drhoS_dmuS;

}

void Csampler::CalcRvisc(Chyper *hyper){
	double xB,xI,xS,xxB,xxI,xxS;
	double epsilon,dedt,P,p4overE3,eEbar,m2dens,Rshear,Rbulk,RTbulk,A;
	if(!forMU0_calculated)
		GetNHMu0();
	if(!hyper->epsilon_calculated)
		CalcNHadronsEpsilonP(hyper);
	epsilon=hyper->epsilon;
	P=hyper->P;

	xB=exp(hyper->muB);
	xI=exp(0.5*hyper->muII);
	xS=exp(hyper->muS);
	xxB=1.0/xB;
	xxI=1.0/xI;
	xxS=1.0/xS;

	dedt=dedth0_b0i0s0+0.5*dedth0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*dedth0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*dedth0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*dedth0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*dedth0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*dedth0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*dedth0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*dedth0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
									+0.5*dedth0_b2i0s0*(xB*xB+xxB*xxB);
	
	p4overE3=p4overE3h0_b0i0s0+0.5*p4overE3h0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*p4overE3h0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*p4overE3h0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*p4overE3h0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*p4overE3h0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*p4overE3h0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*p4overE3h0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*p4overE3h0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
									+0.5*p4overE3h0_b2i0s0*(xB*xB+xxB*xxB);

	eEbar=eEbarh0_b0i0s0+0.5*eEbarh0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*eEbarh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*eEbarh0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*eEbarh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*eEbarh0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*eEbarh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*eEbarh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*eEbarh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
									+0.5*eEbarh0_b2i0s0*(xB*xB+xxB*xxB);
	
	m2dens=m2densh0_b0i0s0+0.5*m2densh0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*m2densh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*m2densh0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*m2densh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*m2densh0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*m2densh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*m2densh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*m2densh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
									+0.5*m2densh0_b2i0s0*(xB*xB+xxB*xxB);
	
	A=P/(Tf*Tf*dedt-eEbar);
	Rshear=(2.0/15.0)*p4overE3-2*P;
	Rbulk=-A*epsilon*Tf+(A/3.0)*(dedt*Tf*Tf-m2dens);
	Rbulk=Rbulk-(2.0/3.0)*P+p4overE3/9.0;
	RTbulk=-Rbulk/(A*Tf*Tf);
	hyper->RTbulk=RTbulk; hyper->Rshear=Rshear; hyper->Rbulk=Rbulk;
	hyper->dedT=dedt;
	
	
}

// For sampling, generates mass of resonances (uses GetDensPMax below)
double Csampler::GenerateThermalMass(CresInfo *resinfo){
	map<double,double>::iterator it1,it2;
	double E = 0.0, E1 = 0.0, E2 = 0.0;
	int ires=resinfo->ires;
	double p1 = 0.0, p2 = 0.0, netprob=randy->ran();
	if(!resinfo->decay)
		E=resinfo->mass;
	else{
		it1=sfdens0imap[ires].lower_bound(netprob);
		if(it1==sfdens0imap[ires].end()){
			snprintf(message,CLog::CHARLENGTH,"it1 already at end of map\n");
			CLog::Fatal(message);
		}
		it2=it1;
		--it1;
		p1=it1->first;
		p2=it2->first;
		E1=it1->second;
		E2=it2->second;
		E=((netprob-p1)*E2+(p2-netprob)*E1)/(p2-p1);
	}
	if(E<resinfo->minmass || E1>E || E2<E){
		resinfo->PrintBranchInfo();
		resinfo->PrintSpectralFunction();
		map<double,double>::iterator it;
		for(it=sfdens0imap[ires].begin();it!=sfdens0imap[ires].end();++it){
			snprintf(message,CLog::CHARLENGTH,"Emap=%g, netdens=%g\n",it->second,it->first);
		}
		CLog::Info(message);
		snprintf(message,CLog::CHARLENGTH,"something odd in MSU Sampler, GenerateThermalMass, E=%g, (E1,E2)=(%g,%g), (p1,p2)=(%g,%g), minmass=%g, netprob=%g\n",E,E1, E2,p1,p2,resinfo->minmass,netprob);
		CLog::Fatal(message);
	}
	return E; //returns a random mass proportional to dens*SF'
}

// gets nhadrons, epsilon and P
void Csampler::CalcNHadronsEpsilonP(Chyper *hyper){
	CalcNHadronsEpsilonP(hyper->muB*hyper->T0/Tf,hyper->muII*hyper->T0/Tf,hyper->muS*hyper->T0/Tf,hyper->nhadrons,hyper->epsilon,hyper->P);
	hyper->epsilon_calculated=true;
	printf("?? nhadrons=%g\n",hyper->nhadrons);
}

void Csampler::CalcNHadronsEpsilonP(double muB,double muII,double muS,double &nhadronsf,double &epsilonf,double &Pf){
	double xB,xI,xS,xxB,xxI,xxS,xbose;
	int nbose;
	if(SETMU0){
		nhadronsf=nhadrons0;
		epsilonf=epsilon0;
		Pf=P0;
	}
	else{
		xB=exp(muB);
		xI=exp(muII);
		xS=exp(muS);
		xxB=1.0/xB;
		xxI=1.0/xI;
		xxS=1.0/xS;
		xbose=exp(muII);
		
		nhadronsf=nh0_b0i0s0+0.5*nh0_b0i2s0*(xI*xI+xxI*xxI)
			+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
					+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
						+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
							+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
								+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
									+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
										+0.5*nh0_b2i0s0*(xB*xB+xxB*xxB);
		
		Pf=nhadronsf*Tf;
		epsilonf=eh0_b0i0s0+0.5*eh0_b0i2s0*(xI*xI+xxI*xxI)
			+0.25*eh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.5*eh0_b1i0s1*(xB*xxS+xxB*xS)
					+0.5*eh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
						+0.25*eh0_b1i1s0*(xB+xxB)*(xI+xxI)
							+0.25*eh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
								+0.25*eh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
									+0.25*eh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
										+0.5*eh0_b2i0s0*(xB*xB+xxB*xxB);
		if (bose_corr){
			for(nbose=2;nbose<=n_bose_corr;nbose++){
				xbose=(xbose+1.0+1.0/xbose);
				nhadronsf+=pibose_dens0[nbose]*xbose;
				epsilonf+=pibose_epsilon0[nbose]*xbose;
				Pf+=pibose_P0[nbose]*(pow(xbose,nbose)+pow(xbose,-nbose)+1.0);
			}
		}
	}
}

void Csampler::CalcSFDensMap(CresInfo *resinfo,double T,map<double,double> &sfdensmap){
	int iE,nE;
	double z,E,E2,sfnorm,netprob,k1,k0;
	vector<double> dens;
	nE=resinfo->SpectVec.size();
	dens.resize(nE);
	sfnorm=0.0;
	for(iE=0;iE<nE;iE++){
		E=resinfo->SpectEVec[iE];
		z=E/T;
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		dens[iE]=E*E*k0+2.0*E*T*k1;
		sfnorm+=resinfo->SpectVec[iE]*dens[iE];
	}
	netprob=0.0;
	sfdensmap.insert(pair<double,double>(netprob,resinfo->minmass));
	for(iE=0;iE<nE;iE++){
		E=resinfo->SpectEVec[iE];
		netprob+=resinfo->SpectVec[iE]*dens[iE]/sfnorm;
		if(iE!=nE-1)
			E2=0.5*(resinfo->SpectEVec[iE]+resinfo->SpectEVec[iE+1]);
		else
			E2=2.0*resinfo->SpectEVec[iE]-resinfo->SpectEVec[iE-1];
		sfdensmap.insert(pair<double,double>(netprob,E2));
	}
	dens.clear();
}

