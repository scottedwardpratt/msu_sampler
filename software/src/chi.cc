#include "msu_sampler/sampler.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/constants.h"
using namespace NMSUPratt;

void Csampler::CalcChiWithFugacity(Chyper *hyper){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos;
	double Pi,epsiloni,densi,dedti,p4overE3i,Ji;
	double x;
	int B,S,II,fugacity;
	int a,b;
	double f_u=hyper->fugacity_u;
	double f_d=hyper->fugacity_d;
	double f_s=hyper->fugacity_s;
	for(a=0;a<3;a++){
		for(b=0;b<3;b++){
			chi(a,b)=0.0;
		}
	}
	
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			fugacity=f_u*resinfo->Nu+f_d*resinfo->Nd+f_s*resinfo->Ns;
			if(resinfo->baryon==0 && (resinfo->Nu+resinfo->Nd+resinfo->Ns!=2)){
				CLog::Info("quark numbers do not add up\n");
				resinfo->Print();
				exit(1);
			}
			if(abs(resinfo->baryon)==1 && (resinfo->Nu+resinfo->Nd+resinfo->Ns!=3)){
				CLog::Info("quark numbers do not add up\n");
				resinfo->Print();
				exit(1);
			}
			
			MSU_EOS::GetEpsilonPDens_OneSpecies(Tf,resinfo,epsiloni,Pi,densi,dedti,p4overE3i,Ji,true);
			B=resinfo->baryon;
			S=resinfo->strange;
			II=resinfo->q[0]-resinfo->q[1];
			x=fugacity*exp(hyper->muB*B+hyper->muII*II+hyper->muS*S);
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					chi(a,b)+=x*densi*resinfo->q[a]*resinfo->q[b];
				}
			}
		}
	}
	chiinv=chi.inverse();
}

/*
void Csampler::CalcChi(Chyper *hyper){
	GetEpsilonRhoChi(hyper->muB,hyper->muII,hyper->muS,hyper->epsilon,hyper->rhoB,hyper->rhoII,hyper->rhoS,hyper->chi4);
		hyper->chi4inv=(hyper->chi4).inverse();
		hyper->epsilon_calculated=true;
	}
*/

/*
void Csampler::CalcChiSlow(Chyper *hyper){
	GetEpsilonRhoChiSlow(hyper->muB,hyper->muII,hyper->muS,hyper->epsilon,hyper->P,hyper->rhoB,hyper->rhoII,hyper->rhoS,hyper->nhadrons,hyper->chi4);
	hyper->chi4inv=(hyper->chi4).inverse();
	hyper->epsilon_calculated=true;
}
*/

/*
void Csampler::CalcChi4BQS(Chyper *hyper){
	GetEpsilonRhoChi4BQS(hyper->muB,hyper->muII,hyper->muS,hyper->epsilon,hyper->P,hyper->rhoB,hyper->rhoII,hyper->rhoS,hyper->nhadrons,hyper->chi4BQS,hyper->dedT);
	hyper->chi4BQSinv=(hyper->chi4BQS).inverse();
	hyper->epsilon_calculated=true;
}
*/

/*
void Csampler::GetEpsilonRhoChiSlow(double muB,double muII,double muS,double &epsilon,double &P,double &rhoB,double &rhoII,double &rhoS,double &nhadrons,Eigen::Matrix<double,4,4> &chi){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos;
	double Pi,epsiloni,densi,dedti,p4overE3i,Ji;
	double x;
	int B,S,II;
	int a,b;
	rhoII=rhoB=rhoS=0.0;
	epsilon=P=rhoB=rhoII=rhoS=nhadrons=0.0;
	for(a=0;a<4;a++){
		for(b=0;b<4;b++)
			chi(a,b)=0.0;
	}
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			MSU_EOS::GetEpsilonPDens_OneSpecies(Tf,resinfo,epsiloni,Pi,densi,dedti,p4overE3i,Ji,true);
			B=resinfo->baryon;
			S=resinfo->strange;
			II=resinfo->q[0]-resinfo->q[1];
			x=exp(muB*B+muII*II+muS*S);
			rhoB+=B*densi*x;
			rhoII+=II*densi*x;
			rhoS+=S*densi*x;
			nhadrons+=densi*x;
			epsilon+=epsiloni*x;
			P+=Pi*x;
			
			chi(0,0)+=Tf*Tf*dedti*x;
			chi(1,1)+=densi*B*B*x;
			chi(2,2)+=densi*II*II*x;
			chi(3,3)+=densi*S*S*x;
			chi(0,1)+=epsiloni*B*x;
			chi(0,2)+=epsiloni*II*x;
			chi(0,3)+=epsiloni*S*x;
			chi(1,2)+=densi*B*II*x;
			chi(1,3)+=densi*B*S*x;
			chi(2,3)+=densi*II*S*x;
		}
		chi(1,0)=chi(0,1);
		chi(2,0)=chi(0,2);
		chi(3,0)=chi(0,3);
		chi(2,1)=chi(1,2);
		chi(3,1)=chi(1,3);
		chi(3,2)=chi(2,3);
	}			
}
*/

/*
void Csampler::GetEpsilonRhoChi(double muB,double muII,double muS,double &epsilon,double &rhoB,double &rhoII,double &rhoS,Eigen::Matrix<double,4,4> &chi){
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
	
	chi(0,0)=Tf*Tf*de_dT;
	chi(0,1)=de_dmuB;
	chi(0,2)=de_dmuII;
	chi(0,3)=de_dmuS;

	chi(1,0)=Tf*Tf*drhoB_dT;
	chi(1,1)=drhoB_dmuB;
	chi(1,2)=drhoB_dmuII;
	chi(1,3)=drhoB_dmuS;

	chi(2,0)=Tf*Tf*drhoII_dT;
	chi(2,1)=drhoII_dmuB;
	chi(2,2)=drhoII_dmuII;
	chi(2,3)=drhoII_dmuS;

	chi(3,0)=Tf*Tf*drhoS_dT;
	chi(3,1)=drhoS_dmuB;
	chi(3,2)=drhoS_dmuII;
	chi(3,3)=drhoS_dmuS;

}
*/
/*
void Csampler::GetEpsilonRhoChi4BQS(double muB,double muII,double muS,double &epsilon,double &P,double &rhoB,double &rhoII,double &rhoS,double &nhadrons,Eigen::Matrix<double,4,4> &chi4BQS,double &dedT){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos;
	double Pi,epsiloni,densi,dedti,p4overE3i,Ji;
	double x;
	int B,S,II,Q;
	int a,b;
	rhoII=rhoB=rhoS=0.0;
	epsilon=P=rhoB=rhoII=rhoS=nhadrons=dedT=0.0;
	for(a=0;a<4;a++){
		for(b=0;b<4;b++)
			chi4BQS(a,b)=0.0;
	}
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			MSU_EOS::GetEpsilonPDens_OneSpecies(Tf,resinfo,epsiloni,Pi,densi,dedti,p4overE3i,Ji,true);
			B=resinfo->baryon;
			S=resinfo->strange;
			Q=resinfo->charge;
			II=resinfo->q[0]-resinfo->q[1];
			x=exp(muB*B+muII*II+muS*S);
			rhoB+=B*densi*x;
			rhoII+=II*densi*x;
			rhoS+=S*densi*x;
			nhadrons+=densi*x;
			epsilon+=epsiloni*x;
			dedT+=dedti*x;
			P+=Pi*x;
			
			chi4BQS(0,0)+=Tf*Tf*dedti*x;
			chi4BQS(1,1)+=densi*B*B*x;
			chi4BQS(2,2)+=densi*Q*Q*x;
			chi4BQS(3,3)+=densi*S*S*x;
			chi4BQS(0,1)+=epsiloni*B*x;
			chi4BQS(0,2)+=epsiloni*Q*x;
			chi4BQS(0,3)+=epsiloni*S*x;
			chi4BQS(1,2)+=densi*B*Q*x;
			chi4BQS(1,3)+=densi*B*S*x;
			chi4BQS(2,3)+=densi*Q*S*x;
		}
		chi4BQS(1,0)=chi4BQS(0,1);
		chi4BQS(2,0)=chi4BQS(0,2);
		chi4BQS(3,0)=chi4BQS(0,3);
		chi4BQS(2,1)=chi4BQS(1,2);
		chi4BQS(3,1)=chi4BQS(1,3);
		chi4BQS(3,2)=chi4BQS(2,3);
	}
}
*/