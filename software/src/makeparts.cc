#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "msu_sampler/sampler.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/misc.h"
#include <iostream>

using namespace std;
using namespace NMSUPratt;

int Csampler::MakeParts(Chyper *hyper){
	int nparts=0,dnparts,ires,nbose,II3;
	CresInfo *resinfo;
	double udotdOmega=hyper->udotdOmega;
	double dN,dNtot=0,dNtotprime=0,mutot=0,fugacity;
	double dNcheck=0.0;
	CresMassMap::iterator iter;
	double fugacity_u=hyper->fugacity_u;
	double fugacity_d=hyper->fugacity_d;
	double fugacity_s=hyper->fugacity_s;
	
	
	if(SETMU0)
		dNtot=dNtotprime=udotdOmega*nhadrons0;
	else
		dNtot=dNtotprime=udotdOmega*hyper->nhadrons;
	
	dNtot*=NSAMPLE;
	dNtotprime*=NSAMPLE;
	if(randy->test_threshold(dNtot)){
		dNcheck=0.0;
		for(iter=reslist->massmap.begin();iter!=reslist->massmap.end();++iter){
			resinfo=iter->second;
			if(resinfo->pid!=22){
				ires=resinfo->ires;
				II3=2.0*resinfo->charge-resinfo->baryon-resinfo->strange;
				mutot=hyper->muB*resinfo->baryon+hyper->muII*II3+hyper->muS*resinfo->strange;
				mutot=mutot*hyper->T0/Tf;
				fugacity=pow(fugacity_u,abs(resinfo->Nu))*pow(fugacity_d,abs(resinfo->Nd))*pow(fugacity_s,abs(resinfo->Ns));
				dN=NSAMPLE*fugacity*exp(mutot)*density0i[ires]*udotdOmega;
				dNcheck+=dN;
				dNtotprime-=dN;
				if(dNtotprime<-0.0001){
					snprintf(message,CLog::CHARLENGTH,"res=%d, dNtotprime=%g, should not be negative, mutot=%g, dNcheck=%g, dNtot=%g\n",
					ires,dNtotprime,mutot,dNcheck,dNtot);
					CLog::Info(message);
					resinfo->Print();
					snprintf(message,CLog::CHARLENGTH,"nhadrons0=%g, hyper->nhadrons=%g\n",nhadrons0,hyper->nhadrons);
					CLog::Fatal(message);
				}
				dnparts=CheckResInVolume(dN,Tf,resinfo,hyper);
				nparts+=dnparts;
				if(!(randy->test_threshold(dNtotprime))){
					randy->increment_netprob(dNtotprime);
					goto NoMoreParts;
				}
			}
		}
		if(bose_corr){
			fugacity=fugacity_u*fugacity_d;
			
			for(nbose=2;nbose<=n_bose_corr;nbose++){
				resinfo=reslist->GetResInfoPtr(211);
				ires=resinfo->ires;
				mutot=2.0*nbose*hyper->muII*hyper->T0/Tf;
				dN=NSAMPLE*pow(fugacity,nbose)*exp(mutot)*pibose_dens0[nbose]*udotdOmega;
				dNcheck+=dN;
				dNtotprime-=dN;
				dnparts=CheckResInVolume(dN,Tf/double(nbose),resinfo,hyper);
				nparts+=dnparts;
				if(!(randy->test_threshold(dNtotprime))){
					randy->increment_netprob(dNtotprime);
					goto NoMoreParts;
				}
			}
			
			for(nbose=2;nbose<=n_bose_corr;nbose++){	
				resinfo=reslist->GetResInfoPtr(111);
				ires=resinfo->ires;
				mutot=0.0;
				dN=NSAMPLE*pow(fugacity,nbose)*exp(mutot)*pibose_dens0[nbose]*udotdOmega;
				dNtotprime-=dN;
				dnparts=CheckResInVolume(dN,Tf/double(nbose),resinfo,hyper);
				nparts+=dnparts;
				if(!(randy->test_threshold(dNtotprime))){
					randy->increment_netprob(dNtotprime);
					goto NoMoreParts;
				}
			}
				
			for(nbose=2;nbose<=n_bose_corr;nbose++){
				resinfo=reslist->GetResInfoPtr(-211);
				ires=resinfo->ires;
				mutot=-2.0*nbose*hyper->muII*hyper->T0/Tf;
				dN=NSAMPLE*pow(fugacity,nbose)*exp(mutot)*pibose_dens0[nbose]*udotdOmega;
				dNcheck+=dN;
				dNtotprime-=dN;
				dnparts=CheckResInVolume(dN,Tf/double(nbose),resinfo,hyper);
				nparts+=dnparts;
				if(!(randy->test_threshold(dNtotprime))){
					randy->increment_netprob(dNtotprime);
					goto NoMoreParts;
				}
			}
		}
		NoMoreParts:
		
		if(dNcheck>dNtot*1.01){
			snprintf(message,CLog::CHARLENGTH,"Inside Csampler::MakeParts dNcheck=%g > dNtot=%g, dNtotprime=%g, T=%g\n",
			dNcheck,dNtot,dNtotprime,Tf);
			CLog::Fatal(message);
		}
		return nparts;
	}
	else{
		randy->increment_netprob(dNtot);
		return 0;
	}
}

int Csampler::CheckResInVolume(double dN,double T,CresInfo *resinfo,Chyper *hyper){
	int dnparts=0,alpha;
	double eta;
	FourVector p,r,ubj,ptilde,rtilde;
	randy->increment_netprob(dN);
	while(randy->test_threshold(0.0)){
		GetP(hyper,T,resinfo,p);
		if(BJORKEN_2D){
			if(fabs(hyper->u[3])>0.01){
				snprintf(message,CLog::CHARLENGTH,"BJORKEN_2D set, but u_z=%g\n",hyper->u[3]);
				CLog::Fatal(message);
			}
			eta=BJORKEN_YMAX*(1.0-randy->ran());
			ubj[1]=ubj[2]=0.0;
			ubj[0]=cosh(eta);
			ubj[3]=sinh(eta);
			Misc::Boost(ubj,p,ptilde);
			Misc::Boost(ubj,hyper->r,rtilde);
			for(alpha=0;alpha<4;alpha++){
				p[alpha]=ptilde[alpha];
				r[alpha]=rtilde[alpha];
			}
		}
		for(alpha=0;alpha<4;alpha++)
			r[alpha]=hyper->r[alpha];
		partlist->AddPart(resinfo->pid,p,r);
		dnparts+=1;
		randy->increase_threshold();
	}
	return dnparts;
}

void Csampler::GetPInFluidFrame(double m,Chyper *hyper,double T,FourVector &p){
	FourVector pnoshear,pnobulk;
	double Tbulk;
	int alpha;
	if(hyper->Rvisc_calculated==false){
		CalcRvisc(hyper);
		hyper->Rvisc_calculated=true;
	}
	if(INCLUDE_BULK_VISCOSITY){
		Tbulk=T-hyper->PItilde/hyper->RTbulk;
		randy->generate_boltzmann(m,Tbulk,pnobulk);
		BulkScale(hyper,m,pnobulk,pnoshear);
	}
	else
		randy->generate_boltzmann(m,T,pnoshear);
	if(INCLUDE_SHEAR_VISCOSITY){
		ShearScale(hyper,m,pnoshear,p);
	}
	else{
		for(alpha=0;alpha<4;alpha++)
			p[alpha]=pnoshear[alpha];
	}
}

void Csampler::BulkScale(Chyper *hyper,double mass,FourVector &pnobulk,FourVector &p){
	int alpha;
	for(alpha=1;alpha<4;alpha++)
		p[alpha]=pnobulk[alpha]/(1.0+hyper->PItilde/(3.0*hyper->Rbulk));
	p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
}

void Csampler::ShearScale(Chyper *hyper,double mass,FourVector &pnoshear,FourVector &p){
	int alpha,beta;
	double R,Rmax,pmag,ctheta,stheta,phi;
	for(alpha=0;alpha<4;alpha++)
		p[alpha]=pnoshear[alpha];
	pmag=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
	Rmax=1.0-hyper->biggestpitilde*pmag*pmag/(p[0]*Tf*hyper->Rshear);

	R=1.0;
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			R-=hyper->pitilde[alpha][beta]*p[alpha]*p[beta]/(p[0]*Tf*hyper->Rshear);
		}
	}
	while(randy->ran()>R/Rmax){
		ctheta=1.0-2.0*randy->ran();
		stheta=sqrt(1.0-ctheta*ctheta);
		phi=2.0*M_PI*randy->ran();
		p[1]=pmag*stheta*cos(phi);
		p[2]=pmag*stheta*sin(phi);
		p[3]=pmag*ctheta;
		p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+mass*mass);
		R=1.0;
		for(alpha=1;alpha<4;alpha++){
			for(beta=1;beta<4;beta++){
				R-=hyper->pitilde[alpha][beta]*p[alpha]*p[beta]/(p[0]*Tf*hyper->Rshear);
			}
		}
		if(R>Rmax){
			snprintf(message,CLog::CHARLENGTH,"R=%g, pmag=%g, Rmax=%g, biggestpitilde=%g, Rshear=%g\n",R,pmag,Rmax,hyper->biggestpitilde,hyper->Rshear);
			CLog::Info(message);
		}
	};
	
	/*
	for(alpha=1;alpha<4;alpha++){
	p[alpha]=pnoshear[alpha];
	for(beta=1;beta<4;beta++){
	p[alpha]-=hyper->pitilde[alpha][beta]*pnoshear[beta]/hyper->Rshear;
	}
	}
	p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+mass*mass);
	*/
	
}

void Csampler::GetP(Chyper *hyper,double T,CresInfo *resinfo,FourVector &p){
	bool reflect;
	double pdotdOmega,nhatnorm,nhatdotp,wreflect;
	FourVector dOmegaTilde,ptilde;
	double m,nhat[4]={0.0};
	if((!resinfo->decay) || resinfo->width<0.0001 || USE_POLE_MASS){
		m=resinfo->mass;
	}
	else{
		m=GenerateThermalMass(resinfo);
	}
	GetPInFluidFrame(m,hyper,T,ptilde);

	Misc::BoostToCM(hyper->u,hyper->dOmega,dOmegaTilde);  //dOmegaTilde is dOmega in fluid (u=0) frame
	pdotdOmega=ptilde[0]*dOmegaTilde[0]-ptilde[1]*dOmegaTilde[1]-ptilde[2]*dOmegaTilde[2];
	wreflect=pdotdOmega/(ptilde[0]*dOmegaTilde[0]);
	reflect=false;
	if(wreflect<0.0)
		reflect=true;
	if(wreflect<1.0){
		if(wreflect<randy->ran())
			reflect=true;
	}
	if(reflect){
		nhatnorm=sqrt(dOmegaTilde[1]*dOmegaTilde[1]+dOmegaTilde[2]*dOmegaTilde[2]);
		nhat[1]=dOmegaTilde[1]/nhatnorm;
		nhat[2]=dOmegaTilde[2]/nhatnorm;
		nhatdotp=nhat[1]*ptilde[1]+nhat[2]*ptilde[2];
		ptilde[1]-=2.0*nhat[1]*nhatdotp;
		ptilde[2]-=2.0*nhat[2]*nhatdotp;
	}
	Misc::Boost(hyper->u,ptilde,p);
}
