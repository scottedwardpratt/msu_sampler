#include "msu_sampler/hyper.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "msu_sampler/sampler.h"
char *Chyper::message=new char[CLog::CHARLENGTH];

Chyper::Chyper(){
	muB=muII=muS=0.0;
	Rvisc_calculated=false;
	firstcall=true;
	int alpha,beta;
	sigma=0.093;
	sampler=nullptr;
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++)
			pitilde[alpha][beta]=0.0;
	}
	chi4.resize(4,4);
	chi4inv.resize(4,4);
}

void Chyper::SetSampler(Csampler *samplerptr){
	sampler=samplerptr;
}

void Chyper::Print(){
	CLog::Info("----- hyper info -------\n");
	CLog::Info("T0="+to_string(T0)+"\n");
	
	//CLog::Info("T0="+to_string(T0)+",sampler->Tf="+to_string(sampler->Tf)+", sigma="+to_string(sigma)+" rhoB="+to_string(rhoB)+",rhoII="+to_string(rhoII)+"\n");
	
	CLog::Info("muB="+to_string(muB)+",muII="+to_string(muII)+",muS="+to_string(muS)+"\n");
	CLog::Info("epsilon="+to_string(epsilon)+", P="+to_string(P)+", s="+to_string(GetEntropyDensity())+"\n");
	
	snprintf(message,CLog::CHARLENGTH,"dOmega=(%g,%g,%g,%g), udotdOmega=%g\n",dOmega[0],dOmega[1],dOmega[2],dOmega[3],udotdOmega);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"u=(%g,%g,%g,%g).    tau=%g,eta=%g\n",u[0],u[1],u[2],u[3],tau,eta);
	CLog::Info(message);
	for(int alpha=0;alpha<4;alpha++){
			snprintf(message,CLog::CHARLENGTH,"%10.3e %10.3e %10.3e %10.3e\n",
			pitilde[alpha][0],pitilde[alpha][1],pitilde[alpha][2],pitilde[alpha][3]);
			CLog::Info(message);
	}
	snprintf(message,CLog::CHARLENGTH,"-----------------------\n");
	CLog::Info(message);
}

double Chyper::GetEntropyDensity(){
	double s=((epsilon+P)/T0)-muB*rhoB-muII*rhoII-muS*rhoS;
	return s;
}

void Chyper::CalcBiggestpitilde(){
	int alpha,beta;
	Eigen::MatrixXd A(3,3);
	Eigen::VectorXd V(3);
	A(0,0)=0.0;
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++)
			A(alpha-1,beta-1)=pitilde[alpha][beta];
	}
	//Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
	cout << A << endl;
	cout << es.eigenvalues() << endl;
	V=es.eigenvalues();
	snprintf(message,CLog::CHARLENGTH,"eigenvalues: %g,%g,%g\n",V(0),V(1),V(2));
	CLog::Info(message);
	biggestpitilde=0.0;
	for(alpha=0;alpha<3;alpha++){
		if(V(alpha)>biggestpitilde)
			biggestpitilde=V(alpha);
	}
}

void Chyper::Copy(Chyper *oldhyper){
	sigma=oldhyper->sigma;
	rhoB=oldhyper->rhoB;
	rhoII=oldhyper->rhoII;
	rhoS=oldhyper->rhoS;
	tau=oldhyper->tau;
	eta=oldhyper->eta;
	muB=oldhyper->muB;
	muII=oldhyper->muII;
	muS=oldhyper->muS;
	nhadrons=oldhyper->nhadrons;
	Rshear=oldhyper->Rshear;
	Rbulk=oldhyper->Rbulk;
	epsilon=oldhyper->epsilon;
	P=oldhyper->P;
	dedT=oldhyper->dedT;
	biggestpitilde=oldhyper->biggestpitilde;
	udotdOmega=oldhyper->udotdOmega;
	PItilde=oldhyper->PItilde;
	for(int alpha=0;alpha<4;alpha++){
		u[alpha]=oldhyper->u[alpha];
		r[alpha]=oldhyper->r[alpha];
		qmu[alpha]=oldhyper->qmu[alpha];
		dOmega[alpha]=oldhyper->dOmega[alpha];
		for(int beta=0;beta<4;beta++){
			pitilde[alpha][beta]=oldhyper->pitilde[alpha][beta];
		}
	}
}
