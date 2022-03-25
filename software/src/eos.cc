#include "msu_sampler/eos.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/constants.h"
#include "msu_sampler/resonances.h"
using namespace MSU_EOS;

void MSU_EOS::freegascalc_onespecies_finitewidth(double T,CresInfo *resinfo,double &epsilon,double &P,double &dens,double &dedt){
	int iE,nE;
	double E,rho;
	double rhosum=0.0,esum=0.0,psum=0.0,dsum=0.0,dedtsum=0.0;
	nE=resinfo->SpectVec.size();
	for(iE=0;iE<nE;iE++){
		E=resinfo->SpectEVec[iE];
		rho=resinfo->SpectVec[iE];
		freegascalc_onespecies(T,E,epsilon,P,dens,dedt);
		esum+=epsilon*rho;
		psum+=P*rho;
		dsum+=dens*rho;
		dedtsum+=dedt*rho;
		rhosum+=rho;
	}
	epsilon=esum/rhosum;
	P=psum/rhosum;
	dens=dsum/rhosum;
	dedt=dedtsum/rhosum;
}

void MSU_EOS::freegascalc_onespecies_finitewidth(double T,CresInfo *resinfo,double &epsilon,double &P,double &dens,double &dedt,double &p4overE3,double &Ji){
	int iE,nE;
	double E,rho;
	double rhosum=0.0,esum=0.0,psum=0.0,dsum=0.0,p4overE3sum=0.0,dedtsum=0.0,Jisum=0.0;
	nE=resinfo->SpectVec.size();
	for(iE=0;iE<nE;iE++){
		E=resinfo->SpectEVec[iE];
		rho=resinfo->SpectVec[iE];
		freegascalc_onespecies(T,E,epsilon,P,dens,dedt);
		p4overE3=Getp4overE3(T,E,dens);
		Ji=GetJi(T,E,dens);
		esum+=epsilon*rho;
		psum+=P*rho;
		dsum+=dens*rho;
		dedtsum+=dedt*rho;
		p4overE3sum+=p4overE3;
		Jisum+=Ji*rho;
		rhosum+=rho;

	}
	epsilon=esum/rhosum;
	P=psum/rhosum;
	dens=dsum/rhosum;
	dedt=dedtsum/rhosum;
	p4overE3=p4overE3sum/rhosum;
	Ji=Jisum/rhosum;
}

void MSU_EOS::freegascalc_onespecies(double T,double m,double &epsilon,double &P,double &dens,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC_GEV,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	z=m/T;
	if(z>40){ // non-relativistic approximations
		dens=exp(-z)*pow(m*T/(2.0*M_PI*HBARC_GEV*HBARC_GEV),1.5);
		P=dens*T;
		epsilon=(m+1.5*T)*dens;
		dedt=z*z*dens+3.0*z*dens;
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,T);
			exit(1);
		}
		t2=T*T;
		t3=t2*T;
		z=m/T;
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		P=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		dens=P/T;
		epsilon=prefactor*(3.0*m2*t2*k0+(m3*T+6.0*m*t3)*k1);
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*T*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/T)+6.0*m2*T)*k1prime);
	}
}

double MSU_EOS::Getp4overE3(double T,double m,double dens){
	double K2,z,J=0.0,nfact,sign;
	double p4overE3=0.0;
	const int nmax=50;
	double G[nmax+5];
	int n;
	z=m/T;
	if(z<0 || z!=z){
		printf("K0(z<0)\n, z=%g\n",z);
		exit(1);
	}
	if(z<50.0){
		G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
		for(int j=1;j<nmax+5;j++){
			n=5-2*j;
			if(n!=-1)
				G[j]=(-exp(-z)/n)+(G[j-1]*z*z-z*exp(-z))/((n+1.0)*n);
			else G[j]=gsl_sf_gamma_inc(-1,z)*z;
		}
		J=0.0;
		nfact=1.0;
		sign=1.0;
		for(n=0;n<nmax;n+=1){
			if(n>0)
				sign=-1.0;
			J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
			nfact=nfact*0.5/(n+1.0);
			if(n>0)
				nfact*=(2.0*n-1.0);
		}
		//p4overE3=pow(m,4)*(-z*J+15.0*Bessel::Kn(2,z)/(z*z));
		K2=dens*2.0*PI*PI*HBARC_GEV*HBARC_GEV*HBARC_GEV/(m*m*T);
		p4overE3=pow(m,4)*(-z*J+15.0*K2);
		p4overE3=2.0*p4overE3/(4.0*M_PI*M_PI*HBARC_GEV*HBARC_GEV*HBARC_GEV);
	}
	else
		p4overE3=0.0;
	
	/*
	double delp=0.002,pmag,e;
	double p4overE3test;
	printf("__ TESTING, m=%g, p4overE3=%g, J=%g, T=%g, z=%g\n",m,p4overE3,J,T,z);
	p4overE3test=0.0;
	for(pmag=0.5*delp;pmag<10.5;pmag+=delp){
		e=sqrt(m*m+pmag*pmag);
		//p4overE3test+=degen*(4.0*M_PI/(8.0*M_PI*M_PI*M_PI*HBARC_GEV*HBARC_GEV*HBARC_GEV))*pmag*pmag*dp*exp(-e/T)*( (2.0/3.0)*(pmag*pmag/e) - (2.0/15.0)*pow(pmag,4)/pow(e,3) );
		p4overE3test+=(4.0*M_PI/(8.0*M_PI*M_PI*M_PI*HBARC_GEV*HBARC_GEV*HBARC_GEV))*pmag*pmag*delp*exp(-e/T)*pow(pmag,4)/pow(e,3);
	}
	p4overE3=p4overE3test;
	//printf("XXXXXXX p4overE3test=%g, p4overE3test/p4overE3=%g, z=%g\n",p4overE3test,p4overE3test/p4overE3,z);
	*/
		
	return p4overE3;
}

double MSU_EOS::GetJi(double T,double mass,double dens){
	/*
	Ji=\frac{1}{2\pi^2} \int d^3p~e^{-E/T} \frac{p^2}{3*T*E^2}
	*/
	int n;
	double Ji,eta=mass/T;
	double J[100]={0.0};
	double expeta=exp(-eta),kappa;
	J[1]=-gsl_sf_expint_Ei(-eta)/expeta;
	for(n=2;n<100;n++){
		J[n]=(-eta*J[n-1]+1.0)/double(n-1);
	}
	Ji=1.0/eta;
	kappa=1.0;
	for(n=1;n<50;n++){
		kappa*=-(1.5-n)/double(n);
		Ji+=J[2*n]*kappa;
	}
	Ji*=expeta;
	//freegascalc_onespecies(T,mass,epsiloni,Pi,densi,dedti);
	double Jcon=1.0/(2.0*PI*PI*pow(HBARC_GEV,3));
	Ji=dens-mass*mass*mass*Ji*Jcon;
	Ji=Ji/(3.0*T);
	return Ji;
	
	/*
	double p,E,Jtest=0.0,dp=1.0;
	for(p=0.5*dp;p<2000;p+=dp){
		E=sqrt(p*p+mass*mass);
		Jtest+=Jcon*exp(-E/T)*dp*p*p*p*p/(3.0*E*E*T);
	}
	*/
}

void MSU_EOS::GetEpsilonPDens_OneSpecies(double T,CresInfo *resinfo,double &epsiloni,double &Pi,double &densi,double &dedti,double &p4overE3i,double &Ji){
	double degen,m;
	//decay=false;
	m=resinfo->mass;
	degen=resinfo->degen;
	if(resinfo->width>MIN_WIDTH && resinfo->decay && !USE_POLE_MASS){
		freegascalc_onespecies_finitewidth(T,resinfo,epsiloni,Pi,densi,dedti,p4overE3i,Ji);
	}
	else{
		freegascalc_onespecies(T,m,epsiloni,Pi,densi,dedti);
		p4overE3i=Getp4overE3(T,m,densi);
		Ji=GetJi(T,m,densi);
	}
	epsiloni*=degen;
	Pi*=degen;
	dedti*=degen;
	densi*=degen;
	p4overE3i*=degen;
	Ji*=degen;
}

void MSU_EOS::CalcEoSandTransportCoefficients(double T,CresList *reslist,double &epsilon,double &P,double &nh,vector<double> &density,Eigen::Matrix3d &chi,Eigen::Matrix3d &sigma){
	CresInfo *resinfoptr;
	CresInfoMap::iterator rpos;
	double s;
	double pi,epsiloni,densi,dedti,p4overE3i,Ji;
	Eigen::Matrix3d sigmai(3,3);
	int a,b,ires;
	chi.setZero();
	sigma.setZero();
	P=epsilon=s=nh=0.0;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfoptr=rpos->second;
		ires=resinfoptr->ires;
		if(resinfoptr->pid!=22){
			GetEpsilonPDens_OneSpecies(T,resinfoptr,epsiloni,pi,densi,dedti,p4overE3i,Ji);
			P+=pi;
			epsilon+=epsiloni;
			s+=(pi+epsiloni)/T;
			nh+=densi;
			density[ires]=densi;
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					chi(a,b)+=densi*resinfoptr->q[a]*resinfoptr->q[b];
					sigma(a,b)+=Ji*resinfoptr->q[a]*resinfoptr->q[b];
				}
			}
		}
		else{
			density[ires]=0.0;
		}
	}
}

double MSU_EOS::GetLambda(double T,CresList *reslist,double P,double epsilon){
	int i,n;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,dIpp=0.0;
	//double Ipptest=0.0,Ptest=0.0,dIpptest=0.0,dp=4.0,p,e;
	double J,nfact,sign,alpha;
	double lambdafact;
	CresInfo *resinfo;
	CresInfoMap::iterator rpos;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			m=resinfo->mass;
			degen=resinfo->degen;
			z=m/T;
			alpha=0.0;

			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
			for(i=1;i<nmax+5;i++){
				n=5-2*i;
				if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
				else G[i]=gsl_sf_gamma_inc(-1,z)*z;
			}
			J=0.0;
			nfact=1.0;
			sign=1.0;
			for(n=0;n<nmax;n+=1){
				if(n>0) sign=-1.0;
				J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
				nfact=nfact*0.5/(n+1.0);
				if(n>0) nfact*=(2.0*n-1.0);
			}
			dIpp=degen*exp(alpha)*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
			dIpp=dIpp/(60.0*PI*PI*HBARC_GEV*HBARC_GEV*HBARC_GEV);
			/*
			dIpptest=0.0;
			for(p=0.5*dp;p<3000;p+=dp){
			e=sqrt(m*m+p*p);
			dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC_GEV*HBARC_GEV*HBARC_GEV))*p*p*dp*exp(-e/T)*( (2.0/3.0)*(p*p/e) - (2.0/15.0)*pow(p,4)/pow(e,3) );
			//dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC_GEV*HBARC_GEV*HBARC_GEV))*p*p*dp*exp(-e/T)*((2.0/15.0)*pow(p,4)/pow(e,3) );
			}
			Ipptest+=dIpptest;
			*/
			
			Ipp+=dIpp;
			//Ptest+=Ipptest;
		}
	}
	lambdafact=(2.0*P-4.0*Ipp)/(P+epsilon);

	return lambdafact;
}

double MSU_EOS::GetSigma2(double T,double mass){
	double sigma2,Iomega,I1,I2,z;
	z=mass/T;
	Iomega=exp(-z)/(30.0*PI*PI*HBARC*HBARC*HBARC);
	I1=pow(mass,1.5)*pow(T,3.5)*7.5*sqrt(2.0*PI);
	I2=24.0*pow(T,5);
	sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	return sigma2;
}

