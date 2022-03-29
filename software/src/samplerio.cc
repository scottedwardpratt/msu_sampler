#include "msu_sampler/master.h"
#include "msu_commonutils/constants.h"

//#define __TEST_PITILDE_TRACE__

using namespace std;


void CmasterSampler::ReadHyper_BEST_Binary3D(){
	string filename;
	Chyper *elem;
	int ielement=0;
	double u0,ux,uy,uz,utau,ueta;
	double t,x,y,z,tau,eta;
	double udotdOmega,udotdOmega_music;
	double dOmega0,dOmegaX,dOmegaY,dOmegaZ,dOmegaTau,dOmegaEta;
	double epsilonf;
	double Tdec,muB,muS,muC;
	double PIbulk __attribute__((unused)), Pdec __attribute__((unused));
	double qmu_tau, qmu_eta, qmu0,qmu1,qmu2,qmu3;
	double rhoB;
	FourTensor pivisc;

	nelements=0;
	filename=parmap->getS("HYPER_INFO_FILE",string("hyperdata/BEST/surface_2D.dat"));
	printf("opening %s\n",filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"rb");
	if (fptr==NULL) {
		fprintf(stderr,"Can't open hyper info file\n");
		printf("Error %d \n", errno);
	}

	while(!feof(fptr)){
		elem=new Chyper();
		// read from binary file
		float array[34];
		fread(array, sizeof(array),1,fptr);
		tau = array[0];
		x = array[1];
		y = array[2];
		eta = array[3];
		dOmegaTau = array[4];
		dOmegaX = array[5];
		dOmegaY = array[6];
		dOmegaEta = array[7];
		utau = array[8];
		ux = array[9];
		uy = array[10];
		ueta = array[11];
		utau=sqrt(1.0+ux*ux+uy*uy+ueta*ueta);
		
		const double u_milne_sqr = utau * utau - ux * ux - uy * uy - ueta * ueta;
		if (std::abs(u_milne_sqr - 1.0) > 1.e-6) {
			printf("Warning at reading from MUSIC output: "
				"u_Milne (u_eta multiplied by tau) = %9.6f %9.6f %9.6f %9.6f"
					", u^2 == 1 is not fulfilled with error %12.8f.\n",
			utau, ux, uy, ueta, std::abs(u_milne_sqr - 1.0));
		}
		udotdOmega_music = tau * (dOmegaTau * utau +
			dOmegaX * ux + dOmegaY * uy + dOmegaEta * ueta / tau);

		// Transforming from Milne to Cartesian
		const double ch_eta = std::cosh(eta), sh_eta = std::sinh(eta);
		t = tau * ch_eta;
		z = tau * sh_eta;
		u0 = utau * ch_eta + ueta * sh_eta;
		uz = utau * sh_eta + ueta * ch_eta;

		const double usqr = u0 * u0 - ux * ux - uy * uy - uz * uz;
		if (std::abs(usqr - 1.0) > 1.e-3) {
			printf("u*u should be 1, u*u = %12.8f\n", usqr);
		}

		dOmega0 = tau * ch_eta * dOmegaTau - sh_eta * dOmegaEta;
		dOmegaX = -tau * dOmegaX;
		dOmegaY = -tau * dOmegaY;
		dOmegaZ = tau * sh_eta * dOmegaTau - ch_eta * dOmegaEta;

		udotdOmega = dOmega0 * u0 - dOmegaX * ux - dOmegaY * uy - dOmegaZ * uz;
		if (std::abs(udotdOmega - udotdOmega_music) > 1.e-4) {
			printf("u^mu * dsigma_mu should be invariant: %12.9f == %12.9f\n",
			udotdOmega, udotdOmega_music);
		}
 
		epsilonf = array[12]*HBARC_GEV; //was labeled Edec--guessed this was epsilon
		Tdec = array[13]*HBARC_GEV;
		muB  = array[14]*HBARC_GEV/Tdec;
		muS  = array[15]*HBARC_GEV/Tdec;
		muC  = array[16]*HBARC_GEV/Tdec;
		Pdec = array[17]*Tdec - epsilonf;
		
		pivisc[0][0]=array[18]*HBARC_GEV;
		pivisc[0][1]=pivisc[1][0]=array[19]*HBARC_GEV;
		pivisc[0][2]=pivisc[2][0]=array[20]*HBARC_GEV;
		pivisc[0][3]=pivisc[3][0]=array[21]*HBARC_GEV; // /tau;
		pivisc[1][1]=array[22]*HBARC_GEV;
		pivisc[1][2]=pivisc[2][1]=array[23]*HBARC_GEV;
		pivisc[1][3]=pivisc[3][1]=array[24]*HBARC_GEV;  // /tau;
		pivisc[2][2]=array[25]*HBARC_GEV;
		pivisc[2][3]=pivisc[3][2]=array[26]*HBARC_GEV; // /tau;
		pivisc[3][3]=array[27]*HBARC_GEV; // /(tau*tau);
        TransformPiTotz(pivisc, ch_eta, sh_eta);
		
		//printf("reading, trace=%g =? 0\n",pivisc[0][0]-pivisc[1][1]-pivisc[2][2]-pivisc[3][3]);

		PIbulk = array[28]*HBARC_GEV;   // GeV/fm^3
		rhoB = array[29];  // 1/fm^3

		qmu_tau = array[30];
		qmu1 = array[31];
		qmu2 = array[32];
		qmu_eta = array[33];
        qmu0 = qmu_tau*ch_eta + qmu_eta*sh_eta;
        qmu3 = qmu_tau*sh_eta + qmu_eta*ch_eta;

		if(parmap->getB("SAMPLER_BJORKEN_2D",false)){
			double YMAX_ratio=parmap->getD("SAMPLER_BJORKEN_YMAX",1.0)/parmap->getD("HYDRO_BJORKEN_YMAX",1.0);
			dOmega0*=YMAX_ratio;
			dOmegaX*=YMAX_ratio;
			dOmegaY*=YMAX_ratio;
			dOmegaZ*=YMAX_ratio;
		}
		if(udotdOmega >= 0.0) {
			elem->tau=tau;
			elem->dOmega[0]=dOmega0; 
			elem->dOmega[1]=dOmegaX; 
			elem->dOmega[2]=dOmegaY; 
			elem->dOmega[3]=dOmegaZ;

			elem->udotdOmega=udotdOmega;

			elem->r[0]=t;
			elem->r[1]=x;
			elem->r[2]=y;
			elem->r[3]=z;

			elem->u[0]=u0;
			elem->u[1]=ux;
			elem->u[2]=uy;
			elem->u[3]=uz;
			
			GetPitilde(pivisc,elem->pitilde,elem->u);

			elem->muB=muB+0.5*muC;
			elem->muS=muS+0.5*muC;
			elem->muI=muC;
			
			elem->epsilon=epsilonf;
			elem->T0=Tdec;
			elem->rhoB=rhoB;
			elem->rhoS=0.0;
			elem->rhoI=0.0;

			elem->qmu[0]=qmu0;
			elem->qmu[1]=qmu1;
			elem->qmu[2]=qmu2;
			elem->qmu[3]=qmu3;
			elem->rhoB=rhoB;

			hyperlist.push_back(elem);
			ielement+=1;

		}
	}
	nelements=ielement;
}

#include "msu_sampler/master.h"
#include "msu_commonutils/constants.h"

//#define __TEST_PITILDE_TRACE__

using namespace std;


void CmasterSampler::ReadHyper_OSU_2D(){
	string filename;
	Chyper *elem;
	int ielement=0;
	double u0,ux,uy;
	double x,y,tau;
	double udotdOmega,ETAMAX_ratio;
	double dOmega0,dOmegaX,dOmegaY;
	double epsilonf,pitildexx,pitildexy,pitildeyy;
	double Tdec,muB,muS,muC;
	double PIbulk __attribute__((unused)), Pdec __attribute__((unused));
	double qmu0,qmu1,qmu2,qmu3;
	double rhoB;
	char dummy[200];
	double netvolume=0.0;

	if(parmap->getB("SAMPLER_BJORKEN_2D",true)){
		ETAMAX_ratio=2.0*parmap->getD("SAMPLER_BJORKEN_ETAMAX",1.0)/parmap->getD("HYDRO_BJORKEN_ETAMAX",1.0);
	}
	else{
		ETAMAX_ratio=1.0;
	}

	nelements=0;
	filename=parmap->getS("HYPER_INFO_FILE",string("hyperdata/OSU/alice_cent0_5/hyper.txt"));
	printf("opening %s\n",filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"r");
	if (fptr==NULL) {
		fprintf(stderr,"Can't open hyper info file\n");
		printf("Error %d \n", errno);
	}
	fscanf(fptr,"%lf",&Tdec);
	fgets(dummy,200,fptr);	fgets(dummy,200,fptr);
	while(!feof(fptr)){
		elem=new Chyper();
		// read from binary file
		double readstuff[11];
		for(int iread=0;iread<11;iread++)
			fscanf(fptr,"%lf",&readstuff[iread]);
		tau=readstuff[0];
		x=readstuff[1];
		y=readstuff[2];
		ux=readstuff[3];
		uy=readstuff[4];
		u0=sqrt(1.0+ux*ux+uy*uy);
		dOmega0=readstuff[5]*ETAMAX_ratio;
		dOmegaX=readstuff[6]*ETAMAX_ratio;
		dOmegaY=readstuff[7]*ETAMAX_ratio;
		pitildexx=readstuff[8];
		pitildexy=readstuff[9];
		pitildeyy=readstuff[10];

		muB=muS=muC=0.0;
		PIbulk=0.0;
		qmu0=qmu1=qmu2=qmu3=0.0;
		rhoB=0.0;
		epsilonf=-1.0; // not used

		udotdOmega=dOmega0*u0-dOmegaX*ux-dOmegaY*uy;
		if(udotdOmega >= 0.0) {
			netvolume+=udotdOmega;
			elem->tau=tau;
			elem->dOmega[0]=dOmega0; 
			elem->dOmega[1]=dOmegaX; 
			elem->dOmega[2]=dOmegaY; 
			elem->dOmega[3]=0.0;

			elem->udotdOmega=udotdOmega;

			elem->r[0]=tau;
			elem->r[1]=x;
			elem->r[2]=y;
			elem->r[3]=0.0;

			elem->u[0]=u0;
			elem->u[1]=ux;
			elem->u[2]=uy;
			elem->u[3]=0.0;

			elem->pitilde[0][0]=elem->pitilde[0][1]=elem->pitilde[0][2]=elem->pitilde[0][3]=0.0;
			elem->pitilde[3][0]=elem->pitilde[3][1]=elem->pitilde[3][2]=0.0;
			elem->pitilde[1][0]=elem->pitilde[1][3]=0.0;
			elem->pitilde[2][0]=elem->pitilde[2][3]=0.0;
			elem->pitilde[1][1]=pitildexx;
			elem->pitilde[2][2]=pitildeyy;
			elem->pitilde[1][2]=elem->pitilde[2][1]=pitildexy;
			elem->pitilde[3][3]=-pitildexx-pitildeyy;
			
			elem->muB=muB+0.5*muC;
			elem->muS=muS+0.5*muC;
			elem->muI=muC;
			
			elem->epsilon=epsilonf;
			elem->T0=Tdec;
			elem->rhoB=0.0;
			elem->rhoS=0.0;
			elem->rhoI=0.0;

			elem->qmu[0]=qmu0;
			elem->qmu[1]=qmu1;
			elem->qmu[2]=qmu2;
			elem->qmu[3]=qmu3;
			elem->rhoB=rhoB;

			hyperlist.push_back(elem);
			ielement+=1;

		}
	}
	nelements=ielement;
	printf("netvolume=%g\n",netvolume);
}


