#define __RESINFO_CC__
#include "msu_sampler/resonances.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/log.h"


Crandy *CresInfo::randy=NULL;
int CresInfo::NSPECTRAL=100;
string CresInfo::SFDIRNAME="../local/resinfo/spectralfunctions";

CresInfo::CresInfo(){
	minmass=1.0E20;
	SFcalculated=false;
	int b=branchlist.size()-1;
	while(b>=0){
		delete branchlist[b];
		b=branchlist.size()-1;
	}
	branchlist.clear();
	count=0;
}

CresInfo::~CresInfo(){
	int b=branchlist.size()-1;
	for(b=0;b<int(branchlist.size());b++){
		delete branchlist[b];
	}
	branchlist.clear();
	SpectVec.clear();
	SpectEVec.clear();
	GammaVec.clear();
}

CbranchInfo::CbranchInfo(){
}

CbranchInfo::~CbranchInfo(){
	resinfo.clear();
}

void CresInfo::PrintBranchInfo(){
	Print();
	if(decay){
		printf(" ------  branches -------\n");
		for(int ib=0;ib<int(branchlist.size());ib++){
			printf("%2d, branching=%5.3f, ",ib,branchlist[ib]->branching);
			for(int ir=0;ir<int(branchlist[ib]->resinfo.size());ir++)
				printf("%6d ",branchlist[ib]->resinfo[ir]->pid);
			printf("Gamma_i=%g\n",width*branchlist[ib]->branching);
		}
	}
}

void CresInfo::Print(){
	printf("+++++++ ID=%d, M=%g, M_min=%g, %s +++++++++\n",pid,mass,minmass,name.c_str());
	printf("Gamma=%g, Degen=%g, Decay=%d\n",width,degen,int(decay));
	printf("Q=%d, B=%d, S=%d, G_parity=%d\n",charge,baryon,strange,G_Parity);
}

void CresInfo::CalcMinMass(){
	bptr_minmass=NULL;
	if(decay){
		minmass=SpectEVec[0];
		double mbranch;
		int ibranch,nbodies,ibody;
		CbranchInfo *bptr;
		for(ibranch=0;ibranch<int(branchlist.size());ibranch++){
			bptr=branchlist[ibranch];
			mbranch=0.0;
			nbodies=bptr->resinfo.size();
			for(ibody=0;ibody<nbodies;ibody++){
				if(bptr->resinfo[ibody]->minmass>1.0E8){
					bptr->resinfo[ibody]->CalcMinMass();
				}
				mbranch+=bptr->resinfo[ibody]->minmass;
			}
			if(mbranch<minmass){
				bptr_minmass=bptr;
			}
		}
	}
	else
		minmass=mass;
}

void CresInfo::ReadSpectralFunction(){
	char filename[200];
	double E,Gamma,SF,netprob;
	if(abs(baryon)!=0)
		sprintf(filename,"%s/%d.txt",SFDIRNAME.c_str(),abs(pid));
	else
		sprintf(filename,"%s/%d.txt",SFDIRNAME.c_str(),pid);
	FILE *fptr=fopen(filename,"r");
	if (fptr==NULL) {
		fprintf(stderr,"Can't open spectral function file\n");
		printf("Error %d \n", errno);
		exit(1);
	}
	fscanf(fptr,"%lf",&E);
	do{
		fscanf(fptr,"%lf %lf %lf",&Gamma,&SF,&netprob);
		SpectEVec.push_back(E);
		GammaVec.push_back(Gamma);
		fscanf(fptr,"%lf",&E);
		SpectVec.push_back(SF);
	}while(!feof(fptr));
	fclose(fptr);
	SFcalculated=true;
}

void CresInfo::CalcSpectralFunction(){
	if(!decay){
		printf("Calling CalcSpectralFunction() for stable particle\n");
		exit(1);
	}
	int n,ibranch;
	double E,M0=mass,A,rhoratio;
	double rho_ab,rho_ab0,Gamma,Gamma0,dGamma;
	CresInfo *resinfo_a,*resinfo_b;
	CresInfo *resinfoswitch;
	if(SpectVec.size()==0)
		SpectVec.resize(NSPECTRAL);
	for(n=0;n<NSPECTRAL;n++){
		E=GetEofN(n);
		Gamma=0;
		for(ibranch=0;ibranch<int(branchlist.size());ibranch++){
			dGamma=0.0;
			Gamma0=width*branchlist[ibranch]->branching;
			resinfo_a=branchlist[ibranch]->resinfo[0];
			resinfo_b=branchlist[ibranch]->resinfo[1];
			if(E>resinfo_a->minmass+resinfo_b->minmass){
				if(resinfo_b->decay && !resinfo_a->decay){
					resinfoswitch=resinfo_a;
					resinfo_a=resinfo_b;
					resinfo_b=resinfoswitch;
				}
				rho_ab0=GetRhoAB(mass,resinfo_a,resinfo_b,branchlist[ibranch]->L);
				if(E>minmass && rho_ab0>1.0E-5){
					rho_ab=GetRhoAB(E,resinfo_a,resinfo_b,branchlist[ibranch]->L);
					rhoratio=rho_ab/rho_ab0;
					dGamma=Gamma0*rhoratio;
					if(dGamma!=dGamma){
						printf("Disaster, dGamma=%g\n",dGamma);
						Print();
						printf("---------------------------------\n");
						printf("ibranch=%d, E=%g, Gamma_ab=%g\n",ibranch,E,Gamma0*rho_ab/rho_ab0);
						printf("Gamma0=%g, rho_ab=%g, rho_ab0=%g\n",Gamma0,rho_ab,rho_ab0);
						printf("mass-m_amin-m_bmin=%g\n",mass-resinfo_a->minmass-resinfo_b->minmass);
						resinfo_a->Print();
						resinfo_b->Print();
						exit(1);
					}
				}
			}
			Gamma+=dGamma;
		}
		//if(pid==229 && E>minmass)
		//	printf("\n");
		A=GetBW(E,M0,Gamma);
		SpectVec[n]=A/GetBW_base(E,M0,width);
		GammaVec[n]=Gamma;
	}
	NormalizeSF();
	SFcalculated=true;
}

double CresInfo::GetSpectralFunction(double E){  // FIX THIS
	int n;
	double delE,E0;
	E0=SpectEVec[0];
	if(SpectVec.size()>1){
		delE=SpectEVec[1]-E0;
		n=lrint((E-E0)/delE);
		if(n<int(SpectVec.size())){
			return SpectVec[n]/delE;
		}
		else
			return 0.0;
	}
	else{
		return 0.0;
	}
}

double CresInfo::GetEofN(int n){
	return mass+0.5*width*tan(M_PI*(double(n)+0.5-0.5*NSPECTRAL)/double(NSPECTRAL));
}

double CresInfo::GetMeshE(double E){
	int n=floorl(NSPECTRAL*(0.5+atan2(E-mass,0.5*width)/M_PI));
	return GetEofN(n);
}

double CresInfo::GetRhoAB(double E,CresInfo *resinfo_a,CresInfo *resinfo_b,int L){
	//printf("CHECK, starting CresInfo::GetRhoAB()\n");
	double pf=0.0,rho_ab=0.0,drho=0.0;
	double Ea,Aa,Eb,Ab;
	int na,nb;
	if(!resinfo_a->decay){  // both a & b stable
		pf=GetDecayMomentum(E,resinfo_a->mass,resinfo_b->mass);
		rho_ab=(pf/E)*GetBL2(pf,L);
	}
	else if(!resinfo_b->decay){ // a is unstable, b is stable
		if(!(resinfo_a->SFcalculated))
			resinfo_a->CalcSpectralFunction();
		for(na=0;na<NSPECTRAL;na++){
			Ea=resinfo_a->GetEofN(na);
			if(Ea>=resinfo_a->minmass){
				Aa=resinfo_a->SpectVec[na];
				pf=GetDecayMomentum(E,Ea,resinfo_b->mass);
				drho=Aa*(pf/E)*GetBL2(pf,L)*GetFF(E,Ea,resinfo_b->mass,resinfo_a,resinfo_b);
				rho_ab+=drho;
			}
		}
	}
	else{   // both a and b are unstable
		if(!(resinfo_a->SFcalculated))
			resinfo_a->CalcSpectralFunction();
		if(!(resinfo_b->SFcalculated))
			resinfo_b->CalcSpectralFunction();
		for(na=0;na<NSPECTRAL;na++){
			Ea=resinfo_a->GetEofN(na);
			if(Ea>=resinfo_a->minmass){
				Aa=resinfo_a->SpectVec[na];
				for(nb=0;nb<NSPECTRAL;nb++){
					Eb=resinfo_b->GetEofN(nb);
					if(Ea+Eb<E && Eb>=resinfo_b->minmass){
						Ab=resinfo_b->SpectVec[nb];
						pf=GetDecayMomentum(E,Ea,Eb);
						drho=Aa*Ab*(pf/E)*GetBL2(pf,L)*GetFF(E,Ea,Eb,resinfo_a,resinfo_b);
						rho_ab+=drho;
					}
				}
			}
		}
	}
	return rho_ab;
}

double CresInfo::GetFF(double E,double Ea,double Eb,CresInfo *resinfo_a,CresInfo *resinfo_b){
	// For Factor From Eq. (36), arXiv:1606.06642v2 -- only if at least one daughter product is unstable
	double lambda=0.0,s0,FF;
	s0=(Ea+Eb)*(Ea+Eb);
	if(!resinfo_a->decay && !resinfo_b->decay)
		return 1.0;

	if(resinfo_a->decay && resinfo_b->decay){
		lambda=0.6;
		//lambda=0.8;
	}
	else if(resinfo_a->baryon!=0){
		lambda=2.0;
		//lambda=0.5;
	}
	else{
		lambda=1.6;
		if( (resinfo_a->pid==113 || abs(resinfo_a->pid)==213) && (resinfo_b->pid==111 || abs(resinfo_b->pid)==211))
			lambda=0.8;
	}
	FF=(pow(lambda,4)+0.25*pow(s0-mass*mass,2))
		/(pow(lambda,4)+pow( E*E-0.5*(s0+mass*mass),2));
	//if(pid==229 && resinfo_a->pid==223 && resinfo_b->pid==223 && FF>18){
	/* if(FF>25 && E<mass){
	printf("pid=%d, mass=%g, E=%g, Ea=%g, Eb=%g, lambda=%g, FF=%g\n",pid,mass,E,Ea,Eb,lambda,FF);
	resinfo_a->Print();
	}*/
	return FF*FF;
}

double CresInfo::GetBL2(double k,int L){
	const double R=1.0;
	double BL2=1.0;
	double x2;
	if(L>0){
		x2=pow(k*R/HBARC_GEV,2);
		BL2=pow(x2/(1.0+x2),L);
	}
	return BL2;
}

double CresInfo::GetBW(double E,double M0,double Gamma){
	return (2.0*E*E*Gamma/M_PI)/(pow(E*E-M0*M0,2)+E*E*Gamma*Gamma);
}

double CresInfo::GetBW_base(double E,double M0,double Gamma0){ // simple lorentzian used for MC weights
	return (0.5*Gamma0/M_PI)/((E-M0)*(E-M0)+0.25*Gamma0*Gamma0);
	//return (2.0*E*E*Gamma0/M_PI)/(pow(E*E-M0*M0,2)+E*E*Gamma0*Gamma0);
}

double CresInfo::GenerateMass_BW(){
	double r1=randy->ran();
	return ((width/2)*tan(M_PI*(r1-0.5))) + mass;  // mass according to BW distribution
}

double CresInfo::GenerateMass_T0(){
	map<double,double>::iterator it1,it2;
	double E = 0.0, E1 = 0.0, E2 = 0.0;
	double p1 = 0.0, p2 = 0.0, netprob=randy->ran();
	if(!decay)
		E=mass;
	else{
		it1=sfmassmap.lower_bound(netprob);
		if(it1==sfmassmap.end()){
			printf("GenerateMass_T0: it1 at end of map, netprob=%g\n",netprob);
			Print();
			it2=sfmassmap.begin();
			printf(" -----SF massmap----- \n");
			while(it2!=sfmassmap.end()){
				printf("%g   %g\n",it2->first,it2->second);
				it2++;
			}
			exit(1);
		}
		it2=it1;
		--it1;
		p1=it1->first;
		p2=it2->first;
		E1=it1->second;
		E2=it2->second;
		E=((netprob-p1)*E2+(p2-netprob)*E1)/(p2-p1);
	}
	if(E<minmass || E1>E || E2<E){
		printf("something odd in CresInfo::GenerateMass_T0, E=%g, (E1,E2)=(%g,%g), (p1,p2)=(%g,%g), minmass=%g, netprob=%g\n",E,E1, E2,p1,p2,minmass,netprob);
		PrintBranchInfo();
		PrintSpectralFunction();
		for(it1=sfmassmap.begin();it1!=sfmassmap.end();++it1){
			printf("Emap=%g, netdens=%g\n",it1->second,it1->first);
		}
		exit(1);
	}
	return E; //returns a random mass proportional to SF'
}

void CresInfo::NormalizeSF(){
	int n,N=SpectVec.size();
	double netnorm,norm=0.0;
	for(n=0;n<N;n++){
		norm+=SpectVec[n];
	}
	netnorm=0.0;
	for(n=0;n<N;n++){
		SpectVec[n]=SpectVec[n]/norm;
		netnorm+=SpectVec[n];
		sfmassmap.insert(pair<double,double>(netnorm,SpectEVec[n]));
	}
}

void CresInfo::PrintSpectralFunction(){
	int n;
	printf("- Spectral Function for pid=%d -\n",pid);
	printf("__ E ___ Gamma ____  A/A_BW ___ A _ (A/A_BW)*dens(M)/dens(M0) _ A*dens(M)/dens(M0) \n");
	for(n=0;n<int(SpectVec.size());n++){
		printf("%8.4f %8.4f %8.4f\n",
		SpectEVec[n],GammaVec[n],SpectVec[n]);
	}
}

double CresInfo::GetDecayMomentum(double M,double ma,double mb){ // Gives relative momentum
	double pf=0.0;
	if(ma+mb<M)
		pf=sqrt(fabs(pow((M*M-ma*ma-mb*mb),2.0)-(4.0*ma*ma*mb*mb)))/(2.0*M);
	return pf;
}

void CresInfo::DecayGetResInfoPtr(int &nbodies,array<CresInfo *,5> &daughterresinfo){
	char message[100];
	double r,bsum;
	int ibody,ibranch;
	CbranchInfo *bptr;
	bptr=NULL;
	bsum=0.0;
	r=randy->ran();
	ibranch=0;
	do{
		bptr=branchlist[ibranch];
		bsum+=bptr->branching;
		ibranch++;
		if(bsum>1.00000001){
			sprintf(message,"In DecayGetResInfo: bsum too large, = %g\n",bsum);
			CLog::Fatal(message);
		}
	}while(bsum<r);
	nbodies=bptr->resinfo.size();
	for(ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr->resinfo[ibody];
	}
}

void CresInfo::DecayGetResInfoPtr(double mothermass,int &nbodies,array<CresInfo *,5> &daughterresinfo){
	char message[100];
	double r,bsum;
	int ibody,ibranch;
	CbranchInfo *bptr;
	bptr=NULL;
	bsum=0.0;
	r=randy->ran();
	ibranch=0;
	do{
		bptr=branchlist[ibranch];
		bsum+=bptr->branching;
		ibranch++;
		if(bsum>1.00000001){
			sprintf(message,"In DecayGetResInfo: bsum too large, = %g\n",bsum);
			CLog::Fatal(message);
		}
	}while(bsum<r);
	nbodies=bptr->resinfo.size();
	for(ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr->resinfo[ibody];
	}
}

void CresInfo::DecayGetResInfoPtr_minmass(int &nbodies,array<CresInfo *,5> &daughterresinfo){
	nbodies=bptr_minmass->resinfo.size();
	for(int ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr_minmass->resinfo[ibody];
	}
}
