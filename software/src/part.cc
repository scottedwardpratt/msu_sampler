#include "msu_sampler/part.h"

CresList *CpartList::reslist=NULL;
char *Cpart::message=new char[CLog::CHARLENGTH];
char *CpartList::message=new char[CLog::CHARLENGTH];

Cpart::Cpart(){
	msquared=0.0;
	pid=0;
	for(int alpha=0;alpha<4;alpha++){
		r[alpha]=p[alpha]=0.0;
	}
}
Cpart::~Cpart(){
}

void Cpart::Print(){
	snprintf(message,CLog::CHARLENGTH,"________________ PART INFO FOR PART, pid=%d _____________________________\n",pid);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"m^2=%g\\p=(%g,%g,%g,%g)\n",msquared,p[0],p[1],p[2],p[3]);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"r=(%g,%g,%g,%g)\n",r[0],r[1],r[2],r[3]);
	CLog::Info(message);
}

double Cpart::GetMass(){
	if(pid==22)
		return 0.0;
	else
		return sqrt(msquared);
}

void Cpart::SetMsquared(){
	msquared=p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3];
}

void Cpart::Setp0(){
	p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+msquared);
}

void Cpart::Boost(FourVector &u){
	BoostP(u);
	BoostR(u);
}

void Cpart::BoostP(FourVector &u){
	int alpha;
	FourVector pprime;
	Misc::Boost(u,p,pprime);
	for(alpha=0;alpha<4;alpha++)
		p[alpha]=pprime[alpha];
}

void Cpart::BoostR(FourVector &u){
	int alpha;
	FourVector rprime;
	Misc::Boost(u,r,rprime);
	for(alpha=0;alpha<4;alpha++)
		r[alpha]=rprime[alpha];
}

void Cpart::SetEQWeight(Chyper *hyper,Eigen::VectorXd &EQTarget){
	int B=resinfo->baryon,S=resinfo->strange,II=resinfo->q[0]-resinfo->q[1];
	double chipinv=1.0/((hyper->P+hyper->epsilon)*hyper->T0);
	double nhadrons=hyper->nhadrons;
	EQWeightVec.resize(7);
	EQWeightVec[0]=(hyper->chiinv(0,0)*p[0] +hyper->chiinv(0,1)*B
		 +hyper->chiinv(0,2)*II +hyper->chiinv(0,3)*S)*nhadrons;
	
	EQWeightVec[1]=chipinv*p[1]*nhadrons;
	EQWeightVec[2]=chipinv*p[2]*nhadrons;
	EQWeightVec[3]=chipinv*p[3]*nhadrons;
	
	EQWeightVec[4]=(hyper->chiinv(1,0)*p[0] +hyper->chiinv(1,1)*B
		+hyper->chiinv(1,2)*II +hyper->chiinv(1,3)*S)*nhadrons;
	EQWeightVec[5]=(hyper->chiinv(2,0)*p[0] +hyper->chiinv(2,1)*B
		 +hyper->chiinv(2,2)*II +hyper->chiinv(2,3)*S)*hyper->nhadrons;
	EQWeightVec[6]=(hyper->chiinv(3,0)*p[0] +hyper->chiinv(3,1)*B
		 +hyper->chiinv(3,2)*II +hyper->chiinv(3,3)*S)*hyper->nhadrons;
	EQWeight=0.0;
	for(int i=0;i<7;i++)
		EQWeight+=EQWeightVec(i)*EQTarget(i);
}


/////////////////
/////////////////
/////////////////



CpartList::CpartList(CparameterMap *parmap,CresList *reslist_in){
	nparts_blocksize=parmap->getI("MSU_SAMPLER_NPARTS_BLOCKSIZE",2000);
	partvec.resize(nparts_blocksize);
	reslist=reslist_in;
	Reset();
}

CpartList::~CpartList(){
	partvec.clear();
}

Cpart* CpartList::GetPart(){
	if(int(partvec.size())==nparts){
		partvec.resize(partvec.size()+nparts);
	}
	nparts+=1;
	return &partvec[nparts];
}

void CpartList::Clear(){
	partvec.clear();
	nparts=0;
}

void CpartList::Reset(){
	nparts=0;
	for(int alpha=0;alpha<4;alpha++){
		for(int beta=0;beta<4;beta++)
			SE[alpha][beta]=0.0;
	}
}

void CpartList::WriteParts(string filename){
	FILE *fptr=fopen(filename.c_str(),"w");
	if (fptr==NULL) {
		snprintf(message,CLog::CHARLENGTH,"Can't open file to write parts\n");
		CLog::Info(message);
	}
	for(int ipart=0;ipart<nparts;ipart++){
		fprintf(fptr,"%5d %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e\n",
		partvec[ipart].pid,partvec[ipart].msquared,partvec[ipart].p[0],partvec[ipart].p[1],partvec[ipart].p[2],partvec[ipart].p[3],
		partvec[ipart].r[0],partvec[ipart].r[1],partvec[ipart].r[2],partvec[ipart].r[3]);
	}
	fclose(fptr);
}

void CpartList::CountResonances(){
	CresInfo *resinfo;
	for(int ipart=0;ipart<nparts;ipart++){
		resinfo=reslist->GetResInfoPtr(partvec[ipart].pid);
		resinfo->count+=1;
	}
}

long long int CpartList::CountResonances(int pid){
	long long int count=0;
	for(int ipart=0;ipart<nparts;ipart++){
		if(partvec[ipart].pid==pid)
		count+=1;
	}
	return count;
}

void CpartList::IncrementSpectra(int pid,double dp,vector<double> &spectra){
	int ip,ipart,np=spectra.size();
	double pt;
	for(ipart=0;ipart<nparts;ipart++){
		if(abs(partvec[ipart].pid)==pid){
			//pt=sqrt(partvec[ipart].p[1]*partvec[ipart].p[1]+partvec[ipart].p[2]*partvec[ipart].p[2]);
			pt=sqrt(partvec[ipart].p[1]*partvec[ipart].p[1]+partvec[ipart].p[2]*partvec[ipart].p[2]
				+partvec[ipart].p[3]*partvec[ipart].p[3]);
			ip=lrint(floor(pt/dp));
			if(ip<np)
				spectra[ip]+=1;
		}
	}
}

void CpartList::IncrementMassDist(int pid,double dm,vector<double> &massdist){
	int im,ipart,nm=massdist.size();
	double m;
	for(ipart=0;ipart<nparts;ipart++){
		if(abs(partvec[ipart].pid)==pid){
			partvec[ipart].SetMsquared();
			m=sqrt(partvec[ipart].msquared);
			im=lrint(floor(m/dm));
			if(im<nm)
				massdist[im]+=1;
		}
	}
}

double CpartList::SumEnergy(){
	double energy=0.0,ipart;
	for(ipart=0;ipart<nparts;ipart++){
		energy+=partvec[ipart].p[0];
	}
	return energy;
}

double CpartList::SumEnergy(int pid){
	double energy=0.0;
	int ipart;
	for(ipart=0;ipart<nparts;ipart++){
		if(partvec[ipart].pid==pid)
			energy+=partvec[ipart].p[0];
	}
	return energy;
}

void CpartList::AddPart(int pidset,FourVector &pset,FourVector &rset){
	int alpha;
	if(int(partvec.size())==nparts){
		partvec.resize(partvec.size()+nparts_blocksize);
	}
	partvec[nparts].pid=pidset;
	partvec[nparts].resinfo=reslist->GetResInfoPtr(pidset);
	for(alpha=0;alpha<4;alpha++){
		partvec[nparts].p[alpha]=pset[alpha];
		partvec[nparts].r[alpha]=rset[alpha];
	}
	partvec[nparts].SetMsquared();
	nparts+=1;
}

void CpartList::SumSETensor(){
	int ipart,alpha,beta;
	//nparts=partvec.size();
	for(ipart=0;ipart<nparts;ipart++){
		for(alpha=0;alpha<4;alpha++){
			for(beta=0;beta<4;beta++){
				SE[alpha][beta]+=partvec[ipart].p[alpha]*partvec[ipart].p[beta]/partvec[ipart].p[0];
			}
		}
	}
}


void CpartList::SetEQWeight(Chyper *hyper,Eigen::VectorXd &EQTarget){
	//int II,nparts=partvec.size();
	int i;
	for(int ipart=0;ipart<nparts;ipart++){
		for(i=0;i<7;i++)
		partvec[ipart].SetEQWeight(hyper,EQTarget);
	}
}

void CpartList::IncrementEQTot(Eigen::VectorXd &EQtot){
	CresInfo *resinfo;
	//int II,nparts=partvec.size();
	for(int ipart=0;ipart<nparts;ipart++){
		resinfo=partvec[ipart].resinfo;
		int II=resinfo->q[0]-resinfo->q[1];
		EQtot[0]+=partvec[ipart].EQWeight*partvec[ipart].p[0];
		EQtot[1]+=partvec[ipart].EQWeight*partvec[ipart].p[1];
		EQtot[2]+=partvec[ipart].EQWeight*partvec[ipart].p[2];
		EQtot[3]+=partvec[ipart].EQWeight*partvec[ipart].p[3];
		EQtot[4]+=partvec[ipart].EQWeight*partvec[ipart].resinfo->baryon;
		EQtot[5]+=partvec[ipart].EQWeight*II;
		EQtot[6]+=partvec[ipart].EQWeight*resinfo->strange;
	}
}
