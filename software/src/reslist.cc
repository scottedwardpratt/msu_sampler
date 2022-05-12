#include "msu_sampler/resonances.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/log.h"

double CresList::MIN_DECAY_WIDTH=0.0001;
char *CresList::message=new char[200];

CresList::CresList(){
}

CresList::~CresList(){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos=resmap.begin();
	while(rpos!=resmap.end()){
		resinfo=rpos->second;
		delete resinfo;
		resmap.erase(rpos);
		rpos=resmap.begin();
	}
	resmap.clear();
	massmap.clear();
}

CresInfo* CresList::GetResInfoPtr(int pid){
	CresInfoMap::iterator rpos;
	rpos=resmap.find(pid);
	if(rpos!=resmap.end())
		return rpos->second;
	else{
		sprintf(message,"GetResInfoPtr() can't find match for PID=%d\n",pid);
		CLog::Fatal(message);
		return NULL;
	}
}

void CresList::CalcMinMasses(){
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		resinfo->CalcMinMass();
	}
}

CresList::CresList(CparameterMap* parmap_in){
	parmap=parmap_in;
	CresInfo::NSPECTRAL=parmap->getI("MSU_SAMPLER_NSPECTRAL",100);
	CresInfo::SFDIRNAME=parmap->getS("MSU_SAMPLER_SFDIRNAME","../local/resinfo/spectralfunctions");
	//RESONANCE_DECAYS=parmap->getB("RESONANCE_DECAYS",true);
	ReadResInfo();
	//CalcSpectralFunctions();
	ReadSpectralFunctions();
	CalcMinMasses();
}

void CresList::CalcSpectralFunctions(){
	CresMassMap::iterator rpos;
	CresInfo *resinfo;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->decay && !resinfo->SFcalculated){
			resinfo->CalcSpectralFunction();
		}
		else
			resinfo->SFcalculated=true;
	}
}

void CresList::PrintMassMaps(){
	CresMassMap::iterator rpos;
	map<double,double>::iterator it;
	CresInfo *resinfo;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->decay){
			it=resinfo->sfmassmap.begin();
			sprintf(message," ----- SF massmap for pid=%d ----- \n",resinfo->pid);
			CLog::Info(message);
			while(it!=resinfo->sfmassmap.end()){
				sprintf(message,"%g   %g\n",it->first,it->second);
				CLog::Info(message);
				it++;
			}
		}
	}
}
