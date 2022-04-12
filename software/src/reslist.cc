#include "msu_sampler/resonances.h"
#include "msu_commonutils/constants.h"

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
		printf("Warning GetResInfoPtr() can't find match for PID=%d\n",pid);
		exit(1);
		return NULL;
	}
}

void CresList::CalcMinMasses(){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos;
	for(rpos=resmap.begin();rpos!=resmap.end();++rpos){
		resinfo=rpos->second;
		resinfo->CalcMinMass();
	}
}

CresList::CresList(CparameterMap* parmap_in){
	parmap=parmap_in;
	CresInfo::NSPECTRAL=parmap->getI("SAMPLER_NSPECTRAL",100);
	CresInfo::SFDIRNAME=parmap->getS("SAMPLER_SFDIRNAME","../local/resinfo/spectralfunctions");
	//RESONANCE_DECAYS=parmap->getB("RESONANCE_DECAYS",true);
	ReadResInfo();
	CalcMinMasses();
	//CalcSpectralFunctions();
	ReadSpectralFunctions();
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
