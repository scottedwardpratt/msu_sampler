#include "msu_sampler/meanfield.h"

CmeanField::CmeanField(){
	//
}

CmeanField_Simple::CmeanField_Simple(CparameterMap *parmap) : CmeanField(){
	char message[200];
	if(parmap->getS("MF_DESCRIPTION","NO_MF_DESCRIPTION")!="NO_MF_DESCRIPTION"){
		sprintf(message,"MF_DESCRIPTION=%s\n",(parmap->getS("MF_DESCRIPTION","NO_MF_DESCRIPTION")).c_str());
		CLog::Info(message);
	}
}

double CmeanField_Simple::GetMass(CresInfo *resinfo,double sigma){
	if(resinfo->baryon!=0)
		return (sigma/0.093)*resinfo->mass;
	else
		return resinfo->mass;
}

double CmeanField::GetMass(CresInfo *resinfo,double sigma){
	// Dummy function
	return (sigma/0.093)*resinfo->mass;
}
