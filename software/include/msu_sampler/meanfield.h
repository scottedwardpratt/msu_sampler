#ifndef __SAMPLER_MEANFIELD_H__
#define __SAMPLER_MEANFIELD_H__

#include "msu_sampler/classdefs.h"
#include "msu_eos/resonances.h"
#include "msu_commonutils/parametermap.h"

class CmeanField{
public:
	bool crap;
	CmeanField();
	virtual double GetMass(CresInfo* resinfo,double sigma);
};

class CmeanField_Simple : public CmeanField {
public:
	CmeanField_Simple(CparameterMap *parmap);
	double GetMass(CresInfo* resinfo,double sigma);
};

#endif
