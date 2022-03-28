#include "msu_sampler/master.h"
#include <cstring>
using namespace std;

int main(){
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters_sampler.txt");
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	long long int deltacount=0;
	ms.randy->reset(time(NULL));
	
	ms.ReadHyper_BEST_Binary3D();
	
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS_TOT",10);
	for(int ievent=0;ievent<nevents;ievent++){
		nparts+=ms.MakeEvent();
		deltacount+=ms.partlist->CountResonances(2224)+ms.partlist->CountResonances(-2224)
			+ms.partlist->CountResonances(2214)+ms.partlist->CountResonances(-2214)
				+ms.partlist->CountResonances(2224)+ms.partlist->CountResonances(-2224)
					+ms.partlist->CountResonances(1114)+ms.partlist->CountResonances(-1114);
		if((10*(ievent+1))%nevents==0)
			printf("ievent=%lld nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
	}
	printf("Ndelta=%g\n",double(deltacount)/double(nevents));
	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
