#include "msu_sampler/resonances.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/log.h"
using namespace std;

void CresList::ReadResInfo(){
	//Cmerge *merge;
	int motherpid,pid;
	bool nophotons;
	double bmax,bsum;
	int ires,ichannelres,ichannel,ibody,nbodies,NResonances,LDecay=1;
	int ires1,ires2,iresflip;
	//int netq,netb,nets;
	string name, filename;
	CresInfo *resinfo=NULL,*aresinfo=NULL,*temp=NULL;
	CdecayInfo *decayinfo, *adecayinfo;
	CbranchInfo *bptr=NULL,*firstbptr=NULL;
	FILE *resinfofile;
	char cname[200],dummy[200];
	CdecayInfoMap decaymap;
	double gisospin;
	int dummy_int;
	CresInfoMap::iterator iter;
	CdecayInfoMap::iterator diter;
	Cmerge *merge;
	filename=parmap->getS("RESONANCES_INFO_FILE",string("../software/resinfo/pdg-SMASH.dat"));
	sprintf(message,"will read resonance info from %s\n",filename.c_str());
	resinfofile=fopen(filename.c_str(),"r");
	CLog::Info("will read resonance info from "+filename+"\n");
	if (resinfofile==NULL) {
		sprintf(message,"Can't open resinfofile\n");
		CLog::Fatal(message);
	}

	ires=0;
	NResonances=0;
	while(fscanf(resinfofile," %d",&pid)!=EOF && NResonances<20000){
		resinfo=new CresInfo();
		NResonances+=1;
		resinfo->pid=pid;

		//main reading
		fscanf(resinfofile, " %s %lf %lf %lf %d %d %d %d %lf %d %d", cname,&resinfo->mass,&resinfo->width,&resinfo->degen,&resinfo->baryon,&resinfo->strange,&resinfo->charm,&resinfo->bottom,&gisospin,&resinfo->charge,&resinfo->nchannels);
		resinfo->minmass=resinfo->mass;
		if(resinfo->width<MIN_DECAY_WIDTH){
			resinfo->decay=false;
		}
		else{
			resinfo->decay=true;
			decayinfo=new CdecayInfo();
		}
		cname[int(strlen(cname))-1]='\0';
		resinfo->name=cname;

		//decay reading
		//reads into map values: will access for decays when done creating resonances
		for (ichannel=0; ichannel<resinfo->nchannels; ichannel++){
			if(resinfo->decay){
				fscanf(resinfofile, " %d %d %lf %d %d %d %d %d %d", 
					&dummy_int,&decayinfo->Nparts[ichannel],&decayinfo->branchratio[ichannel],&decayinfo->products[ichannel][0],
					&decayinfo->products[ichannel][1],&decayinfo->products[ichannel][2],&decayinfo->products[ichannel][3],
					&decayinfo->products[ichannel][4],&decayinfo->d_L[ichannel]);
			}
			else{
				fscanf(resinfofile,"%d",&dummy_int);
				fgets(dummy,200,resinfofile);
			}
		}
		if(resinfo->decay){
			bsum=0.0;
			for(ichannel=0;ichannel<resinfo->nchannels;ichannel++)
				bsum+=decayinfo->branchratio[ichannel];
			for(ichannel=0;ichannel<resinfo->nchannels;ichannel++)
				decayinfo->branchratio[ichannel]=decayinfo->branchratio[ichannel]/bsum;
		}

		if(resinfo->pid!=22){ //copied from old pid
			resinfo->ires=ires;
			ires+=1;
		}

		resinfo->branchlist.clear();

		resmap.insert(CresInfoPair(resinfo->pid,resinfo));
		massmap.insert(CresMassPair(resinfo->mass,resinfo));
		if(resinfo->decay)
			decaymap.insert(CdecayInfoPair(resinfo->pid,decayinfo));

		//antiparticle creation

		if(resinfo->baryon!=0){
			aresinfo=new CresInfo();
			adecayinfo=new CdecayInfo();
			NResonances+=1;

			aresinfo->pid=-resinfo->pid;
			aresinfo->mass=resinfo->mass;
			aresinfo->minmass=resinfo->minmass;
			aresinfo->width=resinfo->width;
			aresinfo->degen=resinfo->degen;
			aresinfo->baryon=-resinfo->baryon;
			aresinfo->strange=-resinfo->strange;
			aresinfo->charm=-resinfo->charm;
			aresinfo->bottom=-resinfo->bottom;
			aresinfo->charge=-resinfo->charge;
			aresinfo->decay=resinfo->decay;
			aresinfo->nchannels=resinfo->nchannels;
			cname[int(strlen(cname))-1]='\0';
			string s(cname);
			aresinfo->name="Anti-"+s;
			aresinfo->ires=ires;
			ires+=1;

			aresinfo->branchlist.clear();
			for(ichannel=0; ichannel<resinfo->nchannels; ichannel++) { //reads into map values: will access for decays when done creating resonances
				for(int i=0; i<decayinfo->Nparts[ichannel]; i++) {
					pid=decayinfo->products[ichannel][i];
					if(pid!=0){
						temp=GetResInfoPtr(pid);
						if(temp->baryon==0 && temp->charge==0 && temp->strange==0){
							adecayinfo->products[ichannel][i]=decayinfo->products[ichannel][i];
						}
						else{
							adecayinfo->products[ichannel][i]=-decayinfo->products[ichannel][i];
						}
					}
					else adecayinfo->products[ichannel][i]=0;
				}
				adecayinfo->Nparts[ichannel]=decayinfo->Nparts[ichannel];
				adecayinfo->branchratio[ichannel]=decayinfo->branchratio[ichannel];
				adecayinfo->d_L[ichannel]=decayinfo->d_L[ichannel];
			}
			resmap.insert(CresInfoPair(aresinfo->pid,aresinfo));
			massmap.insert(CresMassPair(aresinfo->mass,aresinfo));
			if(aresinfo->decay)
				decaymap.insert(CdecayInfoPair(aresinfo->pid,adecayinfo));
		}
	}
	fclose(resinfofile);
	CLog::Info("NResonances:"+to_string(NResonances)+"\n");
	//------------------------------------------
	MergeArray=new Cmerge **[NResonances];
	//SigmaMaxArray=new double *[NResonances];
	for(ires=0;ires<NResonances;ires++){
		MergeArray[ires]=new Cmerge *[NResonances];
		//SigmaMaxArray[ires]=new double[NResonances];
		for(ichannelres=0;ichannelres<NResonances;ichannelres++){
			MergeArray[ires][ichannelres]=NULL;
			//SigmaMaxArray[ires][ichannelres]=0.0;
		}
	}

	//now, use the stored decay information to create branchlists
	for(iter=resmap.begin();iter!=resmap.end();++iter){
		resinfo=iter->second;
		if(resinfo->decay){
			motherpid=iter->first;
			decayinfo=(decaymap.find(motherpid))->second; //decaymap[motherpid];
			bmax=0.0;
			for(ichannel=0; ichannel<resinfo->nchannels; ichannel++) {
				resinfo->Print();
				nbodies=decayinfo->Nparts[ichannel];
				bptr=new CbranchInfo();
				bptr->resinfo.clear();
				resinfo->branchlist.push_back(bptr);
				bptr->branching=decayinfo->branchratio[ichannel];

				//netq=-resinfo->charge;
				//netb=-resinfo->baryon;
				//nets=-resinfo->strange;

				nophotons=true;
				for(ibody=0; ibody<nbodies; ibody++) {
					pid=decayinfo->products[ichannel][ibody];
					if(pid==22)
						nophotons=false;
					bptr->resinfo.push_back(GetResInfoPtr(pid));
					//netq+=bptr->resinfo[ibody]->charge;
					//netb+=bptr->resinfo[ibody]->baryon;
					//nets+=bptr->resinfo[ibody]->strange;
				}
				bptr->L=decayinfo->d_L[ichannel];

			//total charge and baryon number should be conserved, and shouldn't be larger than single strangeness
				/*
				if(netq!=0 || netb!=0 || abs(nets)>1){
					sprintf(message,"Charge conservation failure while reading decay info,\nnetq=%d, netb=%d, nets=%d\n",netq,netb,nets);
					CLog::Info(message);
					resinfo->Print();
					sprintf(message,"nchannels=%d, ichannel=%d\n",resinfo->nchannels,ichannel);
					CLog::Info(message);
					sprintf(message,"DAUGHTERS:\n");
					CLog::Info(message);
					for(ibody=0;ibody<nbodies;ibody++)
						bptr->resinfo[ibody]->Print();
					exit(1);
				}*/
				if(nophotons){
					ires1=bptr->resinfo[0]->ires;
					ires2=bptr->resinfo[1]->ires;
					bptr->resinfo[0]->Print();
					if(ires1>ires2){
						iresflip=ires1; ires1=ires2; ires2=iresflip;
					}
					merge=MergeArray[ires1][ires2];
					if(merge==NULL){
						MergeArray[ires1][ires2]=new Cmerge(resinfo,bptr->branching, LDecay);
					}
					else{
						while(merge->next!=NULL){
							merge=merge->next;
						}
						merge->next=new Cmerge(resinfo,bptr->branching, LDecay);
					}
				}

				// switch places to make sure first branch has largest
				if(bptr->branching>bmax){
					bmax=bptr->branching;
					if(ichannel>0){
						firstbptr=resinfo->branchlist[0];
						resinfo->branchlist[0]=bptr;
						resinfo->branchlist[ichannel]=firstbptr;
					}
				}
			}  //out of channel loops
		}
	}
	for(iter=resmap.begin();iter!=resmap.end();++iter){
		resinfo=iter->second;
		motherpid=iter->first;
		decayinfo=decaymap[motherpid];
		delete decayinfo;
	}
	decaymap.clear();
}

void CresList::ReadSpectralFunctions(){
	CresMassMap::iterator rpos;
	CresInfo *resinfo;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->decay && !resinfo->SFcalculated){
			resinfo->ReadSpectralFunction();
			resinfo->NormalizeSF();
		}
		else
			resinfo->SFcalculated=true;
	}
}


