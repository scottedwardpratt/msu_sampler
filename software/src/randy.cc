#include "msu_sampler/randy.h"

using namespace std;
using namespace msu_sampler;

Crandy::Crandy(int iseed){
    seed=iseed;
    mt.seed(seed);
    netprob=0.0;
    threshold=-log(ran());
}

void Crandy::reset(int iseed){
    seed=iseed;
    mt.seed(seed);
}

double Crandy::ran(){
    return ranu(mt);
}

double Crandy::ran_exp(){
    return -log(ran());
}

void Crandy::generate_boltzmann_alt(double mass,double T,FourVector &p){
    double r1,r2,r3,r0,I1,I2,I3,Itot;
    double pmag,E,ctheta,stheta,phi,K;
    do {
        r0=ran();
        I1=mass*mass;
        I2=2.0*mass*T;
        I3=2.0*T*T;
        Itot=I1+I2+I3;
        if(r0<I1/Itot){
            r1=ran();
            K=-T*log(r1);
        }
        else if(r0<(I1+I2)/Itot){
            r1=ran();
            r2=ran();
            K=-T*log(r1*r2);
        }
        else{
            r1=ran();
            r2=ran();
            r3=ran();
            K=-T*log(r1*r2*r3);
        }
        E=K+mass;
        pmag=sqrt(E*E-mass*mass);
        r0=ran();
    } while (r0>pmag/E);
    phi=2.0*M_PI*ran();
    ctheta=1.0-2.0*ran();
    stheta=sqrt(1.0-ctheta*ctheta);
    p[3]=pmag*ctheta;
    p[1]=pmag*stheta*cos(phi);
    p[2]=pmag*stheta*sin(phi);
    p[0]=E;
    //printf("p0=%lf\t",p[0]);
}

void Crandy::generate_boltzmann(double mass,double T,FourVector &p){
    double r1,r2,r3,a,b,c;
    double pmag,ctheta,stheta,phi;
    if(T/mass>0.6){
        do {
            r1=ran();
            r2=ran();
            r3=ran();
            a=-log(r1); b=-log(r2); c=-log(r3);
            pmag=T*(a+b+c);
            p[0]=sqrt(pmag*pmag+mass*mass);
        } while (ran()>exp((pmag-p[0])/T));
        ctheta=(a-b)/(a+b);
        stheta=sqrt(1.0-ctheta*ctheta);
        phi=T*T*pow(a+b,2)/(pmag*pmag);
        phi=2.0*M_PI*phi;
        p[3]=pmag*ctheta;
        p[1]=pmag*stheta*cos(phi);
        p[2]=pmag*stheta*sin(phi);
    }
    else generate_boltzmann_alt(mass,T,p);
}

bool Crandy::test_threshold(double delprob){
    if((delprob+netprob)<threshold){
        //netprob+=delprob;
        return false;
    }
    else
        return true;
}

void Crandy::increase_threshold(){
    threshold-=log(ran());
}

void Crandy::increment_netprob(double delN){
    netprob+=delN;
}
