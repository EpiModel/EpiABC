/************************************************************************************
    Massive simulations of trait-based filtered communities for ABC estimation
    Franck Jabot, 2nd September 2008.

Last modified on 31th december 2008.
	windows-compatible compilation: g++-3 -O3 trait_model.cpp -lm -lgsl -o trait_model -L/usr/bin -mno-cygwin
*************************************************************************************/

// Libraries
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <string>
#include <math.h>
#include <list>
#include <R.h>
#include <Rmath.h>

using namespace std; // peut-êtr à enlever, je ne sais pas.

//Random number generator
double genrand2(void);
void sgenrand2(unsigned long);
unsigned long genrand2i(void);
void sgenrand2i(unsigned long);

double calculfitness(double *tabfitness,int S,double **trait,double *h,double *sig,double SS,double ntrait){
    double min=SS;
    double max=0.0;
    for (int j=0;j<ntrait;j++){
        for (int i=0;i<S;i++){
            tabfitness[i]*=exp(log(SS)/ntrait)*exp(-(trait[i][j]-h[j])*(trait[i][j]-h[j])/(2*sig[j]*sig[j]));
        	//cerr<<trait[i]<<" "<<h<<" "<<sig<<" "<<SS<<" "<<(tabfitness[i])<<endl;
        }
    }
    for (int i=0;i<S;i++){
        if (tabfitness[i]>max){
            max=tabfitness[i];
        }
        if (tabfitness[i]<min){
            min=tabfitness[i];
        }
    }
    for (int i=0;i<S;i++){
        tabfitness[i]=((1+tabfitness[i])/(1+min)-1);
    }
    
return (((1.0+max)/(1.0+min))-1.0);
}

int draw(double *abondreg){
	double u=genrand2();
	int k=0;
	while (abondreg[k]<u){
		u-=abondreg[k];
		k++;
	}
return k;
}

int draw2(double *abondreg,double *sumabondfit,double *tabfitness){
	double u=sumabondfit[1]*genrand2();
	int k=0;
	while ((abondreg[k]*(1+tabfitness[k]))<u){
		u-=(abondreg[k]*(1+tabfitness[k]));
		k++;
	}
return k;
}

int drawint (int *abondloc,double *sumabondfit,double *tabfitness){
	double u=sumabondfit[0]*genrand2();
	int k=0;
	while ((abondloc[k]*(1+tabfitness[k]))<u){
		u-=(abondloc[k]*(1+tabfitness[k]));
		k++;
	}
return k;
}

int drawint2 (int *abondloc,int J){
	int u=genrand2i()%J;
	int k=0;
	while (abondloc[k]<(1+u)){
		u-=abondloc[k];
		k++;
	}
return k;
}


void initialiser(int *abondloc,double *abondreg,int J,double *sumabondfit,double *tabfitness){
	for (int i=0;i<J;i++){
		int k=draw2(abondreg,sumabondfit,tabfitness);
		abondloc[k]++;
		sumabondfit[0]+=(1+tabfitness[k]);
	}
}

void stepdyn(int *abondloc,double *tabfitness,double *abondreg,double maxSS,double I,int J,double *sumabondfit){
	int k=drawint2(abondloc,J);
	abondloc[k]--;
	//cerr<<k<<" "<<abondloc[k]<<" "<<sumabondfit[0]<<" ";
	sumabondfit[0]-=(1+tabfitness[k]);
	//cerr<<sumabondfit[0]<<endl;
	//int toto;
	//cin >> toto;
	double u=genrand2();
	if (u<(I/(I+J-1))){
		//int kk=draw(abondreg);
		int kk=draw2(abondreg,sumabondfit,tabfitness);
		abondloc[kk]++;
		sumabondfit[0]+=(1+tabfitness[kk]);
	}
	else{
		int kk=drawint(abondloc,sumabondfit,tabfitness);
		abondloc[kk]++;
		sumabondfit[0]+=(1+tabfitness[kk]);
	}
}

void forwarddynamics(int *abondloc,double *tabfitness,double *abondreg,int S,int J,double I,double maxSS,int lsimul){
	for (int i=0;i<S;i++){
		abondloc[i]=0;
	}
	double *sumabondfit;
	sumabondfit=new double[2];
	sumabondfit[0]=0.0;
	sumabondfit[1]=0.0;
	for (int i=0;i<S;i++){
		sumabondfit[1]+=(abondreg[i]*(1+tabfitness[i]));
	}
	initialiser(abondloc,abondreg,J,sumabondfit,tabfitness);
	for (int i=0;i<(J*lsimul);i++){
		stepdyn(abondloc,tabfitness,abondreg,maxSS,I,J,sumabondfit);
	}
	//for (int i=0;i<S;i++){
	//	cerr<<abondloc[i]<<" ";
	//}
	//int toto;
	//cin >> toto;
	
}

void calculstat(double *stat,int *abondloc, int l, double **trait,double ntrait){
    //STATISTIQUES CALCULEES:
    //J 
    //S
    //Shan 
    //Var(Ni) - deleted
    //Mean(trait)
    //Var(trait) - deleted
    //Skewness(trait)
    //Mean(trait)_espece - deleted
    //Var(trait)_espece - deleted
    //Skewness(trait)_espece - deleted
    //Fourth Moment - deleted
    //Fourth Moment espece - deleted
    
    for (int i=0;i<(4+ntrait);i++){
        stat[i]=0;
    }
    
    for (int i=0;i<l;i++){
        stat[0]+=abondloc[i];
        if (abondloc[i]>0){
            stat[1]+=1;
            stat[2]-=abondloc[i]*log(double(0.0+abondloc[i]));
        }
    }
    stat[2]+=stat[0]*log(stat[0]);
    stat[2]/=stat[0];
    /*double meanN=stat[0]/stat[1];
    for (int i=0;i<l;i++){
        if (abondloc[i]>0){
            stat[3]+=(abondloc[i]-meanN)*(abondloc[i]-meanN);
        }
    }
    stat[3]/=stat[1];
    if (stat[1]>1){
        stat[3]*=(stat[1]/(stat[1]-1));
    }*/

 for (int k=0;k<ntrait;k++){

    for (int i=0;i<l;i++){
        if (abondloc[i]>0){
            stat[(3+2*k)]+=abondloc[i]*trait[i][k];
            //stat[(7+8*k)]+=trait[i][k];
        }
    }
    stat[(3+2*k)]/=stat[0];
    //stat[(7+8*k)]/=stat[1];
    
    for (int i=0;i<l;i++){
        if (abondloc[i]>0){
            //stat[(4+2*k)]+=abondloc[i]*(trait[i][k]-stat[(3+2*k)])*(trait[i][k]-stat[(3+2*k)]);
            //stat[(8+8*k)]+=(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)]);
            stat[(4+2*k)]+=abondloc[i]*(trait[i][k]-stat[(3+2*k)])*(trait[i][k]-stat[(3+2*k)])*(trait[i][k]-stat[(3+2*k)]);
            //stat[(9+8*k)]+=(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)]);
	    //stat[(10+8*k)]+=abondloc[i]*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)]);
	    //stat[(11+8*k)]+=(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)]);
        }
    }
    stat[(4+2*k)]/=stat[0];
    //stat[(6+8*k)]/=stat[0];
    //stat[(8+8*k)]/=stat[1];
    //stat[(9+8*k)]/=stat[1];
    //stat[(10+8*k)]/=stat[0];
    //stat[(11+8*k)]/=stat[1];
    
    //rendre les estimateurs non biais�s :
    //if (stat[1]>1){
    //    stat[(8+8*k)]*=(stat[1]/(stat[1]-1));
    //}
    //stat[(4+8*k)]*=(stat[0]/(stat[0]-1));
  }
}



extern "C" {
void trait_model(double *input,double *stat_to_return){

//LECTURE DES FICHIERS D'ENTREE 
    char buffer[256];
    long double seed1p = input[0];
    unsigned long seed2 = 0;
    unsigned long seed1;
    if (seed1p>10000.0){
	    seed2=long(floor(seed1p/1000.0));
	    seed1=long(long(floor(seed1p))%10000);
    }
    else{
	    seed1=long(floor(seed1p));
    }
    sgenrand2(seed1*100000+seed2*124);
    sgenrand2i(1027+seed1*10000+seed2*127);
    int J = input[1];
    double I = input[2];
    double SS = input[3];
    double ntrait = input[4];
    double *h = new double[int(ntrait)];
    for (int i=0;i<ntrait;i++){
        h[i] = input[(5+i)];
    }
	double *sig; sig= new double[int(ntrait)];
    for (int i=0;i<int(ntrait);i++){
        sig[i]=input[(5+int(ntrait)+i)];
    }
    int lsimul;
    lsimul=J;
    
    // fixed species traits for EasyABC sample
    if (ntrait != 1) {
      return;
    }
    int S = 1000;
    double **trait = new double*[S];
    double *abondreg = new double[S];
    double sumabondreg=1.0;
    for (int i=0;i<S;i++){
        abondreg[i] = 0.001;
	trait[i]=new double[1];
	trait[i][0] = (i+1)/10.0;
    }

    double maxSS;
    double *tabfitness;
    tabfitness=new double[S];
    for (int i=0;i<S;i++){
	tabfitness[i]=1;
    }
	
    int *abondloc;
    abondloc= new int[S];
    for (int i=0;i<S;i++){
	    abondloc[i]=0;
    }
    double *stat;
    int nstat=3+ntrait*2;
    stat= new double[nstat];

    int i=0;
    //DRAWING OF PARAMETERS
    I=exp(I);
    SS=exp(SS);
    for (int j=0;j<ntrait;j++){   
	sig[j]=exp(sig[j]);
    }
    maxSS=calculfitness(tabfitness,S,trait,h,sig,SS,ntrait);

    //SIMULATION
    forwarddynamics(abondloc,tabfitness,abondreg,S,J,I,maxSS,lsimul);
    
    //COMPUTATION OF STATISTICS
    calculstat(stat,abondloc,S,trait,ntrait);

    for (int j=1;j<nstat;j++){
	  stat_to_return[(j-1)]=stat[j];
    }
}
}

/**** GENERATEUR DE NOMBRES ALEATOIRES ****/
/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

//#include<stdio.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializing the array with a NONZERO seed */
void sgenrand2(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

/* generating reals */
/* unsigned long */ /* for integer generation */
double genrand2()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand2(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y / (unsigned long)0xffffffff ); /* reals */
    /* return y; */ /* for integer generation */
}


/* initializing the array with a NONZERO seed */
void sgenrand2i(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

unsigned long genrand2i()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand2i(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return y; 
}
