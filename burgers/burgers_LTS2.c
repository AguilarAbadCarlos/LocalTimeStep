#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define L 100.0			// Longitude [m]
//#define D 0.5			// Diameter [m]
//#define c 1.0		// wave speed [m/s]
#define simTime 10.0	// Simulation time [s]
//#define nCells 100		    // Number of nodes
//#define dx (L/(imax-1))     // Spatial increment
#define CFL 1.0		// CFL number
#define g 9.81

#define INITIAL_CONDITION 1
#define MESH_REGULARITY 2

#define v1 -1.			// Initial velocity
#define v2 3.0
#define x1 30.
#define x2 55.
#define PI (4*atan(1))

#define EPS 1E-5
#define INFTY 1E8
#define M_LIMITE 7
#define TOL12 1e-12

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define NormRANu (2.3283063671E-10F)
#define Frec_Max 100
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
int frec[Frec_Max];
extern float Random(void);
extern void ini_ran(int SEMILLA);

int PrintToFile(int *PrintCounter, double *x, double *v, double *exact, double t);
double GlobalTimeStep (double *v, double *dx, double *dtLocal);
double FluxFunction (int i, double v, double vL, double vR, int *FluxCounter);
double SquareWave (double x);
double Gaussian (double x);
void ErrorNorms (double *v, double *dx, double *exact);
int exactBurgers(double t, double *x, double *exact);

int power(int base, unsigned int exp);

//double L, simTime, dx, v1, v2, x1, x2;
int nCells;

int main () {

/*

FALTA INCORPORAR GESTIÓN DE FRONTERAS DE REGIONES DE DIFERENTE M-VALUE.
AHORA REPRODUCE LA INFORMACIÓN PERO NO ESTIMA BIEN LA VELOCIDAD DE PRO-
PAGACIÓN DE LA ONDA CUADRADA. ¡¡¡CHECK!!!

*/

    nCells = 10000;

    srand(time(NULL));
//    ini_ran(rand());

    int semilla = 32207;
    ini_ran(semilla);

//    FILE *inParams;
//    inParams=fopen("parametrosBurgers.txt","r");
////    f'{L:.1f} {simTime:.1f} {nCells:d} {dx:f} {v1:.1f} {v2:.1f} {x1:.1f} {x2:.1f}'
//    fscanf(inParams, "%lf %lf %i %lf %lf %lf %lf %lf", &L, &simTime, &nCells, &dx, &v1, &v2, &x1, &x2);
//    fclose(inParams);

# if MESH_REGULARITY == 1
    int irregularity = 5;
    nCells = 40 * (1<<irregularity) + 10 * irregularity + 60;
# endif // MESH_REGULARITY

    int i,iter=0,k=0,printFreq=1,j=0;
	double *v,vLeft,vRight;		// Velocity arrays (conserved variables)
	double *vL,*dv;      // Interpolation points
	double *Q, *x, *exact, *dx;			// Discharge, distance
	double t=0.0, dt, dtGlobal, *dtLocal, tLocal=0., *tFrozen;                // Time
	double *flux;
	double eps=1e-5;                       // Area
	int PrintCounter=0, FluxCounter=0, IterCounter=0;             // Print counter
	int *mLocal, mMax=0, *mPow2Local, mAux=0, mDif=0, mPow2Dif=0, domainLength=0;
	int p=0, pTotal=0;
//	bool allowsFrozenFlux[nCells];
	double log2 = log(2.);

    printf("\n   Linear convection");
	printf("\n   --------------------------------------------------");


// Reading simulation parameters

	printf("\n>> Simulation setup loaded");
	printf("\n   Simulation time: %.5lf s",simTime);
	printf("\n   CFL: %.2lf",CFL);
	printf("\n   Number of cells: %d",nCells);
//	printf("\n   dx: %.2lf m",dx);
	printf("\n   Pipe length: %.2lf m\n",L);


	// Variable initialization
//	for (i=0;i<nCells;i++){
//        v[i] = 0.0;
//        vL[i] = 0.0;
//        dv[i] = 0.0;
//        Q[i] = 0.0;
//        x[i] = 0.0;
//        flux[i] = 0.0;
//        tFrozen[i] = 0.0;
//	}
	v = calloc(nCells,sizeof(double));
	vL = calloc(nCells,sizeof(double));
	dv = calloc(nCells,sizeof(double));
	Q = calloc(nCells,sizeof(double));
	x = calloc(nCells,sizeof(double));
	dx = calloc(nCells,sizeof(double));
	flux = calloc(nCells,sizeof(double));
	tFrozen = calloc(nCells,sizeof(double));
	exact = calloc(nCells,sizeof(double));
	dtLocal = calloc(nCells,sizeof(double));
	mLocal = calloc(nCells,sizeof(int));
	mPow2Local = calloc(nCells,sizeof(int));

//	Irregular mesh

#if MESH_REGULARITY == 1

	dx[0]=1.0;
	x[0]=dx[0]/2;
	for (i=1; i<25; i++){
        dx[i] = 1.0;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (k=1; k<irregularity; k++){
        for (j=0;j<5;j++){
            dx[i] = 1.0/(1<<k);
            x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
            i++;
        }
	}
	for (j=0; j<20+40*(1<<k); j++){
        dx[i] = 1.0/(1<<k);
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
        i++;
	}
	for (k=irregularity-1; k>0; k--){
        for (j=0;j<5;j++){
            dx[i] = 1.0/(1<<k);
            x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
            i++;
        }
	}
	for (; i<nCells; i++){
        dx[i] = 1.0;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}

#elif MESH_REGULARITY == 2

    int i_celda_fina = nCells/2;
    double toLenght = L/nCells;
    dx[0] = 1.0*toLenght;
    x[0] = dx[0]/2.;
    for (i=1; i<i_celda_fina; i++){
        dx[i] = 1.*toLenght;
        x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i]);
    }
    dx[i] = 0.001*toLenght;
    x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i]);
    i++;
    for (; i<nCells; i++){
        dx[i] = 1.*toLenght;
        x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i]);
    }

#elif MESH_REGULARITY == 3

    int l = nCells/100;
    k = (nCells-l)/2-1;
    int n_factor = 128;
    double f = pow(n_factor,1./(k));

    double d = L/2./(n_factor+f*(1-pow(f,k))/(1-f)+(l+1)/2);

    dx[0] = n_factor*d;
    x[0] = dx[0]/2.;

    for (i=1; i<=k; i++)
    {
        dx[i] = dx[i-1]/f;
        x[i] = x[i-1] + 0.5*(dx[i]+dx[i-1]);
    }

    for (   ; i<=k+1+l; i++)
    {
        dx [i] = d;
        x[i] = x[i-1] + 0.5*(dx[i]+dx[i-1]);
    }

    for (   ; i<nCells-1; i++)
    {
        dx[i] = dx[i-1]*f;
        x[i] = x[i-1] + 0.5*(dx[i]+dx[i-1]);
    }

    dx[i] = n_factor*d;
    x[i] = x[i-1] + 0.5*(dx[i]+dx[i-1]);

#else // REGULAR MESH

//    dx[0] = dx;
    x[0] = dx/2.;
    for (i=1; i<nCells; i++){
//        dx[i] = dx;
        x[i] = x[i-1] + dx;
    }

#endif // MESH_REGULARITY

// INITIAL CONDITION SETUP

#if INITIAL_CONDITION == 1

        for (i=0; i<nCells; i++)
        {
            vL[i] = SquareWave(x[i]);
            v[i] = SquareWave(x[i]);
        }

#elif INITIAL_CONDITION == 2

        for (i=0; x[i]<x1; i++) v[i] = v2;
        for (; x[i]<x2; i++)    v[i] = v2 - 2* v2 *(i*dx-x1) / (x2-x1);
        for (; i<nCells; i++)   v[i] = -v2;

#else

        for (i=0; i<nCells; i++) v[i] = Gaussian(x[i]);

#endif // INITIAL_CONDITION

//  Testing mesh initialization
//    for (i=0;i<nCells;i++){
//        printf("%d\t%.6f\t%.6f\n",i,dx[i],x[i]);
//    }

//  Testing velocity initialization
//    for (i=0;i<nCells;i++){
//        printf("%d\t%.3f\n",i,v[i]);
//    }
//  getchar();

    exactBurgers(t,x,exact);
    PrintToFile(&PrintCounter, x, v, exact, t);

    while (t < simTime){

//      Initiating Global Time Step

//        printf("Iniciando Global Time Step...\n");

        dt = GlobalTimeStep(v,dx,dtLocal);

        if (t+dt > simTime){
            dt = simTime - t;
            if ((simTime-t)<TOL12) break;
        }


//        if (IterCounter % 100 == 0) printf("t=%lfdt=%lf\n",t,dt);
//        getchar();

//        mLocal[0] = (int) MIN((log(dtLocal[0]/dt)/log2+eps),M_LIMITE);
//        mMax = mLocal[0];
//        mPow2Local[0] = (int) pow(2,mLocal[0]);
//        printf("mLocal[%d]=%d\n",0,mLocal[0]);

        for (i=0; i<nCells; i++){
            mLocal[i] = (int) MIN((log(dtLocal[i]/dt)/log2+eps),M_LIMITE);

//        printf("mLocal[%d]=%d\n",i,mLocal[i]);
            mPow2Local[i] = power(2,mLocal[i]);
            mMax = MAX(mMax,mLocal[i]);
        }
        pTotal = power(2,mMax);
//        printf("pTotal = %d\n",pTotal);
//        getchar();

        if (t+pTotal*dt > simTime)
        {
            dt = (simTime-t)/(double)pTotal;

        }

        /// Actualización de m-values
        for (i=1; i<nCells; i++){

            mAux = mLocal[i-1];

            if (mAux != mLocal[i]){

                mDif = mLocal[i] - mAux;

                if (mDif>0)
                {
                    for (j=0; mLocal[i+j] == mLocal[i] && j<nCells-1-i; j++) domainLength++;

                    mPow2Dif = MIN( power(2,mDif-1)+3,domainLength);
                    for (j=0; j<=3; j++){
                        mLocal[i+j] = mAux;
                        mPow2Local[i+j] = power(2,mLocal[i+j]);
                    }
                    for (; j<=mPow2Dif; j++) {
                        mLocal[i+j] = (int) (mAux + 1.*log(j-2) / log2 + eps);
                        mPow2Local[i+j] = power(2,mLocal[i+j]);
                    }
                        i = i + domainLength - 1;
                        domainLength = 0;
                }
                else
                {
                    mDif = -mDif;
                    mAux = mLocal[i];
                    for (j=1; mLocal[i-j] == mLocal[i-1] && i-j>=0; j++) domainLength++;

                    mPow2Dif = MIN( power(2,mDif-1)+3,domainLength);

                    for (j=0; j<=3; j++){
                        mLocal[i-j-1] = mAux;
                        mPow2Local[i-j-1] = power(2,mLocal[i-j-1]);
                    }
                    for (; j<=mPow2Dif; j++) {
                        mLocal[i-j-1] = (int) (mAux + 1.*log(j-2) / log2 + eps);
                        mPow2Local[i-j-1] = power(2,mLocal[i-j-1]);
                    }
                        //i = i + domainLength - 1;
                        domainLength = 0;
                }



            }

        }

//        for (i=0; i<nCells; i++){
//            printf("i=%d m[i]=%d 2^m[i]=%d v[i]=%lf vL[i]=%lf\n",i,mLocal[i],mPow2Local[i],v[i],vL[i]);
//        }
//
//        getchar();

        for (i=0; i<nCells; i++) vL[i] = v[i];

        for (p=0; p<pTotal; p++){

            for (i=0; i<nCells; i++){

                if (p % mPow2Local[i] == 0){

                    if (mLocal[i] >= mLocal[i-1]) {
                        vLeft = v[i-1];
                    } else {
                        vLeft = 0.5*(v[i-1]+vL[i-1]);
//                        dv[i] = dt * mPow2Local[i] / dx[i] * FluxFunction(i, v[i], 0.5*(v[i-1]+vL[i-1]), &FluxCounter);
                    }

                    if (mLocal[i] >= mLocal[i+1]) {
                        vRight = v[i+1];
                    } else {
                        vRight = 0.5*(v[i+1]+vL[i+1]);
//                        dv[i] = dt * mPow2Local[i] / dx[i] * FluxFunction(i, v[i], 0.5*(v[i-1]+vL[i-1]), &FluxCounter);
                    }

                    dv[i] = dt * mPow2Local[i] / dx[i] * FluxFunction(i, v[i], vLeft, vRight, &FluxCounter);




                }

            }

            for (i=0; i<nCells; i++){

                if (p % mPow2Local[i] == 0){

//                if (fabs(dv[i])>1e-12) printf("t=%f\tOld => i=%d\tv[i]=%lf\tvL[i]=%lf\tdv[i]=%lf\n",t,i,v[i],vL[i],dv[i]);
                    vL[i] = v[i];
                    v[i] = v[i] - dv[i];
                    dv[i] = 0.0;

                }

//                printf("New => i=%d\tv[i]=%lf\tvL[i]=%lf\tdv[i]=%lf\n",i,v[i],vL[i],dv[i]);

            }
//            getchar();

            t = t + dt;

        }

        IterCounter++;

//            getchar();

        }


//    t = t + dtGlobal;



    exactBurgers(t,x,exact);
    PrintToFile(&PrintCounter, x, v, exact, t);
    printf("Flux function evaluations:%d, iterations=%d, COUNTER=%d",FluxCounter,IterCounter,0);
    ErrorNorms(v,dx,exact);

    free(v);
	free(vL);
	free(dv);
	free(Q);
	free(x);
	free(dx);
	free(flux);
	free(tFrozen);
	free(exact);
	free(dtLocal);
	free(mLocal);
	free(mPow2Local);

    return 0;
}

int PrintToFile(int *PrintCounter, double *x, double *v, double *exact, double t){

    int i;
    char filename[80];

    FILE *output;

//	Write the numbering in the file name with the corresponding extension
    sprintf(filename,"OutputFiles/burgersLocalTimeStep2_%d.dat",*PrintCounter);
    output=fopen(filename,"w");

    fprintf(output,"#t=%lf\n",t);
    for(i=0;i<nCells;i++){
       fprintf(output,"%lf %lf %lf\n",x[i],v[i],exact[i]);
    }

    fclose(output);
    *PrintCounter=*PrintCounter+1;

    return 0;

}

double GlobalTimeStep(double *v, double *dx, double *dtLocal){

    int i = 0;
    double dt = 0.0;

    dtLocal[0] = fabs(dx[0]/v[0]);
    dt = dtLocal[0];

    for (i = 1; i < nCells; i++){
        dtLocal[i] = fabs(dx[i]/(v[i]+TOL12));
        dt = MIN(dt,dtLocal[i]);
    }

    return dt;
}

double FluxFunction(int i, double v, double vL, double vR, int *FluxCounter){

    *FluxCounter = *FluxCounter + 1;

    if (i==0 || i==nCells-1){
        return 0.;
    }

    if (v*vL<0. && vL<0.){
//        return v[i-1]*(v[i]-0.5*(v[i]+v[i-1]));
        return v*(0.5*(v+vL)-vL);
    }

    if (vR*v<0. && v<0.){
//        return v[i+1]*(0.5*(v[i+1]+v[i])-v[i]);
        return v*(vR-0.5*(vR+v));
    }

    double cPlus,cMinus;
    cPlus = MAX(0.5*(v+vL),0.);
    cMinus = MIN(0.5*(v+vR),0.);

//    printf("i=%d, v=%lf, vL=%lf, vR=%lf, cPlus=%lf, cMinus=%lf\n",i,v,vL,vR,cPlus,cMinus);

    return cPlus*(v-vL)+cMinus*(vR-v);

}

double SquareWave (double x){

    if (x < x1 || x > x2) return v1;

    return v2;
}

double Gaussian (double x){

    double sigma = 5.0;
    double v0 = 1.;

    return v0*exp(-(x-40.)*(x-40.)/2./sigma/sigma);

}

void ErrorNorms (double *v, double *dx, double *exact){ // CHANGE FOR EXACT SOLUTION GIVEN BELOW

    double L1=0., L2=0., Linf=0.;
    int i;


    for (i=0; i<nCells; i++){
        L1 = L1 + fabs(v[i]-exact[i])*dx[i];
        L2 = L2 + pow(fabs(v[i]-exact[i]),2)*dx[i];
        Linf = MAX(Linf, fabs(v[i]-exact[i]));
    }

    printf("\n\nError norms:\n");
    printf("\tL1=%lf",L1);
    printf("\tL2=%lf",sqrt(L2));
    printf("\tLinf=%lf",Linf);

}

int exactBurgers(double t, double *x, double *exact){      // Exact solution for the Burgers equation with square initial condition
    int i;

#if INITIAL_CONDITION == 1
//    int i1=(int)round((x1+v1*t)/dx);
//    int i2=(int)round((x1+v2*t)/dx);
//    int i3=(int)round((x2+(v2+v1)/2.*t)/dx);

    for (i=0;x[i]<(x1+v1*t);i++){
        exact[i]=v1;
    }
    for (;x[i]<(x1+v2*t);i++){
        exact[i]=(x[i]-x1)/t;
    }
    for (;x[i]<(x2+(v2+v1)/2.*t);i++){
        exact[i]=v2;
    }
    for (;i<nCells;i++){
        exact[i]=v1;
    }
    return 0;

#elif INITIAL_CONDITION == 2

    int i1=(int) MIN(round((x1+v2*t)/dx),0.5*(x1+x2)/dx);
    int i2=(int) MAX(round((x2-v2*t)/dx),0.5*(x1+x2)/dx);

    for (i=0; i<i1; i++)        exact[i] = v2;
    for (i=i1; i<i2; i++)       exact[i] = v2 - 2.* v2 * (i*dx-x1-v2*t)/(x2-x1+2.*v2*t);
    for (i=i2; i<nCells; i++)   exact[i] = -v2;

#else

// UNDEFINED

#endif // INITIAL_CONDITION
}



float Random(void)
{
float r;
ig1=ind_ran-24;
ig2=ind_ran-55;
ig3=ind_ran-61;
irr[ind_ran]=irr[ig1]+irr[ig2];
ir1=(irr[ind_ran]^irr[ig3]);
ind_ran++;
r=ir1*NormRANu;
//printf("r=%f\n",r);
return r;
}

void ini_ran(int SEMILLA)
{
printf("%d",SEMILLA);
int INI,FACTOR,SUM,i;
srand(SEMILLA);
INI=SEMILLA;
FACTOR=67397;
SUM=7364893;
for(i=0;i<256;i++)
{
INI=(INI*FACTOR+SUM);
irr[i]=INI;
}
ind_ran=ig1=ig2=ig3=0;
}

int power(int base, unsigned int exp){

    if (exp == 0)
        return 1;
    int temp = power(base, exp/2);
    if (exp%2 == 0)
        return temp*temp;
    else
        return base*temp*temp;

}







