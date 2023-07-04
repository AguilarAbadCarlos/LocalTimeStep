#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define L 100.0			// Longitude [m]
#define D 0.5			// Diameter [m]
#define c 1.0		// wave speed [m/s]
#define simTime 10.0	// Simulation time [s]
#define nCells 10000		    // Number of nodes
//#define dx (L/(imax-1))     // Spatial increment
#define CFL 1.0		// CFL number
#define g 9.81

#define mLimite 8

#define INITIAL_CONDITION 2
#define MESH_REGULARITY 2

#define v0 5.0			// Initial velocity
#define PI (4*atan(1))

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

int PrintToFile (int *PrintCounter, double x[], double v[], double t);
double GlobalTimeStep (double dx[], double t, double *dtLocal);
double FluxFunction (int i, double v, double vL, long *FluxCounter);
double SquareWave (double x);
double Gaussian (double x);
void ErrorNorms (double v[], double dx[], double x[]);
int power(int base, unsigned int exp);

int main () {

    srand(time(NULL));
//    ini_ran(rand());

    int semilla = 32207;
    ini_ran(semilla);

    int i,iter=0,k=0,printFreq=1,j=0,l=0;
	double v[nCells];		// Velocity arrays (conserved variables)
	double vL[nCells],dv[nCells];      // Interpolation points
	double Q[nCells], x[nCells], dx[nCells];			// Discharge, distance
	double t=0.0, dt, dtGlobal, dtLocal[nCells], tLocal=0., tFrozen[nCells];                // Time
	double flux[nCells];
	double eps=1e-5;                       // Area
	int PrintCounter=0, IterCounter=0;             // Print counter
	long FluxCounter=0;
	int mLocal[nCells], mMax=0, mPow2Local[nCells], mAux=0, mDif=0, mPow2Dif=0, domainLength=0;
	int p=0, pTotal=0;
//	bool allowsFrozenFlux[nCells];
	double log2 = log(2.);

    printf("\n   Linear convection");
	printf("\n   --------------------------------------------------\n");

//    printf("2^5=%d\n",1<<5);
//    printf("2^7=%d\n",1<<7);
//    printf("2^0=%d\n",1<<0);
//    getchar();

// Reading simulation parameters

	printf("\n>> Simulation setup loaded");
	printf("\n   Simulation time: %.5lf s",simTime);
	printf("\n   CFL: %.2lf",CFL);
	printf("\n   Number of cells: %d",nCells);
	printf("\n   dx: %.2lf m",dx);
	printf("\n   Pipe length: %.2lf m\n",L);


	// Variable initialization
	for (i=0;i<nCells;i++){
        v[i] = 0.0;
        vL[i] = 0.0;
        dv[i] = 0.0;
        Q[i] = 0.0;
        x[i] = 0.0;
        flux[i] = 0.0;
        tFrozen[i] = 0.0;
	}

//	Irregular mesh

#if MESH_REGULARITY == 1 // Irregular

	dx[0]=1.0;
	x[0]=dx[0]/2;
	for (i=1;i<35;i++){
        dx[i] = 1.0;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (i=35;i<40;i++){
        dx[i] = 0.5;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (i=40;i<60;i++){
        dx[i] = 0.25;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (i=60;i<65;i++){
        dx[i] = 0.5;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}
	for (i=65;i<nCells;i++){
        dx[i] = 1.0;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}

#elif MESH_REGULARITY == 2 // Random

    double to100m = L/nCells;
    dx[0] = 2.*to100m*Random();
    x[0] = dx[0]/2.;
    for (i=1; i<nCells; i++){
        dx[i] = 2.*to100m*Random();
        x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i]);
    }

#else // Regular

    dx[0] = L/nCells;
    x[0] = dx[0]/2.;
    for (i=1;i<nCells;i++){
        dx[i] = dx[0];
        x[i] = x[i-1] + dx[i];
    }

#endif // MESH_REGULARITY

// INITIAL CONDITION SETUP

#if INITIAL_CONDITION == 1

        for (i=0; i<nCells; i++) v[i] = SquareWave(x[i]);

#else

        for (i=0; i<nCells; i++) v[i] = Gaussian(x[i]);

#endif // INITIAL_CONDITION

//  Testing mesh initialization
//    for (i=0;i<nCells;i++){
//        printf("%d\t%.3f\t%.3f\n",i,dx[i],x[i]);
//    }

//  Testing velocity initialization
//    for (i=0;i<nCells;i++){
//        printf("%d\t%.3f\n",i,v[i]);
//    }
//  getchar();

    PrintToFile(&PrintCounter, x, v, t);

    while (t < simTime){

//      Initiating Global Time Step

//        printf("Iniciando Global Time Step...\n");

        dt = GlobalTimeStep(dx,t,dtLocal);

        if (t+dt > simTime){
            dt = simTime - t;
        }


//        printf("dt=%lf\n",dt);
//        getchar();

//      Asignación de índices m[i]=int log2(dt[i]/dt), eps para evitar redondeos, mLimite para evitar inestabilidad por demasiada diferencia

        mLocal[0] = (int) MIN((log(dtLocal[0]/dt)/log2+eps),mLimite);
        mMax = mLocal[0];
        mPow2Local[0] = (int) power(2,mLocal[0]);
//        printf("mLocal[%d]=%d\n",0,mLocal[0]);
        for (i=1; i<nCells; i++){
            mLocal[i] = (int) MIN((log(dtLocal[i]/dt)/log2+eps),mLimite);

//        printf("mLocal[%d]=%d expected=%lf\n",i,mLocal[i],log(dtLocal[i]/dt)/log2);
            mPow2Local[i] = (int) power(2,mLocal[i]);
            mMax = MAX(mMax,mLocal[i]);
        }
//        getchar();
        pTotal = (int) power(2,mMax);
//        printf("pTotal = %d\n",pTotal);
        if (dt == simTime-t) pTotal = 1;

//        printf("dtGlobal=%lf\n",dtGlobal);
//        printf("kMax=%d\n",kMax);
//        getchar();

//        if (t+dtGlobal > simTime){
//            dtGlobal =  simTime - t;
//        }


//      Bucle de reasignación en fronteras de m-regiones para asegurar correcta propagación de información

        vL[0] = v[0];
        mAux = mLocal[0];

        for (i=1; i<nCells; i++){

            vL[i] = v[i];
            mAux = mLocal[i-1];

            if (mAux < mLocal[i]){ // solo reasignar si la información de propagación rápida se encuentra a izda porque es el sentido de propagación

                mDif = mLocal[i] - mAux;

                for (j=0; mLocal[i+j] == mLocal[i] && j<nCells-1-i; j++) domainLength++; // Contar número de celdas hacia dcha con mismo m

                mPow2Dif = (int) MIN( power(2,mDif)-1,domainLength);    // evitar reasignar fuera de región
//                printf("mPow2Dif=%d, domainLneght=%d\n",mPow2Dif,domainLength);
                for (j=0; j<mPow2Dif; j++) {
                    mLocal[i+j] = (int) (mAux + log(j+1) / log2 + eps); // reasignar. si mDif = 3 (p. ej.), log2(j+1) asigna 0 1 1 2 2 2 2, para las 2^3-1=7 celdas
                    mPow2Local[i+j] = (int) power(2,mLocal[i+j]);       // en las que se modifica el valor para propagar la información
                }

//                {
//
//                    mPow2Dif = (int) pow(2,j);
//                    for (k=0; k<mPow2Dif; k++){
//                        mLocal[i+k+mPow2Dif-1] = mAux + j;
//                    }
//                }

                    i = i + domainLength - 1; // saltar hasta la siguiente celda en la que hay una diferencia de m's

            }

        }

//        for (i=0; i<nCells; i++){
//            printf("i=%d m[i]=%d 2^m[i]=%d\n",i,mLocal[i],mPow2Local[i]);
//        }

//        getchar();

        for (p=0; p<pTotal; p++){ // Bucle de actualización

            for (i=0; i<nCells; i++){

                if (p % mPow2Local[i] == 0){ // cada celda actualiza cada 2^m[i] pasadas, con m=0 actualizan en todas las pasadas

                    if (mLocal[i] == mLocal[i-1] && i>0) {
                        dv[i] = dt * mPow2Local[i] / dx[i] * FluxFunction(i, v[i], v[i-1], &FluxCounter);
                    } else {
                        dv[i] = dt * mPow2Local[i] / dx[i] * FluxFunction(i, v[i], 0.5*(v[i-1]+vL[i-1]), &FluxCounter);
                    }




                }

            }

            for (i=0; i<nCells; i++){

                if (p % mPow2Local[i] == 0){ // actualizar aquellas que se han integrado. vL se asigna al valor en el nivel temporal anterior para poder interpolar si necesario

                    vL[i] = v[i];
                    v[i] = v[i] - dv[i];
                    dv[i] = 0.0;

                }

//                printf("i=%d,v[i]=%lf\n",i,v[i]);

            }

            t = t + dt;
//            printf("dt=%lf, t=%lf\n",dt,t);
//            getchar();

        }

        IterCounter++;

//            getchar();

        }


//    t = t + dtGlobal;



    PrintToFile(&PrintCounter, x, v, t);
    printf("Flux function evaluations:%d, iterations=%d",FluxCounter,IterCounter);
    ErrorNorms(v,dx,x);
        return 0;
}

int PrintToFile(int *PrintCounter, double x[], double v[], double t){

    int i;
    char filename[40];

    FILE *output;

//	Write the numbering in the file name with the corresponding extension
    sprintf(filename,"OutputFiles/conveccionLinealLocalTimeStep%d.dat",*PrintCounter);
    output=fopen(filename,"w");

    fprintf(output,"#t=%lf\n",t);
    for(i=0;i<nCells;i++){
       fprintf(output,"%f %lf\n",x[i],v[i]);
    }

    fclose(output);
    *PrintCounter=*PrintCounter+1;

    return 0;

}

double GlobalTimeStep(double dx[], double t, double *dtLocal){

    int i = 0;
    double dt = 0.0;

    dtLocal[0] = dx[0]/c;
    dt = dtLocal[0];

    for (i = 1; i < nCells; i++){
        dtLocal[i] = dx[i]/c;
        dt = MIN(dt,dtLocal[i]);
    }

    return dt;
}

double FluxFunction(int i, double v, double vLeft, long *FluxCounter){

    *FluxCounter = *FluxCounter + 1;

    if (i==0){
        return 0.;
    }

    return c * (v - vLeft);

}

double SquareWave (double x){

    if (x < 25. || x > 55.) return 0.;

    return v0;
}

double Gaussian (double x){

    double sigma = 5.0;

    return v0*exp(-(x-40.)*(x-40.)/2./sigma/sigma);

}

void ErrorNorms (double v[], double dx[], double x[]){

    double L1=0., L2=0., Linf=0.;
    int i;

#if INITIAL_CONDITION == 1

    for (i=0; i<nCells; i++){
        L1 = L1 + fabs(v[i]-SquareWave(x[i]-c*simTime))*dx[i];
        L2 = L2 + pow(fabs(v[i]-SquareWave(x[i]-c*simTime)),2)*dx[i];
        Linf = MAX(Linf, fabs(v[i]-SquareWave(x[i]-c*simTime)));
    }

#else

    for (i=0; i<nCells; i++){
            L1 = L1 + fabs(v[i]-Gaussian(x[i]-c*simTime))*dx[i];
            L2 = L2 + pow(fabs(v[i]-Gaussian(x[i]-c*simTime)),2)*dx[i];
            Linf = MAX(Linf, fabs(v[i]-Gaussian(x[i]-c*simTime)));
    }

#endif

    printf("\n\nError norms:\n");
    printf("\tL1=%lf",L1);
    printf("\tL2=%lf",sqrt(L2));
    printf("\tLinf=%lf",Linf);

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









