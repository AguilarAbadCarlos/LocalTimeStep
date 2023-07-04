#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 100.0			    // Longitude [m]
#define simTime 10.0	    // Simulation time [s]
//#define imax 100		    // Number of nodes
//#define dx (L/(imax))     // Spatial increment
#define CFL 1.0		        // CFL number
#define CI 1                // Condición inicial: 1 CUADRADA, 2 RAMPA

#define INITIAL_CONDITION 1
#define MESH_REGULARITY 2

#define v1 -1.			// Initial velocity
#define v2 3.0
#define x1 30.
#define x2 55.
#define PI (4*atan(1))

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)


int PrintToFile(int *PrintCounter, double *x, double *v, double *exact, double t);
double TimeStep (double *v, double *dx);
double FluxFunction (int i, double *v, long *FluxCounter);
double SquareWave (double x);
void ErrorNorms (double *v, double *dx, double *exact);
int exactBurgers(double t, double *x, double *exact);

//double L, simTime, dx, v1, v2, x1, x2;
int imax, nCells;

int main () {

//    FILE *inParams;
//    inParams=fopen("parametrosBurgers.txt","r");
////    f'{L:.1f} {simTime:.1f} {nCells:d} {dx:.5f} {v1:.1f} {v2:.1f} {x1:.1f} {x2:.1f}'
//    fscanf(inParams, "%lf %lf %i %lf %lf %lf %lf %lf", &L, &simTime, &imax, &dx, &v1, &v2, &x1, &x2);
//    fclose(inParams);

    imax = 10000;
    nCells = imax;

# if MESH_REGULARITY == 1
    int irregularity = 7;
    imax = 40 * (1<<irregularity) + 10 * irregularity + 60;
# endif // MESH_REGULARITY

    int i,iter=0,k=0,printFreq=10, j=0;
	double *v, *exact, *dx;		// Velocity arrays (conserved variables)
	double *vL;                // Interpolation points
	double *x;			// Discharge, distance
	double t=0.0,dt,tPrint=1.0;                // Time
	double eps = 1e-7;                       // Area
	double cPlus=0.0,cMinus=0.0, c=0.0;
	int PrintCounter=0, IterCounter=0;             // Print counter
    long FluxCounter=0;

    printf("\n   Burgers convection");
	printf("\n   --------------------------------------------------");


// Reading simulation parameters

	printf("\n>> Simulation setup loaded");
	printf("\n   Simulation time: %.5lf s",simTime);
	printf("\n   CFL: %.2lf",CFL);
	printf("\n   Number of cells: %d",imax);
//	printf("\n   dx: %.2lf m",dx);
	printf("\n   Pipe length: %.2lf m\n\n",L);


// Variable initialization
//	for (i=0;i<imax;i++){
//        v[i]=0.0;
//        vL[i]=0.0;
//        vR[i]=0.0;
//        Q[i]=0.0;
////        x[i]=i*dx;
//        exact[i]=0.0;
//	}
	//A=PI*D*D/4.0;
	v = calloc(imax,sizeof(double));
	vL = calloc(imax,sizeof(double));
//	vR = calloc(imax,sizeof(double));
	dx = calloc(imax,sizeof(double));
	x = calloc(imax,sizeof(double));
	exact = calloc(imax,sizeof(double));

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
	for (; i<imax; i++){
        dx[i] = 1.0;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
	}

#elif MESH_REGULARITY == 2

    int i_celda_fina = imax/2;
    double toLenght = L/imax;
    dx[0] = 1.0*toLenght;
    x[0] = dx[0]/2.;
    for (i=1; i<i_celda_fina; i++){
        dx[i] = 1.*toLenght;
        x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i]);
    }
    dx[i] = 0.001*toLenght;
    x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i]);
    i++;
    for (; i<imax; i++){
        dx[i] = 1.*toLenght;
        x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i]);
    }

#elif MESH_REGULARITY == 3

    int l = imax/100;
    k = (imax-l)/2-1;
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

    for (   ; i<imax-1; i++)
    {
        dx[i] = dx[i-1]*f;
        x[i] = x[i-1] + 0.5*(dx[i]+dx[i-1]);
    }

    dx[i] = n_factor*d;
    x[i] = x[i-1] + 0.5*(dx[i]+dx[i-1]);

#else // REGULAR MESH

//    dx[0] = dx;
    x[0] = dx/2.;
    for (i=1; i<imax; i++){
//        dx[i] = dx;
        x[i] = x[i-1] + dx;
    }

#endif // MESH_REGULARITY

// INITIAL CONDITION SETUP

#if INITIAL_CONDITION == 1

        for (i=0; i<imax; i++) v[i] = SquareWave(x[i]);

#elif INITIAL_CONDITION == 2

        for (i=0; x[i]<x1; i++) v[i] = v2;
        for (; x[i]<x2; i++)    v[i] = v2 - 2* v2 *(i*dx-x1) / (x2-x1);
        for (; i<imax; i++)   v[i] = -v2;

#else

        for (i=0; i<imax; i++) v[i] = Gaussian(x[i]);

#endif // INITIAL_CONDITION

//  Testing mesh initialization
//    for (i=0;i<imax;i++){
//        printf("%d\t%.3f\t%.3f\n",i,dx[i],x[i]);
//    }

//  Testing velocity initialization
//    for (i=0;i<imax;i++){
//        printf("%d\t%.3f\n",i,v[i]);
//    }
//  getchar();

    exactBurgers(t,x,exact);
    PrintToFile(&PrintCounter, x, v, exact, t);

// -------------------------------- TIME LOOP --------------------------------
    while (t<simTime){

        // Obtención del paso temporal a utilizar en cada momento
        dt=MIN(CFL*TimeStep(v,dx),simTime-t);

//        if (IterCounter % 100 == 0) printf("t=%lfdt=%lf\n",t,dt);
//        printf("dt=%lf\n",dt);

//        if(t + dt > simTime && simTime - t >= eps){
//                dt = simTime - t;
//                printf("Patata dt=%lf\n",dt);
//        }

        // Solver for each cell using the explicit first-order upwind scheme
        for (i=0;i<imax;i++){

            vL[i]=v[i]-FluxFunction(i,v,&FluxCounter)*(dt/dx[i]);

        }
//        getchar();

        // Actualización con la variable conservada
        for (i=0;i<imax;i++){
//            printf("i=%d\tv[i]=%lf\tvL[i]=%lf\n",i,v[i],vL[i]);
            v[i]=vL[i];
        }
//        getchar();

        // Contour conditions for each case
//        switch (CI){
//            case 1:
//                v[0]=-v1; v[imax-1]=-v1;
//                break;
//            case 2:
//                v[0]=v2; v[imax-1]=v[imax-2];
//                break;
//            case 3:
//                v[0]=-v1; v[imax-1]=v[imax-2];
//                break;
//        }

        // Printing to file every second
//        if((t/tPrint)>=k){
//            exactBurgers(t,exact);
//            PrintToFile(&PrintCounter,v,exact,t);
//            printf(">> Volcado %d en t=%.3lf s\n",PrintCounter-1,t);
//            k++;
//        }

        t=t+dt;
//        printf("%lf  %lf\n", dt, t);

        IterCounter++;

        /*printf("\nIteración:\n");
        for (i=320;i<331;i++){
            printf("i=%d\tu=%lf\n",i,v[i]);
        }*/

        //getchar();

    }
// -------------------------------- END OF TIME LOOP --------------------------------

// Print final state to file 'dataBurgers999.dat' and plotting of the results
    exactBurgers(t,x,exact);
    PrintToFile(&PrintCounter, x, v, exact, t);
    printf("Flux function evaluations:%ld, iterations=%d",FluxCounter,IterCounter);
    ErrorNorms(v,dx,exact);

    free(v);
    free(vL);
    free(exact);
    free(x);
    free(dx);


    return 0;
}


// ------------------------------ FUNCTIONS ------------------------------
int PrintToFile(int *PrintCounter, double *x, double *v, double *exact, double t){

    int i;
    char filename[80];

    FILE *output;

//	Write the numbering in the file name with the corresponding extension
    sprintf(filename,"OutputFiles/burgersUpwind_%d.dat",*PrintCounter);
    output=fopen(filename,"w");

    fprintf(output,"#t=%lf\n",t);
    for(i=0;i<imax;i++){
       fprintf(output,"%lf %lf %lf\n",x[i],v[i],exact[i]);
    }

    fclose(output);
    *PrintCounter=*PrintCounter+1;

    return 0;

}

double TimeStep (double *v, double *dx) {                  // Calculated as the minimum of the time step for each point

    int i;
    double aux;
    double dt=dx[0]/fabs(v[0]);

    for (i=1;i<imax;i++){
        aux=dx[i]/fabs(v[i]);
        if (aux<dt) dt=aux;
    }

    return dt;
}

double FluxFunction(int i, double *v, long *FluxCounter){

    *FluxCounter = *FluxCounter + 1;

    if (i==0 || i==imax-1){
        return 0.;
    }

    if (v[i]*v[i-1]<0. && v[i-1]<0.){
//        return v[i-1]*(v[i]-0.5*(v[i]+v[i-1]));
        return v[i]*(0.5*(v[i]+v[i-1])-v[i-1]);
    }

    if (v[i+1]*v[i]<0. && v[i]<0.){
//        return v[i+1]*(0.5*(v[i+1]+v[i])-v[i]);
        return v[i]*(v[i+1]-0.5*(v[i+1]+v[i]));
    }

    double cPlus,cMinus;
    cPlus = MAX(0.5*(v[i]+v[i-1]),0.);
    cMinus = MIN(0.5*(v[i]+v[i+1]),0.);

//    printf("i=%d, v=%lf, vL=%lf, vR=%lf, cPlus=%lf, cMinus=%lf\n",i,v[i],v[i-1],v[i+1],cPlus,cMinus);

    return cPlus*(v[i]-v[i-1])+cMinus*(v[i+1]-v[i]);

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


    for (i=0; i<imax; i++){
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
    for (;i<imax;i++){
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

