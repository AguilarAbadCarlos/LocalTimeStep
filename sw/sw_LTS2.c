#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define L 100.0			// Longitude [m]
#define simTime 5.0	// Simulation time [s]
#define CFL 1.0		// CFL number
#define g 9.81

#define INITIAL_CONDITION 4
#define MESH_REGULARITY 2

#define v1 1.0			// Initial velocity
#define v2 3.0
#define x1 30.
#define x2 50.
#define PI (4*atan(1))
#define n_manning 0.

#define EPS 1E-5
#define INFTY 1E8
#define M_LIMITE 2
#define TOL14 1e-14

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

//int PrintToFile(int *PrintCounter, double x[], double v[], double exact[], double t);
double GlobalTimeStep (double *Q, double *H, double *dx, double *dtLocal);
double FluxFunction (int i, double v, double vL, double vR, int *FluxCounter);
double SquareWave (double x);
void ErrorNorms (double v[], double dx[], double exact[]);
void WaveSpeeds_wall (double Q1, double H1, double Q2, double H2, double *u, double *c);
void WaveSpeeds_cell (double Q, double H, double *u, double *c);
double Entropy_fix (double lambda, double lambdaL, double lambdaR, int LR);
int PrintToFile (int *PrintCounter, double *H, double *Q, double *x, double *z, double t);
double SourceTerms(double QL, double HL, double zL, double QR, double HR, double zR, double dx);
double RadioHidraulico(double H);
double signo(double x);

int power(int base, unsigned int exp);

//double L, simTime, dx, v1, v2, x1, x2;
int nCells = 1000;

int main()
{

    int i,j,k;
    int PrintCounter=0, IterCounter=0;

    double *Q, *H, *dQ, *dH, *QOld, *HOld;
    double *u, *c, *uOld, *cOld;

    double *x, *dx, dt=0., t=0., dtGlobal=0., *dtLocal;
    double *z;
    double lambdaL1=0., lambdaL2=0., lambdaR1=0., lambdaR2=0., lambda1=0., lambda2=0., uAux=0., cAux=0.;
    double *lambda_plus,*lambda_minus;
    double alpha1=0., alpha2=0.;
    double source=0.;

    int *mLocal, mMax=0, *mPow2Local, mAux=0, mDif=0, mPow2Dif=0, domainLength=0;
    int p=0, pTotal=0;

    double log2 = log(2.), eps = 1e-6;

    /// Inicialización de variables

    Q = calloc(nCells, sizeof(double));
    QOld = calloc(nCells, sizeof(double));
    dQ = calloc(nCells, sizeof(double));
    H = calloc(nCells, sizeof(double));
    HOld = calloc(nCells, sizeof(double));
    dH = calloc(nCells, sizeof(double));
    u = calloc(nCells, sizeof(double));
    c = calloc(nCells, sizeof(double));

    x = calloc(nCells, sizeof(double));
    dx = calloc(nCells, sizeof(double));
    dtLocal = calloc(nCells, sizeof(double));
    z = calloc(nCells, sizeof(double));

    mLocal = calloc(nCells, sizeof(int));
    mPow2Local = calloc(nCells, sizeof(int));

    lambda_plus = calloc(nCells-1,sizeof(double));
    lambda_minus = calloc(nCells-1,sizeof(double));

/// Mesh setup
#if MESH_REGULARITY == 0
    // Undefined
#elif MESH_REGULARITY == 2

    int i_celda_fina = nCells/4;
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

    k = nCells/4-1;
    int n_factor = 128;
    double f = pow(n_factor,1./(k));

    double d = L/2./(n_factor+f*(1-pow(f,k))/(1-f)+(nCells+2)/4);

    dx[0] = n_factor*d;
    x[0] = dx[0]/2.;

    for (i=1; i<=k; i++)
    {
        dx[i] = dx[i-1]/f;
        x[i] = x[i-1] + 0.5*(dx[i]+dx[i-1]);
    }

    for (   ; i<=k+1+nCells/2; i++)
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
#else
    // Regular mesh
    dt = L/nCells;
    x[0] = dt/2.;
    dx[0] = dt;
    for (i=1; i<nCells; i++){
        dx[i] = dt;
        x[i] = x[i-1] + 0.5 * (dx[i-1]+dx[i]);
    }
#endif // MESH_REGULARITY

/// Initial condition setup
#if INITIAL_CONDITION == 1
    for (i=0; i<nCells; i++){
        Q[i] = (x[i] < L/2. ? 0. : 0.);
        H[i] = (x[i] < L/2. ? 1. : 0.30179953); //1 : 0.30179953
        z[i] = (x[i] < L/2. ? 0. : 0.05); // 0. : 0.05
    }
#elif INITIAL_CONDITION == 2
    for (i=0; i<nCells; i++){
        Q[i] = (x[i] < L/2. ? 0.4 : 0.);
        H[i] = (x[i] < L/2. ? 4. : 0.50537954); //1 : 0.30179953
        z[i] = (x[i] < L/2. ? 0. : 1.5); // 0. : 0.05
    }
#elif INITIAL_CONDITION == 3
    for (i=0; i<nCells; i++){
        Q[i] = (x[i] < L/2. ? 3.75 : 0.);
        H[i] = (x[i] < L/2. ? 2.5 : 2.49977381); //1 : 0.30179953
        z[i] = (x[i] < L/2. ? 0. : 0.25); // 0. : 0.05
    }
#elif INITIAL_CONDITION == 4
    for (i=0; i<nCells; i++){
        Q[i] = (x[i] < L/2. ? 3.0 : 0.);
        H[i] = (x[i] < L/2. ? 1.5 : 0.16664757); //1 : 0.30179953
        z[i] = (x[i] < L/2. ? 0. : 2.0); // 0. : 0.05
    }
#elif INITIAL_CONDITION == 5
    for (i=0; i<nCells; i++){
        Q[i] = (x[i] < L/2. ? 0.2 : 0.);
        H[i] = (x[i] < L/2. ? 1. : 0.04112267); //1 : 0.30179953
        z[i] = (x[i] < L/2. ? 0.25 : 0.0); // 0. : 0.05
    }
#elif INITIAL_CONDITION == 6
    for (i=0; i<nCells; i++){
        Q[i] = (x[i] < L/2. ? 0.21 : 0.);
        H[i] = (x[i] < L/2. ? 0.6 : 0.02599708); //1 : 0.30179953
        z[i] = (x[i] < L/2. ? 1.2 : 0.0); // 0. : 0.05
    }
#else
    for (i=0; i<nCells; i++){
        Q[i] = (x[i] < L/2. ? 0. : 0.);
        H[i] = (x[i] < L/2. ? 1. : 0.30179953); //1 : 0.30179953
        z[i] = (x[i] < L/2. ? 0. : 0.05); // 0. : 0.05
    }
#endif // INITIAL_CONDITION

//  Testing mesh initialization
//    for (i=0;i<nCells;i++){
//        printf("%d\t%.3f\t%.3f\n",i,dx[i],x[i]);
//    }

//  Testing velocity initialization
//    for (i=0;i<nCells;i++){
//        printf("%d\t%.3f\t%.3f\n",i,H[i],Q[i]);
//    }
//  getchar();

    PrintToFile(&PrintCounter,H,Q,x,z,t);

    while (t<simTime){

        for (i=0; i<nCells; i++)
        {
            WaveSpeeds_cell(Q[i], H[i], &u[i], &c[i]);
//            printf("\ni=%d\tH=%f\tQ=%f\tu=%f\tc=%f",i,H[i],Q[i],u[i],c[i]);
        }
//        getchar();

        for (j=1; j<nCells; j++)
        {
            WaveSpeeds_wall(Q[j-1],H[j-1],Q[j],H[j],&uAux,&cAux);
            lambda_plus[j-1]=uAux+cAux;
            lambda_minus[j-1]=uAux-cAux;
        }

        dt = GlobalTimeStep(Q,H,dx,dtLocal);

        if (t+dt>simTime)
        {
            dt = simTime-t;
            if ((simTime-t)<1e-12) break;
        }

//        if (IterCounter % 1000 ==0) printf("t=%.15lf\tdt=%.15lf\n",t,dt);
//        getchar();

        /// Inicialización de LTS2
        mMax = 0;
        for (i=0; i<nCells; i++){
            mLocal[i] = (int) MIN((log(dtLocal[i]/dt)/log2+eps),M_LIMITE);

//            printf("mLocal[%d]=%d\n",i,mLocal[i]);
            mPow2Local[i] = power(2,mLocal[i]);
            mMax = MAX(mMax,mLocal[i]);
        }
        pTotal = power(2,mMax);
//        printf("pTotal = %d\n",pTotal);
        if (t+pTotal*dt > simTime)
        {
            dt = (simTime-t)/pTotal;
        }
//        getchar();

        /// Actualización de m-values

        for (i=1; i<nCells; i++){

            mAux = mLocal[i-1];

            if (mAux != mLocal[i]){

                mDif = mLocal[i] - mAux;

                if (mDif>0)
                {
                    for (j=0; mLocal[i+j] == mLocal[i] && j<nCells-1-i; j++) domainLength++;

                    mPow2Dif = (int) MIN( pow(2,mDif-1)+3,domainLength);
                    for (j=0; j<=3; j++){
                        mLocal[i+j] = mAux;
                        mPow2Local[i+j] = (int) pow(2,mLocal[i+j]);
                    }
                    for (; j<=mPow2Dif; j++) {
                        mLocal[i+j] = (int) (mAux + 1.*log(j-2) / log2 + eps);
                        mPow2Local[i+j] = (int) pow(2,mLocal[i+j]);
                    }
                        i = i + domainLength - 1;
                        domainLength = 0;
                }
                else
                {
                    mDif = -mDif;
                    mAux = mLocal[i];
                    for (j=0; mLocal[i-j] == mLocal[i] && i-j>=0; j++) domainLength++;

                    mPow2Dif = (int) MIN( pow(2,mDif-1)+3,domainLength);

                    for (j=0; j<=3; j++){
                        mLocal[i-j-1] = mAux;
                        mPow2Local[i-j-1] = (int) pow(2,mLocal[i]);
                    }
                    for (; j<=mPow2Dif; j++) {
                        mLocal[i-j-1] = (int) (mAux + 1.*log(j-2) / log2 + eps);
                        mPow2Local[i-j-1] = (int) pow(2,mLocal[i-j-1]);
                    }
                        //i = i + domainLength - 1;
                        domainLength = 0;
                }



            }

        }

//        for (i=0; i<nCells; i++)
//        {
//            printf("mLocal[%d]=%d\n",i,mLocal[i]);
//        }
//        getchar();

        /// Algoritmo de integración

        for (i=0; i<nCells; i++)
        {
            QOld[i] = Q[i];
            HOld[i] = H[i];
//            WaveSpeeds_cell(QOld[i], HOld[i], &uOld[i], &cOld[i]);
        }

        for (p=0; p<pTotal; p++)
        {
            for (i=1; i<nCells-1; i++)
            {
                if (p % mPow2Local[i] == 0)
                {
                    if (mLocal[i] >= mLocal[i-1]){
                        WaveSpeeds_wall(Q[i-1],H[i-1],Q[i],H[i],&uAux,&cAux);
                        lambda1 = MAX(uAux+cAux,0.);
                        lambda2 = MAX(uAux-cAux,0.);
                        source = SourceTerms(Q[i-1], H[i-1], z[i-1], Q[i], H[i], z[i], dx[i]);
                        alpha1 = (cAux-uAux)*(H[i]-H[i-1])/2./cAux + (Q[i]-Q[i-1])/2./cAux-source/fabs(uAux+cAux);
                        alpha2 = (cAux+uAux)*(H[i]-H[i-1])/2./cAux - (Q[i]-Q[i-1])/2./cAux+source/fabs(uAux-cAux);
                        WaveSpeeds_cell(Q[i-1],H[i-1],&uAux,&cAux);
                        lambdaL1 = uAux+cAux;
                        lambdaL2 = uAux-cAux;
                        WaveSpeeds_cell(Q[i],H[i],&uAux,&cAux);
                        lambdaR1 = uAux+cAux;
                        lambdaR2 = uAux-cAux;
                    } else {
                        WaveSpeeds_wall(0.5*(Q[i-1]+QOld[i-1]),0.5*(H[i-1]+HOld[i-1]),Q[i],H[i],&uAux,&cAux);
                        lambda1 = MAX(uAux+cAux,0.);
                        lambda2 = MAX(uAux-cAux,0.);
                        source = SourceTerms(0.5*(Q[i-1]+QOld[i-1]), 0.5*(H[i-1]+HOld[i-1]), z[i-1], Q[i], H[i], z[i], dx[i]);
                        alpha1 = (cAux-uAux)*(H[i]-0.5*(H[i-1]+HOld[i-1]))/2./cAux + (Q[i]-0.5*(Q[i-1]+QOld[i-1]))/2./cAux -source/fabs(uAux+cAux);
                        alpha2 = (cAux+uAux)*(H[i]-0.5*(H[i-1]+HOld[i-1]))/2./cAux - (Q[i]-0.5*(Q[i-1]+QOld[i-1]))/2./cAux +source/fabs(uAux-cAux);
                        WaveSpeeds_cell(0.5*(Q[i-1]+QOld[i-1]),0.5*(H[i-1]+HOld[i-1]),&uAux,&cAux);
                        lambdaL1 = uAux+cAux;
                        lambdaL2 = uAux-cAux;
                        WaveSpeeds_cell(Q[i],H[i],&uAux,&cAux);
                        lambdaR1 = uAux+cAux;
                        lambdaR2 = uAux-cAux;
                    }

                    dH[i] += dt*mPow2Local[i]/dx[i] * (Entropy_fix(fabs(lambda1),lambdaL1,lambdaR1,1)*alpha1+Entropy_fix(fabs(lambda2),lambdaL2,lambdaR2,1)*alpha2);
                    dQ[i] += dt*mPow2Local[i]/dx[i] * ((Entropy_fix(fabs(lambda1),lambdaL1,lambdaR1,1)*alpha1)*lambda1+(Entropy_fix(fabs(lambda2),lambdaL2,lambdaR2,1)*alpha2)*lambda2);
//                    dQ[i] -= dt*mPow2Local[i]/dx[i] * (lambda1-lambda2) * source;

                    if (mLocal[i] >= mLocal[i+1]){
                        WaveSpeeds_wall(Q[i+1],H[i+1],Q[i],H[i],&uAux,&cAux);
                        lambda1 = MIN(uAux+cAux,0.);
                        lambda2 = MIN(uAux-cAux,0.);
                        source = SourceTerms(Q[i], H[i], z[i], Q[i+1], H[i+1], z[i+1], dx[i]);
                        alpha1 = (cAux-uAux)*(H[i]-H[i+1])/2./cAux + (Q[i]-Q[i+1])/2./cAux-source/fabs(uAux+cAux);
                        alpha2 = (cAux+uAux)*(H[i]-H[i+1])/2./cAux - (Q[i]-Q[i+1])/2./cAux+source/fabs(uAux-cAux);
                        WaveSpeeds_cell(Q[i+1],H[i+1],&uAux,&cAux);
                        lambdaL1 = uAux+cAux;
                        lambdaL2 = uAux-cAux;
                        WaveSpeeds_cell(Q[i],H[i],&uAux,&cAux);
                        lambdaR1 = uAux+cAux;
                        lambdaR2 = uAux-cAux;
                    } else {
                        WaveSpeeds_wall(0.5*(Q[i+1]+QOld[i+1]),0.5*(H[i+1]+HOld[i+1]),Q[i],H[i],&uAux,&cAux);
                        lambda1 = MIN(uAux+cAux,0.);
                        lambda2 = MIN(uAux-cAux,0.);
                        source = SourceTerms(Q[i], H[i], z[i], 0.5*(Q[i+1]+QOld[i+1]), 0.5*(H[i+1]+HOld[i+1]), z[i+1], dx[i]);
                        alpha1 = (cAux-uAux)*(H[i]-0.5*(H[i+1]+HOld[i+1]))/2./cAux + (Q[i]-0.5*(Q[i+1]+QOld[i+1]))/2./cAux-source/fabs(uAux+cAux);
                        alpha2 = (cAux+uAux)*(H[i]-0.5*(H[i+1]+HOld[i+1]))/2./cAux - (Q[i]-0.5*(Q[i+1]+QOld[i+1]))/2./cAux+source/fabs(uAux+cAux);
                        WaveSpeeds_cell(0.5*(Q[i+1]+QOld[i+1]),0.5*(H[i+1]+HOld[i+1]),&uAux,&cAux);
                        lambdaL1 = uAux+cAux;
                        lambdaL2 = uAux-cAux;
                        WaveSpeeds_cell(Q[i],H[i],&uAux,&cAux);
                        lambdaR1 = uAux+cAux;
                        lambdaR2 = uAux-cAux;
                    }

                    dH[i] += dt*mPow2Local[i]/dx[i] * (Entropy_fix(fabs(lambda1),lambdaL1,lambdaR1,0)*alpha1+Entropy_fix(fabs(lambda2),lambdaL2,lambdaR2,0)*alpha2);
                    dQ[i] += dt*mPow2Local[i]/dx[i] * ((Entropy_fix(fabs(lambda1),lambdaL1,lambdaR1,0)*alpha1)*lambda1+(Entropy_fix(fabs(lambda2),lambdaL2,lambdaR2,0)*alpha2)*lambda2);
//                    dQ[i] -= dt*mPow2Local[i]/dx[i] * (lambda1-lambda2) * source;

                }

            }

            for (i=0; i<nCells; i++)
            {
                if (p % mPow2Local[i] == 0)
                {
                    HOld[i] = H[i];
                    QOld[i] = Q[i];

                    H[i] -= dH[i];
                    Q[i] -= dQ[i];

                    dH[i] = 0.;
                    dQ[i] = 0.;
                }
            }

            IterCounter++;

            t = t + dt;

        }


//        PrintToFile(&PrintCounter,H,Q,x,z,t);

    }

    PrintCounter = 999;
    PrintToFile(&PrintCounter,H,Q,x,z,t);

    free(Q);
    free(QOld);
    free(dQ);
    free(H);
    free(HOld);
    free(dH);
    free(u);
    free(c);
    free(x);
    free(dx);
    free(dtLocal);
    free(z);
    free(mLocal);
    free(mPow2Local);
    free(lambda_plus);
    free(lambda_minus);

    return 0;
}

void WaveSpeeds_wall (double Q1, double H1, double Q2, double H2, double *u, double *c)
{
    *u = (Q1/sqrt(H1)+Q2/sqrt(H2)) / (sqrt(H1)+sqrt(H2));
    *c = sqrt(0.5*g*(H1+H2));
}
void WaveSpeeds_cell (double Q, double H, double *u, double *c)
{
    *u = Q/H;
    *c = sqrt(g*H);
}

double GlobalTimeStep (double *Q, double *H, double *dx, double *dtLocal) {                  // Calculated as the minimum of the time step for each cell boundary

    //printf("Calculo del paso temporal para la siguiente iteración:\n");
    int i;
    double c, u, dt;
    dt=INFTY;
    //printf("i=0,\tdt=%.3lf\n",dt);

    for (i=0;i<nCells-1;i++){
        WaveSpeeds_wall(Q[i],H[i],Q[i+1],H[i+1],&u,&c);   // Wall-based velocidad de propagación
        dtLocal[i] = dx[i]/(fabs(u)+c);
        dt=MIN(dt,dtLocal[i]);
        WaveSpeeds_cell(Q[i],H[i],&u,&c);                   // Cell-based velocidad de propagación
        dtLocal[i] = dx[i]/(fabs(u)+c);
//        dt=MIN(dt,dtLocal[i]);
        //printf("i=%d,\tdt=%.3lf\n",i,dt);
    }
    WaveSpeeds_cell(Q[i],H[i],&u,&c);                   // Cell-based velocidad de propagación
    dtLocal[i] = dx[i]/(fabs(u)+c);

    return dt;

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

double Entropy_fix (double lambda, double lambdaL, double lambdaR, int LR)
{
    if (lambda==0.) return 0.;
    // LR = 0 ==> updating left cell
    // LR = 1 ==> updating right cell
    if (lambdaL*lambdaR<0. && lambdaL < 0.)
    {
        if (LR == 0) // updating L
        {
            return lambdaL * (lambdaR-lambda) / (lambdaR-lambdaL);
        }
        if (LR == 1) // updating R
        {
            return lambdaR * (lambda-lambdaL) / (lambdaR-lambdaL);
        }
    }
//    return MAX(lambda,lambda_fix);
    return lambda;
}

int PrintToFile (int *PrintCounter, double *H, double *Q, double *x, double *z, double t){

    int i=0;
    char filename[200];
    double c=0.0;

    FILE *output;

//	Write the numbering in the file name with the corresponding extension
    sprintf(filename,"OutputFiles/dataSW_LTS%d.dat",*PrintCounter);
    output=fopen(filename,"w");

    fprintf(output,"#t=%lf\n",t);
    fprintf(output,"#x\th\tQ/b\tz\n");
    for(i=0;i<nCells;i++){
//        c=sqrt(g*H[i]);
        fprintf(output,"%f\t%.15lf\t%.15lf\t%lf\n",x[i],H[i],Q[i],z[i]);
    }

    fclose(output);
    *PrintCounter=*PrintCounter+1;

    return 0;

}

double SourceTerms(double QL, double HL, double zL, double QR, double HR, double zR, double dx)
{

    double c,u,H,s,dz,hj,dzprime;
    int j;

    WaveSpeeds_wall(QL,HL,QR,HR,&u,&c);
    H = 0.5*(HL+HR);

    dz=zR-zL;

    if (dz>=0){
        hj=HL;
    } else {
        hj=HR;
    }

    if (dz>=0 && (HL+zL)<zR){
        dz = HL;
    } else {
        if (dz<0 && (HR+zR)<zL){
            dz = -HR;
        }
    }

    s = g/2./c*(-(hj-fabs(dz)/2.)*dz-H*dx*n_manning*n_manning*u*fabs(u)/pow(MAX(HL,HR),4./3.));

//    s = g/2./c*(-H*(zR-zL)-H*dx*n_manning*n_manning*u*fabs(u)/pow(MAX(HL,HR),4./3.));

    return s;
}

double RadioHidraulico(double H)
{
    return H/(1.+2.*H);
}

double signo(double x)
{
    if (x>=0.) return 1.;
    return -1.;
}
