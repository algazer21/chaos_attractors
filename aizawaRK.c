/**********************************************************
Alan García Zermeño
5/12/22
Resuelve el sistema dinámico de Aizawa usando Runge Kutta 4
***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DELTA 3.5
#define ALPHA 0.95 //0.5
#define BETA 0.7
#define EPSILON 0.25
#define ZETA 0.1
#define GAMMA 0.6

void createvalues(size_t,double*,double*,double*);
void RKutta4(size_t,double,double*,double*,double*);
double dx(double t, double x, double y, double z){ return((z-BETA)*x - DELTA*y); }
double dy(double t, double x, double y, double z){ return((z-BETA)*y + DELTA*x); }
double dz(double t, double x, double y, double z){ return(GAMMA + ALPHA*z - (z*z*z/3.0) - (x*x + y*y)*(1.0 + EPSILON*z) + ZETA*z*x*x*x); }

int main(){
	unsigned int n = 20000;
	double dt = 0.02;

	double *x = (double *)malloc(n*sizeof(double));
   double *y = (double *)malloc(n*sizeof(double));
   double *z = (double *)malloc(n*sizeof(double));

	x[0] = 0.1; y[0] = 0.1; z[0] = 0.08;
	RKutta4(n,dt,x,y,z);

	createvalues(n,x,y,z);
	free(x); free(y); free(z);
	printf("Memoria liberada, datos guardados en archivo dat\n");
	return 0;
}

void createvalues(size_t n, double *x, double *y, double *z){  
   FILE *pf;
   pf = fopen("data/aiz.dat", "w");
   if (pf==NULL){
      printf("No se encontro archivo\n");
      exit(1);}
   for (int i = 0; i < n; ++i){
      fprintf(pf,"%lf %lf %.10lf\n",x[i],y[i],z[i]);
   }
   fclose(pf);
}

void RKutta4(size_t n, double dt, double *x, double *y, double *z){
	double t = 0.0,k[4],l[4],m[4];
	for (int i = 0; i < (n-1); ++i){
		t += dt;
		k[0] = dx(t,x[i],y[i],z[i]);
		l[0] = dy(t,x[i],y[i],z[i]);
		m[0] = dz(t,x[i],y[i],z[i]);

		k[1] = dx(t + 0.5*dt, x[i] + 0.5*k[0]*dt, y[i] + 0.5*l[0]*dt, z[i] + 0.5*m[0]*dt);
		l[1] = dy(t + 0.5*dt, x[i] + 0.5*k[0]*dt, y[i] + 0.5*l[0]*dt, z[i] + 0.5*m[0]*dt);
		m[1] = dz(t + 0.5*dt, x[i] + 0.5*k[0]*dt, y[i] + 0.5*l[0]*dt, z[i] + 0.5*m[0]*dt);

		k[2] = dx(t + 0.5*dt, x[i] + 0.5*k[1]*dt, y[i] + 0.5*l[1]*dt, z[i] + 0.5*m[1]*dt);
		l[2] = dy(t + 0.5*dt, x[i] + 0.5*k[1]*dt, y[i] + 0.5*l[1]*dt, z[i] + 0.5*m[1]*dt);
		m[2] = dz(t + 0.5*dt, x[i] + 0.5*k[1]*dt, y[i] + 0.5*l[1]*dt, z[i] + 0.5*m[1]*dt);

		k[3] = dx(t + dt, x[i] + k[2]*dt, y[i] + l[2]*dt, z[i] + m[2]*dt);
		l[3] = dy(t + dt, x[i] + k[2]*dt, y[i] + l[2]*dt, z[i] + m[2]*dt);
		m[3] = dz(t + dt, x[i] + k[2]*dt, y[i] + l[2]*dt, z[i] + m[2]*dt);

		x[i+1] = x[i] + dt*(k[0] +2.0*k[1] + 2.0*k[2] + k[3])/6.0;
		y[i+1] = y[i] + dt*(l[0] +2.0*l[1] + 2.0*l[2] + l[3])/6.0;
		z[i+1] = z[i] + dt*(m[0] +2.0*m[1] + 2.0*m[2] + m[3])/6.0;
	}
}
