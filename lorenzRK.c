/**********************************************************
Alan García Zermeño
5/12/22
Resuelve el sistema dinámico de Lorenz usando Runge Kutta 4
***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIGMA 10.0
#define RHO 28.0  //17.0
#define BETA 8.0/3.0

void createvalues(size_t,double*,double*,double*);
void RKutta4(size_t,double,double*,double*,double*);
double dx(double t, double x, double y, double z){ return(SIGMA*(y-x)); }
double dy(double t, double x, double y, double z){ return(RHO*x - y - x*z); }
double dz(double t, double x, double y, double z){ return(x*y- BETA*z); }

int main(){
	unsigned int n = 30000;
	double dt = 0.005;

	double *x = (double *)malloc(n*sizeof(double));
   double *y = (double *)malloc(n*sizeof(double));
   double *z = (double *)malloc(n*sizeof(double));

	x[0] = 1.0; y[0] = 2.0; z[0] = 7.0;
	RKutta4(n,dt,x,y,z);

	createvalues(n,x,y,z);
	free(x); free(y); free(z);
	printf("Memoria liberada, datos guardados en archivo dat\n");
	return 0;
}

void createvalues(size_t n, double *x, double *y, double *z){  
   FILE *pf;
   pf = fopen("data/lor.dat", "w");
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
