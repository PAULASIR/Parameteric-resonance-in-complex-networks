/********* duffing oscillator with parametric perturbation**************/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<limits.h>

#define N 30
double om,om1, f, p ,mu,y[20],x, eps;
void RK4(int,double,double,double[],void (*DGL)(double,double[],double[]));
void DGL(double, double[],double[]);

double sinc(const double x)
{
if (x==0)
return 1;
return sin(x)/x;
}

void main()
{
    int i,j;
    int nn=2;
    double t,h,pi;
    double xmax,xmin,z[50001];
   

 //clrscr();
 
 FILE *fp1;
 fp1=fopen("pend.dat","w");


   y[1]= (float) rand()/(double)RAND_MAX*-4.0+2.0; 
   y[2]= (float) rand()/(double)RAND_MAX*-4.0+2.0;

 mu = 0.05;
//step size
 h=0.01;
 f = 1.5;
 om = 1.0;
 om1 = 1.0;
 eps =0.2;

xmin=INT_MAX;  xmax=-INT_MAX;
    for(i=1;i<=50000;i++)
	  { 
	    t=h*(double)(i);
	    RK4(nn,h,t,y,DGL);
	    if(i>=10000)
          {
             fprintf(fp1,"%f  %f   %f\n", sinc(y[1]), y[1], y[2]);
             
            z[i]=y[1];
             if(z[i]>xmax)
              xmax=z[i];
             if(z[i]<xmin)
              xmin=z[i];}       
	      }

	printf("process over!!");
}

//************************RK4 SUBROUTINE*********************************//
void RK4(int nn,double h,double t,double y[20],
	   void (*DGL)(double,double[],double[]))
{
	   int i;
	   double k1[10],k2[10],k3[10],k4[10];
	   double yaux[20];

	   DGL(t,y,k1);
	   for(i=1;i<=nn;i++)
	   {
	   yaux[i]=y[i]+h*k1[i]/2.0;
	   }
	   DGL(t+h/2.0,yaux,k2);
	   for(i=1;i<=nn;i++)
	   {
	   yaux[i]=y[i]+h*k2[i]/2.0;
	   }
	   DGL(t+h/2.0,yaux,k3);
	   for(i=1;i<=nn;i++)
	   {
	   yaux[i]=y[i]+h*k3[i];
	   }
	   DGL(t+h,yaux,k4);
	   for(i=1;i<=nn;i++)
	   {
	      y[i]=y[i]+h*((k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0);
	   }
}
//*********************FUNCTION SUBROUTINE********************************//
void DGL(double t,double y[20],double F[5])
{ 

  
     p =3*(pow(om,4) - 4*om*om + 8*mu*mu)/2*(3*om*om - 4*mu*mu);

 F[1]= y[2];
 F[2]= -mu*y[1] - eps*( 1 + p*sinc(y[1]))*sin(y[1]) + f*sin(om1*t) ; 
}


