/* Local coupling -parametric resonance in a Pendula */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>
#include<limits.h>

int i,j,nn=20,n=2;
double y[101][501],k;
int p;  //# of nearest neighbours
double k_c,om,om1, f, p_coef ,mu,x, eps;
double pwr;
double complex pwr_i;
void RK4(int,int,double,double,double[101][501],
                   void (*DGL)(double,double[101][501],double[101][501]));
void DGL(double, double[101][501],double[101][501]);


double sinc(const double x)
{
if (x==0)
return 1;
else
return sin(x)/x;
}


void main()
{
//nn= no 0f oscillators, n=dimension of the model//
double t,h,pi;
double z[100000];
double amplitude,vmin,vmax,c_max,c_min;
double x_max,x_min;
int cnt;
srand(time(0));

FILE *fp1,*fp2,*fp3;
fp1=fopen("om_bif1.dat","w");
fp2 = fopen("tim.dat", "w");
 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-2.0+1.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-2.0-1.0;
 }

for(om = 6.0   ; om >= 0.0; om -= 0.05 )
{
//parameters
//scanf("k_c = %lf",&k_c);
 printf("%lf  \n", om);

 mu = 0.05;
 f = 0.25;
 om1 = 0.2;
 eps = 0.2;
 k_c = 0.1;
 
//***time step***//
h=0.01; t=0.0;
x_max=  -LONG_MAX ; x_min=LONG_MAX;
cnt = 0;
 
 
for(k=0;k<= 5*M_PI; k+= 0.25 )
	{               
        t=h*(double)(k);
    	RK4(nn,n,h,t,y,DGL);    
   for(j=1;j<=nn;j++)
	{  
/*    leaving transients*/               
//	if(k>=4000)
//	{          
/*----------------max & min of the series--------*/
       if(y[j][1]>x_max){x_max=y[j][1];} 
       if(y[j][1]<x_min){x_min=y[j][1];}   
                                                       
      // pwr_i = cexp(fabs(z[k]));

             fprintf(fp1,"%lf     %lf  %lf \n", om, x_max, x_min);
          }	
	}  
   }
	     
  printf("process over!!\n");
}
//************************RK4 SUBROUTINE*********************************//
void RK4(int nn,int n,double h,double t,double y[101][501],
	   void (*DGL)(double,double[101][501],double[101][501]))
{
	   
	   double k1[101][501],k2[101][501],k3[101][501],k4[101][501];
	   double yaux[101][501];

	   DGL(t,y,k1);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k1[j][i]/2.0;
	   }
	   
	   DGL(t+h/2.0,yaux,k2);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k2[j][i]/2.0;
	   }
	   
	   DGL(t+h/2.0,yaux,k3);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k3[j][i];
	   }
	   
	   DGL(t+h,yaux,k4);
	   for(j=1;j<=nn;j++)
	   {
             for(i=1;i<=n;i++)
	      y[j][i]=y[j][i]+h*((k1[j][i]+2*k2[j][i]+2*k3[j][i]+k4[j][i])/6.0);
	   }
}
//*********************FUNCTION SUBROUTINE********************************//
void DGL(double t,double y[101][501],double F[101][501])
{

    int aa[300][300];  
      
      p=15;
       
   for(j=1;j<=nn;j++)
     { 
     
     for(i=1;i<=nn;i++)
     {
      
      if(i != j && abs(i-j)%(nn-p)<=p)    
         aa[j][i]=1;
          else
            aa[j][i]=0; 
   
     //  printf("%d\n" ,aa[j][i]);
 
       p_coef =3*(pow(om,4) - 4*om*om + 8*mu*mu)/(2*(3*om*om - 4*mu*mu)); 
       
       
 F[j][1]= y[j][2] -(k_c/((float)2.0*p))*aa[j][i]*(y[i][1] - y[j][1]);
 F[j][2]= -mu*y[j][1] - eps*( 1 + p_coef*sinc(y[j][1]))*cos(y[j][1]) + f*sin(om1*t) ;
     }
}

}

