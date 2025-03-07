/* Local coupling -parametric resonance in a Pendula */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>
#include<limits.h>

int i,j,k,nn=30,n=2;
double y[101][501];
int p;  //# of nearest neighbours
double k_c,om,om1, f, p_coef ,mu,x, eps;
double pwr, sum;
double complex pwr_i;
double p_c;
double del, ptr[200];
double x_max;
int num = 30;
int l;

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


/**************************************
double getmx( double pwr_spec[], int num)
  {
     int i; 
     double max;
       for(i =1;i<=num; i++){
         if(pwr_spec[i] > pwr_spec[i+1] &&  pwr_spec[i] > pwr_spec[i-1] )
               max = pwr_spec[i];
       }
              
    return max;
  }
/*********************************/

void main()
{
//nn= no 0f oscillators, n=dimension of the model//
double t,h,pi;
double z[100000];
srand(time(0));

FILE *fp1,*fp2,*fp3;
fp1=fopen("ress1_larc_ompc08.dat","w");

 
 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-4.0+2.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-4.0-2.0;
 }


x_max = -LONG_MAX; 
for (om  = 4.0;  om <= 6.0; om +=0.1)
 {
 
 for(i =1; i<= 1; i++)
  {
 //  printf("%f \n", p_coef);
   
//parameters
//scanf("k_c = %lf",&k_c);

p_c  =0.8;
p = p_c *nn;
 mu = 0.05;
 f = 0.2;
 om1 = 5.0;
 eps =0.2;
 k_c = 0.3;

 
//***time step***//
h=0.01; t=0.0;

for(k=1;k<=50000;k++)
	{               
        t=h*(double)(k);
    	RK4(nn,n,h,t,y,DGL);    
   for(j=1;j<=nn;j++)
	{  
/*    leaving transients*/               
	if(k>=4000)
	{                 
       sum =0.0;
       pwr = fabs(y[j][1])*fabs(y[j][1])*om/ (2.0*M_PI*(float) nn);
       if (pwr == NAN)
          pwr = 1.0;
       sum+= pwr;            
      // pwr_i = cexp(fabs(z[k]));
   
          }
	}    
}

     z[l] = pwr;
   // printf("%lf \n", z[l]);
    if(z[l] >= x_max){
      x_max  = z[l];}
      fprintf(fp1,"%20.12lf   %20.12lf\n", om, pwr);
      }    
}
      printf("x_max = %lf \n", x_max);
       
      del = x_max;

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


    double complex phase;
    double pert, del;
    int aa[300][300];  
      
       
   for(j=1;j<=nn;j++)
     { 
     
     for(i=1;i<=nn;i++)
     {
         if(i != j && abs(i-j)%(nn-p)<=p) 
           aa[j][i]=1;
          else
            aa[j][i]=0; 
 

       phase = cexp(sinc(y[j][1])*I);
       //  phase = cproj(sinc(y[j][1])*I);

       pert  = sinc(y[j][1]);


       if(cimag(phase) >= M_PI  || cimag(phase) <= -M_PI)
         pert    = sinc(M_PI);


       if (pwr >=  x_max/(float) 16.0)
         del  =   -del;
       else if( pwr < x_max/ (float) 16.0)
         del = +del;

 
       p_coef =3*(pow(om,4) - 4*om*om + 8*mu*mu)/(2*(3*om*om - 4*mu*mu)); 

 F[j][1]= y[j][2] -(k_c/(float)2*p)* aa[j][i]*(y[i][1] - y[j][1]);
 F[j][2]= -mu*y[j][1] - eps*( 1 + p_coef* sinc(y[j][1]))*sin(y[j][1]) + f*sin(om1*t) ;
     }
}

}

