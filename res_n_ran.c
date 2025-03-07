/* Local coupling -parametric resonance in a Pendula */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>

int i,j,k,n=2;
int nn;
double y[101][501];
int p;  //# of nearest neighbours
double k_c,om,om1, f, p_coef ,mu,x, eps;
double pwr, sum;
double arr[200][200];
double p_c;
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

/***********************/
int shuffle(int arr[nn][nn])
{
  int temp;
  int l,m;

  for (int j = 0; j < nn; j++)
   {
      for (int i = 0; i < nn; i++)
     {
       
        if((j < p) &&(i < p))
       {  arr[j][i] =1;}
        else
        {arr[j][i] = 0;} 
     
        l  = rand() % (nn-1);
        m = rand() % (nn-1);

      // Swap the elements
      temp = arr[j][i];
      arr[j][i] = arr[l][m];
      arr[l][m] = temp;
    }
  }
}
/***************************/





void main()
{
//nn= no 0f oscillators, n=dimension of the model//
double t,h,pi;
double z[100000];
double amplitude,vmin,vmax,c_max,c_min;
double x_max,x_min;
int ll;
srand(time(0));

FILE *fp1,*fp2,*fp3;
fp1=fopen("ran_n_pwr1.dat","w");

 
 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-4.0+2.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-4.0-2.0;
 }


for(nn = 2; nn <= 100; nn += 2){
printf("%d\n",nn);
for(p_c =0.0;  p_c <=1.0; p_c +=0.1){
p = p_c*nn;
double sm1 = 0;  
for(ll =0; ll <= 10; ll++){
x_max=x_min=0.0; 
for(om = 4.0; om <= 6.0; om += 0.1)
  {
//parameters
//scanf("k_c = %lf",&k_c);

 mu = 0.05;
 f = 0.2;
// om = 1.0;
 om1 = 5.0;
 eps =0.2;
 k_c = 0.4;
 
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
                    // printf("%f   %20.30lf\n", y[j][1], pwr);                                  
      // pwr_i = cexp(fabs(z[k]));
   
          }
	}    
}

     z[i] = sum;
	  if (z[i] > x_max)
	      x_max = z[i];
	
     }
     
     sm1 += x_max;
       printf("%f \n", sm1);
        }
       double  mean = sm1/  ll;
         printf("%f  %d %d\n", mean, p, ll);
        fprintf(fp1,"%d     %lf  %20.12lf\n", nn , p_c, mean);
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
       
   for(j=1;j<=nn;j++)
     { 
     
     for(i=1;i<=nn;i++)
     {
        
       p_coef =3*(pow(om,4) - 4*om*om + 8*mu*mu)/(2*(3*om*om - 4*mu*mu)); 

 F[j][1]= y[j][2] -k_c* arr[j][i]*(y[i][1] - y[j][1]);
 F[j][2]= -mu*y[j][1] - eps*( 1 + p_coef*sinc(y[j][1]))*sin(y[j][1]) + f*sin(om1*t) ;
     
     }
}

}

