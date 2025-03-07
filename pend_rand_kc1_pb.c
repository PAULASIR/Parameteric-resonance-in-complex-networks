/* Local coupling -parametric resonance in a Pendula  --  Resonance Probability*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>
#include<limits.h>

int i,j,k,nn=20,n=2;
double y[101][501];
int arr[101][101];
int p;  //# of nearest neighbours
double k_c,om,om1, f, p_coef ,mu,x, eps;
double pwr, sum, del;
double complex pwr_i;
double pc;
double z[10000];

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

double fpeaks(double z[])
{
    double peak;
    for( int i = 1; i <= 21; i++)
       {
           if(z[i] >= z[i-1] && z[i] >= z[i+1])
           peak = z[i];     
      }
return peak;
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
double amplitude,vmin,vmax,c_max,c_min;
double x_max,x_min;
int ll, cnt, mean, sum_cnt;
double avg;
srand(time(0));

FILE *fp1,*fp2,*fp3;
fp1=fopen("ran_pb_rand_kc1.dat","w");

 
 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-4.0+2.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-4.0-2.0;
 }

for(k_c = 0.0;  k_c <= 2.0; k_c += 0.05)
   {
//parameters
 //scanf("k_c = %lf",&k_c);
 printf("%lf \n", k_c);
 sum_cnt =0;
for(ll =1; ll <= 10; ll++){
     cnt = 0;
       x_max = LONG_MIN;  x_min = LONG_MAX;
for(om = 4.0; om <= 6.0; om += 0.1)
  {
   // printf("omega =%lf \n" , om );
    pc = 0.2;

   p  = pc * (float)(nn);

 mu = 0.05;
 f = 0.2;
// om = 5.0;
  om1 = 5.0;
// k_c = 1.5;
 eps = 0.2;
 
//***time step***//
h=0.01; t=0.0;

for(k=1;k<=20000;k++)
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
          }
	}    
}
          z[i] = pwr;
          
//        printf("%lf \n", z[i]);
        
        if( z[i] > x_max)
            x_max = z[i];
            if( z[i] < x_min)
            x_min = z[i];
         
          avg = (x_max + x_min)/ (float)2.0;  
        
        if(z[i] >= avg) 
            {
		cnt = cnt + 1;  }
         
        }
             sum_cnt +=cnt;
      }
  
       mean = sum_cnt/  ll;
     
      printf("%d \n", cnt);
      fprintf(fp1,"%20.12lf    %d \n",p_coef, mean);
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
 
 F[j][1]= y[j][2] - k_c*arr[j][i]*(y[i][1] - y[j][1]);
 F[j][2]= -mu*y[j][1] - eps*( 1 + p_coef*sinc(y[j][1]))*sin(y[j][1]) + f*sin(om1*t) ;    
     }
}

}

