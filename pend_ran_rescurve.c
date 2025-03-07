/* random coupling -parametric resonance in a Pendula */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>

int i,j,k,nn= 30,n=2;
double y[101][501];
double arr[101][101];
int p;  //# of degerees
double p_c;
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


/*****************************/
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

/******************************/
void printarr(int arr[nn][nn]){
int i,j;
for(j =0 ; j<nn; j++){
     for(i =0; i< nn; i++)
     {
         printf("%2d  ", arr[j][i]);     
     }
        printf("\n");
 }
}
/*******************************/
void main()
{
//nn= no 0f oscillators, n=dimension of the model//
double t,h,pi;
double z[100000];
double amplitude,vmin,vmax,c_max,c_min;
double x_max,x_min;
double sum;
srand(time(0));

FILE *fp1,*fp2,*fp3;
fp1=fopen("ran_ares_curve1.dat","w");

 
 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-4.0+2.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-4.0-2.0;
 }


for(om = 0.5; om <= 2.1; om += 0.1)
 {
   printf("%f\n",om);
/*printf("p_c =");
scanf("%lf", &p_c);
printf("%d   %d\n", nn, p); */


//parameters
 mu = 0.05;
 f = 0.2;
 om1 = 1.5;
 eps =0.2;
 k_c = 1.5;
 p_c = 0.8;
 
 p =  nn* p_c;

//***time step***//
h=0.01; t=0.0;
x_max=x_min=0.0;

for(k=1;k<= 30000;k++)
	{               
        t=h*(double)(k);
    	RK4(nn,n,h,t,y,DGL);    
   for(j=1;j<=nn;j++)
	{  
/*    leaving transients*/               
	if(k>=4000)
	 {                 
           sum =0.0;
           pwr = (fabs(y[j][1])*fabs(y[j][1])*om)/ (2.0*M_PI*(float) nn);
           if (pwr == NAN)
           pwr = 1.0;
           sum+= pwr;  
          }
	       }
	} 
	fprintf(fp1,"%20.12lf   %20.12lf\n", om , pwr);
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

   double complex phase;
   double pert, del;

   for(j=1;j<=nn;j++)
     { 
     for(i=1;i<=nn;i++)
       {
       p_coef =3*(pow(om,4) - 4*om*om + 8*mu*mu)/(2*(3*om*om - 4*mu*mu)); 
       
          phase = cexp(sinc(y[j][1])*I);      
    //  phase = cproj(sinc(y[j][1])*I);

     pert  = sinc(y[j][1]);
       
       if(cimag(phase) >= M_PI  || cimag(phase) <= -M_PI)
               pert    = sinc(M_PI);
      

       /*  if (pwr >=  x_max/(float) 16.0)
           del  =   -del;
         else if( pwr < x_max/ (float) 16.0)
           del = +del; */

 F[j][1]= y[j][2] - k_c*arr[j][i]*(y[i][1] - y[j][1]);
 F[j][2]= -mu*y[j][1] - eps*( 1 + p_coef)*sin(y[j][1]) + f*sin(om1*t) ;
     }
 }
}

