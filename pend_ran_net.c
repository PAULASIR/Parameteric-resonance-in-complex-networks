/* random coupling -parametric resonance in a Pendula */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>

int i,j,k,nn= 20,n=2;
double y[101][201];
int arr[101][101];
int p;  //# of degerees
double p_c;
double k_c,om,om1, f, p_coef ,mu,x, eps;
double pwr;
double complex pwr_i;
void RK4(int,int,double,double,double[101][201],
                   void (*DGL)(double,double[101][201],double[101][201]));
void DGL(double, double[101][201],double[101][201]);


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

/*****************************
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


srand(time(0));

FILE *fp1,*fp2,*fp3;
fp1=fopen("ts1.dat","w");
fp2=fopen("chm.dat","w");
fp3=fopen("snp.dat","w");
 
 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-4.0+2.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-4.0-2.0;
 }

//printf("p_c =");
//scanf("%lf", &p_c);
//p =  nn* p_c;
//printf("%d   %d\n", nn, p);

 
//parameters
 mu = 0.05;
 f = 0.2;
 om1 = 1.5;
 om=1.4;
 eps =0.2;
 k_c = 0.25;
 p_c = 0.8;
// shuffle(arr);
// printarr(arr); 
 
 
//***time step***//
h=0.01; t=0.0;
x_max=x_min=0.0;

for(k=1;k<= 20000;k++)
	{               
        t=h*(double)(k);
    	RK4(nn,n,h,t,y,DGL);    
   for(j=1;j<=nn;j++)
	{  
/*    leaving transients*/               
	if(k>=4000)
	{          
/*----------------max & min of the series--------
       if(y[j][1]>x_max){x_max=y[j][1];} 
       if(y[j][1]<x_min){x_min=y[j][1];}   */
                                                       
      // pwr_i = cexp(fabs(z[k]));


 // printf("%f   %f  %f\n", pwr, cimag(pwr_i), creal(pwr_i));

/* -----------------Printing time series-----------------*/        
       if( k == 4950)
          fprintf(fp3,"%d %f \n",j,y[j][1]);
   
      for(i =1;i<=nn;i++)
        fprintf(fp2,"  %f  %f  \t", y[i][1], y[i][2]);
      fprintf(fp2," \n");

       fprintf(fp1,"%f\t",t);
         for(i=1;i<=nn;i++)
           fprintf(fp1,"%f\t",y[i][1]);    
           fprintf(fp1,"\n");
/**************************************************/    
          }
	}    
}

  printf("process over!!\n");
}
//************************RK4 SUBROUTINE*********************************//
void RK4(int nn,int n,double h,double t,double y[101][201],
	   void (*DGL)(double,double[101][201],double[101][201]))
{
	   
	   double k1[101][201],k2[101][201],k3[101][201],k4[101][201];
	   double yaux[101][201];

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
void DGL(double t,double y[101][201],double F[101][201])
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

