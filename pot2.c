#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>



void main()
{

complex double y;
double theta , mu, A,beta, omega, P, k;

mu = 0.05;
beta = 0.2;
A = 0.2;
k =0.2;
P = 0.1;
omega = 0.1;

  FILE *fp1;
  fp1 = fopen("2.dat", "w" );


for(theta = -4*M_PI;  theta <= 4*M_PI; theta += 0.05){
      
    //  y =  - mu * theta*theta/ 2 - (1 - (beta*cos(theta)*theta*log(theta))*sin(theta) - A * (cos(omega*theta))*omega + (k *theta*theta) /4*P) ;
      
     y =  3.0*(pow(x,2))/(2*x) - (0.6*pow(x,4))/(12*x)+3.0*(pow(x,2)/2) - 0.6*(pow(x,4)/12) + 2*cos(2*x); 
     fprintf(fp1,"%f  %20.20lf  %20.20lf\n",theta, cimag(y), creal(y));
    }
}
