#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*solves the Poisson Equation for a grounded NxNxN cube with a point charge in the center (only works with even N obviously)*/
int main(void){

int N=10;
int i,j,k;
int count=0;

double h=1;
double V_old[N][N][N],V_new[N][N][N],rho[N][N][N],V_real[N][N][N];
double r[N][N][N];
double E_x[N][N], E_y[N][N];

double w=1.8;
double E_max=0;
double r_max=2;
double conv_crit=0.000001;

/*initialize the grids*/
for(i=0;i<N;i++){
	for(j=0;j<N;j++){
		for(k=0;k<N;k++){
			V_old[i][j][k]=0;
			V_new[i][j][k]=0;
			rho[i][j][k]=0;
			V_real[i][j][k]=0;
		}
	}
}

/*set up the charge density*/
rho[N/2][N/2][N/2]=1.0;

/*relaxation loop*/
do{
for(i=1;i<N-1;i++){
         for(j=1;j<N-1;j++){
                 for(k=1;k<N-1;k++){
			
			//Gauss-Seidel Method
			V_new[i][j][k]=(V_new[i-1][j][k]+V_old[i+1][j][k]+V_new[i][j-1][k]+V_old[i][j+1][k]+V_new[i][j][k-1]+V_old[i][j][k+1])/6+rho[i][j][k]*h*h/6;    
			r[i][j][k]=V_new[i][j][k]-V_old[i][j][k];
			
			//Successive Over-relaxation (SOR)
			V_new[i][j][k]=V_old[i][j][k]+w*r[i][j][k];		
			}
		}
	}
	//reset r_max and find new r_max
//printf("%f\n",r_max); /*debugging*/
r_max=0;
for(i=1;i<N-1;i++){
	for(j=1;j<N-1;j++){
		for(k=1;k<N-1;k++){
			r[i][j][k]=V_new[i][j][k]-V_old[i][j][k];

			if(r[i][j][k]>r_max){
				r_max=r[i][j][k];
				}
			
			V_old[i][j][k]=V_new[i][j][k];	
			}
		}
	}
count++;
}while(r_max>conv_crit);

/*calculate electric field in x-y plane (z=N/2)*/
for(i=1;i<N-1;i++){
	for(j=1;j<N-1;j++){
		E_x[i][j]=(V_new[i-1][j][N/2]-V_new[i+1][j][N/2])/2/h;
		E_y[i][j]=(V_new[i][j-1][N/2]-V_new[i][j+1][N/2])/2/h;
		if(E_x[i][j]>E_max){
			E_max=E_x[i][j];
		}
	}
}

/*calculate the actual (known) potential*/
for(i=0;i<N;i++){
	for(j=0;j<N;j++){
		for(k=0;k<N;k++){
			V_real[i][j][k]=(1/(h*sqrt((i-N/2)*(i-N/2)+(j-N/2)*(j-N/2)+(k-N/2)*(k-N/2))))/12.566;
		}
	}	
}


/*print results (x,y,z,V)*/
for(i=0;i<N;i++){
	for(j=0;j<N;j++){
		for(k=0;k<N;k++){
	//		printf("%d %d %d %f\n",i,j,k,V_new[i][j][k]);
		}
	}
}

/*prints results (x,y,z=0,V)*/

for(i=0;i<N;i++){
	for(j=0;j<N;j++){
//		printf("%d %d %f\n",i,j,V_new[i][j][N/2]);
	}
}


/*print results (x[0:N],V(x,N/2,N/2),V_real(x,N/2,N/2))*/
for(i=N/2;i<N;i++){
//	printf("%d %f %f\n",i,V_new[i][N/2][N/2],V_real[i][N/2][N/2]);
}


/*print results (x,y,E_x,E_y)*/
for(i=0;i<N;i++){
	for(j=0;j<N;j++){
		printf("%d %d %f %f\n",i,j,E_x[i][j]/E_max,E_y[i][j]/E_max);
	}
}


/*print iteration data*/
//printf("%d\n",count);


return 0;
}
