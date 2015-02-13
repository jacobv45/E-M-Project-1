#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(void){

int N=100;
double V_old[N][N], V_new[N][N], r[N][N];
double E_x[N][N],E_y[N][N];
int i,j;
int count=0;
//int n=300;
//
//over-relaxation parameter
double w=1.94; 
double r_max=0.;
double E_max=0.;
double conv_crit=0.00001;


/*initialize grids*/
for(i=0;i<N;i++){
	for(j=0;j<N;j++){
		V_old[i][j]=0;
		V_new[i][j]=0;
		E_x[i][j]=0;
		E_y[i][j]=0;
	}
}

/*initialize geometry*/
for(i=N/4;i<3*N/4;i++){
	V_old[i][N/4]=-100;
	V_new[i][N/4]=-100;
	V_old[i][3*N/4]=100;
	V_new[i][3*N/4]=100;
}


/*carry out relaxation until convergence is satisfied and avoiding the source points*/
do{	
	for(i=1;i<N-1;i++){
		for(j=1;j<N-1;j++){
			//Jacobi
			//V_new[i][j]=(V_old[i+1][j]+V_old[i-1][j]+V_old[i][j+1]+V_old[i][j-1])/4;		
	
			//Gauss-Seidel
			V_new[i][j]=(V_old[i+1][j]+V_new[i-1][j]+V_old[i][j+1]+V_new[i][j-1])/4;		
		
			
			r[i][j]=V_new[i][j]-V_old[i][j];
			
			//Successive Over-Relaxation
	                V_new[i][j]=V_old[i][j]+w*r[i][j];
		}
	}
	//reset geometry
	for(i=N/4;i<3*N/4;i++){
        	V_old[i][N/4]=-100;
        	V_new[i][N/4]=-100;
		V_old[i][3*N/4]=100;
		V_new[i][3*N/4]=100;
	}
	
	//V_new[i][j]=V_old[i][j]+w*r[i][j];


	//resetting and finding r_max
	r_max=0;
	for(i=1;i<N-1;i++){
		for(j=1;j<N-1;j++){
			r[i][j]=V_new[i][j]-V_old[i][j];
	
			if(r[i][j]>r_max){
				r_max=r[i][j];
			}
			
			//copy grid over for new loop
			V_old[i][j]=V_new[i][j];
		}
	}
count++;
}while(r_max>conv_crit);

/*Calculate Electric field components and find biggest*/
for(i=1;i<N-1;i++){
	for(j=1;j<N-1;j++){
		E_x[i][j]=(V_new[i-1][j]-V_new[i+1][j])/2;
		E_y[i][j]=(V_new[i][j-1]-V_new[i][j+1])/2;
		if(E_x[i][j]>E_max){
			E_max=E_x[i][j];
		}	
		if(E_y[i][j]>E_max){
			E_max=E_y[i][j];
		}
	}
}
for(i=N/4;i<3*N/4;i++){
    		E_x[i][N/4]=0;
        	E_y[i][N/4]=0;
		E_x[i][3*N/4]=0;
		E_y[i][3*N/4]=0;
}


/*print voltage results in grid form*/
/*
for(i=0;i<N;i++){
	for(j=0;j<N;j++){
		printf("%f\t",V_new[i][j]);
	}
	puts("");
}
*/

/*print results in (x,y,V,Ex,Ey) list form*/
/*
for(i=0;i<N;i++){
	for(j=0;j<N;j++){
		printf("%d %d %f %f %f\n",i,j,V_new[i][j],E_x[i][j]/E_max,E_y[i][j]/E_max);
	}
}
*/

/*print iteration data*/

printf("%d\n",count);

return 0;
}
