#include <stdio.h>
#include <stdlib.h>
#include <math.h>

# define ROWS 40
# define COLS 40

int main(void){

double V_old[ROWS][COLS], V_new[ROWS][COLS], r[ROWS][COLS];
double E_x[ROWS][COLS],E_y[ROWS][COLS];
int i,j;
int count=0;
//int k;
//int n=300;
//
//over-relaxation parameter
double w=2/(1+3.14/ROWS); 
double r_max=2.0;
double E_max=0;
double conv_crit=0.1;


/*initialize grids*/
for(i=0;i<ROWS;i++){
	for(j=0;j<COLS;j++){
		V_old[i][j]=0;
		V_new[i][j]=0;
		E_x[i][j]=0;
		E_y[i][j]=0;
	}
}

/*initialize geometry*/
for(i=15;i<25;i++){
	for(j=20;j<30;j++){
		V_old[i][j]=100;
		V_new[i][j]=100;
		E_x[i][j]=0;
		E_y[i][j]=0;
	}
}

for(i=10;i<30;i++){
	V_old[i][10]=-100;
	V_new[i][10]=-100;
	E_x[i][10]=0;
	E_y[i][10]=0;
}

/*carry out relaxation until convergence is satisfied and avoiding the source points*/
while(r_max>conv_crit){
for(i=1;i<ROWS-1;i++){
	for(j=1;j<COLS-1;j++){
		//Jacobi
		V_new[i][j]=(V_old[i+1][j]+V_old[i-1][j]+V_old[i][j+1]+V_old[i][j-1])/4;		
	
		//Gauss-Seidel
		//V_new[i][j]=(V_old[i+1][j]+V_new[i-1][j]+V_old[i][j+1]+V_new[i][j-1])/4;		
	
		r[i][j]=V_new[i][j]-V_old[i][j];
	
		//Successive Over-Relaxation
		//V_new[i][j]=V_old[i][j]+w*r[i][j];
					
		V_old[i][j]=V_new[i][j];
		}
	}
for(i=15;i<25;i++){
	for(j=20;j<30;j++){
        	V_old[i][j]=100;
        	V_new[i][j]=100;
	}
}
for(i=10;i<30;i++){
	V_old[i][10]=-100;
	V_new[i][10]=-100;
}
//resetting and finding r_max
r_max=0;
for(i=1;i<ROWS-1;i++){
	for(j=1;j<COLS-1;j++){
		if(r[i][j]>r_max){
			r_max=r[i][j];
		}
	}
}
//count++;
}

/*Calculate Electric field components and find biggest*/
for(i=1;i<ROWS-1;i++){
	for(j=1;j<COLS;j++){
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
for(i=15;i<25;i++){
	for(j=20;j<30;j++){
    		E_x[i][j]=0;
        	E_y[i][j]=0;
	}
}
for(i=10;i<30;i++){
	E_x[i][10]=0;
	E_y[i][10]=0;
}

/*print voltage results in grid form*/
/*
for(i=0;i<ROWS;i++){
	for(j=0;j<COLS;j++){
		printf("%f\t",V_new[i][j]);
	}
	puts("");
}
*/

/*print results in (x,y,V,Ex,Ey) list form*/

for(i=0;i<ROWS;i++){
	for(j=0;j<COLS;j++){
		printf("%d %d %f %f %f\n",i,j,V_new[i][j],E_x[i][j]/E_max,E_y[i][j]/E_max);
	}
}

/*print iteration data*/
/*
printf("%d",count);
*/
return 0;
}
