#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "linear_algebra.h"

#define PI 4*atan(1) //a precise calculation of Pi


void vector_add(double a, struct vector *x, double b, struct vector *y, struct vector *z){
	//calculates z = a * x + b * y, for the vectors x, y, z and scalars a and b
	//only takes two vectors at a time, so several commands may appear together
	assert(x->n == y->n); // cute little test: can only add vectors of equal length
	assert(x->n == z->n);
	for (int i = 0; i < x->n; i++) {
		VEC(z, i) = a * VEC(x, i) + b * VEC(y, i);
	}
}

void rhs(struct vector *x, struct vector *v, double dx){
	//calculates the spatial finite central differences for uu_x
	//takes a vector x and applies the finite differences to it
	//these are returned in v, while x is preserved
	//hard codes the periodic points--self explainatory 
	int n=v->n;
	int i=0;
	VEC(v, i) = -VEC(x,i) * (VEC(x,i+1) - VEC(x,n-1)) / 2 / dx;
	i=1;
	VEC(v, i) = -VEC(x,i) * (VEC(x,i+1) - VEC(x,i-1)) / 2 / dx;
	for (i = 2; i < n - 2; i++){
		VEC(v, i) = -VEC(x,i) * (VEC(x,i+1) - VEC(x,i-1)) / 2 / dx;
	}
	i = n-2;	
	VEC(v, i) = -VEC(x,i) * (VEC(x,i+1) - VEC(x,i-1)) / 2 / dx;
	i = n-1;
  	VEC(v, i) = -VEC(x,i) * (VEC(x,i+1-n) - VEC(x,i-1)) / 2 / dx;
}










//One Step Adams-Bashforth, i.e., Backward Euler
void ab_1(struct vector *ic, double dx, double dt, int N, int J){
	int k;
	struct vector *next  = vector_create(J);
	struct vector *dummy = vector_create(J);
	
	FILE *fout3 = NULL; 
	fout3 = fopen("results_ab1.txt", "w");	

	//Print IC
	for(k=0; k<J-1; k++){
		fprintf(fout3, "%.17le\t", ic->vals[k]);
	}
	fprintf(fout3,"\n");
	
	//STEP ONE: INITIALISE AB1/EULER
	rhs(ic, dummy, dx);
	vector_add(1., ic, dt, dummy, next);

	//Print first step
	for(k=0; k<J-1; k++){
		fprintf(fout3, "%.17le\t", next->vals[k]);
	}
	fprintf(fout3,"\n");


	//STEP TWO-->ONWARDS
	for(int j=0; j < N-1; j++){
		rhs(next, dummy, dx);
		vector_add(1., next, dt, dummy, next);
		//Print results:
		if (j%20==0){ 
		for(k=0; k<J-1; k++){
			fprintf(fout3, "%.17le\t", next->vals[k]);
		}
		fprintf(fout3,"\n");
		}
	}
	fclose(fout3);
}





//Two Step Adams-Bashforth
void ab_2(struct vector *ic, double dx, double dt, int N, int J){
	int k;
	struct vector *next     = vector_create(J);
	struct vector *current  = vector_create(J);
	struct vector *previous = vector_create(J);
	struct vector *dummy    = vector_create(J);	

	FILE *fout4 = NULL; 
	fout4 = fopen("results_ab2.txt", "w");	

	//Print IC		
	for(k=0; k < J-1; k++){
		fprintf(fout4, "%.17le\t", ic->vals[k]);
	}
	fprintf(fout4,"\n");
	

	//STEP ONE: WE MUST FIND `PREVIOUS' -- USE AB1/EULER for this
	rhs(ic, dummy, dx);
	vector_add(1., ic, dt, dummy, next);

	//Print first result		
	for(k=0; k < J-1; k++){
		fprintf(fout4, "%.17le\t", next->vals[k]);
	}
	fprintf(fout4,"\n");


	//STEP TWO: INITIALISE AB2
	rhs(next, current, dx);
	vector_add(3., current, -1., next, dummy);
	vector_add(1., next, 0.5 * dt, dummy, next);
	
	//Print SECOND result		
	for(k=0; k < J-1; k++){
		fprintf(fout4, "%.17le\t", next->vals[k]);
	}
	fprintf(fout4,"\n");


	//STEP THREE-->ONWARDS TO STEP N-2
	for(int j=0; j < N-3; j++){

		for(int i=0; i < J; i++){
			previous->vals[i] = current->vals[i];
		}

		rhs(next, current, dx);
		vector_add(3., current, -1., previous, dummy);
		vector_add(1., next, 0.5 * dt, dummy, next);
		
		//Print results
		if (j%20==0){ 
		for(k=0; k<J-1; k++){
			fprintf(fout4, "%.17le\t", next->vals[k]);
		}
		fprintf(fout4,"\n");
		}
	}

	fclose(fout4);
}





//Three Step Adams-Bashforth
void ab_3(struct vector *ic, double dx, double dt, int N, int J){
	int k;
	int i;
	struct vector *next         = vector_create(J);
	struct vector *current      = vector_create(J);
	struct vector *previous     = vector_create(J);
	struct vector *uberprevious = vector_create(J);
	struct vector *dummy        = vector_create(J);	

	FILE *fout5 = NULL; 
	fout5 = fopen("results_ab3.txt", "w");	

	
	//Print IC
	for(k=0; k<J-1; k++){
		fprintf(fout5, "%.17le\t", ic->vals[k]);
	}
	fprintf(fout5,"\n");


	//STEP ONE: WE MUST FIND UBERPREVIOUS -- USE AB1/EULER
	rhs(ic, dummy, dx);
	vector_add(1., ic, dt, dummy, uberprevious);

	//Print first result		
	for(k=0; k < J-1; k++){
		fprintf(fout5, "%.17le\t", uberprevious->vals[k]);
	}
	fprintf(fout5,"\n");


	//STEP TWO: WE MUST FIND PREVIOUS -- USE AB2
	rhs(uberprevious, current, dx);
	vector_add(3., current, -1., uberprevious, dummy);
	vector_add(1., uberprevious, 0.5 * dt, dummy, previous);

	//Print SECOND result		
	for(k=0; k < J-1; k++){
		fprintf(fout5, "%.17le\t", previous->vals[k]);
	}
	fprintf(fout5,"\n");


	//STEP THREE: INITIALISE AB3
	rhs(previous, current, dx);
	vector_add(23., current, -16., previous, dummy);
	vector_add(1., dummy, 5., uberprevious, next);
	vector_add(1., previous, dt/12., next, next);

	//Print THIRD result
	for(k=0; k<J-1; k++){
		fprintf(fout5, "%.17le\t", next->vals[k]);
	}
	fprintf(fout5,"\n");


	//STEP FOUR-->ONWARDS 
	for(int j=0; j < N-4; j++){

		for(i=0; i < J; i++){
			uberprevious->vals[i] = previous->vals[i];
			previous->vals[i] = current->vals[i];
		}

		rhs(next,current,dx);
		vector_add(23., current, -16., previous, dummy);	
		vector_add(1., dummy, 5., uberprevious, dummy);
		vector_add(1., next, dt/12., dummy, next);
	
		//Print results
		if (j%20==0){ 
			for(k=0; k<J-1; k++){
				fprintf(fout5, "%.17le\t", next->vals[k]);
			}
			fprintf(fout5,"\n");
		}
	}
	fclose(fout5);
	
}




void print_array(struct vector *t){
	FILE *fout2 = NULL;									 ////
	fout2 = fopen("results_timesteps_xaxis.csv", "w");   ////
	for(i=0; i<t->n;i++){								 ////
		fprintf(fout1, "%.17lg\n", t->vals[i]);			 ////
	}													 ////
	fclose(fout2);	
}



void burgers(double InitialCondition(double x), double x0, double xf, double dt, int N, int J){
	int i;
	//spatial step
	double dx = (xf - x0) / (J - 1.);
	struct vector *x = vector_create(J);
	for (i = 0; i < x->n; i++){
		x->vals[i] = i * dx;
	}

	////This prints spacesteps to be used for x-axis plotting////
		FILE *fout1 = NULL;									 ////
		fout1 = fopen("results_spacesteps_xaxis.txt", "w");	 ////
		for(i=0; i<x->n;i++){								 ////
			fprintf(fout1, "%.17lg\n", x->vals[i]);			 ////
		}													 ////
		fclose(fout1);										 ////
	/////////////////////////////////////////////////////////////

	//temporal step
	struct vector *t = vector_create(N);
	for (i = 0; i < N; i++){
		t->vals[i] = i * dt;
	}
	//print_array(t);


	//initial condition @t=0
	struct vector *ic = vector_create(J);
	for (i = 0; i < J; i++){
		ic->vals[i] = InitialCondition(x->vals[i]);	
	}

	ab_1(ic, dx, dt, N, J);
	ab_2(ic, dx, dt, N, J);
	ab_3(ic, dx, dt, N, J);
}









double InitialCondition(double x){
	return sin(2 * PI * x);
}

int main(){
	double x0 = 0.;			//initial space
	double xf = 1.;			//final space
	double dt = 0.0001;		//time step size
	int     N = 9000;		//time steps
	int 	J = 1000;		//space steps
	burgers(InitialCondition, x0, xf, dt, N, J);
}











