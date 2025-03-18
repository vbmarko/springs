#include <bits/pthreadtypes.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <math.h>

gsl_matrix *csv_to_matrix(char *path,char delim){
	FILE *f = fopen(path, "rb");
	int c, n = -1,m = 0, k = 0;
	double *buff = malloc(0*sizeof(0.0));
	double val = 0;
	while((c = getc(f)) != EOF){
		if(c != delim & c != '\n'){

			val = atof((char*) &c);
			buff = realloc(buff, (k+1)*sizeof(0.0));
			buff[k] = val;
			k = k+1;
		}
		if(c == '\n'){
			if(n == -1){
			n = k;}
			m = m +1;
		}

	}
	printf("\n\n n = %d, m = %d",n,m);
	fclose(f);
	gsl_matrix *M = gsl_matrix_calloc(m,n);
	for(int i = 0;i <n;++i){
		for(int j = 0;j <m;++j){
			gsl_matrix_set(M,j,i,buff[m*i+j]);
		}
	}
	free(buff);
	return M;

}



void compute_eigenvalues(gsl_vector *eval, gsl_matrix *evec,const gsl_matrix *A) {
    gsl_matrix *A_copy = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(A_copy, A);

    gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(A->size1);
    gsl_eigen_symmv(A_copy, eval, evec, workspace);

    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(A_copy);

    for (size_t i = 0; i < eval->size; i++) {
        printf("Eigenvalue %zu: %g\n", i, gsl_vector_get(eval, i));
    }
}


void binarize_matrix(gsl_matrix *out,const gsl_matrix *in,double thresh){
	int in1 = in->size1, in2 = in->size2, on1 = out->size1, on2 = out->size2;
	double val = 0;
	assert(in1 == on1); assert(in2 == on2);
	gsl_matrix_set_zero(out);
	for(int i = 0;i<in1;++i){
		for(int j = 0;j<in2;++j){
			val = gsl_matrix_get(in, i, j);
			if(fabs(val) > thresh){
				gsl_matrix_set(out, i, j,(val > 0) - (val < 0));				
			}
		}
	}
}

void propagator(gsl_matrix *P,double tau,const gsl_vector *eval, const gsl_matrix *evec){
	for(int i = 0;i<P->size1;++i){
		gsl_matrix_set(P,i,i,exp(-tau*gsl_vector_get(eval,i)));
	}
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1,P, evec,0,P);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,evec,P,0,P);

}

void compute_degree_matrix(gsl_matrix *D,const gsl_matrix *A) {
    size_t n = A->size1;

    // Initialize degree matrix to zeros
    gsl_matrix_set_zero(D);

    for (size_t i = 0; i < n; i++) {
        double degree = 0.0;

        for (size_t j = 0; j < n; j++) {
            degree += gsl_matrix_get(A, i, j);
        }

        gsl_matrix_set(D, i, i, degree); // Set diagonal entry to degree
    }
}

void compute_degree_vector(gsl_vector *d, const gsl_matrix *A) {
    size_t n = A->size1;


    for (size_t i = 0; i < n; i++) {
        double degree = 0.0;

        for (size_t j = 0; j < n; j++) {
            degree += gsl_matrix_get(A, i, j);
        }

        gsl_vector_set(d, i,degree); // Set diagonal entry to degree
    }
}


void compute_laplacian_matrix(gsl_matrix *L,const gsl_matrix *A){
	assert(L->size1 == L->size2);
	int n = L->size1;
	gsl_matrix *D = gsl_matrix_calloc(n,n);
	compute_degree_matrix(D,A);
	gsl_matrix_add(L,D);
	gsl_matrix_sub(L,A);
	gsl_matrix_free(D);
}



// Function to print a GSL matrix
void print_matrix(const gsl_matrix *m) {
    for (size_t i = 0; i < m->size1; i++) {
        for (size_t j = 0; j < m->size2; j++) {
            // Using %g to print the double in a compact format.
            printf("%.2g ", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }
}

// Function to print a GSL vector
void print_vector(const gsl_vector *v) {
    for (size_t i = 0; i < v->size; i++) {
        // Print each element followed by a space
        printf("%.2g ", gsl_vector_get(v, i));
    }
    printf("\n");
}



void rand_tree(gsl_matrix *A,int (*depth)(int n), double(*ro)(int d) ){
	int n = A->size1;
	int dep = depth(n);
	double p[dep];
	int n_d[dep];
	double suma = 0;
	int check = 0;
	for(int i = 0;i<dep;++i){
		p[i] = ro(i);
		suma += p[i];
	}
	for(int i = 0;i<dep;++i){
		p[i] = p[i]/suma;
		n_d[i] = p[i]*n;
		check += n_d[i];
	
	}
	if(check != n){
		printf("oops sum(n_depth) != n <-> %d != %d",check,n);
		printf("\nwe add difrence to leaves");
		n_d[dep-1] += n - check;
	}
	int n_above = n;
	for(int i = 0;i<dep;++i){
		if(i == dep-1){
			printf("depth 0");
			continue;
		}
		printf("\n----------------------------------\nhello from depth %d. There are %d of us here, we are: \n",dep-i-1,n_d[dep-1-i]);
		n_above -= n_d[dep-1-i];
		for(int j = 0;j<n_d[dep-1-i];++j){
			int my_idx = n_above+j;
			int pap_idx = n_above -1 -rand()%n_d[dep-2-i] ; 
			gsl_matrix_set(A,my_idx,pap_idx,1);
			gsl_matrix_set(A,pap_idx,my_idx,1);
			printf("node %d my pap is %d\n",my_idx,pap_idx);
		}
		printf("and there are %d above us \n",n_above);

	}

}
int depth(int n){
	return rand()%n + 1 ;
}
double ro(int d){
	return 1;
}

