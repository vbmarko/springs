#pragma once
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void print_matrix(const gsl_matrix *m);
void print_vector(const gsl_vector *v);
void rand_tree(gsl_matrix *A,int (*depth)(int n), double(*ro)(int d));
int depth(int n);
double ro(int d);
void compute_degree_matrix(gsl_matrix *D,gsl_matrix *A) ;
void compute_laplacian_matrix(gsl_matrix *L,const gsl_matrix *A);
void compute_degree_vector(gsl_vector *d,const gsl_matrix *A);
void compute_eigenvalues( gsl_vector *eval, gsl_matrix *evec,const gsl_matrix *A);
void binarize_matrix(gsl_matrix *out,const gsl_matrix *in,double thresh);
void matrix_to_file(const gsl_matrix *m,char *fpath,char delim,char *fmt);
void vector_to_file(const gsl_vector *v,char *fpath,char delim,char *fmt);
gsl_matrix *csv_to_matrix(char *path,char delim);




