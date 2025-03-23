#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "raylib.h"
#include "gsl_funs.h"


#define dt 0.1



typedef struct point point ;
typedef struct edge edge;
typedef struct Cursor Cursor;
typedef struct springs springs;
typedef struct plane plane;
typedef struct pix_homotopy pix_homotopy ;


double dist(double x0 ,double y0, double x1, double y1){
	return fabs(x0-x1) + fabs(y0-y1);

}


int dist_e(int x0 ,int y0, int x1, int y1){
	int dx = x0 -x1, dy = y0 - y1;
	return round(sqrt(dx*dx + dy*dy));
}

struct pix_homotopy{
	int x00,y00, x01, y01;
	int x10,y10, x11, y11;
	Color c;
};

int max(int a, int b){
	return a*(a>b) + b*(b>=a);
}
int min(int a, int b){
	return b*(a>b) + a*(b>=a);
}

double lintrp(double a, double b, double t){
	return (1-t)*a + t*b;
}


void upd_homotopy_pix(pix_homotopy *f, Color *pix,int w ,int h, double (*intrp)(double a, double b,double t)){
	int d = max(dist(f->x00,f->y00,f->x10,f->y10),dist(f->x01,f->y01,f->x11,f->y11));
	//printf("d = %d\n",d);
	double x0,y0,x1,y1,t_d;
	int x,y;
	double t_s, t_t;
	for(int s = 0;s<=d;++s){
		t_s = 1.0*s/d;
		x0 = intrp(f->x00,f->x10,t_s);
		y0 = intrp(f->y00,f->y10,t_s);
		x1 = intrp(f->x01,f->x11,t_s);
		y1 = intrp(f->y01,f->y11,t_s);
		t_d = dist(x0, y0, x1,y1);
		
		//printf("x0 = %d y0 = %d, x1 = %d, y1 = %d \n",x0,y0,x1,y1);
		//printf("s = %d d = %d, t_d = %d \n",s,d,t_d);
		for(int t = 0;t <= t_d;++t){
			t_t = 1.0*t/t_d;
			x = round(intrp(x0, x1, t_t));
			y = round(intrp(y0, y1, t_t));
			//printf("(x,y) = (%d,%d)\n------------------------\n",x,y);
			if(x < 0 | x >= w | y < 0 | y >= h){
					continue;
				}
			pix[w*y +x] = (f->c);

		}
	}
		
}




struct Cursor{
	gsl_vector *x_2d;
	gsl_vector *r;
	int h,w;
};

struct point{
	gsl_vector *x;
	gsl_vector *v;
	gsl_vector *n;
	gsl_vector *x_2d;
	int int_x, int_y;
	double m,r,q;
};
struct edge{
	point *p_i, *p_j;
	double k,l,beta;
	int id, s, next_s;
	Color c;
};

struct springs{
	point *p;
	edge *e;
	gsl_vector *mx;
	int n_p, n_e;
	gsl_vector *scrtch;

};

struct plane{
	gsl_vector *r;
	gsl_vector *P_x;
	gsl_vector *P_y;
	gsl_vector *n;
	gsl_vector *rot_ax;
	gsl_matrix *rot_mat;
	gsl_vector *scrtch;
};


void free_plane(plane *Pi){
	gsl_vector_free(Pi->r);
	gsl_vector_free(Pi->P_x);
	gsl_vector_free(Pi->P_y);
	gsl_vector_free(Pi->n);
	gsl_vector_free(Pi->rot_ax);
	gsl_vector_free(Pi->scrtch);
	gsl_matrix_free(Pi->rot_mat);

}



void init_point(point *p,int n,double m,double r){
	gsl_vector *x = gsl_vector_calloc(n);
	for(int i = 0 ;i<n;++i){
		gsl_vector_set(x,i,1.0*rand()/RAND_MAX);
	}
	gsl_vector *v = gsl_vector_calloc(n);
	gsl_vector *normal = gsl_vector_calloc(n);
	//gsl_vector_set(normal, 0, 1);
	gsl_vector_set(normal, n-1, 1);
	double nn = gsl_blas_dnrm2(normal);
	if(nn){
	gsl_vector_scale(normal, 1/(nn));};
	double xn;
	gsl_blas_ddot(x, normal, &xn);
	gsl_vector_axpby(-xn,normal , 1, x);
	gsl_vector *x_2d = gsl_vector_calloc(2);
	p->x = x; p->x_2d = x_2d; p->v = v; p->n = normal;
	p->m = m; p->r = r, p->q = 10;
}

void free_point(point *p){
	gsl_vector_free(p->x);
	gsl_vector_free(p->x_2d);
	gsl_vector_free(p->v);
	gsl_vector_free(p->n);
}






springs *init_springs_from_adjacency(gsl_matrix *A,gsl_vector *m,gsl_vector *r,int d){
	springs *s = malloc(sizeof(springs));
	int n = A->size1;
	point *p = malloc(n*sizeof(point));
	int n_e = 0;
	for(int i = 0; i<n;++i){
		printf("%d \n",i);
		init_point(p+i,d, gsl_vector_get(m,i), gsl_vector_get(r,i));
		for(int j=i+1;j<n;++j){
			n_e += gsl_matrix_get(A,i,j);
		}
	}
	edge *e = malloc(n_e*sizeof(edge));
	int e_id = 0;
	s->n_p = n; s->n_e = n_e;
	for(int i = 0; i<n;++i){
		for(int j=i+1;j<n;++j){
			double l = gsl_matrix_get(A,i,j);
			if(l > 0){
				(e+e_id)->p_i = p + i;
				(e+e_id)->p_j = p + j;
				(e+e_id)->l = l;
				(e+e_id)->k = 0.1;
				(e+e_id)->beta = 1;
				(e+e_id)->id = e_id;
				(e+e_id)->s = 2;
				(e+e_id)->c = BLACK;
				e_id += 1;

			}
		}
	}
	s->scrtch = gsl_vector_calloc(d);
	s->p = p; s->e = e;
	return s;
}

void free_springs(springs *s){
	for(int i = 0;i<s->n_p;++i){
		free_point(s->p + i);
	}
	free(s->p); 
	free(s->e);
}

void electric(point *p1,point *p2,gsl_vector *scrtch){
	gsl_vector_memcpy(scrtch,p1->x);
	gsl_vector_axpby(-1, p2->x, 1, scrtch);
	for(int i = 0;i<scrtch->size;++i){
		scrtch->data[scrtch->stride*i] += 0.0000000000000000000000000000000000000000*(1.0*rand()/RAND_MAX - 0.5);
	}
	double d;
	d = gsl_blas_dnrm2(scrtch);
	double dqdt = dt*(p1->q)*(p2->q)/(d*d*d);
	gsl_vector_axpby(dqdt,scrtch, 1, p1->v);
	gsl_vector_axpby(-dqdt,scrtch, 1, p2->v);
	double ivn,jvn;
	gsl_blas_ddot(p1->v,p1->n, &ivn);
	gsl_blas_ddot(p2->v,p2->n, &jvn);
	gsl_vector_axpby(-jvn,p2->n,1,p2->v);
	gsl_vector_axpby(-ivn,p1->n,1,p1->v);

}

void upd_cursor_x_2d(Cursor *c,plane *Pi){
	gsl_vector_memcpy(Pi->scrtch, c->r);
	gsl_blas_ddot(Pi->scrtch,Pi->P_x,&c->x_2d->data[c->x_2d->stride*0]);
	gsl_blas_ddot(Pi->scrtch,Pi->P_y,&c->x_2d->data[c->x_2d->stride*1]);
}

void upd_cursor(Cursor *c,plane *Pi, double sc,double xw,double yh){
	double x = GetMouseX()*(c->w/xw) - 1.0*c->w/2;
	double y = GetMouseY()*(c->h/yh) - 1.0*c->h/2;
	x = x*(sc/c->w);
	y = y*(sc/c->h);
	c->x_2d->data[0] = x;
	c->x_2d->data[c->x_2d->stride] = y;
	gsl_vector_memcpy(c->r,Pi->r);
	gsl_vector_axpby(x, Pi->P_x, 1, c->r);
	gsl_vector_axpby(y, Pi->P_y, 1, c->r);
	upd_cursor_x_2d(c,Pi);

};

void cross_mat(const gsl_vector *u,gsl_matrix *ux){
	double u1 = gsl_vector_get(u,0);
	double u2 = gsl_vector_get(u,1);
	double u3 = gsl_vector_get(u,2);
	gsl_matrix_set_zero(ux);
	gsl_matrix_set(ux, 0, 1, -u3);
	gsl_matrix_set(ux, 1, 0, u3);
	gsl_matrix_set(ux, 0, 2, u2);
	gsl_matrix_set(ux, 2, 0, -u2);
	gsl_matrix_set(ux, 1, 2, -u1);
	gsl_matrix_set(ux, 2, 1, u1);

}

void merge_edges(edge *e1,edge *e2){
	double k1 = e1->k, l1 = e1->l;
	double k2 = e2->k, l2 = e2->l;
	double k = k1+k2;
	double l = (k1*l1 + k2*l2)/k;
	e1->k = k; e1->l = l;
}

edge *remove_edge(edge *e, int n_e,int i_r){
	for(int i = i_r;i<(n_e-1);++i){
		*(e +i) = *(e + i +1);
	}
	return realloc(e,(n_e-1)*sizeof(edge));
}

void reduce_edges(springs *s){
	int n_e = s->n_e, b = 0;
	point *p1i, *p1j, *p2i, *p2j;
	edge *e = s->e;
	int k_1, k_2;
	for(k_1 = 0; k_1<n_e;++k_1){
		p1i = (e+k_1)->p_i;
		p1j = (e+k_1)->p_j;
		for(k_2 = k_1+1;k_2<n_e;++k_2){
			p2i = (e+k_2)->p_i;
			p2j = (e+k_2)->p_j;
			b = (p1i == p2i)*(p1j == p2j) + (p1i == p2j)*(p1j == p2i);
			if(b){
				break;
			}

		}
		if(b){
			break;
		}
	}
	if(b){
	printf("merging edges %d and %d\n",k_1,k_2);
	merge_edges(e+k_1, e+k_2);
	s->e = remove_edge(e,n_e,k_2);
	s->n_e = n_e-1;
	reduce_edges(s);
	}
}

void rotate(plane *Pi, double eps){
	gsl_blas_dgemv(CblasNoTrans, eps, Pi->rot_mat,Pi->P_x, 1, Pi->P_x);
	gsl_blas_dgemv(CblasNoTrans, eps, Pi->rot_mat,Pi->P_y, 1, Pi->P_y);
	//double dx = gsl_blas_dnrm2(s->P_x);
	//double dy = gsl_blas_dnrm2(s->P_y);
	//gsl_vector_scale(s->P_x, 1.0/(dx + 0.0000000000000000000000000001));
	//gsl_vector_scale(s->P_y, 1.0/(dy + 0.0000000000000000000000000001));
}

void upd_point_x(point *p,springs *s){
	//gsl_vector_axpby(-1,p->x,1,p->n);
	gsl_vector_axpby(dt,p->v,1,p->x);
	//gsl_vector_axpby(1,p->x,1,p->n);

}

void upd_point_x_2d(point *p,plane *Pi){
	gsl_vector_memcpy(Pi->scrtch, Pi->r);
	gsl_vector_axpby(1, p->x, -1, Pi->scrtch);
	gsl_blas_ddot(Pi->scrtch,Pi->P_x,&p->x_2d->data[p->x_2d->stride*0]);
	gsl_blas_ddot(Pi->scrtch,Pi->P_y,&p->x_2d->data[p->x_2d->stride*1]);
}



void upd_edge(edge *e,gsl_vector *scrtch){
	point *p_i = e->p_i, *p_j = e->p_j;
	gsl_vector *x_ij = scrtch;
	gsl_vector_memcpy(x_ij,p_i->x);
	gsl_vector_sub(x_ij,p_j->x);
	for(int i = 0;i<x_ij->size;++i){
		x_ij->data[i*x_ij->stride] += 0.000000000000000000000000000000000000*(1.0*rand()/RAND_MAX - 0.5);
	}
	double d = gsl_blas_dnrm2(x_ij);
	gsl_vector_scale(x_ij, 1/(d));
	double diff_dt_k = e->k*dt*(e->l-d)/e->beta;
	gsl_vector_axpby(-diff_dt_k/p_j->m,x_ij,1/e->beta,p_j->v);
	gsl_vector_axpby(diff_dt_k/p_i->m,x_ij,1/e->beta,p_i->v);
	double ivn,jvn;
	gsl_blas_ddot(p_i->v,p_i->n, &ivn);
	gsl_blas_ddot(p_j->v,p_j->n, &jvn);
	gsl_vector_axpby(-jvn,p_j->n,1,p_j->v);
	gsl_vector_axpby(-ivn,p_i->n,1,p_i->v);
	e->beta *= 1.00001;
}
void upd_springs(springs *s){
	//double M = 0;
	//gsl_vector_set_zero(s->mx);

	for(int i = 0;i<s->n_p;++i){
		upd_point_x(s->p + i,s);
		for(int j=i+1;j<s->n_p;++j){
			electric(s->p +i, s->p +j, s->scrtch);
		}
	//	gsl_vector_axpby(s->p->m,s->p->x,1,s->mx);
	//	M = M + s->p->m;
	}
	for(int i = 0;i<s->n_e;++i){
		upd_edge(s->e + i,s->scrtch);
	}
	//gsl_vector_scale(s->mx,1/M);

}

void upd_edge_pix(edge *e,Color *pix,int w, int h){
	int ix = e->p_i->int_x, iy = e->p_i->int_y;
	int jx = e->p_j->int_x, jy = e->p_j->int_y;
	int dx = ix-jx, dy = iy-jy;
	double d = sqrt(dx*dx + dy*dy);
	if(((ix<0) &(jx<0)) | (((ix>= w) & (jx >= w))) |  ((iy<0) &(jy<0)) | (((iy>= h) & (jy >= h))) ){
		//printf("we skippin \n");
		return;
	}
	for(int i = 0;i<d;++i){
		for(int r = -1;r<1 ;++r){
			for(int s = -1;s<1;++s){
				int x = round((i/(d-1))*jx + (d-1 -i)/(d-1)*ix) + r;
				int y = round((i/(d-1))*jy + (d-1 -i)/(d-1)*iy) + s;
				//printf("(x,y) =  (%d,%d)    ;",x,y);
				if(x < 0 | x >= w | y < 0 | y >= h){
					continue;
				}
				pix[w*y +x] = (e->c);


			}
		}
	}
	//printf("\n");



}




void update_springs_pix(springs *s,plane *Pi,Color *pix ,int w, int h,double sc){
	for(int i = 0;i<s->n_p;++i){

		double int_r = sqrtf(w*w+h*h)*((s->p+i)->r/sc);
		upd_point_x_2d(s->p+i,Pi);
		gsl_vector *px = (s->p+i)->x_2d;
		double x = gsl_vector_get(px,0);
		double y = gsl_vector_get(px,1);
		x = x/sc;
		y = y/sc;
		int int_x = x*w + 1.0*w/2;
		int int_y = y*h + 1.0*h/2;
		(s->p+i)->int_x = int_x;
		(s->p+i)->int_y = int_y;
		int ro, co, idx;
		if(int_x*int_x+int_y*int_y - int_r*int_r > w*w +h*h){
			//printf("we skippin\n");
			continue;
		}
		for(int k = -int_r/2;k<int_r/2;++k){
 			co = int_x+k;
			if((co >=  w) | (co < 0)){
				continue;
			}
			for(int j = -int_r/2;j<int_r/2;++j){

			    ro = (int_y + j);
			    if(ro >= h | ro < 0){
				    continue;
			    };

			    idx = ro*w + co;
			    pix[idx] = BLUE;
			}
		}
		printf("\n");
	
	}
	for(int i = 0;i<s->n_e;++i){
		edge *e = s->e +i;
		upd_edge_pix(e,pix,  w,  h);
	}

}


int are_incident(edge *e1, edge *e2){
	point *p1i = e1->p_i, *p1j = e1->p_j;
	point *p2i = e1->p_i, *p2j = e1->p_j;
	return 0 < ((p1i == p2i) + (p1i == p2j) + (p1j == p2i) + (p1j == p2j));
}

void upd_edge_s(edge* e,int n_e,double p12,double p20,double p_extern,double p_transmit){
	double t = 1.0*rand()/RAND_MAX;
	if(e->s == 1){
		if(t< p12){
			e->next_s = 2;
			return;
		}
	}else if (e->s == 2) {
		if(t < p20){
			e->next_s = 0;
			return;
		}

	
	}else{
		if(t < p_extern){
			e->next_s = 1;
			return;
		}else {
			edge *n = e - e->id;
			for(int i = 0 ;i <n_e;++i){
				if(i == e->id){
					continue;
				}
				//printf("are incident = %d",are_incident(e,n+i));
				if(are_incident(e, n+i) & (((n+i)->s) == 1)){
					if (1.0*rand()/RAND_MAX < p_transmit){
						printf("we transmitin\n");
						e->next_s = 1;
						return;
					}

				}

			}
		
		}


	}
}

void upd_edges_s(edge* e,int n_e,double p12,double p20,double p_extern,double p_transmit){
	for(int i = 0;i<n_e;++i){
		upd_edge_s(e+i, n_e, p12, p20,p_extern, p_transmit);
	}
	for(int i = 0;i<n_e;++i){
		(e+i)->s = (e+i)->next_s;
		int s = (e+i)->s;
		if(s == 0){
			(e+i)->c = BLUE;
		}
		if(s == 1){
			(e+i)->c = RED;
		}
		if(s == 2){
			(e+i)->c = BLACK;
		}

	}
}

void upd_edges_s_rand(edge* e,int n_e,double p12,double p20,double p_extern,double p_transmit){
	int i = floor(n_e*1.0*rand()/RAND_MAX);
	
		upd_edge_s(e+i, n_e, p12, p20,p_extern, p_transmit);
	

		(e+i)->s = (e+i)->next_s;
		int s = (e+i)->s;
		if(s == 0){
			(e+i)->c = BLUE;
		}
		if(s == 1){
			(e+i)->c = RED;
		}
		if(s == 2){
			(e+i)->c = BLACK;
		}

	
}


int main(){
	int d = 3;

	int w = 1000, h = 1000;
	gsl_vector *mx = gsl_vector_calloc(d);
	gsl_vector *scrtch = gsl_vector_calloc(d);
	gsl_vector *pi_scrtch = gsl_vector_calloc(d);
	gsl_vector *r = gsl_vector_calloc(d);
	
	gsl_vector *P_x = gsl_vector_calloc(d);
	gsl_vector *P_y = gsl_vector_calloc(d);
	gsl_vector *P_z = gsl_vector_calloc(d);
	gsl_matrix *rot_mat = gsl_matrix_calloc(d,d);
	P_x->data[0] = 1;
	P_y->data[1*P_y->stride] = 1;
	P_z->data[0*P_z->stride] = 1;
	P_z->data[1*P_z->stride] = 1;


	gsl_vector *r_c = gsl_vector_calloc(d);
	gsl_vector *x_2d_c = gsl_vector_calloc(2);
	Cursor cursor = {.h = h,.w = w,.r=r_c,.x_2d=x_2d_c};

	plane Pi = {.P_x=P_x,.P_y = P_y,.rot_ax=P_z,.rot_mat=rot_mat,.scrtch = pi_scrtch,.r = r};
	cross_mat(Pi.rot_ax,Pi.rot_mat);




	gsl_matrix *A = csv_to_matrix("sA.csv",',');
	int n = A->size1, g = A->size2;
	printf("%d %d\n",n,g);
	for(int i = 0;i<n;++i){

		printf("\n");
		for(int j = 0;j<g;++j){
			printf("%g,",gsl_matrix_get(A,j,i));

		}
	}
	gsl_vector *m = gsl_vector_calloc(n);
	gsl_vector_set_all(m,1);

	//m->data[0] = 10;
	gsl_vector *ra = gsl_vector_calloc(n);
	gsl_vector_set_all(ra,0.01);

	//ra->data[0] = 0.1;
	springs *s = init_springs_from_adjacency(A, m, ra, d);



	SetConfigFlags(FLAG_WINDOW_RESIZABLE);
	InitWindow(w, h, "");
	Image img = GenImageColor(w, h, BLACK);
	Texture tex = LoadTextureFromImage(img);

	Color rgba_pixels[h * w]; // 4 channels


	int sw = GetScreenWidth(), sh = GetScreenHeight();	
	Rectangle rec = (Rectangle) {.x = 0,.y = 0,.width = w,.height=h};
	Rectangle big_rec =(Rectangle) {.x= 0,.y=0,.width = sw,.height = sh};
	int a = 0;
	int scale = 1;
	int b = 1;
	double x0 =0, y0 = 0,sc= 1;
	double p12 = 1, p20 = 0.5, p_transmit = 0.8;
	double ha = 0.0000000000000000000000000001;
	double hq = 1000000000000000;
	double p_extern = 1 - exp(-ha*dt);
	p12 = 1 - exp(-hq*dt);


	pix_homotopy f = {.x00=w/2+50,.y00=h/2,.x01=w/2+50,.y01=h/2,.x10=w/2 -100,.y10=h/2 + 100,.x11=w/2 + 100,.y11= h/2,.c=BLACK};

	while(!WindowShouldClose()){
		p_extern = 1 - exp(-ha*dt);
		p12 = 1 - exp(-hq*dt);
		sw = GetScreenWidth(); sh = GetScreenHeight();	
		upd_cursor(&cursor, &Pi,sc,sw,sh);
		big_rec.width = sw; big_rec.height = sh;
		if(b | IsKeyPressed(KEY_UP)){
		upd_springs(s);
		upd_edges_s_rand(s->e, s->n_e,p12,p20,p_extern,p_transmit);
				}
		for(int i = 0;i<h*w;++i){
			rgba_pixels[i] = WHITE;
		}

		update_springs_pix(s,&Pi, rgba_pixels, w, h,sc);
			

		if(IsKeyDown(KEY_LEFT)){
			rotate(&Pi,dt);
		}else if (IsKeyDown(KEY_RIGHT)) {

			rotate(&Pi,-dt);
		
		}
		if(IsKeyPressed(KEY_SPACE) & b){
			b = 0;
		}else if (IsKeyPressed(KEY_SPACE)){
			b = 1;
		
		}
		if(IsKeyDown(KEY_W)){
			gsl_vector_axpby(dt, Pi.P_y, 1, Pi.r);
		}if(IsKeyDown(KEY_S)){

			gsl_vector_axpby(-dt, Pi.P_y, 1, Pi.r);
		}	
		if(IsKeyDown(KEY_A)){

			gsl_vector_axpby(-dt, Pi.P_x, 1, Pi.r);
		}if(IsKeyDown(KEY_D)){

			gsl_vector_axpby(dt, Pi.P_x, 1, Pi.r);
		}
		if(IsKeyDown(KEY_T)){
			ha = ha*(1+dt);

		}
		if(IsKeyDown(KEY_G)){
			ha = ha*(1-dt);
		}		
		if(IsKeyDown(KEY_Z)){
			hq = hq*(1+dt);

		}
		if(IsKeyDown(KEY_H)){
			hq = hq*(1-dt);
		}


		upd_homotopy_pix(&f,rgba_pixels,w,h,lintrp);
		
	
		UpdateTexture(tex, rgba_pixels);
		//GenTextureMipmaps(&tex);

		//SetTextureFilter(tex, 	TEXTURE_FILTER_BILINEAR);
		//	DrawCircle(jx,jy, 10, BLUE);
		BeginDrawing();
		//reset window and draw texture
		ClearBackground(RAYWHITE);
		DrawTexturePro(tex,rec,big_rec,(Vector2){0,0},0,WHITE);
		DrawFPS(0,0);
		DrawText(TextFormat("p_e = %g",p_extern), 0, 20,20, BLACK);
		DrawText(TextFormat("p_12 = %g",p12), 0, 40,20, BLACK);
		for(int i = 0;i<s->n_e;++i){
			//printf("%d\n",i);
			edge *e = s->e + i;

			int ix = e->p_i->int_x, iy = e->p_i->int_y;
			int jx = e->p_j->int_x, jy = e->p_j->int_y;
		
		}
		//DrawCircle(w/2,h/2 ,10, BLACK);
		EndDrawing();
		float dsc = GetMouseWheelMove();
		sc = (1+0.01*dsc)*sc;
		}
	free_springs(s);
	free_plane(&Pi);
	free(s);
	UnloadTexture(tex);
	CloseWindow();
	return 0;}
