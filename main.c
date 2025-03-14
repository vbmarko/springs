#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "raylib.h"

#define dt 0.1

typedef struct point point ;
typedef struct edge edge;
typedef struct Cursor Cursor;
typedef struct springs springs;
typedef struct plane plane;

struct Cursor{
	double x,y,vx,vy,s;
	int h,w;
};

struct point{
	gsl_vector *x;
	gsl_vector *v;
	gsl_vector *x_2d;
	double m,r;
};
struct edge{
	point *p_i, *p_j;
	double k,l;
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
	gsl_vector *rot_ax;
	gsl_matrix *rot_mat;
	gsl_vector *scrtch;
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
	gsl_vector_axpby(dt,p->v,1,p->x);

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
	double d = gsl_blas_dnrm2(x_ij);
	gsl_vector_scale(x_ij, 1/(d+0.00000000000001));
	double diff_dt_k = e->k*dt*(e->l-d);
	gsl_vector_axpby(-diff_dt_k/p_j->m,x_ij,1,p_j->v);
	gsl_vector_axpby(diff_dt_k/p_i->m,x_ij,1,p_i->v);
}
void upd_springs(springs *s){
	//double M = 0;
	//gsl_vector_set_zero(s->mx);
	for(int i = 0;i<s->n_e;++i){
		upd_edge(s->e + i,s->scrtch);
	}
	for(int i = 0;i<s->n_p;++i){
		upd_point_x(s->p + i,s);
	//	gsl_vector_axpby(s->p->m,s->p->x,1,s->mx);
	//	M = M + s->p->m;
	}
	//gsl_vector_scale(s->mx,1/M);

}



void update_springs_pix(springs *s,plane *Pi,Color *pix ,int w, int h,double xw,double yh){
	double int_r = sqrtf(w*w+h*h)*(s->p->r/xw);
	for(int i = 0;i<s->n_p;++i){
		upd_point_x_2d(s->p+i,Pi);
		gsl_vector *px = (s->p+i)->x_2d;
		double x = gsl_vector_get(px,0);
		double y = gsl_vector_get(px,1);
		x = x/xw;
		y = y/xw;
		int int_x = x*w;
		int int_y = y*h;
		int ro, co, idx;
		for(int k = -int_r/2;k<int_r/2;++k){
 			co = w/2+int_x+k;
			if((co >=  w) | (co < 0)){
				continue;
			}
			for(int j = -int_r/2;j<int_r/2;++j){

			    ro = (h/2  - int_y + j);
			    if(ro >= h | ro < 0){
				    continue;
			    };

			    idx = ro*w + co;
			    pix[idx] = BLUE;
			}
		}
	}

}

void init_point(point *p,int n,double m,double r){
	gsl_vector *x = gsl_vector_calloc(n);
	gsl_vector *v = gsl_vector_calloc(n);
	gsl_vector *x_2d = gsl_vector_calloc(2);
	p->x = x; p->x_2d = x_2d; p->v = v;
	p->m = m; p->r = r;
}

void free_point(point *p){
	gsl_vector_free(p->x);
	gsl_vector_free(p->x_2d);
	gsl_vector_free(p->v);
}





int main(){
	int n = 3;
	gsl_vector *mx = gsl_vector_calloc(n);
	gsl_vector *scrtch = gsl_vector_calloc(n);
	gsl_vector *pi_scrtch = gsl_vector_calloc(n);
	gsl_vector *r = gsl_vector_calloc(n);
	
	gsl_vector *P_x = gsl_vector_calloc(n);
	gsl_vector *P_y = gsl_vector_calloc(n);
	gsl_vector *P_z = gsl_vector_calloc(n);
	gsl_matrix *rot_mat = gsl_matrix_calloc(n,n);
	P_x->data[0] = 1;
	P_y->data[1*P_y->stride] = 1;
	P_z->data[0*P_z->stride] = 1;
	P_z->data[1*P_z->stride] = 1;



	int np = 3;
	int ne = 5;
	point *p = malloc(np*sizeof(point));
	edge *e = malloc(ne*sizeof(edge));

	for(int i = 0;i < np; ++i){
		init_point(p+i, n, 1, 0.01);
	}

	p->x->data[0] = 1;
	(p+1)->x->data[1] = 0;
	(p+2)->x->data[0] = -1;

	edge e12 = {.k=0.33,.l=1,.p_i = p,.p_j = p+1};
	edge e12p = {.k=0.33,.l=1,.p_i = p,.p_j = p+1};
	edge e12pp = {.k=0.33,.l=1,.p_i = p,.p_j = p+1};
	edge e13 = {.k=1,.l=1,.p_i = p,.p_j =p+2};
	edge e23 = {.k=1,.l=1,.p_i = p+1,.p_j = p+2};
	int k = 0;

	*e = e12;*(e+1) = e13;*(e+2) = e23;; *(e+3) = (e12p), *(e+4) = e12pp;

	

	springs s = {.n_p=np,.n_e=ne,.e=e,.p=p,.mx=mx,.scrtch = scrtch};
	
	plane Pi = {.P_x=P_x,.P_y = P_y,.rot_ax=P_z,.rot_mat=rot_mat,.scrtch = pi_scrtch,.r = r};
	cross_mat(Pi.rot_ax,Pi.rot_mat);





	int w = 1000, h = 1000;
	SetConfigFlags(FLAG_WINDOW_RESIZABLE);
	InitWindow(w, h, "");
	SetTargetFPS(60);
	Image img = GenImageColor(w, h, BLACK);
	// set the image's format so it is guaranteed to be aligned with our pixel buffer format below
	//ImageFormat(&img, PIwELFORMAh_UNCOMPRESSED_R8G8B8A8);
	Texture tex = LoadTextureFromImage(img);

	Color rgba_pixels[h * w]; // 4 channels


	// make sure to set ALL ALPHA CHANNELS of the pixels (every fourth index) to 255, otherwise you'll have transparent pixels
	int sw = GetScreenWidth(), sh = GetScreenHeight();	
	// push the rgba pixel values to the texture
	Rectangle rec = (Rectangle) {.x = 0,.y = 0,.width = w,.height=h};
	Rectangle big_rec =(Rectangle) {.x= 0,.y=0,.width = sw,.height = sh};
	int a = 0;
	int scale = 1;
	int b = 1;
	double x0 =0, y0 = 0,xw = 1,yh = 1;
	Cursor cursor = {.x = 1.0*w*GetMouseX()/sw,.y = 1.0*GetMouseY()/sh,.vx = 0,.vy = 0,.s = 0,.h = h,.w = w};
	while(!WindowShouldClose()){
		
		sw = GetScreenWidth(); sh = GetScreenHeight();	
		big_rec.width = sw; big_rec.height = sh;
		if(b){
		upd_springs(&s);}
		for(int i = 0;i<h*w;++i){
			rgba_pixels[i] = WHITE;
		}
		update_springs_pix(&s,&Pi, rgba_pixels, w, h,xw,yh);
		//printf("\n-------------------------\n");
		edge *e = s.e;
	

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
	
	
		UpdateTexture(tex, rgba_pixels);
		BeginDrawing();
		//reset window and draw texture
		ClearBackground(RAYWHITE);
		DrawTexturePro(tex,rec,big_rec,(Vector2){0,0},0,WHITE);
		EndDrawing();
		float sc = GetMouseWheelMove();
		xw = (1+0.01*sc)*xw;
		yh = (1+0.01*sc)*yh;
		}
	UnloadTexture(tex);
	CloseWindow();
	return 0;}
