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
	gsl_vector *P_x;
	gsl_vector *P_y;
	gsl_vector *rot_ax;
	gsl_matrix *rot_mat;
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

void rotate(springs *s, double eps){
	gsl_blas_dgemv(CblasNoTrans, eps, s->rot_mat,s->P_x, 1, s->P_x);
	gsl_blas_dgemv(CblasNoTrans, eps, s->rot_mat,s->P_y, 1, s->P_y);
	double dx = gsl_blas_dnrm2(s->P_x);
	double dy = gsl_blas_dnrm2(s->P_y);
	//gsl_vector_scale(s->P_x, 1.0/(dx + 0.0000000000000000000000000001));
	//gsl_vector_scale(s->P_y, 1.0/(dy + 0.0000000000000000000000000001));
}

void upd_point_x(point *p,springs *s){
	gsl_vector_axpby(dt,p->v,1,p->x);
	gsl_blas_ddot(p->x,s->P_x,&p->x_2d->data[p->x_2d->stride*0]);
	gsl_blas_ddot(p->x,s->P_y,&p->x_2d->data[p->x_2d->stride*1]);
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
	double M = 0;
	gsl_vector_set_zero(s->mx);
	for(int i = 0;i<s->n_e;++i){
		upd_edge(s->e + i,s->scrtch);
	}
	for(int i = 0;i<s->n_p;++i){
		upd_point_x(s->p + i,s);
		gsl_vector_axpby(s->p->m,s->p->x,1,s->mx);
		M = M + s->p->m;
	}
	gsl_vector_scale(s->mx,1/M);

}



void update_springs_pix(springs *s,Color *pix ,int w, int h){
	double int_r = sqrtf(w*w+h*h)*s->p->r;
	for(int i = 0;i<s->n_p;++i){
		gsl_vector *px = (s->p+i)->x_2d;
		int int_x = gsl_vector_get(px,0)*w/10;
		int int_y = gsl_vector_get(px,1)*h/10;
		for(int k = -int_r/2;k<int_r/2;++k){
			for(int j = -int_r/2;j<int_r/2;++j){
			    int idx = (h/2  - int_y + j)*w + w/2+int_x+k;
			    if(idx > w*h){
				    continue;
			    };
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
	
	gsl_vector *P_x = gsl_vector_calloc(n);
	gsl_vector *P_y = gsl_vector_calloc(n);
	gsl_vector *P_z = gsl_vector_calloc(n);
	gsl_matrix *rot_mat = gsl_matrix_calloc(n,n);
	P_x->data[0] = 1;
	P_y->data[1*P_y->stride] = 1;
	P_z->data[0*P_z->stride] = 1;
	P_z->data[1*P_z->stride] = 1;



	int np = 3;
	int ne = 3;
	point *p = malloc(np*sizeof(point));
	edge *e = malloc(np*sizeof(edge));

	for(int i = 0;i < np; ++i){
		init_point(p+i, n, 1, 0.01);
	}

	p->x->data[0] = 1;
	(p+1)->x->data[1] = -1;
	(p+2)->x->data[2] = 1;

	edge e12 = {.k=1,.l=1,.p_i = p,.p_j = p+1};
	edge e13 = {.k=1,.l=1,.p_i = p,.p_j =p+2};
	edge e23 = {.k=1,.l=1,.p_i = p+1,.p_j = p+2};
	int k = 0;

	*e = e12;*(e+1) = e13;*(e+2) = e23;

	

	springs s = {.n_p=np,.n_e=ne,.e=e,.p=p,.mx=mx,.scrtch = scrtch,.P_x=P_x,.P_y = P_y,.rot_ax=P_z,.rot_mat=rot_mat};
	cross_mat(s.rot_ax,s.rot_mat);





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
	Cursor cursor = {.x = 1.0*w*GetMouseX()/sw,.y = 1.0*GetMouseY()/sh,.vx = 0,.vy = 0,.s = 0,.h = h,.w = w};
	while(!WindowShouldClose()){
		
		sw = GetScreenWidth(); sh = GetScreenHeight();	
		big_rec.width = sw; big_rec.height = sh;
		upd_springs(&s);
		for(int i = 0;i<h*w;++i){
			rgba_pixels[i] = WHITE;
		}
		update_springs_pix(&s, rgba_pixels, w, h);
		printf("\n-------------------------\n");
		edge *e = s.e;
		
		printf("\nx = ");
		for(int i = 0;i<n;++i){
			printf("%g",s.P_x->data[i*s.P_x->stride]);
		}
		printf("\ny = ");
		for(int i = 0;i<n;++i){
			printf("%g",s.P_y->data[i*s.P_x->stride]);
		}

		if(IsKeyDown(KEY_LEFT)){
			rotate(&s,dt);
		}else if (IsKeyDown(KEY_RIGHT)) {

			rotate(&s,-dt);
		
		}
	
		UpdateTexture(tex, rgba_pixels);
		BeginDrawing();
		//reset window and draw texture
		ClearBackground(RAYWHITE);
		DrawTexturePro(tex,rec,big_rec,(Vector2){0,0},0,WHITE);
		EndDrawing();
		int mx = w*GetMouseX()/sw,my = h*GetMouseY()/sh;
		float s = GetMouseWheelMove();
		cursor.vx = mx - cursor.x;
		cursor.vy = my - cursor.y;
		cursor.s = s;
		cursor.x = mx;
		cursor.y = my;
		}
	UnloadTexture(tex);
	CloseWindow();
	return 0;}
