#include "raylib.h"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rl_panes.h"

int main() {
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
	
	// push the rgba pixel values to the texture
	
	int a = 0;
	int scale = 1;	
	Rectangle rec = {.x = 0,.y = 0,.width = w,.height = h };
	Pane w_pane = (Pane) {.r = (Rectangle) {.x=0,.y=0,.width=1,.height=1},.is_drawable = 0,.id = 0};
	int n_panes = 2;
	Pane *p_ptr = malloc(n_panes*sizeof(Pane));
	*p_ptr = w_pane;
	w_pane.is_grabable = 1;
	w_pane.is_movable = 1;
	*(p_ptr +1) = *init_sub_pane(p_ptr, 0, 0,0.5, 0.5, c_fun_from_pix,1,&RED,1,1);
	//*(p_ptr +2) = *init_sub_pane(p_ptr+1, 0, 0,0.5, 0.5, c_fun_grad_rb,2);
	//*(p_ptr +3) = *init_sub_pane(p_ptr+1, 0.1, 0.1,0.5, 0.5, c_fun_rand,3);
	int held_idx = -1;
	int *panes_held = malloc(n_panes);
	int *panes_selected = malloc(n_panes);
	int *cursor_on_panes = malloc(n_panes);
	Color *sub_pix = &BLUE;
	Cursor cursor = {.x = 1.0*w*GetMouseX()/GetScreenWidth(),.y = 1.0*GetMouseY()/GetScreenHeight(),.vx = 0,.vy = 0,.s = 0,.h = h,.w = w};
	while(!WindowShouldClose()){
		for(int i = 0;i<h*w;++i){
			rgba_pixels[i] = WHITE;
		}
		Rectangle big_rec = {.x = 0,.y = 0,.width = GetScreenWidth(),.height = GetScreenHeight() };
		IsCursorOnPanes(p_ptr,&cursor,cursor_on_panes,n_panes,w,h);
		//printf("%d \n",cursor_on_panes[0]);
		//printf("cursor pos (x,y) = (%d,%d),held_idx = %d \n",cursor.x,cursor.y,held_idx);
		int b = IsMouseButtonDown(0);
		
		if(held_idx == -1 && b){
			for(int i = 0;i<n_panes;++i){
				if(cursor_on_panes[i] && (p_ptr+i)->is_grabable){
					held_idx = i;
				}
		};};
		if(!b){held_idx = -1;};
		//update pixels before drawin 
		// update_pixels(rgba_pixels,w,h,...)
		//p_ptr[0].r.x = 1.0*w*rand()/RAND_MAX;
		//p_ptr[0].r.y = 1.0*h*rand()/RAND_MAX;
		update_panes_pixels(p_ptr,n_panes,rgba_pixels, w, h);
		//load new pixs into texture
		UpdateTexture(tex, rgba_pixels);
		BeginDrawing();
		//reset window and draw texture
		ClearBackground(RAYWHITE);
		DrawTexturePro(tex,rec,big_rec,(Vector2){0,0},0,WHITE);
		EndDrawing();
		int mx = w*GetMouseX()/GetScreenWidth(),my = h*GetMouseY()/GetScreenHeight();
		float s = GetMouseWheelMove();
		cursor.vx = mx - cursor.x;
		cursor.vy = my - cursor.y;
		cursor.s = s;
		cursor.x = mx;
		cursor.y = my;
		if(held_idx != -1){
			if(IsMouseButtonPressed(1)){
				n_panes += 1;
				p_ptr = realloc(p_ptr, n_panes*sizeof(Pane));
				*(p_ptr + (n_panes - 1)) = *init_sub_pane(p_ptr+held_idx, 0, 0, 0.5, 0.5,c_fun_from_pix,n_panes-1,sub_pix,1,1);

			}

			Rectangle ol_r = (p_ptr+held_idx)->r;		
			move_pane(p_ptr+held_idx, &cursor);
	//		move_childs(p_pPane *ptr, &cursor,n_panes,held_idx);;
					
			upd_childs(p_ptr, n_panes, held_idx,ol_r.width,ol_r.height,ol_r.x,ol_r.y);}

	}
	UnloadTexture(tex);
	CloseWindow();
	return 0;
}

