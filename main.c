#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>

#include <SDL3/SDL.h>
#include <SDL3/SDL_render.h>


#include "src/dot_isoline.c"
#include "src/pixel_solve_isoline.c"

float f(il_vec2 arg, float * param, size_t param_len){
    return (sinf(arg[0]*param[1]) + sqrtf(arg[1]*param[1]));
}


#define P_LEN 3

#define LINES_COUNT 64
#define LINE_BORDER_START -3.0f
#define LINE_BORDER_DIFF 0.2f

float DRAW_SCALING = 100.0f;
float DRAW_OFFSET = 10.0f;


int main(){
    float * param = (float*) malloc( P_LEN * sizeof(float) );
    param[0] = 1.0f;
    param[1] = 1.0f;
    param[2] = 0.1f;

    il_pixelsolve_config cfg;
    il_pixelsolve_data * data = (il_pixelsolve_data*)malloc(
            LINES_COUNT*sizeof(il_pixelsolve_data)
    );


    for(size_t i=0; i<LINES_COUNT; ++i){
        cfg.offt[0] = 0.0f; cfg.offt[1] = 0.0f;
        cfg.pixel_scale = 0.2f;
        cfg.pixel_len[0] = 50; cfg.pixel_len[1] = 50;
        cfg.f = &f;
        cfg.f_param = param;
        cfg.f_param_len = P_LEN;
        cfg.f_border_value = LINE_BORDER_START + LINE_BORDER_DIFF*i;
        
        cfg.equ_iter = 100;

        il_pixelsolve_isoline(&data[i],cfg);
    }

    



    // WINDOW 
    
    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL_Init failed: %s\n", SDL_GetError());
        return 1;
    }

    // Create a window
    SDL_Window *window = SDL_CreateWindow("Isoline demo",
        1000, 1000, SDL_WINDOW_RESIZABLE);

    if (!window) {
        printf("Window creation failed: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    
    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    if (!renderer) {
        printf("Renderer creation failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }



    int running = 1;
    SDL_Event e;
    while (running) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT)
                running = 0;
        }

        // Clear screen (black)
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        
        for(size_t di=0; di<LINES_COUNT; ++di){
        il_pixelsolve_data d = data[di];

        SDL_SetRenderDrawColor(renderer, 255, di*50, di*5, 255);


        for(size_t i=0; i<d.edges_len; ++i){
            il_connection ci;
            ci[0] = d.edges[i][0];
            ci[1] = d.edges[i][1];

            il_vec2 p0, p1;
            //coords 0
            p0[0] = d.vertex[ci[0]][0];
            p0[1] = d.vertex[ci[0]][1];
            
            //coords 1
            p1[0] = d.vertex[ci[1]][0];
            p1[1] = d.vertex[ci[1]][1];

            SDL_RenderLine(
                    renderer,
                    p0[0]*DRAW_SCALING + DRAW_OFFSET,
                    p0[1]*DRAW_SCALING + DRAW_OFFSET,
                    p1[0]*DRAW_SCALING + DRAW_OFFSET,
                    p1[1]*DRAW_SCALING + DRAW_OFFSET
            );
        }

        }
        //SDL_RenderLine(renderer, 0, 0, 300, 300); 
        // Present the result
        SDL_RenderPresent(renderer);
        

        SDL_Delay(100);
    }


     for(size_t di=0; di<LINES_COUNT; ++di){
        il_pixelsolve_data * d = &data[di];
        free(d->edges);
        free(d->vertex);
     }


    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}






