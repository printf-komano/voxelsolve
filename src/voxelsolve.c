#pragma once


#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>



typedef float vs_vec3[3]; // regular vector that will be used here 
                          // (x,y,z)
                          // Y up

typedef int32_t vs_vec3i[3];
typedef size_t vs_tri[3]; // triangle described by 3 index




typedef struct {
    vs_vec3 offt; 

    float vscale;
    vs_vec3i vox_len;

    size_t solve_steps; 

    // any mathematical function 
    // that returns float value;
    float (*f)(
            vs_vec3,
            float *, // f() additional arguments
            size_t // f() arguments count
            );       
    float * farg;
    size_t fargc;
    float f_edge; // value where line is going

} vs_voxelsolve_con;



typedef struct {
    vs_vec3 * vertex;
    size_t vertex_len;

    vs_tri * triangles;
    size_t triangles_len;
} vs_voxelsolve_data;





#define VS_VEC3_SET(to,from)( to[0] = from[0]; to[1] = from[1]; to[2] = from[2] )




static inline void lerp3(
        vs_vec3 start,
        vs_vec3 end,
        vs_vec3 out,
        float v
        )
{
    out[0] = start[0]  +  (end[0] - start[0]) * v;
    out[1] = start[1]  +  (end[1] - start[1]) * v;
    out[2] = start[2]  +  (end[2] - start[2]) * v;
}





static inline bool is_edge_voxel(
        vs_vec3 start,
        vs_voxelsolve_con con
        )
{
    const size_t CUBE_VLEN = 8;

    vs_vec3 verts[CUBE_VLEN];
    float vals[CUBE_VLEN];
    
    for(size_t i=0; i<CUBE_VLEN; ++i){
        VS_VEC3_SET(verts[i],start);
    }
    
    // unfold cube

    // verts[0] untouched
    verts[1][0] += con.vscale;
    verts[2][1] += con.vscale; 
    verts[3][0] += con.vscale; verts[3][1] += con.vscale;

    
    verts[4][2] += con.vscale;
    verts[5][2] += con.vscale;
    verts[6][2] += con.vscale;
    verts[7][2] += con.vscale;

    // verts[4] untouched
    verts[5][0] += con.vscale;
    verts[6][1] += con.vscale; 
    verts[7][0] += con.vscale; verts[3][1] += con.vscale;

    float first_val = con.f(verts[0], con.farg, con.fargc);

    for(size_t i=1; i<CUBE_VLEN; ++i){
        float vi = con.f(verts[i], con.farg, con.fargc);
        // any value change above or belov
        // f_edge will result in calling 
        // voxel solve method for this point
        if(
            (val > con.f_edge && first_val < con.f_edge) || 
            (val < con.f_edge && first_val > con.f_edge) 
        ) return true;
    }

    return false;
}






static size_t edge_solve(
        vs_vec3 start,
        vs_vec3 dir,
        vs_vec3 out,
        vs_voxelsolve_con con
        )
{
    size_t solutions = 0;
    float sum_offt = 0.0f;  // summ of all solutions.
                            // each solution described
                            // by value in range [0;1]

    float step_size = con.vscale / (float)con.solve_steps;

    vs_vec3 point; // iterable point
    VS_VEC3_SET(point,start);
    
    float val = con.f(point, con.farg, con.fargc);
    float last_val = val;

    for (size_t i=1; i<=con.solve_steps; ++i){
        point[0] += dir[0] * step_size;
        point[1] += dir[1] * step_size;
        point[2] += dir[2] * step_size;

        float val = con.f(point, con.farg, con.fargc);
        
        if(
            (val > con.f_edge && last_val < con.f_edge) || 
            (val < con.f_edge && last_val > con.f_edge) 
        ){
            float last_val = val; // evaluate again 
                                  // (after another solution found)
            sum_offt = (float)i / (float)con.solve_steps;
            ++solutions;
        }
    }

    if(solutions==0) return 0;

    float real_offt = (sum_offt / (float)solutions) * con.vscale;
    out[0] = start[0] + dir[0] * real_offt;
    out[1] = start[1] + dir[1] * real_offt;
    out[2] = start[2] + dir[2] * real_offt;
    
    return solutions;
}




static size_t voxel_solve(
        vs_vec3 start,
        vs_voxelsolve_data * data,
        vs_voxelsolve_con con
        )
{

}




void voxelsolve_isosurface(
        vs_voxelsolve_data * data,
        vs_voxelsolve_con con
        )
{
    vs_vec3 doti;

    for(size_t xi=0; xi<con.vox_len[0];++xi){
        for(size_t yi=0; yi<con.vox_len[1];++yi){
            for(size_t zi=0; zi<con.vox_len[2];++zi){

                //iterable dot
                doti[0] = con.offt[0] + xi*con.vscale;
                doti[1] = con.offt[1] + xi*con.vscale;
                doti[2] = con.offt[2] + xi*con.vscale;


                if(is_edge_voxel(doti,con)){
                    voxel_solve(doti, data, con);
                }

            }
        }
    }
}



