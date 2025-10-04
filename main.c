#include <stdio.h>

#include "src/vs_opengl.h"
#include "src/vs_voxelsolve.h"


float f(vs_vec3 v, float * args, size_t argc){
    return v[0] * v[2] + v[1];
} 
float ARGS[0] = {};
size_t ARGC = 0;

const vs_vec3 OFFT = {0.0f, 0.0f, 0.0f};
const float VSCALE = 1.0f;
const vs_vec3i VOX_LEN = {100,100,100};
const size_t STEPS = 8;

const float THRESHOLD = 1.0f;
const float PRECISION = 0.1f;

int main(){
    return 0;
    vs_voxelsolve_con con;

    con.vscale      =    VSCALE;
    con.solve_steps =    STEPS;
    con.f           =    f;
    con.farg        =    ARGS;
    con.fargc       =    ARGC;
    con.f_threshold =    THRESHOLD;
    con.prec        =    PRECISION;
   
    con.offt[0] = OFFT[0];
    con.offt[1] = OFFT[1];
    con.offt[2] = OFFT[2];
    
    con.vox_len[0] = VOX_LEN[0];
    con.vox_len[1] = VOX_LEN[1];
    con.vox_len[2] = VOX_LEN[2];

    vs_voxelsolve_data data;
    
    voxelsolve(&data,con);
    
    return 0;
}
