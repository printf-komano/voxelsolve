#include <stdio.h>

#include "src/vs_opengl.h"
#include "src/vs_voxelsolve.h"


float f(vs_vec3 v, float * args, size_t argc){
    return v[1];
}

float ARGS[0] = {};
size_t ARGC = 0;

const vs_vec3 OFFT = {0.0f, 0.0f, 0.0f};
const float VSCALE = 1.0f;
const vs_vec3i VOX_LEN = {1,1,1};
const size_t STEPS = 8;

const float THRESHOLD = 0.5f;
const float PRECISION = 0.1f;

vs_voxelsolve_data DATA;

int main(){
    vs_vec3 a = {0.0f,0.5001f,0.0f};
    vs_vec3 b = {0.0f,0.5f,0.0f};
    printf("a=b? %b\n\n\n", vs_vec3_cmp(a,b,PRECISION));

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
 
    voxelsolve(&DATA,con);
    printf("VERTEX [%d]:\n", DATA.vertex_len);
    for(int i=0; i<DATA.vertex_len; ++i){
        printf("(%f %f, %f)\n",
                DATA.vertex[i][0],
                DATA.vertex[i][1],
                DATA.vertex[i][2]
        );
    }

        
    return 0;
}
