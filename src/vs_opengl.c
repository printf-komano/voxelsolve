#pragma once 
#include "vs_voxelsolve.h"


static int8_t tri_dir(float * tri, float * ref_axis){
    vs_vec3 a;
    vs_vec3 b;

    /* find 2 directions (01 and 02) */
    a[0] = tri[1][0] - tri[0][0];
    a[1] = tri[1][1] - tri[0][1];
    a[2] = tri[1][2] - tri[0][2];

    b[0] = tri[2][0] - tri[0][0];
    b[1] = tri[2][1] - tri[0][1];
    b[2] = tri[2][2] - tri[0][2];

    vs_vec3 n;
    norm(a,b, n);

#define TRI_DIR_UNDEF   0
#define TRI_DIR_FORW    1
#define TRI_DIR_BACK   -1

    float sign = dot3(n,ref_axis);
    if(cmpf(sign,0.0f, 0.001f)) return TRI_DIR_UNDEF;
    if(sign<0.0f) return TRI_DIR_BACK;
    return TRI_DIR_FORW;
}


void vs_glmakebuffers(
        vs_voxelsolve_data * data,
        float * vbo,    size_t * vbo_len,
        size_t * ebo,   size_t * ebo_len,
        
        size_t vertex_floats,   // how much floats used
                  // to describe 1 vertex
        
        size_t vertex_cofft  // vertex coordinate offset
){
    /* allocate memory for output */
    size_t vlen, elen; // length (measured in 1-elements)

    vlen = data->vertex_len * vertex_floats;
    *vbo_len = vlen
    *vbo = (float*)malloc((*vbo_len) * sizeof(float));
    
    elen = data->tris_len * 3; // 3 vertices for each triangle
    *ebo_len = elen;
    *ebo = (size_t*)malloc((*ebo_len) * sizeof(size_t));


    /* fill vbo */ 
    for(size_t i=0; i<data->vertex_len; ++i){
        vbo[ (i*vertex_floats)+vertex_cofft+0 ] = data->vertex[i][0];
        vbo[ (i*vertex_floats)+vertex_cofft+1 ] = data->vertex[i][1];
        vbo[ (i*vertex_floats)+vertex_cofft+2 ] = data->vertex[i][2];
    }


    for(size_t i=0; i<data->tris_len; ++i){
        
    }

}




