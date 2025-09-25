#pragma once 
#include "vs_voxelsolve.h"


static int8_t tri_facedir_sign(const vs_vec3 * tri, const vs_vec3 ref_axis){
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

#define VS_TRI_FACEDIR_UNDEF   0
#define VS_TRI_FACEDIR_FORW    1
#define VS_TRI_FACEDIR_BACK   -1

    float sign = dot3(n,ref_axis);
    if(cmpf(sign,0.0f, 0.001f)) return TRI_DIR_UNDEF;
    if(sign<0.0f) return TRI_DIR_BACK;
    return TRI_DIR_FORW;
}




#define VS_TRI_ORDER_NONE 0
#define VS_TRI_ORDER_CW   1
#define VS_TRI_ORDER_CCW  2

void vs_makebuffers(
        vs_voxelsolve_data * data,  // input
        float * vbo,    size_t * vbo_len,
        size_t * ebo,   size_t * ebo_len,
        
        size_t vertex_floats,       // how much floats used
                                    // to describe 1 vertex
        
        size_t vertex_cofft,        // vertex coordinate offset
        
        uint8_t tri_order           // VS_TRI_ORDER_...
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

    /* 
        for each triangle try to form ebo;
        in vs_voxelsolve_data, all vertex are chaotic;
        here we putting them in the right order to describe
        openGL polygon;
    */
    
    /* just a duplicated references */
    float * tris = data->tris;
    size_t tris-len = data->tris_len;
    float vertex = data->vertex;
    size_t vertex_len = data->vertex_len;

    for(size_t i=0; i<tris_len; ++i){
        /* check triange's face direction */
        int8_t tfd_axis = -1; 
        bool tfd_forward = true;
        

        // if order is not specified, algorithm is the simpliest
        if(tri_order == VS_TRI_ORDER_NONE){
            for(size_t vi=0; vi<3; ++vi){
                ebo[ (i*3)+vi ] = data->tris[i][vi];
            }
            return;
        }




        /* try to find the right axis */
        for(size_t ax=0; ax<3; ++ax){
            int8_t tfd_sign = tri_facedir_sign(data->tris[i]);
            // found the right sign
            if(tfd_sign != TRI_DIR_UNDEF){
                tfd_axis = ax; 
                tfd_forward = (tfd_axis >= TRI_FACEDIR_FORW);
                break;
            }
        }
        if(tfd_axis == -1) continue; // invalid triangle
        
        /* 
            get two axes that describe a plain;
            for example, if tfd_axis is Y axis, then
            plain is described by X,Z axes;

            It's ised to order vertex in CW CCW order.
        */
        int8_t plain_axes[2];
        if(tfd_axis != 0) plain_axes[0] = 0;
        else plain_axes[0] = 1;
        if(tfd_axis != 1 && plain_axes[0] != 1) plain_axes[1] = 1;
        else plain_axes[1] = 2;
        
        /* now, fill the element buffer in the certain order */
        size_t ordered_index[3] = {0,0,0};
        float cord_max = -1.0f;

        // by plain axis 0
        for(size_t j=0; j<3; ++j){
        
        }
        

    }

}




