#pragma once


#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>





/* ----------------------------------------
    Main types;
    float 3d vectors should work well 
    with GLM.
---------------------------------------- */

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
    float f_threshold; // value where line is going

} vs_voxelsolve_con;



typedef struct {
    vs_vec3 * vertex;
    size_t vertex_len;

    vs_tri * triangles;
    size_t triangles_len;
} vs_voxelsolve_data;




typedef struct {
    vs_vec3 v_global[8];
    float vscale;

    vs_vec3 v_local[8]; // coordinates from 0 to 1
} vs_cube;



typedef struct vs_solution(){
    vs_vec3 start; 
    vs_vec3 end;

    size_t intersections;
    vs_vec3 value;
}



static const CUBE_VLEN = 8;
static const vs_vec3 VS_VEC3_ZERO = {0.0f, 0.0f, 0.0f};
static const vs_vec3 CUBE1[8] = {
    // bottom face
    {0.0f, 0.0f, 0.0f},
    {1.0f, 0.0f, 0.0f},
    {0.0f, 0.0f, 1.0f},
    {1.0f, 0.0f, 1.0f},

    // top face 
    {0.0f, 1.0f, 0.0f},
    {1.0f, 1.0f, 0.0f},
    {0.0f, 1.0f, 1.0f},
    {1.0f, 1.0f, 1.0f},
};

#define CUBE_START 0.0f
#define CUBE_END 1.0f

/* ----------------------------------------
    Math operations and constants
---------------------------------------- */


inline bool cmpf(float a, float b, float prec){
     return (fabsf(a - b) <= prec);
}


#define VS_VEC3_SET(to,from)( to[0] = from[0]; to[1] = from[1]; to[2] = from[2] )


// compares two vec3 with certain precision
inline bool cmp3(vs_vec3 a, vs_vec3 b, float prec){
    return(
        fabsf(a[0] - b[0]) <= prec &&
        fabsf(a[1] - b[1]) <= prec && 
        fabsf(a[2] - b[2]) <= prec 
    );
}

// create 8 vertex that describe a cube
// out array will have the length of 8
static inline void gen_cube(vs_vec3 start, float vscale, vs_vec3 * out){
    
    // basically, it's just offseted and scaled cube-1
    for(size_t i=0; i< CUBE_VLEN; ++i){
        out[i][0] = start[0] + CUBE1[0][0] * vscale;
        out[i][1] = start[1] + CUBE1[0][1] * vscale;
        out[i][2] = start[2] + CUBE1[0][2] * vscale;

    }
}



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

static inline float dist3(vs_vec3 a, vs_vec3 b){
    return sqrt(
            powf(a[0]-b[0],2.0f) +
            powf(a[1]-b[1],2.0f) +
            powf(a[2]-b[2],2.0f)
    );
}


static inline float scalar_diff3(vs_vec3 a, vs_vec3 b){
    return(
            a[0]-b[0]
            a[1]-b[1]
            a[2]-b[2]
    );
}










/* ----------------------------------------
    Triangulation
    algo to turn a set of dots into triangular mesh
---------------------------------------- */



static inline size_t nearest3(vs_vec3 start, vs_vec3 * v, size_t vlen){
    float min_dist = dist3(start,v[0]);
    size_t ret = 0;

    for(size_t i=1; i<vlen; ++i){
        float new_dist = dist3(start,v[i]);
        if(min_dist < new_dist){
            ret = i;
            min_dist = new_dist;
        }
    }
    return ret;
}



/*
    Find the index of first varied coordinate.
    Because all points located on the edges, 
    there's the only one index possible.

    other 2 values should be equal to 0 or 1.
*/
static inline int32_t first_varied3(vs_vec3 v, float prec){
    const size_t ret_max = 2;
    const size_t ret_offt = ret_max+1;

    int32_t ret = (
        (ret_offt+0) * ( !cmpf(v[0],CUBE_START,prec) 
                && !cmpf(v[0],CUBE_END,prec) ) + // x

        (ret_offt+1) * ( !cmpf(v[1],CUBE_START,prec) 
                && !cmpf(v[1],CUBE_END,prec) ) + // y
        
        (ret_offt+2) * ( !cmpf(v[2],CUBE_START,prec) 
                && !cmpf(v[2],CUBE_END,prec) )   // z
    );

    if( ret-ret_offt > 2 ) return -1;
    else return ret - ret_offt;
}



static inline bool is_neighbour_edge(vs_vec3 a, vs_vec3 b, float prec){
    size_t vara = first_varied3(a,prec);
    size_t varb = first_varied3(b,prec);
    
    if(
        vara < 0 || varb < 0 // dot is not located correctly
    ) return false;

    return(
        vara == varb    // it's definitely on the same edge (not neighbour)
                        // or on the opposite one
    );
}





/*
    Update: new idea for dots connection. 


    If cube has 3 or 4 solutions, triangles are being 
    built in the most classical way.


    If not, follow the algorythm:

    Cube has vertices ("CV") placed on (0-1) coordinates. 
        ----------------------------------------------------------
        3-triangulation 

        1. We pick a CVi and check it's value. If value is greater
           than f_threshold, SKIP CVi (ignore CV inside isosurface).

        2. If not, check all solutions on the neighbour edges (form array).

        3. There's always gonnabe 0-3 dots; 
           If 3 dots connected to CVi, build triangle
           WITHOUT ADDING CONNECTIONS and skip the step.
           Else, if number of dots less then 3, add to queue.

        4. Continue this process in a loop until all CV precessed.
        
        ----------------------------------------------------------
        Additional 4-triangulation

        5. At this point, there're left some additional dots 
        also requiring connection.

        6. Cheoose a dot[i] and the nearest neighbour dots[j] ONLY
           from the neighbour faces (not edges this time). Build a triangle.
           EXCLUDE all triangulated dots[j] from i-iteration 
           (they still can be connected).
           EXCLUDE current dot[i] from j-iteraon (cannot be connected at all)

        7. End of algorythm.
*/


static void dots_triang(
        bool * cv_used,         // f(cv) < threshold (must be bool[8])
        vs_vec3 * dots,         // solutions
        size_t dots_len,        
        vs_tri * out_buffer,    // result
        size_t out_len,
        )
{   
    // unable to build figure
    if(dots_len < 3) return;

    // sipliest variant of the algorythm *(add only 1 trangle)
    else if(dots_len == 3) {
        out_len = 1;
        out_buffer = (vs_tri*) malloc(out_len * sizeof(vs_tri));
    }
    
    
    /*  1  */
    vs_vec3 cv[CUBE_VLEN];
    float cv_values[CUBE_VLEN];
    gen_cube(start,con.vscale,cv);    

    for(size_t i=1; i<CUBE_VLEN; ++i){
        
    }




}












/* ----------------------------------------
    Main algorythm
---------------------------------------- */

static inline bool is_border_voxel(
        vs_vec3 start,
        vs_voxelsolve_con con
        )
{
    vs_vec3 verts[8];
    float vals[CUBE_VLEN];
    
    gen_cube(start,con.vscale,verts);

    float first_val = con.f(verts[0], con.farg, con.fargc);

    for(size_t i=1; i<CUBE_VLEN; ++i){
        float vi = con.f(verts[i], con.farg, con.fargc);
        // any value change above or belov
        // f_edge will result in calling 
        // voxel solve method for this point
        if(
            (val > con.f_threshold && first_val < con.f_threshold) || 
            (val < con.f_threshold && first_val > con.f_threshold) 
        ) return true;
    }

    return false;
}





/*
    Iterative solution. Points move slowly towards certain direction,
    looking for f() value changing between ABOVE and BELOW f_threshold
    
    The result is be present as LOCAL UNSCALED coords
    (x,y,z) where x,y,z in (0;1).
    We can calculate real coordinates of the point, but values (0,1)
    are easier to manipulate.
*/
static size_t edge_solve(
        vs_vec3 start,  // starting point (in real coordinates)
        vs_vec3 dir,    // direction (1-vercor)
        vs_vec3 out,    // RESULT. offset with LOCAL UNSCALED coords
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
            (val > con.f_threshold && last_val < con.f_threshold) || 
            (val < con.f_threshold && last_val > con.f_threshold) 
        ){
            float last_val = val; // evaluate again 
                                  // (after another solution found)
            sum_offt = (float)i / (float)con.solve_steps;
            ++solutions;
        }
    }

    if(solutions==0) return 0;

    float avg_offt = (sum_offt / (float)solutions);
    out[0] = start[0] + dir[0] * avg_offt;
    out[1] = start[1] + dir[1] * avg_offt;
    out[2] = start[2] + dir[2] * avg_offt;
    
    return solutions;
}




static size_t add_vertex( 
        vs_voxelsolve_data * data,
        vs_vec3 v,
        float prec,
        )
{
    // compare every existing vertex with new one;
    // if close enough, merge with old vertex.
    for(size_t i=0; i<data->vertex_len; ++i){
        if(cmp3(v,data->vertex[i]) return i;
    }
    
    // if it's really new, add to buffer
    size_t offt = data->vertex_len;
    ++data->vertex_len;
    VS_VEC3_SET(data->vertex[i],v);

    return offt;
}



static int32_t add_triangle(
    vs_voxelsolve_data * data,
    vs_vec3 a,
    vs_vec3 b,
    vs_vec3 c,
    float prec
    )
{
    if(
        cmp3(a,b, prec) ||
        cmp3(b,c, prec) ||
        cmp3(c,a, prec)
    ) return -1;

    size_t ai = add_vertex(data,a, prec);
    size_t bi = add_vertex(data,b, prec);
    size_t ci = add_vertex(data,c, prec);
    size_t offt = data->triangles_len;
    data->triangles[offt][0] = ai;
    data->triangles[offt][1] = bi;
    data->triangles[offt][2] = ci;

    ++data->triangles_len;
    return offt;
}





// check if 2 points share at least one plane 
// (xy) (xz) (yz)
static inline bool is_atplane3(vs_vec3 a, vs_vec3 b, float prec){
    bool ret = false;
    for(size_t i=0; i<3; ++i){
        ret |=  ( (cmpf(a[i],pl0,prec)) && (cmpf(b[i],pl0,prec)) ) ||
                ( (cmpf(a[i],pl1,prec)) && (cmpf(b[i],pl1,prec)) ) ;
    }
    return ret;
}




static inline size_t voxel_solve(
        vs_vec3 start,
        vs_voxelsolve_data * data,
        vs_voxelsolve_con con
        )
{
    // generated 8 vertices 
    vs_vec3 verts[8];
    gen_cube(start,con.vscale,verts);
    
    vs_vec3 points; 
    
}








void voxelsolve_isosurface(
        vs_voxelsolve_data * data,
        vs_voxelsolve_con con
        )
{
    VS_CUBE1 = gen_cube(VS_VEC3_ZERO, );

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



