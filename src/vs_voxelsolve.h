#pragma once


#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>




/* --------------------------------------------
    Main types;
    float 3d vectors should work well 
    with GLM.
-------------------------------------------- */

typedef float vs_vec3[3]; // regular vector that will be used here 
                          // (x,y,z)
                          // Y up

typedef int32_t vs_vec3i[3];
typedef size_t vs_tri[3]; // triangle described by 3 index



/* --------------------------------------------
    Input-output structures
-------------------------------------------- */

typedef struct {
    vs_vec3 offt; 

    float vscale;
    vs_vec3i vox_len;

    size_t solve_steps;

    // any mathematical function 
    // that returns float value;
    float (*f)(
            vs_vec3,
            float *,    // f() additional arguments
            size_t      // f() arguments count
            );       
    float * farg;
    size_t fargc;
    float f_threshold;  // value where line is going
    float prec;         // precision of vector operations
                        // usually don't affect performance but
                        // takes part in float/vector comparisons

} vs_voxelsolve_con;

typedef struct {
    vs_vec3 * vertex;
    size_t vertex_len;

    vs_tri * tris;
    size_t tris_len;
} vs_voxelsolve_data;






/* --------------------------------------------
    Equation solving (dot on the edge)
-------------------------------------------- */

typedef struct {
    size_t cvi[2];          // cube vertex index; start and end of 
                            // the fragment where equation was solved

    size_t intersections;   // number of intersections of the
                            // edge and real surface

    vs_vec3 dot;            // average calculated value
} vs_solution;

// if at leas one cv is shared, return true
static bool vs_solution_sharedcv(vs_solution a, vs_solution b){
    return(
            a.cvi[0] == b.cvi[0] ||
            a.cvi[0] == b.cvi[1] ||
            a.cvi[1] == b.cvi[0] ||
            a.cvi[1] == b.cvi[1]
    );
}





/* ------------------------------------------------
    Math operations and constants
------------------------------------------------ */

const vs_vec3 VS_VEC3_AXIS[3] = {
    {1.0f,0.0f,0.0f},
    {0.0f,1.0f,0.0f},
    {0.0f,0.0f,1.0f}
};


bool vs_cmpf(float a, float b, float prec){
     return (fabsf(a - b) <= prec);
}


void VS_VEC3_SET(vs_vec3 to, const vs_vec3 from){
        to[0] = from[0];
        to[1] = from[1];
        to[2] = from[2];
}


// compares two vec3 with certain precision
bool vs_vec3_cmp(vs_vec3 a, vs_vec3 b, float prec){
    return(
        fabsf(a[0] - b[0]) <= prec &&
        fabsf(a[1] - b[1]) <= prec && 
        fabsf(a[2] - b[2]) <= prec 
    );
}

// linear interpolation with vs_vec3 output
static inline void lerp3(
        const vs_vec3 start,
        const vs_vec3 end,
        float v,
        vs_vec3 out
){
    out[0] = start[0]  +  (end[0] - start[0]) * v;
    out[1] = start[1]  +  (end[1] - start[1]) * v;
    out[2] = start[2]  +  (end[2] - start[2]) * v;
}


static inline float dist3(const vs_vec3 a, const vs_vec3 b){
    return sqrt(
            powf(a[0]-b[0],2.0f) +
            powf(a[1]-b[1],2.0f) +
            powf(a[2]-b[2],2.0f)
    );
}


// calculate normal by 2 vectors
// (cross-product)
void vs_norm3(
    // two vectors
    const vs_vec3 p,
    const vs_vec3 q,

    vs_vec3 out // out 
){
    /*
        we create 3x3 matrix:
            | x  y  z  |
            | px py pz |
            | qx qy qz |
        and calculating it's det
    */
    out[0] = p[1]*q[2] - p[2]*q[1];
    out[1] = p[0]*q[2] - p[2]*q[0];
    out[2] = p[0]*q[1] - p[1]*q[0];
}

// scalar product of 2 vectors 
float vs_dot3(const vs_vec3 a, const vs_vec3 b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline float scalar_diff3(vs_vec3 a, vs_vec3 b){
    return(
            a[0]-b[0] +
            a[1]-b[1] +
            a[2]-b[2]
    );
}





/* ------------------------------------------------
    Voxel/cube 
------------------------------------------------ */

static const size_t CUBE_VLEN = 8;
const vs_vec3 VS_VEC3_ZERO = {0.0f, 0.0f, 0.0f};
static const vs_vec3 CUBE1[8] = {
    // bottom face
    {0.0f, 0.0f, 0.0f},     // 0
    {1.0f, 0.0f, 0.0f},     // 1
    {0.0f, 0.0f, 1.0f},     // 2
    {1.0f, 0.0f, 1.0f},     // 3

    // top face 
    {0.0f, 1.0f, 0.0f},     // 4
    {1.0f, 1.0f, 0.0f},     // 5
    {0.0f, 1.0f, 1.0f},     // 6
    {1.0f, 1.0f, 1.0f},     // 7
};


/*
    Pairs of vertex index (from CUBE1).
    Start-end combination to to describe each cube
    edge (12 in total). Being used to solve the equations.
*/
static const size_t EDGESOLVE_PAIRS[12][2] = {
    {0,1}, {2,3}, {4,5}, {6,7}, // x 0->1
    {0,4}, {1,5}, {2,6}, {3,7}, // y 0->1
    {0,2}, {1,3}, {4,6}, {5,7}, // z 0->1



    /* OLDER VERISON
    {0, 1}, {0, 2}, {0, 4}, // (0,0,0)
    {7, 3}, {7, 5}, {7, 6}, // (1,1,1)
                            
    {4, 5}, {4, 6},         // (0,1,0)
    {3, 1}, {7, 2},         // (1,0,1)
    
    {1, 5},                 // (1,0,0)
    {6, 2}                  // (0,1,1)
    */

}; 

#define CUBE_START 0.0f
#define CUBE_END 1.0f

// create 8 vertex that describe a cube
// out array will have the length of 8
static inline void gen_cube(
        const vs_vec3 start,
        const float vscale,
        vs_vec3 * out
        ){
    
    // basically, it's just offseted and scaled cube-1
    for(size_t i=0; i<CUBE_VLEN; ++i){
        out[i][0] = start[0] + CUBE1[0][0] * vscale;
        out[i][1] = start[1] + CUBE1[0][1] * vscale;
        out[i][2] = start[2] + CUBE1[0][2] * vscale;
    }
}

static inline void local_to_real3(
        const vs_vec3 local,
        vs_vec3 out,
        vs_voxelsolve_con con
){
    out[0] = local[0]*con.vscale + con.offt[0];
    out[1] = local[1]*con.vscale + con.offt[1];
    out[2] = local[2]*con.vscale + con.offt[2];
}






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
        (ret_offt+0) * ( !vs_cmpf(v[0],CUBE_START,prec) 
                && !vs_cmpf(v[0],CUBE_END,prec) ) + // x

        (ret_offt+1) * ( !vs_cmpf(v[1],CUBE_START,prec) 
                && !vs_cmpf(v[1],CUBE_END,prec) ) + // y
        
        (ret_offt+2) * ( !vs_cmpf(v[2],CUBE_START,prec) 
                && !vs_cmpf(v[2],CUBE_END,prec) )   // z
    );

    if( ret-ret_offt > 2 ) return -1;
    else return ret - ret_offt;
}



static inline bool is_shared_edge(vs_vec3 a, vs_vec3 b, float prec){
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




// check if 2 points are located at the same cube face
// (xy) (xz) (yz)
// (one of the coords matches)
static inline bool is_shared_face(
        const vs_vec3 a,
        const vs_vec3 b,
        float prec
){
    bool ret = false;
    for(size_t i=0; i<3; ++i){
        ret |= ( (vs_cmpf(a[i],CUBE_START,prec)) && (vs_cmpf(b[i],CUBE_START,prec)) ) ||
               ( (vs_cmpf(a[i],CUBE_END,prec)) && (vs_cmpf(b[i],CUBE_END,prec)) ) ;
    }
    return ret;
}






/* --------------------------------------------
    Additional algorythm steps 
-------------------------------------------- */

static inline bool is_border_voxel(
        vs_vec3 start,
        vs_voxelsolve_con con
        )
{

    float firstv = con.f(start, con.farg, con.fargc);

    for(size_t i=1; i<CUBE_VLEN; ++i){
        vs_vec3 p;
        p[0] = CUBE1[i][0] * con.vscale;
        p[1] = CUBE1[i][1] * con.vscale;
        p[2] = CUBE1[i][2] * con.vscale;

        float vi = con.f(p, con.farg, con.fargc);
        // any value change above or belov
        // f_edge will result in calling 
        // voxel solve method for this point
        if(
            (vi > con.f_threshold && firstv < con.f_threshold) || 
            (vi < con.f_threshold && firstv > con.f_threshold) 
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
static void edge_solve(
        size_t starti,          // starting point (in real coordinates)
        size_t endi,            // direction (1-vercor)
        
        vs_solution * out,      // RESULT. offset with LOCAL UNSCALED coords
        
        vs_voxelsolve_con con   // config is needed here for scale, offt and f()
        )
{  

    size_t intersections = 0;
    float sum_offt = 0.0f;      // summ of all solutions.
                                // each solution described
                                // by value in range [0;1]

    float step_size = 1.0f / (float)con.solve_steps; // local scale of the step

    vs_vec3 dot;        // iterable value
    vs_vec3 dot_real;   // f(dot_real)


    // calling the function for the first time
    local_to_real3(CUBE1[starti],dot_real,con);
    float val = con.f(dot_real, con.farg, con.fargc) - con.f_threshold;
    float last_val = val;

    // start iterration
    for (size_t i=1; i<=con.solve_steps+1; ++i){
        // get the coordinates of the point
        float progress = step_size * (float)i;
        lerp3(CUBE1[starti],CUBE1[endi], progress, dot);
        local_to_real3(dot,dot_real,con);

        float val = con.f(dot_real, con.farg, con.fargc) - con.f_threshold;
        bool opposite = (val>=0.0f && last_val<0.0f) ||
                        (val<0.0f && last_val>=0.0f);
        if(opposite){
            ++intersections;
            sum_offt += progress;
        }

        last_val = val; // evaluate again 
    }


    float avg_offt = (sum_offt / (float)intersections);

    // form solution
    lerp3(CUBE1[starti],CUBE1[endi], avg_offt, out->dot);
    out->cvi[0] = starti;
    out->cvi[1] = endi;
    out->intersections = intersections;
    
}



    
static size_t add_vertex( 
        vs_voxelsolve_data * data,
        vs_vec3 v,
        float prec
        )
{
    // compare every existing vertex with new one;
    // if close enough, merge with old vertex.
    for(size_t i=0; i<data->vertex_len; ++i){
        if(vs_vec3_cmp(v,data->vertex[i],prec)){
            return i;
        }
    }
    
    // if it's really new, add to buffer
    size_t offt = data->vertex_len;
    ++data->vertex_len;
    //printf("\tdata->vertex_len=%d\n",++data->vertex_len);
    data->vertex[offt][0] = v[0];
    data->vertex[offt][1] = v[1];
    data->vertex[offt][2] = v[2];

    return offt;
}



static int32_t add_triangle(
    vs_voxelsolve_data * data,
    vs_vec3 start,
    vs_vec3 a,
    vs_vec3 b,
    vs_vec3 c,
    vs_voxelsolve_con con
    )
{
    if(
        vs_vec3_cmp(a,b, con.prec) ||
        vs_vec3_cmp(b,c, con.prec) ||
        vs_vec3_cmp(c,a, con.prec)
    ) {
        printf("unnatural triangle intersection!");
        return -1;
    }
    
    
    vs_vec3 a_real = {
        a[0] * con.vscale + start[0],
        a[1] * con.vscale + start[1],
        a[2] * con.vscale + start[2]
    };
    vs_vec3 b_real = {
        b[0] * con.vscale + start[0],
        b[1] * con.vscale + start[1],
        b[2] * con.vscale + start[2]
    };
    vs_vec3 c_real = {
        c[0] * con.vscale + start[0],
        c[1] * con.vscale + start[1],
        c[2] * con.vscale + start[2]
    };


    size_t ai = add_vertex(data,a_real, con.prec);
    size_t bi = add_vertex(data,b_real, con.prec);
    size_t ci = add_vertex(data,c_real, con.prec);

    
    size_t offt = data->tris_len;
    data->tris[offt][0] = ai;
    data->tris[offt][1] = bi;
    data->tris[offt][2] = ci;

    /*printf("%d\tadded triangle[%d]: %d %d %d\n",
        data->tris[offt],
        offt,
        ai,bi,ci
    );*/

    ++data->tris_len;
    return offt;
}




/* ----------------------------------------
    Triangulation   
---------------------------------------- */

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

        6. Choose a dot[i] and the nearest neighbour dots[j] ONLY
           from the neighbour faces (not edges this time). Build a triangle.
           EXCLUDE all triangulated dots[j] from i-iteration 
           (they still can be connected).
           EXCLUDE current dot[i] from j-iteration (cannot be connected at all)

        7. End of algorythm.
*/
static void dots_triang(
        vs_voxelsolve_data * data,
        vs_vec3 start,
        vs_solution * sol,          // solutions 
        size_t sol_len,       
        vs_voxelsolve_con con
        )
{   
    float prec = con.prec; 

    // unable to build figure
    if(sol_len < 3) return;

    // simpliest variant of the algorythm *(add only 1 trangle)
    else if(sol_len == 3) {
        add_triangle(data, start, sol[0].dot, sol[1].dot, sol[2].dot, con);
        return;
    }
    else if(sol_len == 4) {
        add_triangle(data, start, sol[0].dot, sol[1].dot, sol[2].dot, con);
        add_triangle(data, start, sol[3].dot, sol[2].dot, sol[1].dot, con);
        return;
    }

    
    bool queue [8] = { [0 ... 7] = true }; // unprocessed dots to connect
    
    // 3-TRIANGULATION (seems like correct)
    // find all triangles on the edges
    for(size_t i=0; i<(sol_len-2); ++i){
        size_t nbrs[3] = {i,0,0}; // neighbours of the dot (including itself)
        size_t nbrs_len = 1;

        //look for all dots able to 3-triangulate
        for(size_t j=i+1; j<sol_len; ++j){
            if(nbrs_len == 3 || !queue[i]) continue; //end dot[i] processing

            // add all dots on the same edge to 3-triang queue
            if( 
                    vs_solution_sharedcv(sol[i],sol[j]) &&
                    queue[j]
            ){

                //push index in array
                nbrs[nbrs_len] = j;
                ++nbrs_len;
            }
        }

        // if found all 3 dots, triangulate
        // and disable dots in further algorythm
        if(nbrs_len == 3){

            // no longer avilable
            queue[ nbrs[0] ] = false;
            queue[ nbrs[1] ] = false;
            queue[ nbrs[2] ] = false;

            add_triangle(
                    data, 
                    start,
                    sol[nbrs[0]].dot, sol[nbrs[1]].dot, sol[nbrs[2]].dot,
                    con
                    );
        }
    }
    //return; 
    //printf("\t4-triang -----------\n");
    // 4-TRIANGULATION 
    // find and triangulate the rest
    for(size_t i=0; i<sol_len; ++i){        
        int32_t tri[3] = {i,-1,-1}; // basically the copy of previous step
        if( !queue[i] ) continue; //ignore already tiangulated dots        
        

        // try to find all triangle components
        for(size_t j=i+1; j<sol_len; ++j){
            if( !queue[j] ) continue; // idk what I am doing at this point
                                      // do that dudes really have to be ignored?(((
            
            // look for one neighbour-face dot 
            if( tri[1]<0 && is_shared_face(sol[i].dot,sol[j].dot, prec) ){
                tri[1] = j;
            }
            
            // look for one non-neighbour-face dot 
            if( tri[2]<0 && !is_shared_face(sol[i].dot,sol[j].dot,prec) ){
                tri[2] = j;
            }
        }
        
        // able to build a triangle
        if( tri[1]>=0 && tri[2]>=0 ){
            queue[i] = false;
            add_triangle(
                    data, 
                    start,
                    sol[tri[0]].dot, sol[tri[1]].dot, sol[tri[2]].dot,
                    con
                    );
        }

    }



}






/* ----------------------------------------
    Main algorythm    
---------------------------------------- */

static inline size_t voxel_solve(
        vs_vec3 start,
        vs_voxelsolve_data * data,
        vs_voxelsolve_con con
        )
{
    // generated 8 vertices 
    vs_vec3 verts[8];
    gen_cube(start,con.vscale,verts);
    
    vs_solution dots[12];
    size_t dots_len = 0;
    for(size_t i=0; i<12; ++i){
        //printf("edge solve\n");
        vs_solution sol;
        edge_solve(
                EDGESOLVE_PAIRS[i][0],
                EDGESOLVE_PAIRS[i][1],
                &sol,
                con
        );
        // if found, add to other dots
        if(sol.intersections > 0){
            dots[dots_len] = sol;
            ++dots_len;
        }
    }
    dots_triang(data, start, dots, dots_len, con);
    //printf("triangulated. vertices:%d\n", data->vertex_len);  
    
}



void voxelsolve(
        vs_voxelsolve_data * data,
        vs_voxelsolve_con con
        )
{
    size_t voxelc = con.vox_len[0] * con.vox_len[1] * con.vox_len[2];
    
    data->vertex_len = 0;
    data->vertex = (vs_vec3*) malloc(
            voxelc *        // assume all voxels have some shape
            8 *             // assume 8 edges are taken for reach cube
            sizeof(vs_vec3)
    );
    
    data->tris_len = 0;
    data->tris = (vs_tri*) malloc(
            voxelc *        // assume all voxels have some shape
            4 *             // assume 4 triangles for cube
            sizeof(vs_tri)
    );

    vs_vec3 doti;

    for(size_t xi=0; xi<con.vox_len[0];++xi){
        for(size_t yi=0; yi<con.vox_len[1];++yi){
            for(size_t zi=0; zi<con.vox_len[2];++zi){

                //iterable dot
                doti[0] = con.offt[0] + xi*con.vscale;
                doti[1] = con.offt[1] + yi*con.vscale;
                doti[2] = con.offt[2] + zi*con.vscale;

                if(is_border_voxel(doti,con)){
                    voxel_solve(doti, data, con);
                }
            }
        }
    }
}






