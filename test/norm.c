#include <stdint.h>
#include <stdio.h>
#include <math.h>






typedef float vs_vec3[3];



// calculate normal for triangle 
void norm3(
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



float dot3(const vs_vec3 a, const vs_vec3 b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


static const vs_vec3 VS_XUP = {1.0f,0.0f,0.0f};
static const vs_vec3 VS_YUP = {0.0f,1.0f,0.0f};
static const vs_vec3 VS_ZUP = {0.0f,0.0f,1.0f};



int main(){
    vs_vec3 a; 
    vs_vec3 b;

    a[0]=-1.0f; a[1]=0.0f; a[2]=2.0f;
    b[0]=0.0f; b[1]=0.0f; b[2]=1.0f;

    

    vs_vec3 n;
    norm3(a,b, n);
    
    printf("a: %f %f %f\n",    a[0], a[1], a[2]);
    printf("b: %f %f %f\n",    b[0], b[1], b[2]);

    printf("\ndot(a,b): %f\n", dot3(a,b));
    printf("normal: %f %f %f\n",    n[0], n[1], n[2]);

    printf("\ndot(n,a): %f\n", dot3(n,a));

    printf("\ndot(n,X): %f\n", dot3(n,VS_XUP));
    printf("dot(n,Y): %f\n", dot3(n,VS_YUP));
    printf("dot(n,Z): %f\n", dot3(n,VS_ZUP));

}
