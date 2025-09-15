#include <stdint.h>
#include <stdio.h>
#include <math.h>






typedef float vs_vec3[3];



// calculate normal for triangle 
void norm3(
    // two vectors
    vs_vec3 p,
    vs_vec3 q,

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



int main(){
    vs_vec3 a; 
    vs_vec3 b;

    a[0]=-1.0f; a[1]=0.0f; a[2]=0.0f;
    b[0]=0.0f; b[1]=0.0f; b[2]=1.0f;

    vs_vec3 n;
    norm3(a,b, n);

    printf("normal: %f %f %f\n",    n[0], n[1], n[2]);

}
