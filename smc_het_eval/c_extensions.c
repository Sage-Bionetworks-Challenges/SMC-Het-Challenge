#include <omp.h>
//This is the first loop optimization in def filterFPs in file ___ 
//This loop does: 
// the elements at the indicies specified by mask are "picked" and assembled into a matrix
// that grows out of the upper-left corner of the original matrix (the original matrix will
// always be bigger than the eventual masked matrix)
void filterFPs_remove_fps(double *x,int dimx1, int dimx2, double  *mask, int lenmask){
    #pragma omp parallel for simd
    for (int i = 0; i < lenmask; i++){
        for (int j = 0; j < lenmask; j++){
            int mi = mask[i];
            int mj = mask[j];
            int index = i * dimx2 + j;
            int mindex = mi * dimx2 + mj;

            x[index] = x[mindex];
            
        }
    }
}


// int main(){
//     double  x[4][4] = {{0,0,0,0},{1,1,1,1},{2,2,2,2,},{3,3,3,3}};
//     double mask[4] = {1,3,3,2};
//     int lenmask = 4;

//     filterFPs_remove_fps(*x,4,4,mask,lenmask);

//     return(3);
// }

