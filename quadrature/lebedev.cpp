#include <iostream>
//#include <armadillo>
#include "tnt.h"
#include <vector>
#include <math.h>

using namespace TNT;
using namespace std;
//using namespace arma;


//Integrates some function f on the unit sphere which returns a double over the unit sphere in
//spherical coordinates.

//The function f takes two double arguments
//leb_data is n by 3 array.
//phi, theta, weight
double LebedevQuad(double (*f)(double,double), double *leb_data[][3], int n)
{
   double pi, sum;
   pi = 3.1415962;
   sum=0;
   for (int i=0; i<n; i++) {
      sum = sum + *leb_data[i][2] *  f(*leb_data[i][0],*leb_data[i][1]);
      }
   
   sum=sum*4*pi;
   return sum;
}

//Overloaded version
//Lebedev Quad if it needs to be executed numerous times to integrate many times
//The function values need to be evaluated externally and passed by reference
//to the function.
//Parellelized with OpenMP.
//Uses not the function, but an array 
//vector funvals [m number of spatial points][n values (n points of quadrature on sphere)
//compile with option -fopenmp
double LebedevQuad(const vector<vector<double>> &funvals,const double &leb_data[][3],const int m,const int n) {
   
   double sum[m];

   #pragma omp parallel for
   for(int ii=0; ii< m ; ii++) {
      for int(jj=0 ; jj<n ; jj++){
         sum = sum + funvals[ii][jj]*leb_data[ii][jj];
         }
      }
   sum=sum*4*3.1415962;
   return sum;
}

/* Find A(x,y,z)_ij = Int_(S_2)[ outer(a,a) * f(theta, phi) dS]
over spatial problem domain (x,y,z). The integrals are computed with LebedevQuad
Arguments:
X[m][3] coordinates of 'm' points of model domain.
density[m][n] density at 'm' xyz points and 'n' (theta,phi) points of lebedev rule.
leb_data[n][3] (phi, theta, weights) Lebedev rule for 'n' points of quadrature on the sphere.
m <- number of R3 space points.
n <- number of points on sphere for Lebedev quadrature rule.
A[m][6] <- A[[][i]=A_ij for i=j, A[][4]=A12, A[][5]=A13, A[][6]=A23. */

void GetA(const double &X[][3], const vector<vector<double>> &density, const double &leb_data[][3], const int m, const int n, double &A[][6]) {
   
   //First compute the relevant outer products of c-axes at the quadrature points.
   #pragma omp parallel for
   for (ii=0 ii<n; ii++) {
      x=sin(leb_data[ii][0])*cos(leb_data[ii][1]);
      y=sin(leb_data[ii][1])*sin(leb_data[ii][0]);
      z=cos(leb_data[ii][0]);
      cc[ii][0]=x*x;
      cc[ii][1]=y*y;
      cc[ii][2]=z*z;
      cc[ii][3]=x*y;
      cc[ii][4]=x*z;
      cc[ii][5]=y*z;
   }
   //Now use LebedevQuad to compute A[ii][
   #pragma omp parallel for
   for (ii=0; ii<m; ii++) {
      for (jj=0; jj<6; jj++) { 
            for (kk=0;kk<n; kk++) {
               funvals=cc[kk][jj]*density[ii][kk]
               }
            A[ii][ij]= LebedevQuad(&funvals, &leb_data, m,n)
            }
      }
}



void GetB(const double &A[][6], m, n)
   //First diagonalize A
   cube A_diag(m,6,6)
   //get A_diag[m][3]
   for (ii=0,ii<m,i++) {
      A_diag(ii,0,0)=A[ii,0];
      A_diag(ii,1,1)=A[ii,1];
      A_diag(ii,2,2)=A[ii,2];
      A_diag(ii,0,1)=A[ii,3];
      A_diag(ii,0,2)=A[ii,4];
      A_diag(ii,1,2)=A[ii,5];
      }
   
