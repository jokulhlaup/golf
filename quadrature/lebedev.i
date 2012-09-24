%module fec
%include "std_vector.i"
namespace std {
   %template(DoubleVector) vector<double>;
}
%{
extern void LebedevQuad(double *funvals,const vector< double >leb_data,int n, double *sum); 
extern void GetA(const vector< double > density, const vector<double> leb_data, const int m, const int n, double *A);
%}
extern void LebedevQuad(double *funvals,const vector< double >leb_data,int n, double *sum); 
extern void GetA(const vector< double > density, const vector<double> leb_data, const int m, const int n, double *A);

