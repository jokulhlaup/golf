%module lebedev
%{
extern double LebedevQuad(double (*f)(double,double),double *leb_data, int n);
%}
extern double LebedevQuad(double (*f)(double,double),double *leb_data, int n);
