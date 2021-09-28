
//******************************************************************************************
// TestHDind.cpp: Simulate claims for the d=5000 (high deductible) case for the independence model.
//******************************************************************************************
  
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <random>
#include <vinecopulib.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <string> 

using namespace std;
using namespace Eigen;

#define MAXBUFSIZE  ((int) 1e6)
#define BB 1000 // number of replicates.
#define NN 200 // maximum number of claims for a policyholder in a year.
#define LARGENUMBER 1.0e+10

int write_data(const char *filename) {
  std::ofstream file(filename);
  if (file.is_open()) {
    Eigen::Matrix<double, 4, 2, Eigen::DontAlign> m;
    m << 1, 2, 3, 4, 5, 6, 7, 8;
    file << m << std::endl;
  }
  return(0);
}

MatrixXd readMatrix(const char *filename) {
  int cols = 0, rows = 0;
  double * buff = new double[MAXBUFSIZE];
  
  // Read numbers from file into buffer.
  ifstream infile;
  infile.open(filename);
  while (! infile.eof()) {
    string line;
    getline(infile, line);
    
    int temp_cols = 0;
    stringstream stream(line);
    while(! stream.eof())
      stream >> buff[cols*rows+temp_cols++];
    
    if (temp_cols == 0)
      continue;
    
    if (cols == 0)
      cols = temp_cols;
    
    rows++;
  }
  infile.close();
  rows--;
  
  // Populate matrix with numbers.
  MatrixXd result(rows,cols);
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      result(i,j) = buff[ cols*i+j ];
  
  delete [] buff;
  
  return result;
}

MatrixXi readMatrixInt(const char *filename) {
  int cols = 0, rows = 0;
  double * buff = new double[MAXBUFSIZE];
  
  // Read numbers from file into buffer.
  ifstream infile;
  infile.open(filename);
  while (! infile.eof()) {
    string line;
    getline(infile, line);
    
    int temp_cols = 0;
    stringstream stream(line);
    while(! stream.eof())
      stream >> buff[cols*rows+temp_cols++];
    
    if (temp_cols == 0)
      continue;
    
    if (cols == 0)
      cols = temp_cols;
    
    rows++;
  }
  infile.close();
  rows--;
  
  // Populate matrix with numbers.
  MatrixXi result(rows,cols);
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      result(i,j) = buff[ cols*i+j ];
  
  delete [] buff;
  
  return result;
}


double pGB2(double y, double mu, double sigma, double alpha1, double alpha2) {
  double ndf = 2*alpha1;
  double ddf = 2*alpha2;
  double r = (log(y)-mu)/sigma;
  double z = (alpha2/alpha1)*exp(r);
  boost::math::fisher_f f_func(ndf,ddf);
  double ans =  cdf(f_func,z);
  return(ans);
}

double Fcount(double y, double pzero, double pone, double meannb, double size) {
  boost::math::negative_binomial_distribution<double> nb_dist(size, size/(size+meannb));
  double ans = pzero*(y>=0)+pone*(y>=1)+(1-pzero-pone)*cdf(nb_dist,y);
  return(ans);
}


struct obj1dep_params {
  double u;
  double mu1;
  double sigma1;
  double alpha11;
  double alpha21;
  double copDNpar1;
  double copDNpar2;
  double y1;
  double y2dep;
  double pzero;
  double pone;
  double meannb;
  double size;
  double mu3dep;
  double sigma3dep;
  double alpha13dep;
  double alpha23dep;
  double copRYpar1;
  double copRYpar2;
  double copNYpar1;
  double copNYpar2;
};

double obj1dep(double y3, void *params) {
  struct obj1dep_params *p = (struct obj1dep_params *) params;
  double u = p->u;
  double mu1 = p->mu1;
  double sigma1 = p->sigma1;
  double alpha11 = p->alpha11;
  double alpha21 = p->alpha21;
  double y1 = p->y1;
  double y2dep = p->y2dep;
  double pzero = p->pzero;
  double pone = p->pone;
  double meannb = p->meannb;
  double size = p->size;
  double mu3dep = p->mu3dep;
  double sigma3dep = p->sigma3dep;
  double alpha13dep = p->alpha13dep;
  double alpha23dep = p->alpha23dep;
  Matrix<double, 1, 2> utemp; 
  Vector2d copDNparam; copDNparam << p->copDNpar1, p->copDNpar2;
  Vector2d copRYparam; copRYparam << p->copRYpar1, p->copRYpar2;
  Vector2d copNYparam; copNYparam << p->copNYpar1, p->copNYpar2;
  
  //vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::student, 0, copDNparam);
  //vinecopulib::Bicop cop_model_RY(vinecopulib::BicopFamily::frank, 0, copRYparam.head(1));
  //vinecopulib::Bicop cop_model_NY(vinecopulib::BicopFamily::gumbel, 90, -copNYparam.head(1));
  vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::student, 0, copDNparam);
  vinecopulib::Bicop cop_model_RY(vinecopulib::BicopFamily::frank, 0, copRYparam.head(1));
  vinecopulib::Bicop cop_model_NY(vinecopulib::BicopFamily::gumbel, 90, -copNYparam.head(1));
  
  // frequency|deductible
  utemp(0,0) = pGB2(y1,mu1,sigma1,alpha11,alpha21); utemp(0,1) = Fcount(y2dep,pzero,pone,meannb,size);   VectorXd tempvec = cop_model_DN.hfunc1(utemp); double u1 = tempvec(0);
  utemp(0,1) = Fcount(y2dep-1,pzero,pone,meannb,size); tempvec = cop_model_DN.hfunc1(utemp); double u1a = tempvec(0);
  // severity|deductible
  utemp(0,1) = pGB2(y3,mu3dep,sigma3dep,alpha13dep,alpha23dep);
  tempvec = cop_model_RY.hfunc1(utemp);
  double u2 = tempvec(0);
  // numerator
  utemp(0,0) = u1; utemp(0,1) = u2;
  tempvec = cop_model_NY.cdf(utemp);
  double numerator = tempvec(0);
  utemp(0,0) = u1a;
  tempvec = cop_model_NY.cdf(utemp);
  numerator = numerator - tempvec(0);
  // denominator
  utemp(0,0) = pGB2(y1,mu1,sigma1,alpha11,alpha21); utemp(0,1) = Fcount(y2dep,pzero,pone,meannb,size);
  tempvec = cop_model_DN.hfunc1(utemp);
  double denominator = tempvec(0);
  utemp(0,1) = Fcount(y2dep-1,pzero,pone,meannb,size);
  tempvec = cop_model_DN.hfunc1(utemp);
  denominator = denominator - tempvec(0);
  return(numerator/denominator - u);
  // return (pGB2(y3,mu1,sigma1,alpha11,alpha21) - u);
}


double obj1ind(double y3, void *params) {
  struct obj1dep_params *p = (struct obj1dep_params *) params;
  double u = p->u;
  double mu1 = p->mu1;
  double sigma1 = p->sigma1;
  double alpha11 = p->alpha11;
  double alpha21 = p->alpha21;
  double y1 = p->y1;
  double y2dep = p->y2dep;
  double pzero = p->pzero;
  double pone = p->pone;
  double meannb = p->meannb;
  double size = p->size;
  double mu3dep = p->mu3dep;
  double sigma3dep = p->sigma3dep;
  double alpha13dep = p->alpha13dep;
  double alpha23dep = p->alpha23dep;
  Matrix<double, 1, 2> utemp; 
  Vector2d copDNparam; copDNparam << p->copDNpar1, p->copDNpar2;
  Vector2d copRYparam; copRYparam << p->copRYpar1, p->copRYpar2;
  Vector2d copNYparam; copNYparam << p->copNYpar1, p->copNYpar2;
  vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::indep);
  vinecopulib::Bicop cop_model_RY(vinecopulib::BicopFamily::indep);
  vinecopulib::Bicop cop_model_NY(vinecopulib::BicopFamily::indep);
  // frequency|deductible
  utemp(0,0) = pGB2(y1,mu1,sigma1,alpha11,alpha21); utemp(0,1) = Fcount(y2dep,pzero,pone,meannb,size);   VectorXd tempvec = cop_model_DN.hfunc1(utemp); double u1 = tempvec(0);
  utemp(0,1) = Fcount(y2dep-1,pzero,pone,meannb,size); tempvec = cop_model_DN.hfunc1(utemp); double u1a = tempvec(0);
  // severity|deductible
  utemp(0,1) = pGB2(y3,mu3dep,sigma3dep,alpha13dep,alpha23dep);
  tempvec = cop_model_RY.hfunc1(utemp);
  double u2 = tempvec(0);
  // numerator
  utemp(0,0) = u1; utemp(0,1) = u2;
  tempvec = cop_model_NY.cdf(utemp);
  double numerator = tempvec(0);
  utemp(0,0) = u1a;
  tempvec = cop_model_NY.cdf(utemp);
  numerator = numerator - tempvec(0);
  // denominator
  utemp(0,0) = pGB2(y1,mu1,sigma1,alpha11,alpha21); utemp(0,1) = Fcount(y2dep,pzero,pone,meannb,size);
  tempvec = cop_model_DN.hfunc1(utemp);
  double denominator = tempvec(0);
  utemp(0,1) = Fcount(y2dep-1,pzero,pone,meannb,size);
  tempvec = cop_model_DN.hfunc1(utemp);
  denominator = denominator - tempvec(0);
  return(numerator/denominator - u);
  // return (pGB2(y3,mu1,sigma1,alpha11,alpha21) - u);
}



VectorXd SimPolicy(double y1, VectorXd cova, VectorXd cova1, VectorXd gam1, VectorXd gam2, VectorXd gam3,
                   VectorXd copDN, VectorXd copDS, VectorXd copDScind, double coverage) {
  
  // cout << "Simulating a policyholder with parameters: " << y1 << endl;
  
  // Variables related to the solver.
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r, x_lo, x_hi;
  gsl_function F;
  struct obj1dep_params paramsdep;  
  int status, iter, max_iter;
  int BB2 = BB/2;
  
  int k1 = cova1.size();
  int k = cova.size();
  Vector2d cova_zero; 
  Vector2d cova_one;
  cova_zero(0) = cova(0); cova_zero(1) = cova(k-1);
  cova_one(0) = cova(0); cova_one(1) = cova(k-1);
  
  // Dependence parameters
  double rhoRY = copDS[2+k+3+1-1];
  double thetaNY = copDS[2+k+3+2-1];
  double rhoRYcind = copDScind[1+k+3+1-1];
  double fam1 = copDS[0];
  double fam2 = copDS[1];
  double fam1cind = copDScind[0];
  double df1 = 5;
  double df2 = 5;
  double df1cind = 5;
  
  // Deductible parameters
  double mu1 = cova1.dot(gam1.head(k1));
  double sigma1 = gam1(k1);
  double alpha11 = gam1(k1+1);
  double alpha21 = gam1(k1+2);
  
  // Frequency parameters
  VectorXd beta = gam2.head(k);
  VectorXd gamma0 = gam2.segment(k,2);
  VectorXd gamma1 = gam2.segment(k+2,2);
  double size = gam2(k+2+2);
  
  double tempvec0 = exp(cova_zero.dot(gamma0));
  double tempvec1 = exp(cova_one.dot(gamma1)); 
  double pzero = tempvec0; pzero = pzero / (1+tempvec0 + tempvec1);
  double pone = tempvec1; pone = pone / (1+tempvec0 + tempvec1);
  double meannb = exp(cova.dot(beta)); 
  
  // New parameters for severity, after full-likelihood estimation.
  VectorXd gam3dep = copDS.segment(2,k+3);
  VectorXd gam3cind = copDScind.segment(2,k+3);
  
  double mu3dep = cova.dot(gam3dep.head(k));
  double sigma3dep = gam3dep(k);
  double alpha13dep = gam3dep(k+1);
  double alpha23dep = gam3dep(k+2);
  
  double mu3cind = cova.dot(gam3cind.head(k));
  double sigma3cind = gam3cind(k);
  double alpha13cind = gam3cind(k+1);
  double alpha23cind = gam3cind(k+2);
  
  double mu3ind = cova.dot(gam3.head(k));
  double sigma3ind = gam3(k);
  double alpha13ind = gam3(k+1);
  double alpha23ind = gam3(k+2);
  
  // std::default_random_engine generator(time(0));
  // std::uniform_real_distribution<double> uni_dist(0.0,1.0);
  std::uniform_real_distribution<double> uni_dist(0.0,1.0);
  std::random_device rd;
  std::default_random_engine generator(rd());
  
  // cout << "copDN.tail(2): " << copDN.tail(2) << endl;
  vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::student, 0, copDN.tail(2));
  // Matrix<double, 1, 1> tt; tt << 0.4;
  // vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::gaussian, 0, tt);
  
  // Simulating claims for a single policyholder.
  double uR = pGB2(y1,mu1,sigma1,alpha11,alpha21);
  boost::math::negative_binomial_distribution<double> nb_dist(size, size/(size+meannb));
  double tempval1 = pdf(nb_dist,0); 
  double tempval2 = pdf(nb_dist,1);
  double pp0 = pzero + (1-pzero-pone)*tempval1;
  double pp1 = pzero + (1-pzero-pone)*tempval1 + pone + (1-pzero-pone)*tempval2;
  
  MatrixXd y3(BB,NN);
  MatrixXd buff(1,NN);
  VectorXd y3dep_mean;
  VectorXd pureprem_dep;
  pureprem_dep.resize(BB*2);
  
  // simulate a large number of observations.
  for (int b=0; b<BB; b++) {
    
    // 1. Dependent case.
    for (int k=0; k<NN; k++) {
      y3(b,k) = 0.0;
    }
    double uYind = uni_dist(generator);
    Matrix<double, 1, 2> uRuYind; uRuYind << uR, uYind;
    VectorXd tempvec = cop_model_DN.hinv1(uRuYind);
    double uYdep = tempvec(0);
    double y2ind, y2dep, y2cind;
    // sample from a 01-negative binomial distribution.
    if (uYind < pp0) y2ind = 0; else
      if (uYind < pp1) y2ind = 1; else
        y2ind = quantile( nb_dist, pdf(nb_dist,0)+pdf(nb_dist,1) + (uYind-pp1)*(1-pdf(nb_dist,0)-pdf(nb_dist,1))/(1-pp1) );
    // sample from a 01-negative binomial distribution.
    if (uYdep < pp0) y2dep = 0; else
      if (uYdep < pp1) y2dep = 1; else
        y2dep = quantile( nb_dist, pdf(nb_dist,0)+pdf(nb_dist,1) + (uYdep-pp1)*(1-pdf(nb_dist,0)-pdf(nb_dist,1))/(1-pp1) );
    
    y2cind = y2dep;
    
    // Assume frequency doesn't exceed the buffer size.
    if (y2ind>NN) y2ind = NN;
    if (y2dep>NN) y2dep = NN;
    if (y2cind>NN) y2cind = NN;
    
    // cout << "uR: " << uR << ", uYind: " << uYind << ", uYdep: " << uYdep << ", y2dep: " << y2dep << endl;
    if (y2dep>0) {
    
    // Figure out the maximum of the uniform distribution that we can accomodate.
    // paramsdep = {0, mu1, sigma1, alpha11, alpha21, copDN(1), copDN(2), y1, y2dep, pzero, pone, meannb, size, mu3dep, sigma3dep, alpha13dep, alpha23dep, 1, df1, -1, df2}; 
    paramsdep = {0, mu1, sigma1, alpha11, alpha21, copDN(1), copDN(2), y1, y2dep, pzero, pone, meannb, size, mu3dep, sigma3dep, alpha13dep, alpha23dep, rhoRY, df1, thetaNY, df2}; 
    F.function = &obj1dep;
    F.params = &paramsdep;
    double MaxU = F.function(coverage*1000000,&paramsdep);
    for (int j=0; j<y2dep; j++) { buff(0,j) = uni_dist(generator); if(buff(0,j)>MaxU) {buff(0,j)=MaxU;}}
    for (int j=0; j<y2dep; j++) {
    iter = 0; max_iter = 100;
    r = 0.1; 
    x_lo = 0.0;
    x_hi = LARGENUMBER;
    // paramsdep = {buff(0,j), mu1, sigma1, alpha11, alpha21, copDN(1), copDN(2), y1, y2dep, pzero, pone, meannb, size, mu3dep, sigma3dep, alpha13dep, alpha23dep, 1, df1, -1, df2};
    paramsdep = {buff(0,j), mu1, sigma1, alpha11, alpha21, copDN(1), copDN(2), y1, y2dep, pzero, pone, meannb, size, mu3dep, sigma3dep, alpha13dep, alpha23dep, rhoRY, df1, thetaNY, df2};
    F.function = &obj1dep;
    F.params = &paramsdep;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    // printf ("using %s method\n",gsl_root_fsolver_name (s));
    // printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "root", "err(est)");
    do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, 0.01);
    // if (status == GSL_SUCCESS) printf ("Converged:\n");
    // printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, r, x_hi - x_lo);
    // if (status == GSL_SUCCESS) cout << "Converged: " << r << endl;
    
    if (status == GSL_SUCCESS) y3(b,j) = r;
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    }
    }
    
    // Sum g(S) instead of S. In other words, apply the deductible.
    pureprem_dep(b) = 0;
    for (int j=0; j<y2dep; j++) {
    if (y3(b,j)>coverage*1000000)
    pureprem_dep(b) = pureprem_dep(b) + coverage*1000000; // (Error fixed) This was the error. We should multiply 1000,000 to the coverage.
    else {
    pureprem_dep(b) = pureprem_dep(b) + y3(b,j);
    }
    }
    
    // Frequencies.
    pureprem_dep(b+BB) = 0;
    for (int j=0; j<y2dep; j++) {
      pureprem_dep(b+BB) = pureprem_dep(b+BB) + 1;
    }
    
    //pureprem_dep(b) = 0;
    //for (int j=0; j<y2dep; j++) {
    //  if (y3(b,j)>coverage*1000000)
    //    pureprem_dep(b) = pureprem_dep(b) + coverage*1000000; // (Error fixed) This was the error. We should multiply 1000,000 to the coverage.
    //    // pureprem_dep(b) = pureprem_dep(b) + 1; // Test
    //  else
    //    pureprem_dep(b) = pureprem_dep(b) + y3(b,j);
    //    // pureprem_dep(b) = pureprem_dep(b) + 1; // Test
    //}
  }
    
    return(pureprem_dep);
    
}
  
  
  
  VectorXd SimPolicyInd(double y1, VectorXd cova, VectorXd cova1, VectorXd gam1, VectorXd gam2, VectorXd gam3,
  VectorXd copDN, VectorXd copDS, VectorXd copDScind, double coverage) {
  
  // Variables related to the solver.
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r, x_lo, x_hi;
  gsl_function F;
  struct obj1dep_params paramsind;  
  int status, iter, max_iter;
  int BB2 = BB/2;
  
  int k1 = cova1.size();
  int k = cova.size();
  Vector2d cova_zero; 
  Vector2d cova_one;
  cova_zero(0) = cova(0); cova_zero(1) = cova(k-1);
  cova_one(0) = cova(0); cova_one(1) = cova(k-1);
  
  // Dependence parameters
  double rhoRY = copDS[2+k+3+1-1];
  double thetaNY = copDS[2+k+3+2-1];
  double rhoRYcind = copDScind[1+k+3+1-1];
  double fam1 = copDS[0];
  double fam2 = copDS[1];
  double fam1cind = copDScind[0];
  double df1 = 5;
  double df2 = 5;
  double df1cind = 5;
  
  // Deductible parameters
  double mu1 = cova1.dot(gam1.head(k1));
  double sigma1 = gam1(k1);
  double alpha11 = gam1(k1+1);
  double alpha21 = gam1(k1+2);
  
  // Frequency parameters
  VectorXd beta = gam2.head(k);
  VectorXd gamma0 = gam2.segment(k,2);
  VectorXd gamma1 = gam2.segment(k+2,2);
  double size = gam2(k+2+2);
  
  double tempvec0 = exp(cova_zero.dot(gamma0));
  double tempvec1 = exp(cova_one.dot(gamma1)); 
  double pzero = tempvec0; pzero = pzero / (1+tempvec0 + tempvec1);
  double pone = tempvec1; pone = pone / (1+tempvec0 + tempvec1);
  double meannb = exp(cova.dot(beta)); 
  
  // New parameters for severity, after full-likelihood estimation.
  VectorXd gam3dep = copDS.segment(2,k+3);
  VectorXd gam3cind = copDScind.segment(2,k+3);
  
  double mu3dep = cova.dot(gam3dep.head(k));
  double sigma3dep = gam3dep(k);
  double alpha13dep = gam3dep(k+1);
  double alpha23dep = gam3dep(k+2);
  
  double mu3cind = cova.dot(gam3cind.head(k));
  double sigma3cind = gam3cind(k);
  double alpha13cind = gam3cind(k+1);
  double alpha23cind = gam3cind(k+2);
  
  double mu3ind = cova.dot(gam3.head(k));
  double sigma3ind = gam3(k);
  double alpha13ind = gam3(k+1);
  double alpha23ind = gam3(k+2);
  
  // std::default_random_engine generator(time(0));
  // std::uniform_real_distribution<double> uni_dist(0.0,1.0);
  std::uniform_real_distribution<double> uni_dist(0.0,1.0);
  std::random_device rd;
  std::default_random_engine generator(rd());
  
  vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::indep);
  
  // Simulating claims for a single policyholder.
  double uR = pGB2(y1,mu1,sigma1,alpha11,alpha21);
  boost::math::negative_binomial_distribution<double> nb_dist(size, size/(size+meannb));
  double tempval1 = pdf(nb_dist,0); 
  double tempval2 = pdf(nb_dist,1);
  double pp0 = pzero + (1-pzero-pone)*tempval1;
  double pp1 = pzero + (1-pzero-pone)*tempval1 + pone + (1-pzero-pone)*tempval2;
  
  MatrixXd y3(BB,NN);
  MatrixXd buff(1,NN);
  VectorXd y3dep_mean;
  VectorXd pureprem_ind;
  pureprem_ind.resize(BB*2);
  
  // simulate a large number of observations.
  for (int b=0; b<BB; b++) {
  
  // 2. Independent case.
  for (int k=0; k<NN; k++) {
  y3(b,k) = 0.0;
  }
  double uYind = uni_dist(generator);
  Matrix<double, 1, 2> uRuYind; uRuYind << uR, uYind;
  VectorXd tempvec = cop_model_DN.hinv1(uRuYind);
  double uYdep = tempvec(0);
  double y2ind, y2dep, y2cind;
  
  // sample from a 01-negative binomial distribution.
  if (uYind < pp0) y2ind = 0; else
  if (uYind < pp1) y2ind = 1; else
  y2ind = quantile( nb_dist, pdf(nb_dist,0)+pdf(nb_dist,1) + (uYind-pp1)*(1-pdf(nb_dist,0)-pdf(nb_dist,1))/(1-pp1) );
  // sample from a 01-negative binomial distribution.
  if (uYdep < pp0) y2dep = 0; else
  if (uYdep < pp1) y2dep = 1; else
  y2dep = quantile( nb_dist, pdf(nb_dist,0)+pdf(nb_dist,1) + (uYdep-pp1)*(1-pdf(nb_dist,0)-pdf(nb_dist,1))/(1-pp1) );
  
  y2cind = y2dep;
  
  // Assume frequency doesn't exceed the buffer size.
  if (y2ind>NN) y2ind = NN;
  if (y2dep>NN) y2dep = NN;
  if (y2cind>NN) y2cind = NN;
  
  if (y2ind>0) {
    
    // Figure out the maximum of the uniform distribution that we can accomodate.
    paramsind = {0, mu1, sigma1, alpha11, alpha21, copDN(1), copDN(2), y1, y2ind, pzero, pone, meannb, size, mu3ind, sigma3ind, alpha13ind, alpha23ind, rhoRY, df1, thetaNY, df2}; 
    F.function = &obj1ind;
    F.params = &paramsind;
    double MaxU = F.function(coverage*1000000,&paramsind);
    for (int j=0; j<y2ind; j++) { buff(0,j) = uni_dist(generator); if(buff(0,j)>MaxU) {buff(0,j)=MaxU;}}
    for (int j=0; j<y2ind; j++) {
      iter = 0; max_iter = 100;
      r = 0.1; 
      x_lo = 0.0;
      x_hi = LARGENUMBER;
      paramsind = {buff(0,j), mu1, sigma1, alpha11, alpha21, copDN(1), copDN(2), y1, y2ind, pzero, pone, meannb, size, mu3ind, sigma3ind, alpha13ind, alpha23ind, rhoRY, df1, thetaNY, df2};
      F.function = &obj1ind;
      F.params = &paramsind;
      T = gsl_root_fsolver_brent;
      s = gsl_root_fsolver_alloc (T);
      gsl_root_fsolver_set (s, &F, x_lo, x_hi);
      // printf ("using %s method\n",gsl_root_fsolver_name (s));
      // printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "root", "err(est)");
      do {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.01);
        // if (status == GSL_SUCCESS) printf ("Converged:\n");
        // printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, r, x_hi - x_lo);
        // if (status == GSL_SUCCESS) cout << "Converged: " << r << endl;
        if (status == GSL_SUCCESS) y3(b,j) = r;
      }
      while (status == GSL_CONTINUE && iter < max_iter);
      gsl_root_fsolver_free (s);
    }
  }
  
  // Sum g(S) instead of S. In other words, apply the deductible.
  pureprem_ind(b) = 0;
  for (int j=0; j<y2ind; j++) {
    if (y3(b,j)>coverage*1000000)
      pureprem_ind(b) = pureprem_ind(b) + coverage*1000000; // (Error fixed) This was the error. We should multiply 1000,000 to the coverage.
    else {
      pureprem_ind(b) = pureprem_ind(b) + y3(b,j);
    }
  }
  
  // Frequencies.
  pureprem_ind(b+BB) = 0;
  for (int j=0; j<y2ind; j++) {
    pureprem_ind(b+BB) = pureprem_ind(b+BB) + 1;
  }
  
  //pureprem_ind(b) = 0;
  //for (int j=0; j<y2ind; j++) {
    //  if (y3(b,j)>coverage*1000000)
      //    pureprem_ind(b) = pureprem_ind(b) + coverage*1000000; // (Error fixed) This was the error. We should multiply 1000,000 to the coverage.
      //  else
        //    pureprem_ind(b) = pureprem_ind(b) + y3(b,j);
        //}
  
  }


  return(pureprem_ind);

}




int main() {
  
  MatrixXd temp, cova, cova1;
  MatrixXi tempInt;
  VectorXd y1, gam1, gam2, gam3, copDN, copDS, copDScind, coverage;
  VectorXi policynum;
  temp = readMatrix("dat/dat_y1.txt"); y1 = temp.col(0);
  temp = readMatrix("dat/dat_coverage.txt"); coverage = temp.col(0);
  temp = readMatrix("dat/dat_gam1.txt"); gam1 = temp.col(0);
  temp = readMatrix("dat/dat_gam2.txt"); gam2 = temp.col(0);
  temp = readMatrix("dat/dat_gam3.txt"); gam3 = temp.col(0);
  temp = readMatrix("dat/dat_copDN.txt"); copDN = temp.col(0);
  temp = readMatrix("dat/dat_copDS.txt"); copDS = temp.col(0);
  temp = readMatrix("dat/dat_copDScind.txt"); copDScind = temp.col(0);
  cova = readMatrix("dat/dat_cova.txt");
  cova1 = readMatrix("dat/dat_cova1.txt");
  tempInt = readMatrixInt("dat/dat_policynum.txt"); policynum = tempInt.col(0);
  
  VectorXd SimInd; SimInd.resize(BB);   
  
  VectorXd covatemp;  covatemp.resize(7);   
  VectorXd cova1temp; cova1temp.resize(6); 
  
  
  for (int rep=0; rep<250; rep++) {
    
    for (int iter=0; iter<100; iter++) {
      
      cout << "rep:" << rep << " iter: " << iter << endl;
      
      covatemp(0)=1;  covatemp(1)=1;  covatemp(2)=0;  covatemp(3)=0;  covatemp(4)=0;  covatemp(5)=0; covatemp(6)=(cova.row(442))(6);
      cova1temp(0)=1; cova1temp(1)=1; cova1temp(2)=0; cova1temp(3)=0; cova1temp(4)=0; cova1temp(5)=0;
      
      SimInd = SimPolicyInd(5000/coverage(442), covatemp, cova1temp, gam1, gam2, gam3, copDN, copDS, copDScind, coverage(442));
      std::string filenamei_ind = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      std::string filenamei_ind_freq = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      filenamei_ind = "simoutput_ind/" + filenamei_ind;
      filenamei_ind_freq = "simoutput_ind_freq/" + filenamei_ind_freq;
      std::ofstream filei_ind(filenamei_ind);
      if (filei_ind.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind << SimInd(b);
          filei_ind << endl;
        }
      }
      std::ofstream filei_ind_freq(filenamei_ind_freq);
      if (filei_ind_freq.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind_freq << SimInd(b+BB);
          filei_ind_freq << endl;
        }
      }
    }
    
    
    for (int iter=100; iter<200; iter++) {
      
      cout << "rep:" << rep << " iter: " << iter << endl;
      
      covatemp(0)=1;  covatemp(1)=0;  covatemp(2)=1;  covatemp(3)=0;  covatemp(4)=0;  covatemp(5)=0; covatemp(6)=(cova.row(442))(6);
      cova1temp(0)=1; cova1temp(1)=0; cova1temp(2)=1; cova1temp(3)=0; cova1temp(4)=0; cova1temp(5)=0;
      
      SimInd = SimPolicyInd(5000/coverage(442), covatemp, cova1temp, gam1, gam2, gam3, copDN, copDS, copDScind, coverage(442));
      std::string filenamei_ind = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      std::string filenamei_ind_freq = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      filenamei_ind = "simoutput_ind/" + filenamei_ind;
      filenamei_ind_freq = "simoutput_ind_freq/" + filenamei_ind_freq;
      std::ofstream filei_ind(filenamei_ind);
      if (filei_ind.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind << SimInd(b);
          filei_ind << endl;
        }
      }
      std::ofstream filei_ind_freq(filenamei_ind_freq);
      if (filei_ind_freq.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind_freq << SimInd(b+BB);
          filei_ind_freq << endl;
        }
      }
    }
    
    for (int iter=200; iter<300; iter++) {
      
      cout << "rep:" << rep << " iter: " << iter << endl;
      
      covatemp(0)=1;  covatemp(1)=0;  covatemp(2)=0;  covatemp(3)=1;  covatemp(4)=0;  covatemp(5)=0; covatemp(6)=(cova.row(442))(6);
      cova1temp(0)=1; cova1temp(1)=0; cova1temp(2)=0; cova1temp(3)=1; cova1temp(4)=0; cova1temp(5)=0;
      
      SimInd = SimPolicyInd(5000/coverage(442), covatemp, cova1temp, gam1, gam2, gam3, copDN, copDS, copDScind, coverage(442));
      std::string filenamei_ind = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      std::string filenamei_ind_freq = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      filenamei_ind = "simoutput_ind/" + filenamei_ind;
      filenamei_ind_freq = "simoutput_ind_freq/" + filenamei_ind_freq;
      std::ofstream filei_ind(filenamei_ind);
      if (filei_ind.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind << SimInd(b);
          filei_ind << endl;
        }
      }
      std::ofstream filei_ind_freq(filenamei_ind_freq);
      if (filei_ind_freq.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind_freq << SimInd(b+BB);
          filei_ind_freq << endl;
        }
      }
    }
    
    
    for (int iter=300; iter<400; iter++) {
      
      cout << "rep:" << rep << " iter: " << iter << endl;
      
      covatemp(0)=1;  covatemp(1)=0;  covatemp(2)=0;  covatemp(3)=0;  covatemp(4)=1;  covatemp(5)=0; covatemp(6)=(cova.row(442))(6);
      cova1temp(0)=1; cova1temp(1)=0; cova1temp(2)=0; cova1temp(3)=0; cova1temp(4)=1; cova1temp(5)=0;
      
      SimInd = SimPolicyInd(5000/coverage(442), covatemp, cova1temp, gam1, gam2, gam3, copDN, copDS, copDScind, coverage(442));
      std::string filenamei_ind = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      std::string filenamei_ind_freq = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      filenamei_ind = "simoutput_ind/" + filenamei_ind;
      filenamei_ind_freq = "simoutput_ind_freq/" + filenamei_ind_freq;
      std::ofstream filei_ind(filenamei_ind);
      if (filei_ind.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind << SimInd(b);
          filei_ind << endl;
        }
      }
      std::ofstream filei_ind_freq(filenamei_ind_freq);
      if (filei_ind_freq.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind_freq << SimInd(b+BB);
          filei_ind_freq << endl;
        }
      }
    }
    
    
    for (int iter=400; iter<500; iter++) {
      
      cout << "rep:" << rep << " iter: " << iter << endl;
      
      covatemp(0)=1;  covatemp(1)=0;  covatemp(2)=0;  covatemp(3)=0;  covatemp(4)=0;  covatemp(5)=1; covatemp(6)=(cova.row(442))(6);
      cova1temp(0)=1; cova1temp(1)=0; cova1temp(2)=0; cova1temp(3)=0; cova1temp(4)=0; cova1temp(5)=1;
      
      SimInd = SimPolicyInd(5000/coverage(442), covatemp, cova1temp, gam1, gam2, gam3, copDN, copDS, copDScind, coverage(442));
      std::string filenamei_ind = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      std::string filenamei_ind_freq = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      filenamei_ind = "simoutput_ind/" + filenamei_ind;
      filenamei_ind_freq = "simoutput_ind_freq/" + filenamei_ind_freq;
      std::ofstream filei_ind(filenamei_ind);
      if (filei_ind.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind << SimInd(b);
          filei_ind << endl;
        }
      }
      std::ofstream filei_ind_freq(filenamei_ind_freq);
      if (filei_ind_freq.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind_freq << SimInd(b+BB);
          filei_ind_freq << endl;
        }
      }
    }
    
    
    for (int iter=500; iter<600; iter++) {
      
      cout << "rep:" << rep << " iter: " << iter << endl;
      
      covatemp(0)=1;  covatemp(1)=0;  covatemp(2)=0;  covatemp(3)=0;  covatemp(4)=0;  covatemp(5)=0; covatemp(6)=(cova.row(442))(6);
      cova1temp(0)=1; cova1temp(1)=0; cova1temp(2)=0; cova1temp(3)=0; cova1temp(4)=0; cova1temp(5)=0;
      
      SimInd = SimPolicyInd(5000/coverage(442), covatemp, cova1temp, gam1, gam2, gam3, copDN, copDS, copDScind, coverage(442));
      std::string filenamei_ind = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      std::string filenamei_ind_freq = std::to_string(rep+1) + "_" + std::to_string(iter+1) + ".txt";
      filenamei_ind = "simoutput_ind/" + filenamei_ind;
      filenamei_ind_freq = "simoutput_ind_freq/" + filenamei_ind_freq;
      std::ofstream filei_ind(filenamei_ind);
      if (filei_ind.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind << SimInd(b);
          filei_ind << endl;
        }
      }
      std::ofstream filei_ind_freq(filenamei_ind_freq);
      if (filei_ind_freq.is_open()) {
        for (int b=0; b<BB; b++) {
          filei_ind_freq << SimInd(b+BB);
          filei_ind_freq << endl;
        }
      }
    }
    
  }
  
  return(0);
}

