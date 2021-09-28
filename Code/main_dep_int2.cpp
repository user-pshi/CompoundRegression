
//******************************************************************************************
// main_dep_int2.cpp: Predict aggregate claims using the integral approach for the dependence model.
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
#include <gsl/gsl_integration.h>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <string> 
#include <math.h>

using namespace std;
using namespace Eigen;

#define MAXN 200
#define MAXBUFSIZE  ((int) 1e6)
#define LARGE_NUMBER (1e30)

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

double dGB2(double y, double mu, double sigma, double alpha1, double alpha2) {
  double fstterm = exp(alpha1 * (log(y) - mu) / sigma);
  double sndterm = y * abs(sigma)* exp(lgamma(alpha1)) * exp(lgamma(alpha2))/exp(lgamma(alpha1 + alpha2));
  double thdterm = pow(1 + exp((log(y) - mu)/sigma ), alpha1 + alpha2); 
  double ans = fstterm/(sndterm * thdterm);
  return(ans);
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

struct funcY_params {
  double coverage;
  double deduct;
  double mu1;
  double sigma1;
  double alpha11;
  double alpha21;
  double pzero;
  double pone;
  double meannb;
  double size;
  double mu3dep;
  double sigma3dep;
  double alpha13dep;
  double alpha23dep;
  double copDNpar1;
  double copDNpar2;
  double copRYpar1;
  double copRYpar2;
  double copNYpar1;
  double copNYpar2;
};

// Integral approach to mean calculation, for the dependence model.
// The following does not left-censor the losses.
double funcY(double y3, void *params) {
  
  struct funcY_params *p = (struct funcY_params *) params;
  
  double coverage = p->coverage;
  double deduct = p->deduct;
  double mu1 = p->mu1;
  double sigma1 = p->sigma1;
  double alpha11 = p->alpha11;
  double alpha21 = p->alpha21;
  double pzero = p->pzero;
  double pone = p->pone;
  double meannb = p->meannb;
  double size = p->size;
  double mu3dep = p->mu3dep;
  double sigma3dep = p->sigma3dep;
  double alpha13dep = p->alpha13dep;
  double alpha23dep = p->alpha23dep;
  Vector2d copDNparam; copDNparam << p->copDNpar1, p->copDNpar2;
  Vector2d copRYparam; copRYparam << p->copRYpar1, p->copRYpar2;
  Vector2d copNYparam; copNYparam << p->copNYpar1, p->copNYpar2;
  vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::student, 0, copDNparam);
  vinecopulib::Bicop cop_model_RY(vinecopulib::BicopFamily::frank, 0, copRYparam.head(1));
  vinecopulib::Bicop cop_model_NY(vinecopulib::BicopFamily::gumbel, 90, -copNYparam.head(1));
  double y1 = deduct / coverage; // use the given deductible level to calculate dedratio.
  double uR = pGB2(y1,mu1,sigma1,alpha11,alpha21); // Given deductible. 
  Matrix<double, 1, 2> utemp, utemp1, utemp2; 
  VectorXd tempvec, tempvec1, tempvec2;
  double u1, u1a, u2, densNY_R, y2, y3cens;
  double RetVal;

  if(y3>coverage*1000000) { y3cens = coverage*1000000; } else { y3cens = y3; }
  
  RetVal = 0;
  for (int n=1; n<MAXN; n++) {
    y2 = (double)n;
    // frequency|deductible
    utemp(0,0) = pGB2(y1,mu1,sigma1,alpha11,alpha21); utemp(0,1) = Fcount(y2,pzero,pone,meannb,size);  tempvec = cop_model_DN.hfunc1(utemp); u1 = tempvec(0);
    utemp(0,1) = Fcount(y2-1,pzero,pone,meannb,size); tempvec = cop_model_DN.hfunc1(utemp); u1a = tempvec(0);
    // severity|deductible
    utemp(0,1) = pGB2(y3,mu3dep,sigma3dep,alpha13dep,alpha23dep);
    tempvec = cop_model_RY.hfunc1(utemp);
    u2 = tempvec(0);
    utemp1(0,0) = u1; utemp1(0,1) = u2; utemp2(0,0) = u1a; utemp2(0,1) = u2;
    tempvec1 = cop_model_NY.hfunc2(utemp1);
    tempvec2 = cop_model_NY.hfunc2(utemp2);
    densNY_R = dGB2(y3,mu3dep,sigma3dep,alpha13dep,alpha23dep) * (tempvec1(0)-tempvec2(0));
    RetVal = RetVal + y2 * densNY_R;
  }
  return(y3cens*RetVal);
}


int main() {
  MatrixXd temp, cova, cova1;
  MatrixXi tempInt;
  Matrix<double, Dynamic, 2> cova_zero; 
  Matrix<double, Dynamic, 2> cova_one;
  VectorXd mu1, beta, gamma0, gamma1, pzero, pone, meannb, tempvec0, tempvec1, gam3dep, gam3cind;
  VectorXd mu3dep, mu3ind, mu3cind, tempvec;
  int k1, k, nobs;
  double rhoRY, thetaNY, rhoRYcind, fam1, fam2, fam1cind, df1, df2, df1cind;
  double sigma1, alpha11, alpha21, size;
  double sigma3dep, alpha13dep, alpha23dep;
  double sigma3ind, alpha13ind, alpha23ind;
  double sigma3cind, alpha13cind, alpha23cind;
  double uR, uYind, uYdep, pp0, pp1, tempval1, tempval2, y2dep, y2ind, y2cind;
  VectorXd y1, gam1, gam2, gam3, copDN, copDS, copDScind, coverage, deduct;
  VectorXi policynum;
  VectorXd y3dep_mean;
  gsl_function F;
  struct funcY_params paramsdep;  
  double result, error;
  
  temp = readMatrix("dat/dat_y1.txt"); y1 = temp.col(0);
  temp = readMatrix("dat/dat_coverage.txt"); coverage = temp.col(0);
  temp = readMatrix("dat/dat_deduct.txt"); deduct = temp.col(0);
  temp = readMatrix("dat/dat_gam1.txt"); gam1 = temp.col(0);
  temp = readMatrix("dat/dat_gam2.txt"); gam2 = temp.col(0);
  temp = readMatrix("dat/dat_gam3.txt"); gam3 = temp.col(0);
  temp = readMatrix("dat/dat_gam3.txt"); gam3 = temp.col(0);
  temp = readMatrix("dat/dat_copDN.txt"); copDN = temp.col(0);
  temp = readMatrix("dat/dat_copDS.txt"); copDS = temp.col(0);
  temp = readMatrix("dat/dat_copDScind.txt"); copDScind = temp.col(0);
  cova = readMatrix("dat/dat_cova.txt");
  cova1 = readMatrix("dat/dat_cova1.txt");
  tempInt = readMatrixInt("dat/dat_policynum.txt"); policynum = tempInt.col(0);
  k1 = cova1.cols();
  k = cova.cols();
  nobs = cova.rows();
  cova_zero.resize(nobs,2); cova_zero.col(0) = cova.col(0); cova_zero.col(1) = cova.col(k-1);
  cova_one.resize(nobs,2); cova_one.col(0) = cova.col(0); cova_one.col(1) = cova.col(k-1);
  
  cout << "Initializing result matrices" << endl;
  
  // cout << "Maximum double value: " << DBL_MAX << endl;
  
  // Dynamically resize the result matrix, and initialize.
  // y3dep_all.resize(BB*nobs, NN);
  // y3ind_all.resize(BB*nobs, NN);
  // y3cind_all.resize(BB*nobs, NN);
  y3dep_mean.resize(nobs);
  //y3dep_median.resize(nobs);
  for (int i=0; i<nobs; i++) {
    y3dep_mean(i) = 0;
    //y3dep_median(i) = 0;
  }

  cout << "Finished initializing result matrices" << endl;
    
  // Dependence parameters
  rhoRY = copDS[2+k+3+1-1];
  thetaNY = copDS[2+k+3+2-1];
  rhoRYcind = copDScind[1+k+3+1-1];
  fam1 = copDS[0];
  fam2 = copDS[1];
  fam1cind = copDScind[0];
  df1 = 5;
  df2 = 5;
  df1cind = 5;
  
  // Deductible parameters
  mu1 = cova1 * gam1.head(k1);
  sigma1 = gam1(k1);
  alpha11 = gam1(k1+1);
  alpha21 = gam1(k1+2);

  // Frequency parameters
  beta = gam2.head(k);
  gamma0 = gam2.segment(k,2);
  gamma1 = gam2.segment(k+2,2);
  size = gam2(k+2+2);

  tempvec0 = cova_zero * gamma0; for (int i=0; i<tempvec0.size(); i++) tempvec0(i) = exp(tempvec0(i));
  tempvec1 = cova_one * gamma1; for (int i=0; i<tempvec1.size(); i++) tempvec1(i) = exp(tempvec1(i));
  pzero = tempvec0; for (int i=0; i<pzero.size(); i++) pzero(i) = pzero(i) / (1+tempvec0(i) + tempvec1(i));
  pone = tempvec1; for (int i=0; i<pone.size(); i++) pone(i) = pone(i) / (1+tempvec0(i) + tempvec1(i));
  meannb = cova * beta; for (int i=0; i<meannb.size(); i++) meannb(i) = exp(meannb(i));

  // New parameters for severity, after full-likelihood estimation.
  gam3dep = copDS.segment(2,k+3);
  gam3cind = copDScind.segment(2,k+3);
  
  mu3dep = cova * gam3dep.head(k);
  sigma3dep = gam3dep(k);
  alpha13dep = gam3dep(k+1);
  alpha23dep = gam3dep(k+2);
  
  mu3cind = cova * gam3cind.head(k);
  sigma3cind = gam3cind(k);
  alpha13cind = gam3cind(k+1);
  alpha23cind = gam3cind(k+2);
  
  mu3ind = cova * gam3.head(k);
  sigma3ind = gam3(k);
  alpha13ind = gam3(k+1);
  alpha23ind = gam3(k+2);
  
  std::default_random_engine generator;
  std::uniform_real_distribution<double> uni_dist(0.0,1.0);
  vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::student, 0, copDN.tail(2));
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);  
  F.function = &funcY;
  
  cout << "Starting loop." << endl;
  
  // For each policyholder
  for (int i=0; i<nobs; i++) {
    cout << "Simulating claims for policyholder: " << i << endl;
    uR = pGB2(y1(i),mu1(i),sigma1,alpha11,alpha21);
    boost::math::negative_binomial_distribution<double> nb_dist(size, size/(size+meannb(i)));
    tempval1 = pdf(nb_dist,0); 
    tempval2 = pdf(nb_dist,1);
    pp0 = pzero(i) + (1-pzero(i)-pone(i))*tempval1;
    pp1 = pzero(i) + (1-pzero(i)-pone(i))*tempval1 + pone(i) + (1-pzero(i)-pone(i))*tempval2;
    paramsdep = {coverage(i), deduct(i), mu1(i), sigma1, alpha11, alpha21, pzero(i), pone(i), meannb(i), size, mu3dep(i), sigma3dep, alpha13dep, alpha23dep, copDN(1), copDN(2), rhoRY, df1, thetaNY, df2};
    F.params = &paramsdep;
    gsl_integration_qags (&F, 0, coverage(i)*1000000, 0, 1e-7, 1000, w, &result, &error); 
    // gsl_integration_qags (&F, 0, LARGE_NUMBER, 0, 1e-7, 1000, w, &result, &error); 
    y3dep_mean(i) = result;
    
    
    // calculate the adjustment factor u*(1-F(u))...
    double AdjVal = 0;
    for (int n=1; n<MAXN; n++) {
      double y2 = (double)n;
      Vector2d copDNparam; copDNparam << copDN(1), copDN(2);
      Vector2d copRYparam; copRYparam << rhoRY, df1;
      Vector2d copNYparam; copNYparam << thetaNY, df2;
      vinecopulib::Bicop cop_model_DN(vinecopulib::BicopFamily::student, 0, copDNparam);
      vinecopulib::Bicop cop_model_RY(vinecopulib::BicopFamily::frank, 0, copRYparam.head(1));
      vinecopulib::Bicop cop_model_NY(vinecopulib::BicopFamily::gumbel, 90, -copNYparam.head(1));
      Matrix<double, 1, 2> utemp, utemp1, utemp2;
      VectorXd tempvec, tempvec1, tempvec2;
      double u1, u1a, u2;
      // frequency|deductible
      utemp(0,0) = pGB2(y1(i),mu1(i),sigma1,alpha11,alpha21); utemp(0,1) = Fcount(y2,pzero(i),pone(i),meannb(i),size);  tempvec = cop_model_DN.hfunc1(utemp); u1 = tempvec(0);
      utemp(0,1) = Fcount(y2-1,pzero(i),pone(i),meannb(i),size); tempvec = cop_model_DN.hfunc1(utemp); u1a = tempvec(0);
      // severity|deductible
      utemp(0,1) = pGB2(coverage(i)*1000000,mu3dep(i),sigma3dep,alpha13dep,alpha23dep);
      tempvec = cop_model_RY.hfunc1(utemp);
      u2 = tempvec(0);
      utemp1(0,0) = u1; utemp1(0,1) = u2; utemp2(0,0) = u1a; utemp2(0,1) = u2;
      tempvec1 = cop_model_NY.cdf(utemp1);
      tempvec2 = cop_model_NY.cdf(utemp2);
      double InnerIntegral = (u1-u1a) - (tempvec1(0)-tempvec2(0));
      InnerIntegral *= y2;
      AdjVal += InnerIntegral;
    }
    AdjVal *= (coverage(i)*1000000);
    
    y3dep_mean(i) += AdjVal;
    
    // gsl_integration_qags (&F, coverage(i)*1000000, coverage(i)*1000000*10000, 0, 1e-7, 1000, w, &result, &error); 
    cout << "Result: " << y3dep_mean(i) << ", Adjustment: " << AdjVal << endl;
  }
  
  std::ofstream file_dep_mean("mean_dep_int2.txt");
  if (file_dep_mean.is_open()) {
    for (int i=0; i<nobs; i++) {
      file_dep_mean << y3dep_mean(i);
      file_dep_mean << endl;
    }
  }
  
  return(0);
}


