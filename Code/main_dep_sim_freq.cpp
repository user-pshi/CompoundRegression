
//******************************************************************************************
// main_dep_sim_freq.cpp: Predict frequencies using the simulation approach for the dependence model.
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
#define BB 100000 // number of replicates.
#define NN 500 // maximum number of claims for a policyholder in a year.
#define LARGENUMBER 1.0e+20

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
  VectorXd y1, gam1, gam2, gam3, copDN, copDS, copDScind, coverage;
  VectorXi policynum;
  // Matrix<double, BB, NN> y3dep, y3ind, y3cind;
  // double y3[1][NN];
  // double y3ind[BB][NN];
  // double y3cind[BB][NN];
  // Matrix<double, 1, NN> buff;
  MatrixXd y3(BB,NN);
  MatrixXd buff(1,NN);
  VectorXd y3dep_mean;
  // double pureprem_dep[BB], pureprem_ind[BB], pureprem_cind[BB];
  std::vector<double> claims_dep (BB, 0);
  // Matrix<double, Dynamic, NN> y3dep_all, y3ind_all, y3cind_all;
  
  // Variables related to the solver.
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r, x_lo, x_hi;
  gsl_function F;
  struct obj1dep_params paramsdep;  
  int status, iter, max_iter;
  int BB2 = BB/2;
  
  temp = readMatrix("dat/dat_y1.txt"); y1 = temp.col(0);
  temp = readMatrix("dat/dat_coverage.txt"); coverage = temp.col(0);
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
    
    for (int b=0; b<BB; b++) {

            
      // 1. Dependent case.
      for (int k=0; k<NN; k++) {
        y3(b,k) = 0.0;
      }
      uYind = uni_dist(generator);
      Matrix<double, 1, 2> uRuYind; uRuYind << uR, uYind;
      tempvec = cop_model_DN.hfunc1(uRuYind);
      uYdep = tempvec(0);
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
      // if (y2ind>NN) y2ind = NN;
      // if (y2dep>NN) y2dep = NN;
      // if (y2cind>NN) y2cind = NN;

      
      claims_dep[b] = 0;
      for (int j=0; j<y2dep; j++) {
          claims_dep[b] = claims_dep[b] + 1;
      }
      // file_sim_dep << y3[0][0];
      // for (int j=1; j<NN; j++)
      //   file_sim_dep << "," << y3[0][j];
      // file_sim_dep << endl;

      // cout << y2dep << " ";
      // cout << uYdep << endl;
    }
    
    // std::sort (pureprem_dep.begin(), pureprem_dep.begin()+BB);
    // y3dep_median(i) = pureprem_dep[BB2];
    
    for (int b=0; b<BB; b++) {
      y3dep_mean(i) = y3dep_mean(i) + claims_dep[b];
    }    
    y3dep_mean(i) = y3dep_mean(i) / BB;
    
    cout << "Simulated mean: " << y3dep_mean(i) << endl;
    
    std::string filenamei = std::to_string(policynum(i)) + ".txt";
    filenamei = "out_dep_freq/" + filenamei;
    std::ofstream filei(filenamei);
    if (filei.is_open()) {
      for (int b=0; b<BB; b++) {
        filei << y3(b,0);
        for (int j=1; j<NN; j++) 
          filei << "," << y3(b,j);
        filei << endl;
      }
    }
    
    // cout << endl;
    // cout << pp0 << endl;
    // cout << pp1 << endl;
    
    // cout << uYind << endl;
    // cout << uR << endl;
  }
  
  
  
  // cout << mu3cind.head(5) << endl << endl;
  // cout << sigma3cind << endl;
  // cout << alpha13cind << endl;
  // cout << alpha23cind << endl;
  
  // cout << mu3ind.head(5) << endl << endl;
  // cout << sigma3ind << endl;
  // cout << alpha13ind << endl;
  // cout << alpha23ind << endl;
  
  // cout << mu3dep.head(5) << endl << endl;
  // cout << sigma3dep << endl;
  // cout << alpha13dep << endl;
  // cout << alpha23dep << endl;
  
  // cout << gam3dep << endl << endl;
  // cout << gam3cind << endl << endl;
    
  // cout << pzero.head(5) << endl << endl;
  // cout << pone.head(5) << endl << endl;
  // cout << meannb.head(5) << endl << endl;
  
  // cout << beta << endl;
  // cout << gamma0 << endl;
  // cout << gamma1 << endl;
  // cout << size << endl;
  
  // cout << mu1.head(5) << endl;
  // cout << sigma1 << endl;
  // cout << alpha11 << endl;
  // cout << alpha21 << endl;
  
  // cout << rhoRY << endl;
  // cout << thetaNY << endl;
  // cout << rhoRYcind << endl;
  // cout << fam1 << endl;
  // cout << fam2 << endl;
  // cout << fam1cind << endl;
  // cout << df1 << endl;
  // cout << df2 << endl;
  // cout << df1cind << endl;
  
  // std::cout << cova_zero.block(0,0,3,cova_zero.cols()) << endl << endl;
  // std::cout << cova_one.block(0,0,3,cova_one.cols()) << endl << endl;
  // std::cout << cova.block(0,0,3,cova.cols()) << endl;
  
  // Vector2d temptemp; temptemp << 10, 0;
  // vinecopulib::Bicop cop_model_NY(vinecopulib::BicopFamily::gumbel, 90, temptemp.head(1));
  // Matrix<double,5000,2> sim_temp = cop_model_NY.simulate(5000);
  // std::ofstream file("sim_temp.txt");
  // if (file.is_open()) {
  //   for (int j=0; j<sim_temp.rows(); j++) {
  //     file << sim_temp(j,0) << "," << sim_temp(j,1) << endl;
  //   }
  // }
  
  // std::ofstream file("y3dep_all.txt");
  // if (file.is_open()) {
  //   for (int j=0; j<y3dep_all.rows(); j++) {
  //     file << y3dep_all(j,0);
  //     for (int k=1; k<NN; k++) {
  //       file << "," << y3dep_all(j,k);
  //     }
  //     file << endl;
  //   }
  // }

  std::ofstream file_dep_mean("y3dep_mean_freq.txt");
  if (file_dep_mean.is_open()) {
    for (int i=0; i<nobs; i++) {
      file_dep_mean << y3dep_mean(i);
      file_dep_mean << endl;
    }
  }
  
  // std::ofstream file_dep_median("y3dep_median.txt");
  // if (file_dep_median.is_open()) {
  //   for (int i=0; i<nobs; i++) {
  //     file_dep_median << y3dep_median(i);
  //     file_dep_median << endl;
  //   }
  // }

  
  return(0);
}



// class GB2 {
// public:
//   GB2(const double& mu, const double& sigma, const double& alpha1, const double& alpha2) {
//     this->mu = mu;
//     this->sigma = sigma;
//     this->alpha1 = alpha1;
//     this->alpha2 = alpha2;
//     this->ndf = 2*alpha1;
//     this->ddf = 2*alpha2;
//     boost::math::fisher_f temp(ndf,ddf);
//     this->f_func = temp;
//   }
//   pGB2(double& y) {
//     double r = (std::log(y)-this->mu)/this->sigma;
//     double z = (this->alpha2/this->alpha1)*std::exp(r);
//     double ret_val = cdf(this->f_func, z);
//     return(ret_val);
//   }
//   double mu;
//   double sigma;
//   double alpha1;
//   double alpha2;
//   double ndf;
//   double ddf;
//   boost::math::fisher_f f_func;
// };
