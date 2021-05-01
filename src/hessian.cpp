#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List hessian(int n, int p, int kint, int nx, NumericMatrix x, IntegerVector y, IntegerVector delta,
                   NumericVector pr, NumericVector fa, NumericVector fb, NumericVector fpa, NumericVector fpb,
                   NumericVector fppa, NumericVector fppb, int l, int lia){
  
  NumericVector v(l);
  IntegerVector ja(l);
  IntegerVector ia(lia);

  int iv = 0, nkk;
  IntegerVector kk(p);
  for(unsigned int m=0; m < p; m++){
    //  Compute column numbers for nonzero elements: kk
    if(kint > 1){
      if(m == 0){
        nkk = 2;
        kk(0) = 0;
        kk(1) = 1;
      }else if(m > 0 and m < kint-1){
        nkk = 3;
        kk(0) = m - 1;
        kk(1) = m;
        kk(2) = m + 1;
      }else if(m == kint-1){
        nkk = 2;
        kk(0) = m - 1;
        kk(1) = m;
      }else{
        nkk = kint;
        for(unsigned int mm = 0; mm < kint; mm++){
          kk(mm) = mm;
        }
      }
      for(unsigned int mm = kint; mm < p; mm++){
        nkk += 1;
        kk(nkk-1) = mm;
      }
    }else{
      nkk = p;
    }
    
    int k;
    for(unsigned int ik = 1; ik <= nkk; ik++){
      if(kint == 1){
        k = ik-1;
      }else{
        k = kk(ik-1);
      }
      double vmk = 0;
      double a;
      for(unsigned int j = 0; j < n; j++){
        int z = y(j);
        int d = delta(j);
        double pa = fpa(j);
        double pb = fpb(j);
        double ppa = fppa(j);
        double ppb = fppb(j);
        double w = 1 / (pr(j) * pr(j));
        if(m < kint and k < kint){
          if(d == 1 or d == 12 or d == 13){
            a = -w * (pb * (z == m+1) - pa * (z-1 == m+1)) * 
              (pb * (z == k+1) - pa * (z-1 == k+1)) + 
              (ppb * (z == m+1) * (m == k) - ppa * (z-1 == m+1) * (m == k)) / pr(j);
          }else if(d == 2){
            a = (ppb * (z == m+1) * (m == k) - pb * pb * (z == m+1) * (m == k) / fb(j)) / fb(j);
          }else{
            a = -ppa * (z == m+1) * (m == k) - pa * pa * (z == m+1) * (m == k) / (1 - fa(j)) / (1 - fa(j));
          }
        }else if(m >= kint and k >= kint){
          if(d == 1 or d == 12 or d == 13){
            a = x(j, m-kint) * x(j, k-kint) / pr(j) * (-1 / pr(j) * (pb - pa) * (pb - pa) + ppb - ppa);
          }else if(d == 2){
            a = x(j, m-kint) * x(j, k-kint) / fb(j) * (- pb * pb / fb(j) + ppb);
          }else{
            a = x(j, m-kint) * x(j, k-kint) / (1 - fa(j)) * (pa * pa / (1 - fa(j)) + ppa);
          }
        }else if(m >= kint and k < kint){
          if(d == 1 or d == 12 or d == 13){
            a = x(j, m - kint) / pr(j) * (-1 / pr(j) * (pb - pa) * (pb * (z == k+1) - 
              pa * (z-1 == k+1)) + ppb * (z == k+1) - ppa * (z-1 == k+1));
          }else if(d == 2){
            a = -x(j, m - kint) / fb(j) * (ppb * (z == k+1) - pb * pb / fb(j) * (z == k+1));
          }else{
            a = x(j, m - kint) / (1 - fa(j)) * (ppa * (z == k+1) + pa * pa / (1 - fa(j)) * (z == k+1));
          }
        }else if(m < kint and k >= kint){
          if(d == 1 or d == 12 or d == 13){
            a = x(j, k - kint) / pr(j) * (-1 / pr(j) * (pb - pa) * (pb * (z == m+1) - 
              pa * (z-1 == m+1)) + ppb * (z == m+1) - ppa * (z-1 == m+1));
          }else if(d == 2){
            a = -x(j, k - kint) / fb(j) * (ppb * (z == m+1) - pb * pb / fb(j) * (z == m+1));
          }else{
            a = x(j, k - kint) / (1 - fa(j)) * (ppa * (z == m+1)+ pa * pa / (1 - fa(j)) * (z == m+1));
          }
        }
        vmk += a;
      }
      iv += 1;
      v(iv-1) = -vmk;
      if(kint > 1){
        ja(iv-1) = k + 1;
        if(ik == 1){
          ia(m) = iv;
        }
      }
    }
  }
  
  if(kint >1){
    ia(p) = iv + 1;
  }

  // output
  return List::create(Named("v")=v,
                      Named("ja")=ja,
                      Named("ia")=ia);
}
