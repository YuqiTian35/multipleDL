functions {
    real func_link(real x, int link_num) {
      if (link_num == 1) return(inv_logit(x));
      else if (link_num == 2) return(Phi(x));
      else if (link_num == 3) return(gumbel_cdf(x, 0, 1));
      else if (link_num == 4) return(inv_cloglog(x));
      else return x;
    }
    // # 1 = logistic; 2 = probit; 3 = loglog; 4 = cloglog

  vector loglik(vector alpha, vector beta, int[] delta, row_vector[] X, int[] j, int link_num) {
      int N = size(X);
      vector[N] out;
      int J = max(j);
      for (n in 1:N){
        real eta = X[n] * beta;
        int j_n = j[n];
        int delta_n = delta[n];
        if(delta_n == 2) out[n] = log(func_link(alpha[j_n] - eta, link_num));
        else if (delta_n == 3) out[n] = log1m(func_link(alpha[j_n - 1] - eta, link_num));
        else if (delta_n == 12 || j_n == 1) out[n] = log(func_link(alpha[1] - eta, link_num));
        else if (delta_n == 13 || j_n == J) out[n] = log1m(func_link(alpha[J - 1] - eta, link_num));
        else out[n] = log(func_link(alpha[j_n] - eta, link_num) - func_link(alpha[j_n - 1] - eta, link_num));
      }
      return out;
    }
  }



data {
   int<lower = 1> N;   // number of observations
   int<lower = 1> p;   // number of predictors
   row_vector[p] X[N];     // columnwise CENTERED predictors
   int<lower = 2> J;   // number of outcome categories
   int<lower = 1, upper = J> j[N]; // outcome on 1 ... J'
   int delta[N]; //detection indicator
   int link_num; // link function

   // prior standard deviations
   real<lower=0> sds;
 }


parameters {
  vector[p] beta; // coefficients
  simplex[J] pi;  // category probabilities for a person w/ average predictors
}

transformed parameters {
  vector[J - 1] alpha;                               // intercepts
  vector[N] log_lik;                                 // log-likelihood pieces

  if (link_num == 1) for (i in 2:J) alpha[i - 1] = logit(1 - sum(pi[i:J])); // logit link
  if (link_num == 2) for (i in 2:J) alpha[i - 1] = inv_Phi(1 - sum(pi[i:J])); // probit link
  if (link_num == 3) for (i in 2:J) alpha[i - 1] = -log(-log(1 - sum(pi[i:J]))); // loglog link
  if (link_num == 4) for (i in 2:J) alpha[i - 1] = log(-log1m(1 - sum(pi[i:J]))); // cloglog link

  log_lik = loglik(alpha, beta, delta, X, j, link_num);
}

model {
  target += log_lik;
  target += normal_lpdf(beta | 0, rep_vector(sds, p));
  target += dirichlet_lpdf(pi | rep_vector(1, J));
}
