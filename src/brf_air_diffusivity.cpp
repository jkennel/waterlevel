#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double weeks_1979(double lag, 
                  double D,
                  double L, 
                  double precision = 1e-8) {
  
  double d_term = M_PI * M_PI * D * lag / (4 * L * L);
  double term_val = 0.0;
  double ret = 0.0;
  double exp_val = 0.0;
  int m = 1;
  bool more_precise = TRUE;
  
  if (d_term < 0.001) {
    return(1.0);
  }
  
  
  while(more_precise) {
    exp_val  = -m * m * d_term;
    term_val = std::pow(-1, (m - 1) / 2) / m * std::exp(exp_val);
    ret += term_val;
    m   += 2;
    more_precise = std::abs(term_val) > precision;
  }
  
  ret = (4/M_PI * ret);
  if (ret > 1.0) {
    ret = 1.0;
  }
  
  return(ret);
  
}


//==============================================================================
struct weeks_worker : public Worker
{
  // source vector
  
  const arma::vec& input;
  double D;
  double L;
  double precision;
  
  arma::vec& output;
  
  weeks_worker(const arma::vec& input,
               arma::vec& output,
               double D,
               double L,
               double precision)
    : input(input), output(output), D(D), L(L) {}
  
  // calculate the exponential integral
  void operator()(std::size_t begin_row, std::size_t end_row) {
    
    for (std::size_t i = begin_row; i < end_row; i++) {
      output(i) = weeks_1979(input(i), D, L, precision);
    }
    
  }
  
  
};

//==============================================================================
//' @title
//' vadose_response
//'
//' @description
//' weeks_1979 1-D air diffusivity
//'
//' @param time \code{numeric} value of elapsed time
//' @param D \code{numeric} unsaturated zone air diffusivity
//' @param L \code{numeric} unsaturated zone thickness
//'
//' @return brf from weeks 1979 model
//'
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec vadose_response(const arma::vec time, 
                          double D, 
                          double L,
                          double precision = 1e-8) {
  
  int n = time.n_elem;
  
  arma::vec output(n);
  
  weeks_worker weeks(time, output, D, L, precision);
  
  RcppParallel::parallelFor(0, n, weeks);
  
  return(output);
}


/***R

time <- seq(0, 43200, 1)
D <- 0.1
L <- 68
precision <- 1e-8
system.time(
  w <- vadose_response(time, D, L, precision)
)
plot(w, type='l', ylim = c(0, 1))
abline(h = 0.25)
tail(w)

system.time({
  
  for(i in 1:length(time)) {
    weeks_1979(time[i], D, L, precision)
  }
})

*/
