// To be provided to the psd package
// Todo check frequencies

#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;


// For all but short series this should be faster
//' @rdname parabolic_weights_field
//' 
//' @param ntap the maximum number of tapers
//' 
//' @export
// [[Rcpp::export]]
arma::field<arma::vec> parabolic_weights_field(const int ntap) {
  
  //
  // return quadratic spectral weighting factors for a given number of tapers
  // Barbour and Parker (2014) Equation 7
  //
  arma::field<arma::vec> f(ntap, 1);
  arma::vec kseq = arma::pow(arma::regspace<arma::vec>(0, ntap - 1), 2);
  
  double t3;
  for (int i = 1; i < ntap+1; i++) {
    t3 = log(i * (i - 0.25) * (i + 1.0));
    f(i-1,0) = arma::exp(log(1.5) + arma::log(i * i - kseq(arma::span(0, i-1))) - t3);
  }
  
  return(f);
  
}


// Serial version --------------------------------------------------------------



//' @title Resample an fft using varying numbers of sine tapers
//' 
//' @description
//' Produce an un-normalized psd based on an fft and a vector of optimal sine tapers
//' 
//' @details
//' To produce a psd estimate with our adaptive spectrum estimation method, we need only make one 
//' fft calculation initially and then
//' apply the weighting factors given by \code{\link{parabolic_weights_field}}, which this
//' function does.
//' 
//' @param fftz complex; a matrix representing the dual-length \code{\link{fft}}; see also the \code{dbl} argument
//' @param tapers integer; a vector of tapers
//' @param verbose logical; should messages be given?
//' @param dbl logical; should the code assume \code{fftz} is dual-length or single-length?
//' @param tapcap integer; the maximum number of tapers which can be applied; note that the length is
//' automatically limited by the length of the series.
//' 
//' @seealso \code{\link{riedsid}}
//' 
//' @examples
//' fftz <- complex(real=1:8, imaginary = 1:8)
//' taps <- 1:4
//' try(resample_fft_rcpp2(fftz, taps))
//' 
//' @export
// [[Rcpp::export]]
List resample_mvfft( const arma::cx_mat& fftz, 
                     const arma::ivec& tapers,
                     bool verbose = true, 
                     const bool dbl = true, 
                     const int tapcap = 10000 ) {
  
  //
  // resample and reweight an fft estimates for a given number of tapers
  // using quadratic scaling in Barbour and Parker (2014)
  //
  //
  // needs:
  //  - fftz: complex vector -- the FFT of the original signal
  //  - tapers: integer vector -- the number of tapers at each frequency
  //
  // optional args:
  //  - verbose: logical -- should warnings be given?
  //  - dbl: logical -- should the progam assume 'fftz' is for a double-length (padded) series? 
  //                    Otherwise a single-length series is assumed.
  //  - tapcap: integer -- the maximum number of tapers at any frequency
  //
  
  int sc, nf, nt, ne, ne2, nhalf, m2, mleft1, mleft2, Kc, ki, j1, j2;
  int nc = fftz.n_cols;
  double wi;
  
  
  arma::cx_double fdc1, fdc2;
  
  
  if (dbl){
    // double-length fft estimates assumed by default
    sc = 2;
  } else {
    // but could be single-length
    sc = 1;
  }
  
  // even, double, and half lengths
  nf = fftz.n_rows / sc; 
  nt = tapers.n_elem;
  ne = nf - (nf % 2);
  
  if (verbose){
    Function msg("message");
    msg(std::string("\tfft resampling"));
  }
  
  if (ne < nf){
    warning("fft was not done on an even length series");
  }
  
  ne2 = 2 * ne;
  nhalf = ne / 2;
  
  arma::ivec taper_vec(nhalf);
  // arma::vec psd(nhalf);
  // psd.zeros();
  
  
  if (nhalf < 1){
    stop("cannot operate on length-1 series");
  }
  
  
  if (nt == 1){
    warning("forced taper length");
    taper_vec.fill(tapers[0]);
  } else {
    taper_vec = tapers;
  }
  
  
  // set the current number of tapers, limited by a few factors
  arma::uvec wh = arma::find(taper_vec > nhalf);
  taper_vec(wh).fill(nhalf);
  wh = arma::find(taper_vec > tapcap);
  taper_vec(wh).fill(tapcap);
  wh = arma::find(taper_vec <= 0);
  taper_vec(wh).ones();
  
  
  int mm = taper_vec.max();
  
  arma::vec w(mm);
  arma::field<arma::vec> para = parabolic_weights_field(mm);
  
  
  //
  // Calculate the psd by averaging over tapered estimates
  //
  
  
  arma::cx_cube psd(nt, nc, nc);
  psd.fill(arma::cx_double(0.0, 0.0));
  
  // each frequency
  for (int jj = 0; jj < nt; jj++) {
    
    m2 = 2 * jj;
    // number of tapers applied at a given frequency (do not remove m+1 index!)
    
    Kc = taper_vec[jj]; // orig: m+1, but this was leading to an indexing bug
    
    // taper sequence and spectral weights
    w = para(Kc-1, 0);
    
    // scan across ki to get vector of
    // spectral differences based on modulo indices
    for (int ik = 0; ik < Kc; ik++) {
      
      wi = w[ik];
      ki = ik + 1;
      
      mleft1 = m2 + ne2 - ki;
      mleft2 = m2 + ki;
      
      j1 = mleft1 % ne2;
      j2 = mleft2 % ne2;
      
      // sum as loop progresses, auto and cross-spectra
      
      for (int ii = 0; ii < nc; ii++) {
        for (int kk = ii; kk < nc; kk++) {
          
          fdc1 = fftz(j1, ii) - fftz(j2, ii);
          fdc2 = fftz(j1, kk) - fftz(j2, kk);
          
          psd(jj, ii, kk) += (fdc1 * std::conj(fdc2)) * wi;
        }
      }
    }
  }
  
  // Add lower triangle for use in transfer function calculation
  for (int ii = 0; ii < nc; ii++) {
    for (int kk = 0; kk < nc; kk++) {
      if (kk < ii){
        psd(arma::span(), arma::span(ii), arma::span(kk)) = 
          arma::conj(psd(arma::span(), arma::span(kk), arma::span(ii)));
      }
    }
  }
  
  // return list to match previous function definition - jrk
  return Rcpp::List::create(
    Named("freq.inds") = arma::regspace<arma::vec>(1, nt),
    Named("k.capped") = taper_vec(arma::span(0, nt-1)),
    Named("psd") = psd
  );
}

// End serial version ----------------------------------------------------------



// Parallel version ------------------------------------------------------------



// [[Rcpp::export]]
arma::cx_mat calc_psd(const arma::cx_mat& fftz,
                const arma::ivec& taper_vec,
                const arma::field<arma::vec>& para,
                const int j,
                const int ne2) {

  int m2, mleft1, mleft2, Kc, ki, j1, j2;
  int nc = fftz.n_cols;

  double wi = 0.0;
  arma::cx_double fdc1, fdc2;
  
  fdc1 = arma::cx_double(0.0, 0.0);
  fdc2 = arma::cx_double(0.0, 0.0);
  
  arma::cx_mat p(nc, nc);
  p.fill(arma::cx_double(0.0, 0.0));


  m2 = 2 * j;
  // number of tapers applied at a given frequency (do not remove m+1 index!)

  Kc = taper_vec(j); // orig: m+1, but this was leading to an indexing bug

  // taper sequence and spectral weights
  arma::vec w = para(Kc-1, 0);

  // scan across ki to get vector of
  // spectral differences based on modulo indices
  for (int ik = 0; ik < Kc; ik++) {

    wi = w[ik];
    ki = ik + 1;

    mleft1 = m2 + ne2 - ki;
    mleft2 = m2 + ki;

    j1 = mleft1 % ne2;
    j2 = mleft2 % ne2;

    // sum as loop progresses

    for (int ii = 0; ii < nc; ii++) {
      for (int kk = ii; kk < nc; kk++) {

        fdc1 = fftz(j1, ii) - fftz(j2, ii);
        fdc2 = fftz(j1, kk) - fftz(j2, kk);

        p(ii, kk) += (fdc1 * std::conj(fdc2)) * wi;
      }
    }

  }

  return p;

}




struct resample_worker: public Worker {

  const arma::cx_mat& fftz;
  const arma::ivec& taper_vec;
  const arma::field<arma::vec>& para;
  const int ne2;

  arma::cx_cube& psd;


  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  resample_worker(const arma::cx_mat& fftz,
                  const arma::ivec& taper_vec,
                  const arma::field<arma::vec>& para,
                  const int ne2,
                  arma::cx_cube& psd)
    : fftz(fftz), taper_vec(taper_vec), para(para), ne2(ne2), psd(psd) {}

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t ip = begin; ip < end; ip++) {

      psd(arma::span(ip), arma::span(), arma::span()) = calc_psd(fftz, taper_vec, para, ip, ne2);

    }
  };
};


//' @title Resample an fft using varying numbers of sine tapers
//'
//' @description
//' Produce an un-normalized psd based on an fft and a vector of optimal sine tapers
//'
//' @details
//' To produce a psd estimate with our adaptive spectrum estimation method, we need only make one
//' fft calculation initially and then
//' apply the weighting factors given by \code{\link{parabolic_weights_rcpp}}, which this
//' function does.
//'
//' @param fftz complex; a matrix representing the dual-length \code{\link{fft}}; see also the \code{dbl} argument
//' @param tapers integer; a vector of tapers
//' @param verbose logical; should messages be given?
//' @param dbl logical; should the code assume \code{fftz} is dual-length or single-length?
//' @param tapcap integer; the maximum number of tapers which can be applied; note that the length is
//' automatically limited by the length of the series.
//'
//' @seealso \code{\link{riedsid}}
//'
//' @examples
//' fftz <- complex(real=1:8, imaginary = 1:8)
//' taps <- 1:4
//' try(resample_fft_rcpp2(fftz, taps))
//'
//' @export
// [[Rcpp::export]]
List resample_fft_parallel(const arma::cx_mat& fftz,
                           const arma::ivec& tapers,
                           bool verbose = true,
                           const bool dbl = true,
                           const int tapcap = 10000 ) {

  //
  // resample and reweight an fft estimates for a given number of tapers
  // using quadratic scaling in Barbour and Parker (2014)
  //
  //
  // needs:
  //  - fftz: complex vector -- the FFT of the original signal
  //  - tapers: integer vector -- the number of tapers at each frequency
  //
  // optional args:
  //  - verbose: logical -- should warnings be given?
  //  - dbl: logical -- should the progam assume 'fftz' is for a double-length (padded) series?
  //                    Otherwise a single-length series is assumed.
  //  - tapcap: integer -- the maximum number of tapers at any frequency
  //

  int sc, nf, nt, ne, ne2, nhalf, nc;

  if (dbl){
    // double-length fft estimates assumed by default
    sc = 2;
  } else {
    // but could be single-length
    sc = 1;
  }

  // even, double, and half lengths
  nf = fftz.n_rows / sc;
  nc = fftz.n_cols;
  nt = tapers.n_elem;
  ne = nf - (nf % 2);

  if (verbose){
    Function msg("message");
    msg(std::string("\tfft resampling"));
  }

  if (ne < nf){
    warning("fft was not done on an even length series");
  }

  ne2 = 2 * ne;
  nhalf = ne / 2;

  arma::ivec taper_vec(nhalf);

  if (nhalf < 1){
    stop("cannot operate on length-1 series");
  }


  if (nt == 1){
    warning("forced taper length");
    taper_vec.fill(tapers[0]);
  } else {
    taper_vec = tapers;
  }


  // set the current number of tapers, limited by a few factors
  arma::uvec wh = arma::find(taper_vec > nhalf);
  taper_vec(wh).fill(nhalf);
  wh = arma::find(taper_vec > tapcap);
  taper_vec(wh).fill(tapcap);
  wh = arma::find(taper_vec <= 0);
  taper_vec(wh).ones();


  int mm = taper_vec.max();

  arma::field<arma::vec> para = parabolic_weights_field(mm);

  arma::cx_cube psd(nt, nc, nc);
  psd.fill(arma::cx_double(0.0, 0.0));
  
  //
  // Calculate the psd by averaging over tapered estimates
  //


  resample_worker calc_resample_fft(fftz, taper_vec, para, ne2, psd);
  
  RcppParallel::parallelFor(0, nt, calc_resample_fft);

  
  // Add lower triangle for use in transfer function calculation
  for (int ii = 0; ii < nc; ii++) {
    for (int kk = 0; kk < nc; kk++) {
      if (kk < ii){
        psd(arma::span(), arma::span(ii), arma::span(kk)) =
          arma::conj(psd(arma::span(), arma::span(kk), arma::span(ii)));
      }
    }
  }
  
  // return list to match previous function definition - jrk
  return Rcpp::List::create(
    Named("freq.inds") = arma::regspace<arma::vec>(1, nhalf),
    Named("k.capped") = taper_vec,
    Named("psd") = psd
  );
}

// End parallel version --------------------------------------------------------



/*** R
# 
# library(fst)
# library(data.table)
# library(waterlevel)
# library(psd)
# wl <- read_fst('/media/kennel/Data/phd/final/thesis-master (2)/data/dist_lag/hess_2018_kennel_2019.fst', as.data.table = TRUE)
# wl <- na.omit(wl)
# library(microbenchmark)
# library(fftw)
# library(fftwtools)
# 
# n. <- 86400 * 62 - 1
# wl <- wl[(86400*2):(86400*2 + n.)]
# wl[, RD130 := lm(RD130~datetime)$residuals]
# wl[, baro := lm(baro~datetime)$residuals]
# wl[, seismic := lm(seismic~datetime)$residuals]
# wl[, et := lm(et~datetime)$residuals]
# # tmp <- pilot_spec(wl$RD130)
# # taps <- riedsid2(tmp)
# fftz <- matrix(c(FFT(wl$RD130), 
#                  FFT(wl$baro),
#                  FFT(wl$et),
#                  FFT(wl$seismic)),
#                ncol = 4)
# 
# fftz[1,] <- (fftz[1,] + fftz[nrow(fftz),])/2
# fftz[nrow(fftz),] <- (fftz[1,] + fftz[nrow(fftz),])/2
# 
# n <- 100000
# fftz <- fftz[(nrow(fftz)/2):1, ][1:n,]
# taps  <- seq(7, 5000, 30000)
# 
# system.time(
#   psd_p   <- resample_fft_parallel(fftz, taps, verbose = FALSE, tapcap = 12000)$psd
# )
# 
# tf_new <- solve_tf_parallel(psd_p)
# 
# gain  <- Mod(tf_new)
# ph  <- Arg(tf_new)
# 
# nnn <- 200
# n_padded <- nrow(wl)#nrow(pgram)
# df       <- 1 / n_padded
# frequency <- seq.int(from = df, by = df, 
#                      length.out = n_padded) * 86400 * 2
# plot(frequency[1:nnn], as.numeric(gain[1:nnn, 3]), type='l', log = 'xy', ylim = c(1e-7, 0.00001), col = 'black')
# 
# a <- approx(x = frequency[1:nnn], 
#             y = as.numeric(gain[1:nnn, 1]),
#             xout = 10^seq(-1, 3, 0.1))
# plot(a, log = 'x', pch =20)
# 
# 
# all.equal(psd_p, psd_s)
system.time({

taps  <- round(seq(9, 5000, length.out = 50000))
tmp <- transfer_fun(wl, vars = c('RD130', 'baro', 'et'), time = 'datetime', method = 'spec_multitaper', tapers = taps)
}
)
plot(gain_RD130_baro~frequency, tmp, log = 'x', type='l')
pspectrum(wl$RD130[1:86400*3])
plot(RD130~frequency, tmp, log = 'xy', type='l')
grid()
*/
