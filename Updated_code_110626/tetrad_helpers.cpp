// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

// ---------------- Fisher p-values (vectorized) ----------------

// [[Rcpp::export]]
NumericVector rcpp_fisher_pval_vec(NumericVector r, double n) {
  const int m = r.size();
  NumericVector p(m);
  if (n <= 3.0) { std::fill(p.begin(), p.end(), 1.0); return p; }
  const double s = std::sqrt(n - 3.0);
  for (int i = 0; i < m; ++i) {
    double ri = std::max(-0.999999, std::min(0.999999, r[i]));
    double z  = 0.5 * std::log((1.0 + ri) / (1.0 - ri)) * s;
    p[i] = 2.0 * R::pnorm(-std::fabs(z), 0.0, 1.0, true, false);
  }
  return p;
}

// --------------- Vanishing tetrad test (Wishart approx) ----------------

inline bool test_vanishing_tetrad_fast(
    int x, int y, int z, int w,
    const RMatrix<double>& Co, double N, double vtalpha)
{
  // D = (x,z)(y,w) - (x,w)(y,z)
  const double D = Co(x,z) * Co(y,w) - Co(x,w) * Co(y,z);
  
  // 2x2 minors
  const double det_xy = Co(x,x) * Co(y,y) - Co(x,y) * Co(x,y);
  const double det_zw = Co(z,z) * Co(w,w) - Co(z,w) * Co(z,w);
  
  // 4x4 determinant via hard-coded expansion
  const double a11 = Co(x,x), a12 = Co(x,y), a13 = Co(x,z), a14 = Co(x,w);
  const double a22 = Co(y,y), a23 = Co(y,z), a24 = Co(y,w);
  const double a33 = Co(z,z), a34 = Co(z,w);
  const double a44 = Co(w,w);
  
  const double m11 = a22*(a33*a44 - a34*a34) - a23*(a23*a44 - a24*a34) + a24*(a23*a34 - a24*a33);
  const double m12 = a12*(a33*a44 - a34*a34) - a23*(a13*a44 - a14*a34) + a24*(a13*a34 - a14*a33);
  const double m13 = a12*(a23*a44 - a24*a34) - a22*(a13*a44 - a14*a34) + a24*(a13*a24 - a14*a23);
  const double m14 = a12*(a23*a34 - a24*a33) - a22*(a13*a34 - a14*a33) + a23*(a13*a24 - a14*a23);
  const double det4 = a11*m11 - a12*m12 + a13*m13 - a14*m14;
  
  double var = (((N + 1.0) / (N - 1.0)) * det_xy * det_zw - det4) / (N - 2.0);
  if (!std::isfinite(var)) return false;
  if (var <= 0.0) var = -var;          
  if (var < 1e-12) var = 1e-12;       
  
  const double zstat = D / std::sqrt(var);
  const double p = 2.0 * R::pnorm(-std::fabs(zstat), 0.0, 1.0, true, false);
  return p > vtalpha;
}

// ---------------- Pairwise pre-screen (parallel) ----------------

struct PairScreenWorker : public Worker {
  RVector<int>     cidx;     // store views BY VALUE
  int              t_idx, o_idx;
  RMatrix<double>  Co;
  double           N, vtalpha;
  RMatrix<int>     pair_ok;
  
  PairScreenWorker(const IntegerVector& cidx_,
                   int t_idx_, int o_idx_,
                   const NumericMatrix& Co_,
                   double N_, double vtalpha_,
                   IntegerMatrix& pair_ok_)
    : cidx(cidx_), t_idx(t_idx_), o_idx(o_idx_),
      Co(Co_), N(N_), vtalpha(vtalpha_), pair_ok(pair_ok_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t m = static_cast<std::size_t>(cidx.length());
    for (std::size_t i = begin; i < end; ++i) {
      for (std::size_t j = i + 1; j < m; ++j) {
        const bool ok = test_vanishing_tetrad_fast(cidx[i], cidx[j], t_idx, o_idx, Co, N, vtalpha);
        const int  v  = ok ? 1 : 0;
        pair_ok(i, j) = v;
        pair_ok(j, i) = v;   // **make symmetric**
      }
    }
  }
};

// [[Rcpp::export]]
IntegerMatrix rcpp_precompute_pair_ok(IntegerVector cidx,
                                      int t_idx, int o_idx,
                                      NumericMatrix Co,
                                      double N, double vtalpha) {
  const int m = cidx.size();
  IntegerMatrix pair_ok(m, m); // zeros
  if (m >= 2) {
    PairScreenWorker worker(cidx, t_idx, o_idx, Co, N, vtalpha, pair_ok);
    parallelFor(0, m - 1, worker);
  }
  return pair_ok;
}

// --------------- Triplet validity (6 tetrads) ------------------

inline bool is_triplet_valid(int a, int b, int c,
                             int t_idx, int o_idx,
                             const RMatrix<double>& Co, double N, double vtalpha) {
  if (!test_vanishing_tetrad_fast(a, b, c, t_idx, Co, N, vtalpha)) return false;
  if (!test_vanishing_tetrad_fast(a, c, b, t_idx, Co, N, vtalpha)) return false;
  if (!test_vanishing_tetrad_fast(b, c, a, t_idx, Co, N, vtalpha)) return false;
  if (!test_vanishing_tetrad_fast(a, b, c, o_idx, Co, N, vtalpha)) return false;
  if (!test_vanishing_tetrad_fast(a, c, b, o_idx, Co, N, vtalpha)) return false;
  if (!test_vanishing_tetrad_fast(b, c, a, o_idx, Co, N, vtalpha)) return false;
  return true;
}

// ---------------- Count pass triplets (parallel) ---------------

struct CountTripletsWorker : public Worker {
  RVector<int>     cidx;
  int              t_idx, o_idx;
  RMatrix<double>  Co;
  double           N, vtalpha;
  RMatrix<int>     pair_ok;
  bool             use_pair;
  RVector<int>     counts;
  
  CountTripletsWorker(const IntegerVector& cidx_,
                      int t_idx_, int o_idx_,
                      const NumericMatrix& Co_,
                      double N_, double vtalpha_,
                      const IntegerMatrix& pair_ok_,
                      bool use_pair_,
                      IntegerVector& counts_)
    : cidx(cidx_), t_idx(t_idx_), o_idx(o_idx_),
      Co(Co_), N(N_), vtalpha(vtalpha_),
      pair_ok(pair_ok_), use_pair(use_pair_), counts(counts_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    const int m = cidx.length();
    for (std::size_t i = begin; i < end; ++i) {
      int cnt = 0;
      for (int j = static_cast<int>(i) + 1; j < m; ++j) {
        if (use_pair && pair_ok(i, j) == 0) continue;
        for (int k = j + 1; k < m; ++k) {
          if (use_pair && (pair_ok(i, k) == 0 || pair_ok(j, k) == 0)) continue;
          if (is_triplet_valid(cidx[i], cidx[j], cidx[k], t_idx, o_idx, Co, N, vtalpha)) ++cnt;
        }
      }
      counts[i] = cnt;
    }
  }
};

// ---------------- Fill results (parallel) ----------------------

struct FillTripletsWorker : public Worker {
  RVector<int>     cidx;
  int              t_idx, o_idx;
  RMatrix<double>  Co;
  double           N, vtalpha;
  RMatrix<int>     pair_ok;
  bool             use_pair;
  RVector<int>     offsets;
  RMatrix<int>     results;
  
  FillTripletsWorker(const IntegerVector& cidx_,
                     int t_idx_, int o_idx_,
                     const NumericMatrix& Co_,
                     double N_, double vtalpha_,
                     const IntegerMatrix& pair_ok_,
                     bool use_pair_,
                     const IntegerVector& offsets_,
                     IntegerMatrix& results_)
    : cidx(cidx_), t_idx(t_idx_), o_idx(o_idx_),
      Co(Co_), N(N_), vtalpha(vtalpha_),
      pair_ok(pair_ok_), use_pair(use_pair_),
      offsets(offsets_), results(results_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    const int m = cidx.length();
    for (std::size_t i = begin; i < end; ++i) {
      int write = offsets[i];
      for (int j = static_cast<int>(i) + 1; j < m; ++j) {
        if (use_pair && pair_ok(i, j) == 0) continue;
        for (int k = j + 1; k < m; ++k) {
          if (use_pair && (pair_ok(i, k) == 0 || pair_ok(j, k) == 0)) continue;
          if (is_triplet_valid(cidx[i], cidx[j], cidx[k], t_idx, o_idx, Co, N, vtalpha)) {
            results(write, 0) = cidx[i]; // 0-based; map to names in R
            results(write, 1) = cidx[j];
            results(write, 2) = cidx[k];
            ++write;
          }
        }
      }
    }
  }
};

// [[Rcpp::export]]
IntegerMatrix rcpp_screen_triplets(IntegerVector cidx,
                                   int t_idx, int o_idx,
                                   NumericMatrix Co,
                                   double N, double vtalpha,
                                   IntegerMatrix pair_ok, bool use_pair) {
  const int m = cidx.size();
  if (m < 3) return IntegerMatrix(0, 3);
  
  // 1) count
  IntegerVector counts(m, 0);
  CountTripletsWorker counter(cidx, t_idx, o_idx, Co, N, vtalpha, pair_ok, use_pair, counts);
  parallelFor(0, std::max(0, m - 2), counter);
  
  int total = 0;
  for (int i = 0; i < m - 2; ++i) total += counts[i];
  if (total == 0) return IntegerMatrix(0, 3);
  
  // 2) offsets
  IntegerVector offsets(m, 0);
  int s = 0;
  for (int i = 0; i < m - 2; ++i) { offsets[i] = s; s += counts[i]; }
  
  // 3) fill
  IntegerMatrix results(total, 3);
  FillTripletsWorker filler(cidx, t_idx, o_idx, Co, N, vtalpha, pair_ok, use_pair, offsets, results);
  parallelFor(0, std::max(0, m - 2), filler);
  
  return results;
}
