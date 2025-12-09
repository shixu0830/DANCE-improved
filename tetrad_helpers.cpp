#include <Rcpp.h>
#include <algorithm>  // std::min, std::max
#include <cmath> 
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// ---------- Fisher p-values (vectorized) ----------

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

// ---------- Vanishing tetrad helpers ----------

inline double vt_pval(
    int x, int y, int z, int w,
    const RMatrix<double>& Co, double N)
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
  if (!std::isfinite(var)) return 0.0;      // treat as fail
  if (var <= 0.0) var = -var;
  if (var < 1e-12) var = 1e-12;
  
  const double zstat = D / std::sqrt(var);
  return 2.0 * R::pnorm(-std::fabs(zstat), 0.0, 1.0, true, false);
}

inline bool triplet_valid_minp(
    int a, int b, int c,
    int t_idx, int o_idx,
    const RMatrix<double>& Co, double N, double vtalpha,
    double &out_minp)
{
  // six tetrads: (a,b,c;T/O), (a,c,b;T/O), (b,c,a;T/O)
  const double p1 = vt_pval(a,b,c,t_idx,Co,N); if (p1 <= vtalpha) return false;
  const double p2 = vt_pval(a,c,b,t_idx,Co,N); if (p2 <= vtalpha) return false;
  const double p3 = vt_pval(b,c,a,t_idx,Co,N); if (p3 <= vtalpha) return false;
  const double p4 = vt_pval(a,b,c,o_idx,Co,N); if (p4 <= vtalpha) return false;
  const double p5 = vt_pval(a,c,b,o_idx,Co,N); if (p5 <= vtalpha) return false;
  const double p6 = vt_pval(b,c,a,o_idx,Co,N); if (p6 <= vtalpha) return false;
  
  out_minp = std::min(p1, std::min(p2, std::min(p3, std::min(p4, std::min(p5, p6)))));
  return true;
}

// ---------- Pairwise pre-screen (unchanged) ----------

struct PairScreenWorker : public Worker {
  RVector<int>     cidx;
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
        const double p = vt_pval(cidx[i], cidx[j], t_idx, o_idx, Co, N);
        const bool ok = (p > vtalpha);
        const int  v  = ok ? 1 : 0;
        pair_ok(i, j) = v;
        pair_ok(j, i) = v;
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

// ---------- Count (unchanged logic) ----------

inline bool is_triplet_valid_quick(
    int a, int b, int c,
    int t_idx, int o_idx,
    const RMatrix<double>& Co, double N, double vtalpha) {
  double dummy;
  return triplet_valid_minp(a,b,c,t_idx,o_idx,Co,N,vtalpha,dummy);
}

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
          if (is_triplet_valid_quick(cidx[i], cidx[j], cidx[k], t_idx, o_idx, Co, N, vtalpha)) ++cnt;
        }
      }
      counts[i] = cnt;
    }
  }
};

// ---------- Fill + pmins (new) ----------

struct FillTripletsWithPWorker : public Worker {
  RVector<int>     cidx;
  int              t_idx, o_idx;
  RMatrix<double>  Co;
  double           N, vtalpha;
  RMatrix<int>     pair_ok;
  bool             use_pair;
  RVector<int>     offsets;
  RMatrix<int>     results;
  RVector<double>  pmins;
  
  FillTripletsWithPWorker(const IntegerVector& cidx_,
                          int t_idx_, int o_idx_,
                          const NumericMatrix& Co_,
                          double N_, double vtalpha_,
                          const IntegerMatrix& pair_ok_,
                          bool use_pair_,
                          const IntegerVector& offsets_,
                          IntegerMatrix& results_,
                          NumericVector& pmins_)
    : cidx(cidx_), t_idx(t_idx_), o_idx(o_idx_),
      Co(Co_), N(N_), vtalpha(vtalpha_),
      pair_ok(pair_ok_), use_pair(use_pair_),
      offsets(offsets_), results(results_), pmins(pmins_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    const int m = cidx.length();
    for (std::size_t i = begin; i < end; ++i) {
      int write = offsets[i];
      for (int j = static_cast<int>(i) + 1; j < m; ++j) {
        if (use_pair && pair_ok(i, j) == 0) continue;
        for (int k = j + 1; k < m; ++k) {
          if (use_pair && (pair_ok(i, k) == 0 || pair_ok(j, k) == 0)) continue;
          double mp = 0.0;
          if (triplet_valid_minp(cidx[i], cidx[j], cidx[k], t_idx, o_idx, Co, N, vtalpha, mp)) {
            results(write, 0) = cidx[i];
            results(write, 1) = cidx[j];
            results(write, 2) = cidx[k];
            pmins[write] = mp;
            ++write;
          }
        }
      }
    }
  }
};

// [[Rcpp::export]]
List rcpp_screen_triplets_with_p(IntegerVector cidx,
                                 int t_idx, int o_idx,
                                 NumericMatrix Co,
                                 double N, double vtalpha,
                                 IntegerMatrix pair_ok, bool use_pair) {
  const int m = cidx.size();
  if (m < 3) return List::create(_["idx"]=IntegerMatrix(0,3), _["pmin"]=NumericVector(0));
  
  // 1) count
  IntegerVector counts(m, 0);
  CountTripletsWorker counter(cidx, t_idx, o_idx, Co, N, vtalpha, pair_ok, use_pair, counts);
  parallelFor(0, std::max(0, m - 2), counter);
  
  int total = 0;
  for (int i = 0; i < m - 2; ++i) total += counts[i];
  if (total == 0) return List::create(_["idx"]=IntegerMatrix(0,3), _["pmin"]=NumericVector(0));
  
  // 2) offsets
  IntegerVector offsets(m, 0);
  int s = 0;
  for (int i = 0; i < m - 2; ++i) { offsets[i] = s; s += counts[i]; }
  
  // 3) fill (+ pmins)
  IntegerMatrix results(total, 3);
  NumericVector pmins(total);
  FillTripletsWithPWorker filler(cidx, t_idx, o_idx, Co, N, vtalpha, pair_ok, use_pair, offsets, results, pmins);
  parallelFor(0, std::max(0, m - 2), filler);
  
  return List::create(_["idx"]=results, _["pmin"]=pmins);
}