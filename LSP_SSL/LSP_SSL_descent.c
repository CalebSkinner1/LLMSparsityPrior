#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/Memory.h>

double crossprod(double *X, double *y, int n, int j);
double sum(double *x, int n);
double update_sigma2(double *r, int n);
int checkConvergence(double *beta, double *beta_old, double eps, int l, int p);
double SSL(double z, double beta, double lambda0, double lambda1, double theta, double v, double norm, double delta, double sigma2);
double pstar(double x, double theta, double lambda1, double lambda0);
double lambdastar(double x, double theta, double lambda1, double lambda0);
double expectation_approx(double *beta, double a, double b, int p, int l);
double threshold(double theta, double sigma2, double lambda1, double lambda0, double norm);

// Clamp theta to (0, 1) as an open interval.
// theta = 0 : division by zero in pstar -> Inf propagates everywhere.
// theta = 1 : numerically safe but degenerate (collapses to LASSO with lambda1);
//             clamped conservatively to keep behaviour well-defined.
#define THETA_EPS 1e-8
static inline double clamp_theta(double x) {
  if (x < THETA_EPS)       return THETA_EPS;
  if (x > 1.0 - THETA_EPS) return 1.0 - THETA_EPS;
  return x;
}

SEXP cleanupG(double *a, double *r, int *e1, int *e2, double *z,
              SEXP beta, SEXP loss, SEXP iter,
              SEXP thetas_export, SEXP sigmas_export) {
  R_Free(a);
  R_Free(r);
  R_Free(e1);
  R_Free(e2);
  R_Free(z);

  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, thetas_export);
  SET_VECTOR_ELT(res, 4, sigmas_export);
  UNPROTECT(1);
  return(res);
}

double gLoss(double *r, int n) {
  double l = 0;
  for (int i = 0; i < n; i++) l += pow(r[i], 2);
  return(l);
}

SEXP SSL_gaussian(SEXP X_, SEXP y_, SEXP initialbeta_, SEXP penalty_, SEXP variance_,
                  SEXP lambda1_, SEXP lambda0s_, SEXP theta_, SEXP v_vec_,
                  SEXP sigma_, SEXP min_sigma2_, SEXP a_, SEXP b_,
                  SEXP eps_, SEXP max_iter_, SEXP counter_) {

  int n = length(y_);
  int p = length(X_) / n;
  int L = length(lambda0s_);

  double *X           = REAL(X_);
  double *y           = REAL(y_);
  double *initialbeta = REAL(initialbeta_);
  double *v_vec       = REAL(v_vec_);

  const char *penalty  = CHAR(STRING_ELT(penalty_, 0));
  const char *variance = CHAR(STRING_ELT(variance_, 0));

  double lambda1     = REAL(lambda1_)[0];
  double *lambda0s   = REAL(lambda0s_);
  double lambda0;
  double sigma2      = pow(REAL(sigma_)[0], 2);
  double sigma2_init = pow(REAL(sigma_)[0], 2);
  double min_sigma2  = REAL(min_sigma2_)[0];
  double aa          = REAL(a_)[0];
  double bb          = REAL(b_)[0];
  double eps         = REAL(eps_)[0];
  int max_iter       = INTEGER(max_iter_)[0];
  int count_max      = INTEGER(counter_)[0];

    // Read sparsity as a length-p vector (always passed as such from R)
  double *sparsity_in = REAL(theta_);

  // Initialise scalar s for adaptive updates and export.
  // Use the mean of the input vector as the representative scalar.
  double s = 0.0;
  for (int j = 0; j < p; j++) s += sparsity_in[j];
  s /= p;

  // Initialise theta_vec = clamp(s * v_vec[j]).
  // clamp_theta guards against theta=0 (division by zero in pstar)
  // and theta=1 (degenerate collapse to plain LASSO).
  double *theta_vec = R_Calloc(p, double);
  for (int j = 0; j < p; j++) {
    theta_vec[j] = clamp_theta(sparsity_in[j] * v_vec[j]);
  }

  SEXP xnorm_;
  PROTECT(xnorm_ = allocVector(REALSXP, p));
  double *xnorm = REAL(xnorm_);
  for (int j = 0; j < p; j++) {
    xnorm[j] = 0;
    for (int i = 0; i < n; i++) xnorm[j] += pow(X[j * n + i], 2);
  }

  SEXP res, beta, loss, iter, thetas_export, sigmas_export;

  PROTECT(beta = allocVector(REALSXP, L * p));
  double *b = REAL(beta);
  for (int j = 0; j < L * p; j++) b[j] = 0;

  PROTECT(thetas_export = allocVector(REALSXP, L));
  double *thetas = REAL(thetas_export);

  PROTECT(sigmas_export = allocVector(REALSXP, L));
  double *sigmas = REAL(sigmas_export);

  for (int l = 0; l < L; l++) { thetas[l] = NAN; sigmas[l] = NAN; }

  PROTECT(loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i = 0; i < L; i++) INTEGER(iter)[i] = 0;

  double *delta = R_Calloc(p, double);
  double *a     = R_Calloc(p, double);
  double *newa  = R_Calloc(p, double);
  for (int j = 0; j < p; j++) { a[j] = initialbeta[j]; newa[j] = initialbeta[j]; }

  int *e1 = R_Calloc(p, int);
  int *e2 = R_Calloc(p, int);
  for (int j = 0; j < p; j++) { e1[j] = 1 - (a[j] == 0); e2[j] = 1 - (a[j] == 0); }

  double *r   = R_Calloc(n, double);
  double *XTY = R_Calloc(p, double);
  double *XTX = R_Calloc(p * p, double);

  for (int i = 0; i < n; i++) {
    r[i] = y[i];
    for (int j = 0; j < p; j++) r[i] -= X[j * n + i] * a[j];
  }

  double *z = R_Calloc(p, double);
  for (int j = 0; j < p; j++) z[j] = crossprod(X, r, n, j);

  int thres = n;
  if (thres > max_iter) thres = max_iter;

  if (p < thres) {
    for (int j = 0; j < p; j++) {
      XTY[j] = 0;
      for (int i = 1; i < n; i++) XTY[j] += X[j * n + i] * y[i];
    }
    for (int i = 0; i < p; i++)
      for (int j = 0; j < p; j++) {
        XTX[j * p + i] = 0;
        for (int k = 0; k < n; k++) XTX[j * p + i] += X[i * n + k] * X[j * n + k];
      }
  }

  double *thresholds = R_Calloc(p, double);
  double *cutoff     = R_Calloc(p, double);

  int converged      = 0;
  int counter        = 0;
  int violations     = 0;
  int estimate_sigma = 0;

  // ---- Regularisation path ----
  for (int l = 0; l < L; l++) {

    R_CheckUserInterrupt();
    lambda0 = lambda0s[l];

    if (l != 0) {

      if (strcmp(penalty, "adaptive") == 0) {
        s = expectation_approx(b, aa, bb, p, l - 1);
        for (int j = 0; j < p; j++) theta_vec[j] = clamp_theta(s * v_vec[j]);
      }

      if (strcmp(variance, "unknown") == 0) {
        if (INTEGER(iter)[l - 1] < 100) {
          estimate_sigma = 1;
          sigma2 = update_sigma2(r, n);
          if (sigma2 < min_sigma2) { sigma2 = sigma2_init; estimate_sigma = 0; }
        } else {
          estimate_sigma = 0;
          if (INTEGER(iter)[l - 1] == max_iter) sigma2 = sigma2_init;
        }
      }

      for (int j = 0; j < p; j++) {
        thresholds[j] = threshold(theta_vec[j], sigma2, lambda1, lambda0, xnorm[j]);
        cutoff[j]     = thresholds[j];
        delta[j]      = thresholds[j];
      }
      for (int j = 0; j < p; j++) { if (fabs(z[j]) > cutoff[j]) e2[j] = 1; }

    } else {

      for (int j = 0; j < p; j++) {
        thresholds[j] = threshold(theta_vec[j], sigma2, lambda1, lambda0, xnorm[j]);
        cutoff[j]     = thresholds[j];
        delta[j]      = thresholds[j];
      }
      for (int j = 0; j < p; j++) { if (fabs(z[j]) > cutoff[j]) e2[j] = 1; }
    }

    // ---- Inner coordinate descent loops ----
    while (INTEGER(iter)[l] < max_iter) {

      while (INTEGER(iter)[l] < max_iter) {

        while (INTEGER(iter)[l] < max_iter) {

          INTEGER(iter)[l]++;

          for (int j = 0; j < p; j++) {

            if (e1[j]) {

              if (p >= thres)
                z[j] = crossprod(X, r, n, j) + xnorm[j] * a[j];
              else
                z[j] = XTY[j] - crossprod(XTX, newa, p, j) + xnorm[j] * a[j];

              b[l * p + j] = SSL(z[j], a[j], lambda0, lambda1,
                                 theta_vec[j], 1, xnorm[j], delta[j], sigma2);

              if (p >= thres) {
                double shift = b[l * p + j] - a[j];
                if (shift != 0)
                  for (int i = 0; i < n; i++) r[i] -= shift * X[j * n + i];
              } else {
                newa[j] = SSL(z[j], a[j], lambda0, lambda1,
                              theta_vec[j], 1, xnorm[j], delta[j], sigma2);
              }

              counter++;
            }

            if (counter == count_max) {

              if (strcmp(penalty, "adaptive") == 0) {
                s = expectation_approx(b, aa, bb, p, l);
                for (int j = 0; j < p; j++) {
                  theta_vec[j] = clamp_theta(s * v_vec[j]);
                  delta[j]     = threshold(theta_vec[j], sigma2, lambda1, lambda0, xnorm[j]);
                }
              }

              if (strcmp(variance, "unknown") == 0 && estimate_sigma) {
                sigma2 = update_sigma2(r, n);
                if (sigma2 < min_sigma2) sigma2 = sigma2_init;
              }

              counter = 0;
            }
          }

          converged = checkConvergence(b, a, eps, l, p);
          for (int j = 0; j < p; j++) a[j] = b[l * p + j];
          if (converged) break;
        }

        // ---- Check strong-set violations ----
        violations = 0;
        counter    = 0;

        for (int j = 0; j < p; j++) {

          if (e1[j] == 0 && e2[j] == 1) {

            if (p >= thres)
              z[j] = crossprod(X, r, n, j) + xnorm[j] * a[j];
            else
              z[j] = XTY[j] - crossprod(XTX, a, p, j) + xnorm[j] * a[j];

            b[l * p + j] = SSL(z[j], a[j], lambda0, lambda1,
                               theta_vec[j], 1, xnorm[j], delta[j], sigma2);

            if (b[l * p + j] != 0) {
              e1[j] = e2[j] = 1;
              if (p >= thres)
                for (int i = 0; i < n; i++) r[i] -= b[l * p + j] * X[j * n + i];
              a[j] = b[l * p + j];
              violations++;
              counter++;
            }
          }

          if (counter == count_max) {

            if (strcmp(penalty, "adaptive") == 0) {
              s = expectation_approx(b, aa, bb, p, l);
              for (int j = 0; j < p; j++) {
                theta_vec[j] = clamp_theta(s * v_vec[j]);
                delta[j]     = threshold(theta_vec[j], sigma2, lambda1, lambda0, xnorm[j]);
              }
            }

            if (strcmp(variance, "unknown") == 0 && estimate_sigma) {
              sigma2 = update_sigma2(r, n);
              if (sigma2 < min_sigma2) sigma2 = sigma2_init;
            }

            counter = 0;
          }
        }

        if (violations == 0) break;
      }

      // ---- Check full-set violations ----
      int violations = 0;
      counter = 0;

      for (int j = 0; j < p; j++) {

        if (e2[j] == 0) {

          if (p >= thres)
            z[j] = crossprod(X, r, n, j) + xnorm[j] * a[j];
          else
            z[j] = XTY[j] - crossprod(XTX, a, p, j) + xnorm[j] * a[j];

          b[l * p + j] = SSL(z[j], a[j], lambda0, lambda1,
                             theta_vec[j], 1, xnorm[j], delta[j], sigma2);

          if (b[l * p + j] != 0) {
            e1[j] = e2[j] = 1;
            if (p >= thres)
              for (int i = 0; i < n; i++) r[i] -= b[l * p + j] * X[j * n + i];
            a[j] = b[l * p + j];
            violations++;
            counter++;
          }
        }

        if (counter == count_max) {

          if (strcmp(penalty, "adaptive") == 0) {
            s = expectation_approx(b, aa, bb, p, l);
            for (int j = 0; j < p; j++) {
              theta_vec[j] = clamp_theta(s * v_vec[j]);
              delta[j]     = threshold(theta_vec[j], sigma2, lambda1, lambda0, xnorm[j]);
            }
          }

          if (strcmp(variance, "unknown") == 0 && estimate_sigma) {
            sigma2 = update_sigma2(r, n);
            if (sigma2 < min_sigma2) sigma2 = sigma2_init;
          }

          counter = 0;
        }
      }

      if (violations == 0) {
        REAL(loss)[l] = gLoss(r, n);
        thetas[l]     = s;
        sigmas[l]     = sqrt(sigma2);
        break;
      }
    }
  }

  res = cleanupG(a, r, e1, e2, z, beta, loss, iter, thetas_export, sigmas_export);
  UNPROTECT(6);

  R_Free(theta_vec);
  R_Free(delta);
  R_Free(newa);
  R_Free(XTY);
  R_Free(XTX);
  R_Free(thresholds);
  R_Free(cutoff);

  return(res);
}