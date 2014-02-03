/*
 * esf.C: Various functions to calculate the elementary symmetric 
 * functions both for Rasch, partial credit and rating scale models
 * Author: Basil Abou El-Komboz (basil.abou-el-komboz at stat.uni-muenchen.de)
 * Restrictions/conventions used:
 * - oj = number of categories - 1
 * - par does include (first) restricted parameter, but NOT restricted beta_j0
 */


/* includes */
#include <R.h>
#include <Rinternals.h>


/* Main interface function called from R */
SEXP esf(SEXP Par, SEXP Oj, SEXP Order, SEXP Algorithm); 

/* Summation algorithm */
void esf_sum(double eps[], int m, int oj[], int rmax, int rcum[], int eps_position[], 
	     int order, double gamma0[], double gamma1[]); 

/* Difference algorithm */
void esf_diff(double eps[], int npar, int m, int oj[], int rmax, int rcum[], int eps_position[], 
	      double gamma0[], double gamma1[]);


/* 
 * esf(): Main interface function, called from R via .Call() interface. Calls appropriate calculating function.
 */
SEXP esf(SEXP Par, SEXP Oj, SEXP Order, SEXP Algorithm)
{
  /* access variables */
  int npar = length(Par), m = length(Oj), *oj = INTEGER(Oj), order = INTEGER(Order)[0], algo = INTEGER(Algorithm)[0];
  double *par = REAL(Par);

  /* loop variables */
  int  i, j, r, ncol, ocol, pp; 

  /* calculate cumulated scores, maximum score and start positions of parameters */
  int rmax, *rcum = (int *) R_alloc(m, sizeof(int)), *eps_position = (int *) R_alloc(m, sizeof(int));
  rcum[0] = oj[0] + 1;		/* +1 adds score zero */
  eps_position[0] = 0;
  for(i = 1; i < m; i++) {
    rcum[i] = rcum[i-1] + oj[i];
    eps_position[i] = eps_position[i-1] + oj[i-1];
  }
  rmax = rcum[m-1];

  /* prepare parameters: beta -> exp(-beta) =: eps */
  double *eps = (double *) R_alloc(npar, sizeof(double));
  for (i = 0; i < npar; i++) {
    if (ISNA(par[i])) eps[i] = 0;
    else eps[i] = exp(-par[i]);
  }
  
  /* calculate esf matrix cgamma0 with esf functions, then create SEXP gamma0 and copy last column of cgamma0 with esf with all items into it */
  double *cgamma0 = (double *) R_alloc(m * rmax, sizeof(double));
  esf_sum(eps, m, oj, rmax, rcum, eps_position, 0, cgamma0, NULL);
  SEXP gamma0;
  PROTECT(gamma0 = allocVector(REALSXP, rmax));
  double *rgamma0 = REAL(gamma0);
  memcpy(rgamma0, &cgamma0[(m-1)*rmax], rmax * sizeof(double));

  /* create result list, attach SEXP gamma0 */
  SEXP gammas;
  PROTECT(gammas = allocVector(VECSXP, 1 + order));
  SET_VECTOR_ELT(gammas, 0, gamma0);

  /* if first order requested */
  if (order == 1) {
    /* create SEXP gamma1_woi (woi = without inner derivative) */
    SEXP gamma1_woi;
    PROTECT(gamma1_woi = allocMatrix(REALSXP, rmax, m));
    double *rgamma1_woi = REAL(gamma1_woi);

    /* calculate first order with requested algorithm, stored in gamma1_woi */
    switch(algo) {
    case 0: esf_sum(eps, m, oj, rmax, rcum, eps_position, 1, cgamma0, rgamma1_woi); break;     /* sum. alg. */
    case 1: esf_diff(eps, npar, m, oj, rmax, rcum, eps_position, cgamma0, rgamma1_woi); break; /* diff. alg. */
    default: error("Wrong algorithm code.\n"); break;
    }

    /* initialize and clear final SEXP gamma1 */
    SEXP gamma1;
    PROTECT(gamma1 = allocMatrix(REALSXP, rmax, npar));
    double *rgamma1 = REAL(gamma1);
    for (i = 0; i < rmax * npar; i++) 
      rgamma1[i] = 0.0;
    
    /* take rgamma1_woi, expand (copy) each column by the number of parameters per item,
       include inner derivative and shift to the right position */
    for (i = 0, ocol = 0, ncol = 0, pp = 0; i < m; i++, ocol += rmax) {
      for (j = 1; j <= oj[i]; j++, pp++, ncol += rmax) {
	for (r = 0; r < (rmax - j); r++) {
	  rgamma1[r + j + ncol] = rgamma1_woi[r + ocol] * eps[pp];
	}
      }
    }
    
    /* attach to result SEXP gammas */
    SET_VECTOR_ELT(gammas, 1, gamma1);
    UNPROTECT(2);
  }

  UNPROTECT(2);
  return(gammas);
}


/*
 * esf_sum(): Calculates zero and (if requested) first order ESF with summation algorithm.
 * See e.g. Fisher & Ponocny (1995, S.361; 1994, S.183)
 */
void esf_sum(double eps[], int m, int oj[], int rmax, int rcum[], int eps_position[], int order, double gamma0[], double gamma1[])
{
  /* loop variables */
  int i, j, k, r, ncol, ocol;

  /* zero order */
  if (order == 0) {

    /* clear given gamma0 */
    for (k = 0; k < m * rmax; k++) {
      if ((k % rmax) == 0) {	/* gamma_0(X Items) = 1 */
    	gamma0[k] = 1.0;
      } else if (k <= oj[0]) {
    	gamma0[k] = eps[k-1];	/* gamma_r(1 Item) = eps_jr */
      } else {
    	gamma0[k] = 0.0;	/* rest is zero */
      }
    }

    /* summation algorithm, zero order */
    for (i = 1; i < m; i++) {	/* successively add items */
      /* calculate column indices once */
      ncol = i * rmax;
      ocol = (i-1) * rmax;
      /* loop through possible scores with i items */
      for (r = 1; r < rcum[i]; r++) {
	gamma0[r + ncol] = gamma0[r + ocol];	 /* already score r with i-1 items or ... */
	for (k = 0; (k < oj[i]) & (k <= r); k++) /* ..score r-k and now cat k (0 means 1). */
	  gamma0[r + ncol] += gamma0[r - (k + 1) + ocol] * eps[eps_position[i]+k];
      }
    }

    /* first order */
  } else if (order == 1) {
    
    /* create backup element */
    double *gamma1_tmp = (double *) R_alloc(m * rmax, sizeof(double));

    /* clear and initialize gamma1, gamma1_tmp */
    for (k = 0; k < m * rmax; k++) {
      if ((k % rmax) == 0)
     	gamma1[k] = gamma1_tmp[k] = 1.0; /* gamma_0^(j)(x Items) = 1.0 for all j, x */
      else
	gamma1[k] = gamma1_tmp[k] = 0.0;
    }
    
    /* summation algorithm, first order */
    for (i = 1; i < m; i++) {			    /* successively add items */
      for (j = 0; j < i; j++) {			    /* derivatives of item j (given i items) */
	ncol = j * rmax;			    /* col index */
	for (r = 1; r < rcum[i]; r++) {		    /* scores with i items and without item j */
	  gamma1[r + ncol] = gamma1_tmp[r + ncol];  /* elem from round before (because eps_i0=1) ...  */
	  for (k = 0; (k < oj[i]) && (k <= r); k++) /* ... and score r-k and now cat k (0 means 1. */
	    gamma1[r + ncol] += gamma1_tmp[r - (k + 1) + ncol] * eps[eps_position[i]+k];
	}
      }
      memcpy(&gamma1[i* rmax], &gamma0[(i - 1) * rmax], sizeof(double) * rmax); /* cronecker case, j = i */
      memcpy(gamma1_tmp, gamma1, m * rmax * sizeof(double));			/* backup current round for use in next round */
    }
  }
}


/*
 * esf_diff(): Calculates first order ESF with difference algorithm.
 * See e.g. Fisher & Ponocny (1995, S.361; 1994, S.183/184)
 */
void esf_diff(double eps[], int npar, int m, int oj[], int rmax, int rcum[], int eps_position[], double gamma0[], double gamma1[])
{
  /* loop and helper variables */
  int i, r, k, ri_max, ncol, max_col;
  double *eps_tmp = (double *) R_alloc(npar, sizeof(double));

  /* clear and initialize gamma1 */
  for (k = 0; k < m * rmax; k++) {
    if ((k % rmax) == 0)
      gamma1[k] = 1.0;		/* gamma_0^(j)(x Items) = 1.0 for all j, x */
    else
      gamma1[k] = 0.0;
  }

  /* difference algorithm, first order */
  for (i = 0; i < m; i++) {	 /* derivatives for item i */
    /* calculate loop border (adjusted rmax or adjusted rmax/2 if error_handling > 1) */
    ri_max = (rmax - 1) - oj[i]; /* rmax-1 because rmax includes zero */
    /* calculate col index */
    ncol = i * rmax;
    max_col = (m - 1) * rmax;
    /* bottom-up, starting from r = 1 to rmax */
    for (r = 1; r <= ri_max; r++) {	      /* all possible scores until border */
      gamma1[r + ncol] = gamma0[r + max_col]; /* score of r with all items, calculated before */
      for (k = 0; (k < oj[i]) & (k < r); k++) /* subtract possible category/score combinations */
	gamma1[r + ncol] -= gamma1[r - (k + 1) + ncol] * eps[eps_position[i]+k];
    }
  }
}
