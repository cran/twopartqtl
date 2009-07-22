/*******************************************************************************
 * 
 * pmscan_covar.c
 * *****************************************************************************
 * Functions reorg_genoprob, allocate_double, allocate_dmatrix, reorg_errlod,
 *           mstep_cont (slightly modified version of mstep in R/qtl) 
 *
 * copyright (c) 2004-2006, Karl W Broman
 * 
 * last modified Dec, 2006
 * first written Dec, 2004
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 * *****************************************************************************
 * Functions pmscan_covar, R_pmscan_covar 
 * copyright (c) 2009, Sandra L. Taylor
 * 
 * first written Aug, 2009
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the twopartqtl package
 *
 * These functions are for performing a genome scan with a point mass mixture 
 * trait and a single QTL model in the presence of covariates.
 *
 * Contains: reorg_genoprob, allocate_double, allocate_dmatrix, reorg_errlod, 
 * pmscan_covar, R_pmscan_covar, mstep_cont
 ******************************************************************************/
 /**********************************************************************
 *
 * reorg_genoprob
 *
 * Reorganize the genotype probability data so that it is a triply
 * indexed array rather than a single long vector
 *
 * Afterwards, genoprob indexed like Genoprob[gen][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_genoprob(int n_ind, int n_pos, int n_gen,
		    double *genoprob, double ****Genoprob);

/**********************************************************************
 *
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_double(int n, double **vector);

/**********************************************************************
 *
 * allocate_dmatrix
 *
 * Allocate space for a matrix of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_dmatrix(int n_row, int n_col, double ***matrix);

/**********************************************************************
 *
 * reorg_errlod
 *
 * Just like reorg_geno(), only for a matrix of doubles.
 *
 * Afterwards, errlod indexed like Errlod[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_errlod(int n_ind, int n_mar, double *errlod, double ***Errlod);

/**********************************************************************
 * 
 * pmscan_covar
 *
 * Performs genome scan using interval mapping in the presence of
 * covariates.  (The multipoint genotype probabilities have already 
 * been calculated in calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_ind_cont   Number of individuals with continuous observations
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              indexed as Genoprob[gen][pos][ind]
 *
 * Addcov       Matrix of additive covariates indexed as 
 *              Addcov[cov][ind]
 *
 * contAddcov   Matrix of additive covariates for continuous observations only 
 *              indexed as Addcov[cov][ind]
 *
 *
 * n_addcov     Number of columns in Addcov
 *
 * contPheno        Phenotype data of all continuous observations only
 *
 * pheno        Phenotype data of all observations with point mass values 99999
 *
 * pointmass    Phenotype data for binary data (0,1), as a vector
 *
 * bin_start    Starting values for coefficients from logistic regression (gammas) 
 *              vector of length n_gen + n_addcov 
 *
 * result       Result vector of length n_pos; upon return, contains 
 *              the LOD scores.
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * verbose      If 1, print out log likelihood at each iteration
 *
 **********************************************************************/

void R_pmscan_covar(int *n_ind, int *n_ind_cont, int *n_pos, int *n_gen, 
        double *genoprob, double *addcov, double *addcovcont, 
        int *n_addcov, double *pheno, double *contPheno, int *pointmass, 
        double *bin_start, double *result, double *gamma0, double *gamma1,
        double *beta0, double *beta1, int *maxit, double *tol, int *verbose);

/**********************************************************************
 * 
 * pmscan_covar
 *
 **********************************************************************/

void pmscan_covar(int n_ind, int n_ind_cont, int n_pos, int n_gen, 
		      double ***Genoprob, double **Addcov, double **contAddcov,  
          int n_addcov, double *pheno, double *contPheno,
		      int *pointmass, double *bin_start, double *result, double *gamma0, 
          double *gamma1, double *beta0, double *beta1, int maxit,
          double tol, int verbose);
          
/**********************************************************************
 * 
 * mstep_cont:  M-step of the EM algorithm for continuous component of 
 *              point mass mixture with covariates
 *
 * n_ind_cont   Number of individuals with continuous phenotypes
 *
 * n_gen        Number of possible QTL genotypes
 *
 * contAddcov   Additive covariates for continuous observations
 *
 * n_addcov     Number of columns in Addcov
 *
 * contPheno    Continuous Phenotypes only 
 *
 *
 * contwts      Pr(QTL gen | phenotype, model, multipoint marker data),
 *              indexed as wts[gen][ind] for continuous observations only
 *
 * curparcont    On output, the updated parameter estimates for Beta and sigma
 *
 * work1    Workspace of doubles, of length (n_par-1)*(n_par-1)
 *
 * work2    Workspace of doubles, of length (n_par-1)
 *
 * error_flag  Set to 1 if E(X'X) is singular
 *
 **********************************************************************/

void mstep_cont(int n_ind_cont, int n_gen, double **contAddcov, int n_addcov, 
		    double *contPheno, double **contwts, double *curparcont, double *work1, 
		    double *work2, int *error_flag);

