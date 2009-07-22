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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include "pmscan_covar.h"
#define TOL 1e-12

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
		    double *genoprob, double ****Genoprob)
{
  int i, j;
  double **a;

  *Genoprob = (double ***)R_alloc(n_gen, sizeof(double **));

  a = (double **)R_alloc(n_pos*n_gen, sizeof(double *));

  (*Genoprob)[0] = a;
  for(i=1; i< n_gen; i++)
    (*Genoprob)[i] = (*Genoprob)[i-1]+n_pos;

  for(i=0; i<n_gen; i++)
    for(j=0; j<n_pos; j++)
      (*Genoprob)[i][j] = genoprob + i*n_ind*n_pos + j*n_ind;
}

/**********************************************************************
 *
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_double(int n, double **vector)
{
  *vector = (double *)R_alloc(n, sizeof(double));
}

/**********************************************************************
 *
 * allocate_dmatrix
 *
 * Allocate space for a matrix of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_dmatrix(int n_row, int n_col, double ***matrix)
{
  int i;

  *matrix = (double **)R_alloc(n_row, sizeof(double *));

  (*matrix)[0] = (double *)R_alloc(n_col*n_row, sizeof(double));

  for(i=1; i<n_row; i++)
    (*matrix)[i] = (*matrix)[i-1]+n_col;
}

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
void reorg_errlod(int n_ind, int n_mar, double *errlod, double ***Errlod)
{
  int i;

  *Errlod = (double **)R_alloc(n_mar, sizeof(double *));

  (*Errlod)[0] = errlod;
  for(i=1; i< n_mar; i++)
    (*Errlod)[i] = (*Errlod)[i-1] + n_ind;
}

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
 * verbose      If > 1 prints error
 *
 **********************************************************************/

void R_pmscan_covar(int *n_ind, int *n_ind_cont, int *n_pos, int *n_gen, 
        double *genoprob, double *addcov, double *addcovcont, 
        int *n_addcov, double *pheno, double *contPheno, int *pointmass, 
        double *bin_start, double *result, double *gamma0, double *gamma1,
        double *beta0, double *beta1, int *maxit, double *tol, int *verbose)
{
  double ***Genoprob, **Addcov, **contAddcov;

  reorg_genoprob(*n_ind, *n_pos, *n_gen, genoprob, &Genoprob);
 
  reorg_errlod(*n_ind, *n_addcov, addcov, &Addcov);
  reorg_errlod(*n_ind_cont, *n_addcov, addcovcont, &contAddcov);

  pmscan_covar(*n_ind, *n_ind_cont, *n_pos, *n_gen, Genoprob, Addcov, 
         contAddcov, *n_addcov, pheno, contPheno, pointmass, bin_start, 
         result, gamma0, gamma1, beta0, beta1, *maxit, *tol, *verbose);
}

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
          double tol, int verbose)
{

/* Declare variables */
  int g, h, i, j, k, s, ss, m, flag=0, error_flag, n_par_cont, n_par_binary;
  double **wts, **contwts;
  double *newparbin, *curparbin, *newparcont, *curparcont;
  double curloglik, newloglik, *work1, *work2, temp;
  double *jac, **Jac, *grad;
  double fit, sum, mutemp, muHat, logdiff;
  double *temp1s, *temp2s, **temp1, **temp2;
  int info;
  double rcond, *junk;

  n_par_cont = 1 + n_gen + n_addcov;
  n_par_binary = n_gen + n_addcov;

  /* Allocate space */
  allocate_dmatrix(n_gen, n_ind, &wts);
  allocate_dmatrix(n_gen, n_ind_cont, &contwts);
  allocate_double(n_par_binary, &newparbin);
  allocate_double(n_par_binary, &curparbin);
  allocate_double(n_par_cont, &newparcont);
  allocate_double(n_par_cont, &curparcont);
  work1 = (double *)R_alloc((n_par_cont-1)*(n_par_cont-1),sizeof(double));
  work2 = (double *)R_alloc(n_par_cont-1, sizeof(double));

  allocate_double(n_par_binary*n_par_binary, &jac);
  reorg_errlod(n_par_binary, n_par_binary, jac, &Jac);
  allocate_double(n_par_binary, &grad);

  allocate_double(n_ind, &temp1s);
  allocate_double(n_ind, &temp2s);
  allocate_double(n_par_binary, &junk);

  allocate_dmatrix(n_ind, n_gen, &temp1);
  allocate_dmatrix(n_ind, n_gen, &temp2);

 
 /* Begin Genome scane*/
for(g=0; g<n_pos; g++) { /* loop over marker positions */
    /* get initial weights for full data set*/
   
    for(j=0; j<n_ind; j++)  
      for(k=0; k<n_gen; k++) 
	       wts[k][j] = Genoprob[k][g][j];
	    
    /* get initial weights for continuous data set*/
    h=0;
    for(j=0; j<n_ind; j++) {
      if(!pointmass[j]){
        for(k=0; k<n_gen; k++) {
  	       contwts[k][h] = Genoprob[k][g][j];
        }
        h++;
      }
    }  
   
   /* Get starting parameter values for binary and continuous traits*/
  
   for(i=0; i<n_par_binary; i++) curparbin[i] = bin_start[i];    
     
    mstep_cont(n_ind_cont, n_gen, contAddcov, n_addcov, contPheno, contwts,
         curparcont, work1, work2, &error_flag);              

   if(!error_flag) {  /* only proceed if no error */

   /* Calculate starting likelihood */

      /* Calculate starting likelihood */

    curloglik = 0.0;               
    for(j=0; j<n_ind; j++) {
      temp=0.0;
      for(k=0; k<n_gen; k++) {
        fit = curparbin[k];           
        for(s=0; s<n_addcov; s++) 
	       fit += Addcov[s][j] * curparbin[n_gen+s];   
                 
      fit = exp(fit);
      
      if(pointmass[j]) temp += Genoprob[k][g][j]* fit/(1.0+fit);
      else { 
         /* Calculate muhat */
         mutemp=0.0;
         /* calculate fitted values for continuous obs, ie, muhat*/ 
        
        for(s=0, ss=n_gen; s<n_addcov; s++, ss++)
          mutemp += (Addcov[s][j]*curparcont[ss]);
 
          muHat = curparcont[k]+mutemp;
          temp += Genoprob[k][g][j]*dnorm(pheno[j],muHat,curparcont[n_par_cont-1],0)/(1.0+fit);
     }
    }
    curloglik += log10(temp);
   }

   /* Begin EM iterations */

   for (m=0; m<maxit; m++){ 

     R_CheckUserInterrupt(); /* Check for ^C */
     
     /* E step */
     /* calculate w(ij) */

      for(j=0; j<n_ind; j++) {
       sum=0.0;
        for(k=0; k<n_gen; k++) {
          fit = curparbin[k];
 
        for(s=0; s<n_addcov; s++) 
  	       fit += Addcov[s][j] * curparbin[n_gen+s];   

        fit = exp(fit);
        
        if(pointmass[j])
          sum += (wts[k][j] = Genoprob[k][g][j]* fit/(1.0+fit)); 
        else { 
         /* Calculate muhat */
         mutemp=0.0;
         /* calculate fitted values for continuous obs, ie, muhat*/ 
        
            for(s=0, ss=n_gen; s<n_addcov; s++, ss++)
              mutemp += (Addcov[s][j]*curparcont[ss]);
     
              muHat = curparcont[k]+mutemp;
              sum += (wts[k][j] = Genoprob[k][g][j]*dnorm(pheno[j],muHat,curparcont[n_par_cont-1],0)/(1.0+fit));
           }     
        }
        for (k=0; k<n_gen; k++)
            wts[k][j] /= sum;        
      }  
      /* Transfer updated wts for continuous data to contwts */
       h=0;
       for(j=0; j<n_ind; j++) {
        if(!pointmass[j]){
          for(k=0; k<n_gen; k++) {
    	       contwts[k][h] = wts[k][j];
          }
          h++;
        }
      }
    
     /* M step */

    /* Update gammas for binary data */    
    /* 0's in gradient and Jacobian */
    for(i=0; i<n_par_binary; i++) {
      grad[i] = 0.0;
      for(s=0; s<n_par_binary; s++)
	      Jac[i][s] = 0.0;
    }

    /* calculate gradient and Jacobian */
    for(j=0; j<n_ind; j++) {
      temp1s[j] = temp2s[j] = 0.0;
      for(k=0; k<n_gen; k++) {
	     fit = curparbin[k];

        for(s=0; s<n_addcov; s++) 
          fit += Addcov[s][j] * curparbin[n_gen+s];
      
     	fit = exp(fit)/(1.0+exp(fit));

    	temp1s[j] += (temp1[j][k] = wts[k][j]*((double)pointmass[j] - fit));
    	temp2s[j] += (temp2[j][k] = wts[k][j]*fit*(1.0-fit));
          }
    }

    for(k=0; k<n_gen; k++) {
      for(j=0; j<n_ind; j++) {
      	grad[k] += temp1[j][k];
      	Jac[k][k] += temp2[j][k];
      }
    }

    for(s=0; s<n_addcov; s++) {
      for(j=0; j<n_ind; j++) {
	     grad[s + n_gen] += Addcov[s][j] * temp1s[j];      
      	for(ss=s; ss<n_addcov; ss++) 
      	  Jac[ss+n_gen][s+n_gen] += Addcov[s][j]*Addcov[ss][j] * temp2s[j];   	
      	for(k=0; k<n_gen; k++) 
      	  Jac[s+n_gen][k] += Addcov[s][j] * temp2[j][k];      
      }
    }
       
    /* dpoco and dposl from Linpack to calculate Jac^-1 %*% grad */
    F77_CALL(dpoco)(jac, &n_par_binary, &n_par_binary, &rcond, junk, &info);
    if(fabs(rcond) < TOL || info != 0) {
      warning("Jacobian matrix is singular.\n");
      result[g] = curloglik;
      break;
     }
    F77_CALL(dposl)(jac, &n_par_binary, &n_par_binary, grad);

    /* revised Gamma estimates */
    for(i=0; i<n_par_binary; i++) 
      newparbin[i] = curparbin[i] + grad[i];

     /* Update Betas and sigma for continuous data */
      mstep_cont(n_ind_cont, n_gen, contAddcov, n_addcov, contPheno, contwts,
        curparcont, work1, work2, &error_flag);
      
      if(error_flag){ /* error: X'X singular; break out of EM */
         warning("X'X matrix is singular.\n");
         flag=0;
         break;
      }


    newloglik = 0.0;               
    for(j=0; j<n_ind; j++) {
      temp=0.0;
      for(k=0; k<n_gen; k++) {
        fit = newparbin[k];           
        for(s=0; s<n_addcov; s++) 
	       fit += Addcov[s][j] * newparbin[n_gen+s];   
                 
      fit = exp(fit);
      
      if(pointmass[j]) temp += Genoprob[k][g][j]* fit/(1.0+fit);
      else { 
         /* Calculate muhat */
         mutemp=0.0;
         /* calculate fitted values for continuous obs, ie, muhat*/ 
        
        for(s=0, ss=n_gen; s<n_addcov; s++, ss++)
          mutemp += (Addcov[s][j]*curparcont[ss]);
 
          muHat = curparcont[k]+mutemp;
          temp += Genoprob[k][g][j]*dnorm(pheno[j],muHat,curparcont[n_par_cont-1],0)/(1.0+fit);
      }
    }
    newloglik += log10(temp);
    }
    logdiff = newloglik-curloglik;

   for(i=0; i<n_par_binary; i++) curparbin[i] = newparbin[i];
  
    /* Check for convergence */ 
    if(newloglik-curloglik < tol) {
       result[g] = newloglik;
       flag = 1; /* set flag to 1 if converged */
       break;  
    }
    /* if no convergence continue with EM*/
     curloglik = newloglik;  

   } /* end EM iteration*/ 
  } /* end procedure if no error flag */ 
     if(!flag) {
        result[g] = NA_REAL;
        warning("Didn't converge or singular matrix. \n");
     }
        
    if(verbose > 1) {
    	if(error_flag) Rprintf("    %3d NA", g+1);
        	Rprintf("\n\n");
    }
    /* Save parameter estimates */
    gamma0[g] = curparbin[0];
    gamma1[g] = curparbin[1];
    beta0[g] = curparcont[0];
    beta1[g] = curparcont[1];    
  } /*end loop over markers */
}

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
		    double *work2, int *error_flag)
{
  int i, j, k, s, sk, nparm1, info;
  double rcond;

  *error_flag=0;

  nparm1 = n_gen + n_addcov;

  /* calculate {E(X)}' y */
  for(j=0; j<nparm1; j++) work2[j] = 0.0;
  for(i=0; i<n_ind_cont; i++) {
    for(j=0; j<n_gen; j++) /* QTL effects */
      work2[j] += (contwts[j][i]*contPheno[i]);
    for(j=0,k=n_gen; j<n_addcov; j++,k++) /* add covar */
      work2[k] += (contAddcov[j][i]*contPheno[i]);
  }

  /* calculate E{X'X}; only the upper right triangle is needed */
  for(j=0; j<nparm1*nparm1; j++) work1[j] = 0.0;
  for(i=0; i<n_ind_cont; i++) {
    for(j=0; j<n_gen; j++) /* QTLxQTL */
      work1[j+nparm1*j] += contwts[j][i];
    for(j=0, k=n_gen; j<n_addcov; j++, k++) {
      for(s=j, sk=k; s<n_addcov; s++, sk++)  /* add x add */
	       work1[k+nparm1*sk] += (contAddcov[j][i]*contAddcov[s][i]);
      for(s=0; s<n_gen; s++) /* QTL x add */
	       work1[s+nparm1*k] += (contAddcov[j][i]*contwts[s][i]);
    }
  }

  /* solve work1 * beta = work2 for beta */
  F77_CALL(dpoco)(work1, &nparm1, &nparm1, &rcond, curparcont, &info);
  if(fabs(rcond) < TOL || info != 0) { /* error! */
    *error_flag = 1;
  }
  else {
    for(j=0; j<nparm1; j++) curparcont[j] = work2[j];
    F77_CALL(dposl)(work1, &nparm1, &nparm1, curparcont);

    /* calculate residual SD */
    curparcont[nparm1] = 0.0;
    for(i=0; i<n_ind_cont; i++) curparcont[nparm1] += contPheno[i]*contPheno[i];
    for(j=0; j<nparm1; j++) curparcont[nparm1] -= (work2[j]*curparcont[j]);
  
    curparcont[nparm1] = sqrt(curparcont[nparm1] / (double)n_ind_cont);
  }
}
