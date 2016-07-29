# Mediation-Analysis
A Statistical Method for Synthesizing Mediation Analyses Using Product of Coefficient Approach Across Multiple Trials

MultilevelMdiationConfidenceIntervalRandomEffect.R
 Program to calculate a confidence interval for the product of A x B for mediation obtained
 from multiple trials with their own error variances and random effect for both A and B
Copyright (C) Hendricks Brown

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

          Program to calculate a confidence interval for the product of A x B for mediation obtained
          from multiple trials with their own error variances and random effect for both A and B

  Calculates approximate marginal mle's  in a 2 level model and uses this solution to estimate A * B
    then simulates this distribution to obtain confidence intervals

Input: ai and bi are a and b path coefficients from Z to M to Y.
sigma.a and sigma.b are standard errors for each of these

Output: first step is the REML solution. Upper and lower bounds
of a confidence interval for A x B, the EB/ random effects estimate
of this product

 Intermediate Calculations:  the random effect esimates of A x B using R code
     ALSO THE PROPORTION OF THE SIMULATIONS THAT CONVERGE  -- since this simple version is not guaranteed to converge
          the current algorithm used to obtain A * B is based on a bivariate weighted estimate
          If the 2nd level variance is poor



 Part 1:  Define Functions, including minusLogL , MLsolution() to estimate means A, B and 2nd level var-cov
 Part 2:  Input data
 Part 3:  Compute A, B, 2nd level var-cov matrix based on input data
 Part 4:  Use these values as input to Simulate() to generate confidence interval



Input

  ai - vector of regression coefficients for each trial of mediator on covariate (treatment)
  bi - vector of regression coefficients of Y on mediator, adjusted for covariate (or interactions)
  sigma.a - st deviations for each a
  sigma.b -- st deviations for each b

Output
  scalars A, B, and 2 x 2 matrix Sigma
  A is marginal ML for a
  B is marginal ML for b
  Sigma is the trial level var-cov matrix
constants for convergence

 
