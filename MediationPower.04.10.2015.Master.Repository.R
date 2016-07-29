############################################################
##  Author: Hendricks Brown et. al.
##  Date: 5/01/2016

##  Desc: This program estimate marginal means for mediation path a,
##  the relation of the independent variable to the mediator, and path
##  b, the relation of the mediator to the outcome, across multiple
##  trials, and the between-trial level variance-covariance matrix
##  based on a bivariate normal distribution. Monte Carlo Simulation
##  of mediation for one single trial using full pivoting/NOBS subjects.

##  Usage: The program reads from a text file containing one number per line.
##         The program prompts for and reads the name of the file.
############################################################


## MediationPower.R  Monte Carlo Simulation  of mediation for 1 study Using full pivoting
#  chb 2015 04 09
# compute power using Monte Carlo for mediation involving one single trial with NOBS subjects
library (MASS )
NOBS <- 200
BIG <- 4000 # number of replicates for the dataset
R <- 1000     # number of Monte Carlo draws for confidence interval
Type1.Error <- 0.05

# parameter values for indirect, ALPHA * BETA and direct GAMMA
ALPHA <- 0.398 # standardized so Var(M | X ) = 1
BETA  <- 0.20  # standardized so Var (Y | X , M ) = 1
GAMMA <- 0        # no direct effect


# confidence interval for ALPHA * BETA using Monte Carlo
#  code adapted and modified from 05142014_example_1_8_trials.R


MeanAtimesB.star<- rep(NA, R)
AtimesBsigma.hat.star<- rep(NA, R)

generate.CI.singleTrial <- function ( a , b , s.a , s.b , NOBS , R ) {
  # Input
  #   a, b coefficients, standard errors, number obsns, number of Monte Carlo simulations
  mu <- c ( a , b )
  #  Delta is the first level variance-cov matrix defined in the text
  Delta <- matrix ( c ( s.a^2 , 0 , 0 , s.b^2 ), ncol = 2 ) # note: only level 1 variances

  # quantities used in the pivot
  ab <- a * b
  se <- sqrt ( a^2 * s.b^2 + b^2 * s.a^2 ) # Sobel std error

  # generate distribution of pivotal quantity involving means and std errors
  #  first generate a's and b's
  X.star <- mvrnorm (n = R , mu, Delta)


  # average mean of a * b (delta.hat from equation 17)
  MeanAtimesB.star <- mean(X.star[,1] * X.star[,2])
  # average stdev of a * b (se delta.hat from equation 17)

  delta.boot <- X.star[,1] * X.star[,2]  # distinct estimates of a * b

  seDelta.boot <- sqrt( X.star[,1]^2 * var (X.star[,2] ) +
                        X.star[,2]^2 * var (X.star[,1]) )


  critical.values <- quantile (  ( delta.boot - MeanAtimesB.star ) / seDelta.boot , probs = c(.975 , .025))

  CI <- a * b - se * critical.values


  names ( CI) <-c("2.5%", "97.5%")


  return ( CI )

}


total.rejected <- c(0,0)
CI.rejected <- 0
for ( i in 1 : BIG) {
  X <- c(rep(1, NOBS / 2), rep(0,NOBS / 2))
  M <- ALPHA * X + rnorm(NOBS)
  Y <- BETA * M + GAMMA * X + rnorm(NOBS)
  ans.a <- summary(lm(M ~ X))$coeff
  a <- ans.a[ 2,1 ]
  s.a <- ans.a[ 2,2 ]
  pval.a <- ans.a[ 2,4 ]

  ans.b <- summary(lm(Y  ~ X + M ))$coeff
  b <- ans.b[3,1]
  s.b <- ans.b [ 3, 2 ]
  pval.b <- ans.b[ 3,4 ]
  total.rejected <- total.rejected +
    c( pval.a < Type1.Error , pval.b < Type1.Error   )

  ans.ab <- generate.CI.singleTrial (a , b , s.a , s.b , NOBS , R )
  CI.rejected <- CI.rejected +
    ( ans.ab[ 2 ] < 0 | ans.ab [1 ] > 0  )
}
power.a.b.ab <- c( total.rejected, CI.rejected ) / BIG
names(power.a.b.ab) <- c("a" , "b" , "ab")
power.a.b.ab




### check against Selig Preacher http://www.quantpsy.org/medmc/medmc.htm
generate.CI.singleTrial (.398 , .20 , sqrt(.02) , sqrt(.005) , NOBS , R )
# 2.5%      97.5%
#  0.01377829 0.19029537
#  Selig Preacher obtains .0139 .0174 but doesn't pivot on std error





