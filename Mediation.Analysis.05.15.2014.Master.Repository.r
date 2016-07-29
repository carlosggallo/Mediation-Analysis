
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



## MultilevelMdiationConfidenceIntervalRandomEffect.R
## #  Program to calculate a confidence interval for the product of A x B for mediation obtained
## #  from multiple trials with their own error variances and random effect for both A and B
## # Copyright (C) Hendricks Brown

## # This program is free software; you can redistribute it and/or
## # modify it under the terms of the GNU General Public License
## # as published by the Free Software Foundation; either version 2
## # of the License, or (at your option) any later version.

## # This program is distributed in the hope that it will be useful,
## # but WITHOUT ANY WARRANTY; without even the implied warranty of
## # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## # GNU General Public License for more details.

## # You should have received a copy of the GNU General Public License
## # along with this program; if not, write to the Free Software
## # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## #           Program to calculate a confidence interval for the product of A x B for mediation obtained
## #           from multiple trials with their own error variances and random effect for both A and B
## #
## #   Calculates approximate marginal mle's  in a 2 level model and uses this solution to estimate A * B
## #     then simulates this distribution to obtain confidence intervals
## #
## #

## # Input: ai and bi are a and b path coefficients from Z to M to Y.
## # sigma.a and sigma.b are standard errors for each of these

## # Output: first step is the REML solution. Upper and lower bounds
## of a confidence interval for A x B, the EB/ random effects estimate
## of this product

## #  Intermediate Calculations:  the random effect esimates of A x B using R code
## #      ALSO THE PROPORTION OF THE SIMULATIONS THAT CONVERGE  -- since this simple version is not guaranteed to converge
## #           the current algorithm used to obtain A * B is based on a bivariate weighted estimate
## #           If the 2nd level variance is poor
## #



 ## Part 1:  Define Functions, including minusLogL , MLsolution() to estimate means A, B and 2nd level var-cov
 ## Part 2:  Input data
 ## Part 3:  Compute A, B, 2nd level var-cov matrix based on input data
 ## Part 4:  Use these values as input to Simulate() to generate confidence interval



## # Input

## #   ai - vector of regression coefficients for each trial of mediator on covariate (treatment)
## #   bi - vector of regression coefficients of Y on mediator, adjusted for covariate (or interactions)
## #   sigma.a - st deviations for each a
## #   sigma.b -- st deviations for each b
## #
## #  Output
## #   scalars A, B, and 2 x 2 matrix Sigma
## #   A is marginal ML for a
## #   B is marginal ML for b
## #   Sigma is the trial level var-cov matrix
## # constants for convergence







  SMALL <- 0.0001

  Bignegative <- -10  # consider e(-10) as 0


#define functions: - log likelihood value

minusLogL.v <-  function ( v ,  Xvals  ,  sigma.aSE  , sigma.bSE, correction="T"  ) {
# return logL value based on the choleski value, use these instead of Sigma to make range unrestricted
#   input v = 2nd level V11 V12, V22 and A, B, XCurrent the matrix of estimated A, B, sigma.aSE: the st devs of A. sigma.bSE: the st devs of B
#  get var-cov matrix values from a choleski decomposition w = log s11, s12*sqrt(s11), log s22.1
#output   negative restricted likelihood.

# cat(" minusLogL.v: dim(Xvals) = " , dim(Xvals) , v ,  "\n")



	A<-v[4]
	B<-v[5]

# Vs are variances for path coefficients A and B across trials, i.e. 2nd level
	V11 <- v[1]
	V12 <- v[2]
	V22 <- v[3]

	V <- matrix ( c( V11 , V12 , V12, V22), nrow = 2  )


    sigma.asCurrent <- sigma.aSE ^2
    sigma.bsCurrent <- sigma.bSE ^2

    L <- nrow ( Xvals )


  # default settings
#settings =list(correction=="F")
   # if settings is supplied
       if (correction=="T"   )
{  f=(L -1)/L
}

  if (correction=="F")
{  f=1
}



   #- restricted minus log likelihood, formula 7,8.
	rslt2 <- 0.5*L*log(2*pi)
	for ( u in 1 : L ) {

		M <- diag ( c( sigma.asCurrent[ u ] , sigma.bsCurrent [ u] ) ) + V
		if ( det ( M ) > 0 )
			ldetM <- log ( det ( M) )
		else {
		    cat (" ** det(M) = ", det(M) , v  ,  sigma.asCurrent[u] , " " , sigma.bsCurrent[u]  , "\n")
			ldetM <-  1000   #  make exceptionally large to force det to be positive
			}
		rslt2 <- rslt2 +
					1/2 *  f * ldetM +
					1/2 * t((Xvals[ u , ]-c(A,B))) %*% solve ( M ) %*% (Xvals[ u , ]-c(A,B))

	}

	return ( rslt2  )
}


minusLogL.w <-  function ( w ,  Xvals ,  sigma.aSE  , sigma.bSE ) {
# return logL value based on the choleski value, use these instead of Sigma to make range unrestricted
#   input w = log choleski vals and A, B, XCurrent the matrix of estimated A, B, sigma.aSE: the st devs of A. sigma.bSE: the st devs of B
#  get var-cov matrix values from a choleski decomposition w = log s11, s12*sqrt(s11), log s22.1
#output   negative restricted likelihood.

# cat( " w:  minusLogL.w dim(Xvals) = ", dim( Xvals) , "\n")
	v <- WtoV ( w )


	ans <- minusLogL.v ( v , Xvals  , sigma.aSE  , sigma.bSE  )

	return ( ans )
}



# define gradients for optim.
#   gradient.w is derivatives w r t log choleski  i.e. log V11 ,  V11 / sqrt ( V11) , log ( V22.1 ) , meanA1 mean B1
  #   input w = log choleski vals and A, B.
#output   5 derivatives for V11, V12, V22, A, B.


 gradient.v <- function ( v ,  Xvals  ,  sigma.aSE  , sigma.bSE, correction="T" ) {
    #  5 arguments
 	# gradient with respect to v = V11, V12, V22 A and B
 	# returns gradient vector of dimension 5
    # cat(" gradient.v: Check  dim(Xvals) = ", dim(Xvals) , Xvals , "\n")
 	K <- nrow ( Xvals )

 	V <- matrix ( c( v[1], v[2], v[2] , v[3] ), ncol = 2 )
 	A <- v[4]
 	B <- v[5]

 	dsigma.v11 <- matrix ( c( 1 , 0 , 0 , 0 ) , ncol = 2 )
 	dsigma.v12 <- matrix ( c( 0 , 1 , 1 , 0 ) , ncol = 2 )
 	dsigma.v22 <- matrix ( c( 0 , 0 , 0 , 1 ) , ncol = 2 )

 	grv11 <- 0
 	grv12 <- 0
 	grv22 <- 0

 	grA <- 0
 	grB <- 0

 	        if (correction=="T")
{  f=(K -1)/K
}

  if (correction=="F")
{  f=1
}


     for ( i in 1 : K  )  {
         M1<-matrix( c(sigma.aSE [ i ]^2 , 0, 0, sigma.bSE [ i ]^2 ), nr=2) + V
         M1Inv <- solve( M1 )

         ssq <-(Xvals[ i, ] - c(A,B)) %*% t(Xvals [ i , ] - c(A,B))  ## 2 x 2 SSQ matrix


          DiffMatrix <- ssq %*% M1Inv - f* diag( c( 1 , 1 ) )
          grv11 <- grv11 - 1/2 * sum( diag ( DiffMatrix %*% dsigma.v11 %*%  M1Inv  ) )     #formula 11 for d( - R Log L )

          grv12 <- grv12 - 1/2 * sum( diag (DiffMatrix %*% dsigma.v12 %*%  M1Inv  ) )

          grv22 <- grv22 - 1/2 * sum( diag (DiffMatrix %*% dsigma.v22 %*%  M1Inv  ) )

          Value <- M1Inv %*%  ( Xvals[ i, ] - c(A , B) )    #formula 10
          grA <- grA -  Value [1,1]
          grB <- grB -  Value [2,1]
     }

     return( c(grv11 , grv12 , grv22 , grA , grB ))

}

check.gradient.v <- function ( v , h = 0.000001 ,  Xvals  ,  sigma.aSE , sigma.bSE  ) {  ## this gives the same result as gradient.v
  gr <- rep( NA, 5 )
  gr[1] <- ( minusLogL.v( v + c( h , 0 , 0 , 0 , 0 )) - minusLogL.v( v ) )/ h
  gr[2] <- ( minusLogL.v( v + c( 0 , h , 0 , 0 , 0 )) - minusLogL.v( v ) )/ h
   gr[3] <- ( minusLogL.v( v + c( 0 , 0 , h , 0 , 0 )) - minusLogL.v( v ) )/ h
    gr[4] <- ( minusLogL.v( v + c( 0 , 0 , 0 , h , 0 )) - minusLogL.v( v ) )/ h
     gr[5] <- ( minusLogL.v( v + c( 0 , 0 , 0 , 0 , h )) - minusLogL.v( v ) )/ h

     return ( gr)

}



gradient.w <-function(w ,  Xvals  ,  sigma.aSE  , sigma.bSE ) {
   # w is 5 dimensional, transformed Sigma followed by A and B
   # return 5 dimensional vector


  v <- WtoV ( w )
  grad.v <- gradient.v ( v ,  Xvals  ,  sigma.aSE  , sigma.bSE  )

           deriv.w <- t(dVdW(w)) %*% grad.v

          return( deriv.w)
}




dVdW <- function ( w ) {
       # get derivs of V wrt W
          dSigma11dw <- c ( exp ( w[1]) , 0 , 0 , 0 , 0 )  #  = exp ( w[1])
          dSigma12dw <- c( 1/2 * w[2] * exp( 1/2 * w[1] ), exp( 1/2 * w[1] ) , 0 , 0 , 0  ) #  = w[2] * e^ 1/2 w[1]
          dSigma22dw <- c( 0 , 2 * w[2] , exp( w[3]) , 0 , 0      ) # = exp( w[3] ) + w[2]^2

          dvdw <- rbind ( dSigma11dw ,  dSigma12dw,  dSigma22dw ,
         					c( 0 , 0 , 0 , 1 , 0 ) , c( 0 , 0 , 0 , 0 , 1 ) )
		  return( dvdw)
}

dWdV <- function ( v ) {

v22.1 <- v[3] - v[2]^2/v[1]

        dw1dv <- c( 1/v[1] , 0 , 0 , 0 , 0 )
        dw2dv <- c( - 1/2 * v[2] / v[1] ^ (1.5) , v[1]^ (-1/2) , 0 , 0 , 0 )
        dw3dv <- c( 1/v22.1 * v[2]^2/v[1]^2 , - 2 /v22.1  * v[2] / v[1] , 1/v22.1 , 0 , 0 )
        dw4dv <- c( 0 , 0 , 0 , 1 , 0 )
        dw5dv <- c( 0 , 0 , 0 , 0 , 1 )

        dwdv <- rbind ( dw1dv, dw2dv, dw3dv, dw4dv, dw5dv )
        return ( dwdv)
}

WtoV <- function ( w ) {


  # returns Var11, Var12, Var22 of dimension 3 if w is of length 3
  # otherwise returns 5 dimensional vector Var11 Var12 Var22 A B
    #logV11 <- w[1]
	#z <- w[2]
	#logV22.1 <- w[3]


	V11 <- exp ( w[1] )
	V12 <-  w[2] * sqrt ( V11 )
	V22 <- (exp ( w[3] )) + (V12^2/exp(w[1]))

	if ( length ( w ) == 3 )
		return( c( V11 , V12 , V22 ))
	else
		return ( c(  V11 , V12 , V22 , w[4:5] ))
}

VtoW <- function ( v ) {

# returns transformed choleski if V is of dim 3
# otherwise returns transformed choleski and A B
   w1 <- log ( v[1] )
   w2 <- v[2] / sqrt ( v[1] )
   w3 <- log ( v[3] - v[2]^2 / v[1] )
   if ( length (v ) == 3 )
   		return ( c( w1 , w2 , w3 ))

   else
   		return ( c( w1 , w2 , w3 , v[4:5] ))

}

 	# MLsolution.v and MLsolution.w   produce MLE :
 	#		input (aiCurrent , biCurrent , sigma.aVals , sigma.bVals)
 	#  calls optim using minusLogL.v or .w as an argument to get ML for transformed var-cov and A, B.
 	#
    #output : v: var-cov values and A and B,  w: transformed var-cov values and A and B


MLsolution.v <-  function ( Xvals , sigma.aSE , sigma.bSE ) {

#cat(" MLsolution.v dim(Xvals) = " , dim(Xvals) , "\n")

	sigma.asVar <- sigma.aSE ^2
	sigma.bsVar <- sigma.bSE ^2


      # define initial values of the Choleski transformed var-cov matrix
      #     w11s = log sigma11  w12s =  sigma12 / sqrt( sigma11 )   w22s = log sigma22.1
      # set the lower bound for v to exp(-10)

      # initial values

    sV11 <- max (   var ( ai ) -  mean (sigma.asVar ) , SMALL )
    sV12  <-		0
    sV22 <-   max (   var ( bi ) -  mean (sigma.bsVar ) , SMALL )

     a1<-mean(Xvals [ , 1 ] )
     b1<-mean(Xvals [ , 2 ] )

     sV <- matrix ( c( sV11 , sV12, sV12, sV22 ), ncol = 2 )

    # cat ( " MLsolution.v: Check det(V) = ", det(V) , " V = ", c( sV11, sV12, sV22, a1, b1 ), "\n")
     #cat ( " MLsolution.v: Check dim Xvals " , dim(Xvals) , "\n")
     #cat ( " MLsolution.v: -RlogL " , minusLogL.v( c( V11, V12, V22, a1, b1) , Xvals , sigma.aSE , sigma.bSE ), "\n")

 # use optim function for ML estimates of between trial level variance and covariance, A and B

     #  cat( "START No DERIV\n")
               ##try optim without gradient.
      answer2.v <- optim(par=c(sV11, sV12 , sV22, a1 , b1 ) ,
      fn=minusLogL.v,  method =c( "L-BFGS-B"),
      lower = c( 0 , -Inf , 0, -Inf, -Inf ), control=list(trace=6, REPORT=1) , hessian = TRUE,
      Xvals = X ,  sigma.aSE = sigma.a , sigma.bSE = sigma.b )




    #  cat( "START DERIV\n")
             ##try optim with gradient.
        answer.v <- optim(par=c(sV11, sV12 , sV22, a1 , b1 ) ,
        fn=minusLogL.v , gr=gradient.v , method =c( "L-BFGS-B"),
       lower = c( 0 , -Inf , 0, -Inf, -Inf ),
       control=list(trace=6, REPORT=1),
       hessian = TRUE ,
       Xvals = X ,  sigma.aSE = sigma.a , sigma.bSE = sigma.b )





	names ( answer.v$par ) <- c( "v11", "v12" , "v22" , "A", "B" )

  	names ( answer2.v$par ) <- c( "v11", "v12" , "v22" , "A", "B" )

	 	variancematrix<-  solve(   answer.v$hessian )  # for Sigma and means

	 	dimnames( variancematrix ) <- list ( c( "v11", "v12" , "v22" , "A", "B" ), c("v11", "v12" , "v22" , "A", "B" ))


	return (  list( answerDeriv = answer.v ,
					 answerNoDeriv = answer2.v
					)
			)
}



MLsolution.w <-  function ( Xvals , sigma.aSE , sigma.bSE ) {

#cat("Entered MLsolution.w\n")


	L <- nrow ( Xvals )



	sigma.asVar <- sigma.a ^2
	sigma.bsVar <- sigma.b^2


      # define initial values of the Choleski transformed var-cov matrix
      #     w11s = log sigma11  w12s =  sigma12 / sqrt( sigma11 )   w22s = log sigma22.1
      # set the lower bound for v to exp(-10)

      # initial values

    w11s <-   log(max (   var ( ai ) -  mean (sigma.asVar ) , SMALL ))
    w12s <-		0
    w22s <-   log(max (   var ( bi ) -  mean (sigma.bsVar ) , SMALL ))

     a1<-mean(Xvals [ , 1 ] )
     b1<-mean(Xvals [ , 2 ] )

 # use optim function for ML estimates of between trial level variance and covariance, A and B

            #cat( "\tSTART DERIV\n")
            ##try optim with gradient.
        answer <- optim(par=c(w11s,w12s,w22s, a1 , b1 ), fn=minusLogL.w , gr=gradient.w , method =c( "L-BFGS-B"),
       lower = c( Bignegative , -Inf , Bignegative, -Inf, -Inf ),
       control=list(trace=6, REPORT=1),
       hessian = TRUE ,
       Xvals = X ,  sigma.aSE = sigma.a , sigma.bSE = sigma.b )

		   #cat ("\tSTART NO DERIV\n")

     ##try optim without gradient.
      answer2 <- optim(par=c(w11s,w12s,w22s, a1 , b1 ), fn=minusLogL.w ,  method =c( "L-BFGS-B"),
      lower = c( Bignegative , -Inf , Bignegative, -Inf, -Inf ), control=list(trace=6, REPORT=1) , hessian = TRUE,
      Xvals = X ,  sigma.aSE = sigma.a , sigma.bSE = sigma.b )


	return (  list( answerDeriv = answer ,
					 answerNoDeriv = answer2
					)
			)
}

summary.w <- function ( rslt.w ,  Xvals  ,  sigma.aSE  , sigma.bSE) {

   		estimate.w <- rslt.w$par  ## on transformed scale (w)

   		deriv.w <- gradient.w (estimate.w ,  X  ,  sigma.a  , sigma.b)

    # convert choleski parameters to var-cov

	estimate.v <- WtoV ( estimate.w)

	names ( estimate.v ) <- c( "v11", "v12" , "v22" , "A", "B" )

	deriv.v <- gradient.v ( estimate.v ,  X ,  sigma.a , sigma.b)

         dvdw <- dVdW ( estimate.w )

	 	variancematrix.v <-   dvdw %*% solve(   rslt.w$hessian ) %*% t ( dvdw )  # for Sigma and means

	 	dimnames( variancematrix.v ) <- list ( c( "v11", "v12" , "v22" , "A", "B" ), c("v11", "v12" , "v22" , "A", "B" ))

		return( list( convergence = rslt.w$convergence , minusRlogL = rslt.w$value , estimate.w = estimate.w , deriv.w = deriv.w , estimate.v = estimate.v ,
				deriv.v = deriv.v , variancematrix.v = variancematrix.v ) )

}



summary.v <- function ( rslt.v ,  Xvals  ,  sigma.aSE  , sigma.bSE) {

   		estimate.v <- rslt.v$par


    # convert choleski parameters to var-cov

	estimate.w <- VtoW ( estimate.v)    ## on transformed scale (w)

	names ( estimate.w ) <- c( "w11", "z" , "w22" , "A", "B" )

	deriv.v <- gradient.v ( estimate.v ,  X ,  sigma.a  , sigma.b )

	      deriv.w <- gradient.w (estimate.w ,  X  ,  sigma.a  , sigma.b)
	 	variancematrix.v <-    solve(  rslt.v$hessian )  # for Sigma and means

	 	dimnames( variancematrix.v ) <- list ( c( "v11", "v12" , "v22" , "A", "B" ), c("v11", "v12" , "v22" , "A", "B" ))

		return( list( convergence = rslt.v$convergence , minusRlogL = rslt.v$value , estimate.w = estimate.w , deriv.w = deriv.w , estimate.v = estimate.v ,
				deriv.v = deriv.v , variancematrix.v=variancematrix.v ) )

}





  #input data:



ai <- c(3.7456, 5.7483, 2.54526, 3.55833, 3.9242, 3.2381, 5.28333,  3.88333  )
bi <- c( .3890, .31689, .69578, .61406, .36685,  .65524,  .4936,  .7148)


sigma.a <- c(.47772, .60274,  .84179,  .6894,  .77705, .50524,  .79539 , .76949  )

sigma.b <- c( .15988, .35127, .15348, .2261, .19522, .12466, .38804, .20622 )




  X <- cbind ( ai , bi )
    Xvals<-X
  #run function.

rslt.w <- MLsolution.w (X ,  sigma.a , sigma.b) #  minimize - RLogL wrt w , with and w/o derivs

rslt.v <- MLsolution.v (X,  sigma.a , sigma.b)  # minimize - RLogL wrt v , with and w/o derivs

sumryDeriv.w <- summary.w (rslt.w$answerDeriv , X ,  sigma.a , sigma.b )


sumryNoDeriv.w <- summary.w (rslt.w$answerNoDeriv ,  X ,  sigma.a , sigma.b)

sumryDeriv.v <- summary.v ( rslt.v$answerDeriv , X ,  sigma.a , sigma.b)
sumryNoDeriv.v <- summary.v ( rslt.v$answerNoDeriv , X ,  sigma.a , sigma.b)


sumryDeriv.v
sumryNoDeriv.v

sumryDeriv.w
sumryNoDeriv.w

##  Bootstrap method to obtain Confidence Intervals for the Mediated Effect

  library(MASS)
Sigma <- matrix(c(sumryDeriv.w$variancematrix.v[4,4], sumryDeriv.w$variancematrix.v[4,5]
 ,  sumryDeriv.w$variancematrix.v[4,5], sumryDeriv.w$variancematrix.v[5,5]),2,2)
#Sigma <- sumryDeriv.v$variancematrix.v$values[4:5]
Sigma
mu<-sumryDeriv.w$estimate.w[4:5]
mu



R <- 2000 # Use 2,000 bootstrap replicates



 MeanAtimesB.star<- rep(NA, R)
 AtimesBsigma.hat.star<- rep(NA, R)

for (r in 1:R) {
  X.star <- mvrnorm(n=nrow ( Xvals ), mu, Sigma)     #parametrically (bivariate normally) generate sample values of A and B.
  MeanAtimesB.star[r]<-mean(X.star[,1]*X.star[,2])  #compute simulated A*B, ABsimulate$AtimesBvalue
   #compute standard errors ($AtimesBstderrs) of the parametric samples
  AtimesBsigma.hat.star[r] <- sqrt(mean(X.star[,1])^2*var(X.star[,2])+mean(X.star[,2])^2*var(X.star[,1]))
  #Bsigma.hat.star[r] <- sqrt(var(X.star[,2]))
  #mu.star[r]<-mean(X.star)
}


 A<-mu[1]
 B<-mu[2]
# ABse<-0.3   #ABse can be calcualted from hessian matrix.
SeA<-Sigma[1,1]
SeB<-Sigma[2,2]
 ABse<-sqrt(A^2*SeB^2+B^2*SeA^2) #sobel



A*B-  ABse *  quantile (   (MeanAtimesB.star -  A*B) / AtimesBsigma.hat.star, probs = c(.975 , .025))

CL<- c(A*B-  ABse *  quantile (   (MeanAtimesB.star -  A*B) / AtimesBsigma.hat.star, probs = c(.975), names=FALSE) ,
   A*B-  ABse *  quantile (   (MeanAtimesB.star -  A*B) / AtimesBsigma.hat.star, probs = c(.025), names=FALSE)  )

label<-c("2.5%", "97.5%")


rbind(label, CL)








