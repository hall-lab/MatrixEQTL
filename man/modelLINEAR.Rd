\name{modelLINEAR}
\alias{modelLINEAR}
\docType{data}
\title{
Constant for \code{\link{Matrix_eQTL_engine}}.
}
\description{
Use of the constant as a parameter for \code{\link{Matrix_eQTL_engine}} 
to indicates that the effect of genotype on expression is modeled to be linear.
}
\usage{
modelLINEAR
}
\examples{
# Number of columns (samples)
n = 100;
# Number of covariates
nc = 10;

# Generate the standard deviation of the noise
noise.std = 0.1 + rnorm(n)^2;

# Generate the covariates
cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);

# Generate the vectors with single genotype and expression variables
snps.mat = cvrt.mat \%*\% rnorm(nc) + rnorm(n);
gene.mat = 1 + 0.5 * snps.mat + cvrt.mat \%*\% rnorm(nc) + rnorm(n) * noise.std;

# Create 3 SlicedData objects for the analysis
snps1 = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
gene1 = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrt1 = SlicedData$new( t(cvrt.mat) );

# Call the main analysis function
me = Matrix_eQTL_main(
	snps = snps1, 
	gene = gene1, 
	cvrt = cvrt1, 
	'Output_temp.txt', 
	pvOutputThreshold = 1, 
	useModel = modelLINEAR, 
	errorCovariance = diag(noise.std^2), 
	verbose = TRUE,
	pvalue.hist = TRUE );
# remove the output file
file.remove( 'Output_temp.txt' );

# Pull Matrix eQTL results - t-statistic and p-value
tstat = me$all$eqtls[ 1, 3 ];
pvalue = me$all$eqtls[ 1, 4 ];
rez = c( tstat = tstat, pvalue = pvalue);
# And compare to those from linear regression in R
{
	cat('\n\n Matrix eQTL: \n'); 
	print(rez);
	cat('\n R summary(lm()) output: \n')
	lmout = summary(lm(gene.mat ~ snps.mat + cvrt.mat,  weights = 1/noise.std^2))$coefficients[2, 3:4];
	print(lmout)
}
}