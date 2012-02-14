\name{modelANOVA}
\alias{modelANOVA}
\docType{data}
\title{
Constant for \code{\link{Matrix_eQTL_engine}}.
}
\description{
Use of the constant as a parameter for \code{\link{Matrix_eQTL_engine}} 
to indicates that the genotype should be treated as a categorical variable.
}
\usage{
modelANOVA
}
\examples{# Number of columns (samples)
n = 100;
# Number of covariates
nc = 10;

# Generate the standard deviation of the noise
noise.std = 0.1 + rnorm(n)^2;

# Generate the covariates
cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);

# Generate the vectors with single genotype and expression variables
snps.mat = floor(runif(n, min = 0,max = 3));
gene.mat = 1 + (snps.mat==1) + cvrt.mat \%*\% rnorm(nc) + rnorm(n) * noise.std;

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
	useModel = modelANOVA, 
	errorCovariance = diag(noise.std^2), 
	verbose = TRUE,
	pvalue.hist = TRUE );
# remove the output file
file.remove( 'Output_temp.txt' );

# Pull Matrix eQTL results - t-statistic and p-value
fstat = me$all$eqtls[ 1, 3 ];
pvalue = me$all$eqtls[ 1, 4 ];
rez = c( Fstat = fstat, pvalue = pvalue)
# And compare to those from ANOVA in R
{
	cat('\n\n Matrix eQTL: \n'); 
	print(rez);
	cat('\n R anova(lm()) output: \n')
	lmout = anova( lm( gene.mat ~ cvrt.mat + factor(snps.mat), weights = 1/noise.std^2 ) );
	print( c(Fstat = lmout$F[2], pvalue = lmout$"Pr(>F)"[2]) )
}
}