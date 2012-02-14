\name{modelLINEAR_CROSS}
\alias{modelLINEAR_CROSS}
\docType{data}
\title{
Constant for \code{\link{Matrix_eQTL_engine}}.
}
\description{
Use of the constant as a parameter for \code{\link{Matrix_eQTL_main}} 
to indicates that Matrix eQTL should include the interaction 
of SNP and last covariate in the model and
test its significance.
}
\usage{
modelLINEAR_CROSS
}
\examples{
# Number of columns (samples)
n = 25;
# Number of covariates
nc = 10;

# Generate the standard deviation of the noise
noise.std = 0.1 + rnorm(n)^2;

# Generate the covariates
cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);

# Generate the vectors with single genotype and expression variables
snps.mat = cvrt.mat \%*\% rnorm(nc) + rnorm(n);
gene.mat = 1 + 0.5 * snps.mat + cvrt.mat \%*\% rnorm(nc) + 
               snps.mat * cvrt.mat[,nc] + rnorm(n) * noise.std;

# Create 3 SlicedData objects for the analysis
snps1 = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
gene1 = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrt1 = SlicedData$new( t(cvrt.mat) );

# name of temporary output file
filename = tempfile();

# Call the main analysis function
me = Matrix_eQTL_main(
    snps = snps1, 
    gene = gene1, 
    cvrt = cvrt1, 
    filename, 
    pvOutputThreshold = 1, 
    useModel = modelLINEAR_CROSS, 
    errorCovariance = diag(noise.std^2), 
    verbose = TRUE,
    pvalue.hist = TRUE );
# remove the output file
unlink( filename );

# Pull Matrix eQTL results - t-statistic and p-value
tstat = me$all$eqtls[ 1, 3 ];
pvalue = me$all$eqtls[ 1, 4 ];
rez = c( tstat = tstat, pvalue = pvalue)
# And compare to those from linear regression in R
{
    cat(" Matrix eQTL: \n"); 
    print(rez);
    cat("\n R summary(lm()) output: \n")
    lmout = summary(lm(gene.mat ~ snps.mat + cvrt.mat + snps.mat*cvrt.mat[,nc],
        weights = 1/noise.std^2))$coefficients[ , 3:4];
    print(tail(lmout));
}
}