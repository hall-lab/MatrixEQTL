\name{Matrix_eQTL_main}
\alias{Matrix_eQTL_main}
\alias{Matrix_eQTL_engine}
\alias{Matrix_eQTL_engine.cis}
\title{
Perform eQTL analysis.
}
\description{
\code{Matrix_eQTL_main} function tests associations between every row of the \code{snps}
and every row of the \code{gene} using either linear or ANOVA model, as defined by \code{useModel}.
The testing procedure accounts for extra covariates in \code{cvrt}.
To account for heteroskedastic and/or correlated errors,
set the parameter \code{errorCovariance} to the error variance-covariance matrix.
Associations significant at \code{pvOutputThreshold} are saved to \code{output_file_name},
with corresponding test statistics, p-values, and estimated false discovery rate.

Matrix eQTL can perform separate analysis for local (cis) and distant (trans) eQTLs.
A gene-SNP pair is considered local if the distance between them is less than \code{cisDist}.
The genomic location of gene and SNPs is defined by variables \code{snpspos} and {genepos}.
To perform eQTL analysis without regard to gene/SNP location, set \code{pvOutputThreshold.cis = 0} (or do not set it).
To perform eQTL analysis for local gene-SNP pairs only, set \code{pvOutputThreshold = 0} and \code{pvOutputThreshold.cis > 0}.
To perform eQTL analysis with separate treatment of local and distant eQTLs,
set both thresholds to positive values. In this case the false discovery rate is calculated separately for these two groups of eQTLs.

Function \code{Matrix_eQTL_engine} is a wrapper for
\code{Matrix_eQTL_main} provided for compatibility with the previous versions of this package.

There are three linear regression models currently supported by Matrix eQTL (as defined by the \code{useModel} parameter). 
Set \code{useModel} to \code{modelLINEAR} to assume the effect of the genotype to be additive linear.
Use \code{modelANOVA} to treat genotype as a categorical variables and use ANOVA model.
The new special code \code{modelLINEAR_CROSS} adds an new term to the model,
equal to the product of genotype and the last covariate; the significance of this term is then tested using t-statistic.

Extra code examples are provided on the pages for \code{\link{modelLINEAR}}, \code{\link{modelANOVA}}, and \code{\link{modelLINEAR_CROSS}}.

}
\usage{
Matrix_eQTL_main(	
                   snps, 
                   gene, 
                   cvrt = SlicedData$new(), 
                   output_file_name = "", 
                   pvOutputThreshold = 1e-5,
                   useModel = modelLINEAR, 
                   errorCovariance = numeric(), 
                   verbose = TRUE, 
                   output_file_name.cis = "", 
                   pvOutputThreshold.cis = 0,
                   snpspos = NULL, 
                   genepos = NULL,
                   cisDist = 1e6,
                   pvalue.hist = FALSE)

Matrix_eQTL_engine(
                   snps, 
                   gene, 
                   cvrt = SlicedData$new(), 
                   output_file_name, 
                   pvOutputThreshold = 1e-5, 
                   useModel = modelLINEAR, 
                   errorCovariance = numeric(), 
                   verbose = TRUE,
                   pvalue.hist = FALSE)
}
\arguments{
  \item{snps}{
\code{\linkS4class{SlicedData}} object with genotype information. 
Can be real-valued for linear model and 
should take up 2 or 3 distinct values for ANOVA (see \code{useModel} parameter).
}
  \item{gene}{
\code{\linkS4class{SlicedData}} object with gene expression information. 
Should have columns matching those of \code{snps}.
}
  \item{cvrt}{
\code{\linkS4class{SlicedData}} object with additional covariates. 
Can be an empty \code{SlicedData} object in case of no covariates.
}
  \item{output_file_name}{
  character string with the name of the output file. 
Significant (all or distant) associations are saved to this file. 
Is the file with this name exists, it will be overwritten.
}
  \item{output_file_name.cis}{
  character string with the name of the output file. 
Significant local associations are saved to this file. 
Is the file with this name exists, it will be overwritten.
}
  \item{pvOutputThreshold}{
numeric. Only gene-SNP pairs significant at this level will be saved in \code{output_file_name}.	
}
  \item{pvOutputThreshold.cis}{ 
Same as \code{pvOutputThreshold}, but for cis-eQTLs. 
If not thresholds are positive, \code{pvOutputThreshold} determines cut-off for distant (trans) eQTLs.}
  \item{useModel}{
numeric. Set it to \code{modelLINEAR} to model the effect of the genotype to be additive linear,
or \code{modelANOVA} to treat genotype as a categorical variables and use ANOVA model.
The special code \code{modelLINEAR_CROSS} adds an interaction term to the model,
equal to the product of genotype and the last covariate; the significance of this term is then tested. 
}
  \item{errorCovariance}{
numeric. The error covariance matrix, if not multiple of identity matrix. 
Use this parameter to account for heteroscedastic and/or correlated errors.
}
  \item{verbose}{
logical. Set to \code{TRUE} to display detailed report on the progress.
}
  \item{snpspos}{
data.frame with information about SNP locations, with 3 columns - SNP name, chromosome, and position.
}
  \item{genepos}{
data.frame with information about transcript locations, with 4 columns - the name, chromosome, and positions of the left and right ends.
}
  \item{cisDist}{
numeric. SNP-gene pairs within this distance are considered local. The distance is measured from the nearest end of the gene.
}
  \item{pvalue.hist}{
	If \code{pvalue.hist} is not \code{FALSE}, the function returns sufficient information to plot the histogram(s) of (local/distant/all) p-values.
	Set \code{pvalue.hist} to a positive integer to build the histogram with \code{pvalue.hist} bins of equal size.
	Alternatively, a custom set of bin edges can be defined via \code{pvalue.hist}.
}
}
\details{
Note that the the columns of \code{gene}, \code{snps}, and \code{cvrt} must match.
If they do not match in the input files, use \code{ColumnSubsample} method to subset and/or reorder them.
}
\value{
The detected eQTLs are saved in \code{output_file_name} and/or \code{output_file_name.cis}.
The method also returns a list with a summary of the performed analysis.
}
\references{
For more information visit:
\url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\author{
Andrey Shabalin \email{shabalin@email.unc.edu}
}

\seealso{
For more information on the class of the first three arguments see \code{\linkS4class{SlicedData}}.
}



\examples{
# Number of columns (samples)
n = 100;

# Genetate single genotype variable
snps.mat = rnorm(n);

# Generate single expression variable
gene.mat = 0.5*snps.mat + rnorm(n);

# Create 3 SlicedData for the analysis
snps1 = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
gene1 = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrt1 = SlicedData$new();

# Call the main analysis function
me = Matrix_eQTL_main(
	snps = snps1, 
	gene = gene1, 
	cvrt = cvrt1, 
	'Output_temp.txt', 
	pvOutputThreshold = 1, 
	useModel = modelLINEAR, 
	errorCovariance = numeric(), 
	verbose = TRUE,
	pvalue.hist = TRUE );
# remove the output file
file.remove( 'Output_temp.txt' );

# Pull Matrix eQTL results - t-statistic and p-value
tstat = me$all$eqtls[ 1, 3 ];
pvalue = me$all$eqtls[ 1, 4 ];
rez = c( tstat = tstat, pvalue = pvalue)
# And compare to those from linear regression
{
	cat('\n\n Matrix eQTL: \n'); 
	print(rez);
	cat('\n R summary(lm()) output: \n')
	lmout = summary(lm(gene.mat ~ snps.mat))$coefficients[2, 3:4];
	print(lmout)
}
}




% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MatrixeQTL}
\keyword{MatrixEQTL}
\keyword{Shabalin}
\keyword{Matrix eQTL}
