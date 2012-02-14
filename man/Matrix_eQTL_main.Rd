\name{Matrix_eQTL_main}
\alias{Matrix_eQTL_main}
\alias{Matrix_eQTL_engine}
\title{
Perform eQTL analysis.
}
\description{
\code{Matrix_eQTL_engine} function tests associations between every row of the \code{snps} dataset
and every row of the \code{gene} dataset using either linear additive or ANOVA model, as defined by \code{useModel} parameter.
The testing procedure accounts for extra covariates in \code{cvrt} dataset.
The constant is always included in the model and should not be explicitly specified in \code{cvrt}.
The parameter \code{errorCovariance} can be set to the error variance-covariance matrix
to account for heteroskedastic and/or correlated errors.
Associations significant at \code{pvOutputThreshold} are saved to \code{output_file_name},
with corresponding test statistics, p-values, and estimated false discovery rate.

Matrix eQTL can perform separate analysis for local (cis) and distant (trans) eQTLs.
Extra parameters for such analysis are available in \code{Matrix_eQTL_main}.
A gene-SNP pair is considered local if the distance between them is less than \code{cisDist}.
The genomic location of genes and SNPs is defined by variables \code{snpspos} and {genepos}.
The type of analysis is defined by p-value thresholds \code{pvOutputThreshold} and \code{pvOutputThreshold.cis}:
\enumerate{
\item Set \code{pvOutputThreshold > 0} and \code{pvOutputThreshold.cis = 0}, or use \code{Matrix_eQTL_engine} to perform eQTL analysis without regard to gene/SNP location. Associations significant at \code{pvOutputThreshold} level will be recorded in \code{output_file_name}.
\item Set \code{pvOutputThreshold = 0} and \code{pvOutputThreshold.cis > 0} to perform eQTL analysis for local gene-SNP pairs only. Local associations significant at \code{pvOutputThreshold.cis} level will be recorded in \code{output_file_name.cis}.
\item Set \code{pvOutputThreshold > 0} and \code{pvOutputThreshold.cis > 0} to perform eQTL analysis with separate p-value thresholds for local and distant eQTLs. Distant and local associations significant at corresponding thresholds are recorded in \code{output_file_name} and \code{output_file_name.cis} respectively.
In this case the false discovery rate is calculated separately for these two groups of eQTLs.
}
Note that \code{Matrix_eQTL_engine} is a wrapper for \code{Matrix_eQTL_main} provided for easier eQTL analysis without regard to gene/SNP location and for compatibility with the previous versions of this package.

There are three linear regression models currently supported by Matrix eQTL as defined by the \code{useModel} parameter:

       \enumerate{
         \item Set \code{useModel} to \code{\link{modelLINEAR}} to model the effect of the genotype as additive linear and test for its significance using t-statistic.
         \item Use \code{\link{modelANOVA}} to treat genotype as a categorical variables and use ANOVA model and test for its significance using F-test. Note that no more than three distinct values per genotype variable is supported.
         \item The new special code \code{\link{modelLINEAR_CROSS}} adds a new term to the model
equal to the product of genotype and the last covariate; the significance of this term is then tested using t-statistic.
       }

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
should take up to 3 distinct values for ANOVA (see \code{useModel} parameter).
}
  \item{gene}{
\code{\linkS4class{SlicedData}} object with gene expression information. 
Should have columns matching those of \code{snps}.
}
  \item{cvrt}{
\code{\linkS4class{SlicedData}} object with additional covariates. 
Can be an empty \code{SlicedData} object in case of no covariates.
The columns must match those in \code{snps} and \code{gene}.
}
  \item{output_file_name}{
  character string with the name of the output file. 
Significant (all or distant) associations are saved to this file.
If the file with this name exists, it will be overwritten.
}
  \item{output_file_name.cis}{
  character string with the name of the output file. 
Significant local associations are saved to this file. 
If the file with this name exists, it will be overwritten.
}
  \item{pvOutputThreshold}{
numeric. Only gene-SNP pairs significant at this level will be saved in \code{output_file_name}.
}
  \item{pvOutputThreshold.cis}{ 
Same as \code{pvOutputThreshold}, but for local eQTLs.
If both thresholds are positive, \code{pvOutputThreshold} determines cut-off for distant (trans) eQTLs.}
  \item{useModel}{
numeric. Can be \code{modelLINEAR}, \code{modelANOVA}, or \code{modelLINEAR_CROSS}.
See the section above for description.
}
  \item{errorCovariance}{
numeric. The error covariance matrix. Use \code{numeric()} for homoskedastic independent errors. 
}
  \item{verbose}{
logical. Set to \code{TRUE} to display detailed report on the progress.
}
  \item{snpspos}{
\code{data.frame} object with information about SNP locations, must have 3 columns - SNP name, chromosome, and position.
}
  \item{genepos}{
\code{data.frame} with information about transcript locations, must have 4 columns - the name, chromosome, and positions of the left and right ends.
}
  \item{cisDist}{
numeric. SNP-gene pairs within this distance are considered local. The distance is measured from the nearest end of the gene.
}
  \item{pvalue.hist}{
	This parameter defines how the distribution of (all/local/distant) p-values is recorded.
	If \code{pvalue.hist} is \code{FALSE}, the information is not recorded and thus the analysis is performed faster.
	Set \code{pvalue.hist = "qqplot"} to record information sufficient to create a Q-Q plot of the p-values (use \link[=plot.MatrixEQTL]{plot} to create the plot).
	To record information for a histogram set \code{pvalue.hist} to the desired number of bins of equal size.
	Alternatively, a custom set of bin edges can be submitted via \code{pvalue.hist}.
}
}
\details{
Note that the columns of \code{gene}, \code{snps}, and \code{cvrt} must match.
If they do not match in the input files, use \code{ColumnSubsample} method to subset and/or reorder them.
}
\value{
The detected eQTLs are saved in \code{output_file_name} and/or \code{output_file_name.cis}.
The method also returns a list with a summary of the performed analysis.
\item{param}{Keeps all input parameters.}
\item{time.in.sec}{Time required for the performed analysis (in seconds).}
\item{all}{Information about detected eQTLs for the analysis not using genomic locations.}
\item{cis}{Information about detected local eQTLs.}
\item{trans}{Information about detected distant eQTLs.}
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

# Create 3 SlicedData objects for the analysis
snps1 = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
gene1 = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrt1 = SlicedData$new();

# name of temporary output file
filename = tempfile();

# Call the main analysis function
me = Matrix_eQTL_main(
    snps = snps1, 
    gene = gene1, 
    cvrt = cvrt1, 
    filename, 
    pvOutputThreshold = 1, 
    useModel = modelLINEAR, 
    errorCovariance = numeric(), 
    verbose = TRUE,
    pvalue.hist = TRUE );
# remove the output file
unlink( filename );

# Pull Matrix eQTL results - t-statistic and p-value
tstat = me$all$eqtls[ 1, 3 ];
pvalue = me$all$eqtls[ 1, 4 ];
rez = c( tstat = tstat, pvalue = pvalue)
# And compare to those from linear regression
{
    cat("\n\n Matrix eQTL: \n"); 
    print(rez);
    cat("\n R summary(lm()) output: \n")
    lmout = summary(lm(gene.mat ~ snps.mat))$coefficients[2, 3:4];
    print(lmout)
}
}




% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{MatrixEQTL}
