\name{Matrix_eQTL_engine}
\alias{Matrix_eQTL_engine}
\alias{Matrix_eQTL_engine_cis}
\title{
Perform eQTL analysis.
}
\description{
\code{Matrix_eQTL_engine} tests for association between every row of the \code{snps}
and every row of the \code{gene}
using either linear or ANOVA model, as defined by \code{useModel}.
The testing procedure accounts for extra covariates in \code{cvrt}.
To account for heteroskedastic and/or correlated errors,
set the parameter \code{errorCovariance} to the error covariance matrix.
Associations significant at \code{pvOutputThreshold} are saved to \code{output_file_name},
with corresponding test statistics, p-values, and estimates of false discovery rate.




}
\usage{
Matrix_eQTL_engine(snps, 
                   gene, 
                   cvrt = SlicedData$new(), 
                   output_file_name, 
                   pvOutputThreshold = 1e-05, 
                   useModel = modelLINEAR, 
                   errorCovariance = numeric(), 
                   verbose = FALSE)
Matrix_eQTL_engine_cis(snps, 
                       gene, 
                       cvrt = SlicedData$new(), 
                       output_file_name, 
                       pvOutputThreshold_cis = 1e-3,
                       pvOutputThreshold_tra = 1e-6,
                       useModel = modelLINEAR, 
                       errorCovariance = numeric(), 
                       verbose=FALSE, 
                       snpspos, 
                       genepos,
                       cisDist = 1e6 )
}
\arguments{
  \item{snps}{
\code{SlicedData} object with genotype information. 
Can be real-valued for linear model and 
should take up 2 or 3 distinct values for ANOVA (see \code{useModel} parameter).
}
  \item{gene}{
\code{SlicedData} object with gene expression information. 
Should have columns matching those of \code{snps}.
}
  \item{cvrt}{
\code{SlicedData} object with additional covariates. 
Can be an empty \code{SlicedData} object in case of no covariates.
}
  \item{output_file_name}{
  connection or a character string with name of the output file. 
Significant associations will be saved to this file. 
Is the file with this name exists, it will be overwritten.
}
  \item{pvOutputThreshold}{
numeric. Only gene-SNP pairs significant at this level will be saved in \code{output_file_name}.	
}
  \item{pvOutputThreshold_cis}{ Same as \code{pvOutputThreshold}, but for cis-eQTLs.}
  \item{pvOutputThreshold_tra}{ Same as \code{pvOutputThreshold}, but for trans-eQTLs.}
  \item{useModel}{
numeric. Set it to \code{modelLINEAR} to use the linear model 
(additive effect of the SNP on the gene), or
to \code{modelANOVA} to treat genotype as a categorical variables.
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
numeric. SNP-gene pairs within this distance will be considered 'cis-'. The distance is measured from the nearest end of the gene.
}

}
\details{
Note that the the columns of \code{gene}, \code{snps}, and \code{cvrt} must match.
If they do not match in the input files, use \code{ColumnSubsample} method to subset and/or reorder them.
}
\value{
	The method does not return any values.
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
## Settings

# Linear model to use, modelANOVA or modelLINEAR
useModel = modelLINEAR; # modelANOVA or modelLINEAR

# Genotype file name
SNP_file_name = 'Sample_Data/SNP.txt';

# Gene expression file name
expression_file_name = 'Sample_Data/GE.txt';

# Covariates file name
# Set to character() for no covariates
covariates_file_name = 'Sample_Data/Covariates.txt';

# Output file name
output_file_name = 'Sample_Data/eQTL_results_R.txt';

# Only associations significant at this level will be output
pvOutputThreshold = 1e-2;

# Error covariance matrix
# Set to character() for identity.
errorCovariance = character();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 10000; # read file in pieces of 10,000 rows
\dontrun{
snps$LoadFile(SNP_file_name);
}

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t"; # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 10000; # read file in pieces of 10,000 rows
\dontrun{
gene$LoadFile(expression_file_name);
}
## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"; # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
cvrt$fileSliceSize = snps$nCols()+1; # read file in one piece
if(length(covariates_file_name)>0) {
\dontrun{
	cvrt$LoadFile(covariates_file_name);
}
}

## Run the analysis

\dontrun{
Matrix_eQTL_engine(snps,
                   gene,
                   cvrt,
                   output_file_name,
                   pvOutputThreshold,
                   useModel, 
                   errorCovariance, 
                   verbose=TRUE);
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MatrixeQTL}
\keyword{MatrixEQTL}
\keyword{Shabalin}
\keyword{Matrix eQTL}
