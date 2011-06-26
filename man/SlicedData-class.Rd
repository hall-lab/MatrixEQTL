\name{SlicedData-class}
\Rdversion{1.1}
\docType{class}
\alias{SlicedData-class}
\alias{SlicedData}

\title{Class \code{SlicedData} for storing large matrices}
\description{
This class is created for fast and memory efficient manipulations
with large datasets presented in matrix form.
It is used to load, store, and manipulate large datasets, 
e.g. genotype and gene expression matrices.
When a dataset is loaded, it is sliced in blocks of 10,000 rows (default size).
This allows to impute, standardize, and perform other operations with the data with minimal memory overhead.
}
\section{Extends}{

\code{SlicedData} is a reference classes (\code{\linkS4class{envRefClass}}).
Its methods can change the values of the fields of the class.
}
\references{
For more information visit:
\url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\author{
Andrey Shabalin \email{shabalin@email.unc.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
This class is used to load data for eQTL analysis by \code{\link{Matrix_eQTL_engine}}.
}
\note{
To run this example please obtain the sample data from
\url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\examples{

showClass("SlicedData")

SNP_file_name = 'Sample_Data/SNP.txt';
snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 10000; # read file in pieces of 10,000 rows
\dontrun{
snps$LoadFile(SNP_file_name);}

## Display the dimensions and the number of slices

snps$nRows()
snps$nCols()
snps$nSlices()

}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{dataEnv}:}{\code{environment}. Stores the slices of the data matrix.}
    \item{\code{nSlices1}:}{\code{numeric}. Number of slices. }
    \item{\code{rowNameSlices}:}{\code{list}. Slices of row names. }
    \item{\code{columnNames}:}{\code{character}. Column names. }

    \item{\code{fileDelimiter}:}{\code{character}. Delimiter separating values in the input file. }
    \item{\code{fileSkipColumns}:}{\code{numeric}. Number of columns with row labels in the input file. }
    \item{\code{fileSkipRows}:}{\code{numeric}. Number of rows with column labels in the input file. }
    \item{\code{fileSliceSize}:}{\code{numeric}. Maximum number of rows in a slice. }
    \item{\code{fileOmitCharacters}:}{\code{character}. Missing value (NaN) representation in the input file. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{nSlices()}:}{ Return the number of slices. }
    \item{\code{nCols()}:}{ Return the number of columns in the matrix. }
    \item{\code{nRows()}:}{ Return the number of rows in the matrix. }
    \item{\code{Clear()}:}{ Clear the object. Removes the data slices and the row and column names. }
    \item{\code{Clone()}:}{ Make a copy of the object. Changes to the copy do not affect the source object.}

    \item{\code{LoadFile(filename)}:}{ Load a matrix file. Set \code{file*} fields to specify the format. }

    \item{\code{getSlice(sl)}:}{ Get \code{sl}-th slice of the matrix. }
    \item{\code{setSlice(sl, value)}:}{ Set \code{sl}-th slice of the matrix. }

    \item{\code{ColumnSubsample(subset)}:}{ Keep only the columns in the \code{subset}. 
Acts as \code{M = M[ ,subset]}. }
    \item{\code{RowMatrixMultiply(multiplier)}:}{ Multiply each row by the \code{multiplier}.
Acts as \code{M = M*multiplier}. }
    \item{\code{CombineInOneSlice()}:}{ Combines all slices into one. The whole matrix can then be obtained via \code{$getSlice(1)}. }

    \item{\code{RowStandardizeCentered()}:}{ Set the mean of each row to zero and the sum of squares to one. }
    \item{\code{SetNanRowMean()}:}{ Impute rows with row mean. Rows full of NaN's are imputed with zeros. }
  }
}
