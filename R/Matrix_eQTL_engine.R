# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

modelLINEAR = 321;
modelANOVA  = 123;

SlicedData <- setRefClass('SlicedData',
	fields = list( 
		dataEnv = 'environment',
		nSlices1 = 'numeric',
		rowNameSlices = 'list',
		columnNames = 'character',
		fileDelimiter = 'character',
		fileSkipColumns = 'numeric',
		fileSkipRows = 'numeric',
		fileSliceSize = 'numeric',
		fileOmitCharacters = 'character'
	),
	methods = list(
	initialize = function () {
		dataEnv <<- new.env(hash = TRUE, size = 29L);
		nSlices1 <<- 0;
		.self;
	},
	getSlice = function(sl) {
		return(get(paste(sl), dataEnv))	
	},
	setSlice = function(sl,value) {
		assign(paste(sl),value,dataEnv)
		if(nSlices1<sl) {
			nSlices1 <<- sl;
		}
	},	
	LoadFile = function(filename) {
		fid = file(description = filename, open = 'rt', blocking = FALSE, raw = FALSE)
		
		lines = readLines(con = fid, n = max(fileSkipRows,1), ok = TRUE, warn = TRUE)
		line1 = tail(lines,1);
		splt = strsplit(line1,split = fileDelimiter, fixed = TRUE);
		if(fileSkipRows > 0) {
			columnNames <<- splt[[1]][-(1:fileSkipColumns)];
		} else {
			seek(fid,0)
		}		
		
		rm(lines,line1,splt);
		
		rowNameSlices <<- vector('list',15);

		curSliceId = 0;
		repeat
		{
			# preallocate names and data
			if(length(rowNameSlices) < curSliceId) {
				rowNameSlices[[2*curSliceId]] <<- NULL;
			}
			curSliceId = curSliceId + 1;
			
			# read sliceSize rows
			rowtag = vector('character',fileSliceSize);
			rowvals = vector('list',fileSliceSize);
			for(i in 1:fileSliceSize) {
				temp = '';
				if(fileSkipColumns > 0) {
					temp = scan(file = fid, what = character(), n = fileSkipColumns, quiet = TRUE,sep = fileDelimiter);
				}
				rowtag[i] = paste(temp,collapse=' ');
				rowvals[[i]] = scan(file = fid, what = double(), nlines = 1, quiet = TRUE, sep = fileDelimiter, na.strings = fileOmitCharacters);
				if(length(rowvals[[i]])==0) {
					if(i==1) {
						rowtag = matrix('',0,0);
						rowvals = character(0);
					} else 	{
						rowtag = rowtag[1:(i-1)];
						rowvals = rowvals[1:(i-1)];
					}
					break;			
				}
			}
			if(length(rowtag)==0) {
				curSliceId = curSliceId - 1;
				break;
			}
			rowNameSlices[[curSliceId]] <<- rowtag;
			data = c(rowvals, recursive=TRUE);
			dim(data) = c(length(rowvals[[1]]), length(rowvals));
			setSlice(curSliceId, t(data));
			if(length(rowtag) < fileSliceSize) {
				break;
			}
			cat(sprintf('Rows read: %d\n',curSliceId*fileSliceSize));
			flush.console()
		}
		close(fid)
		if(fileSkipRows == 0) {
			columnNames <<- paste('Col_',(1:nCols()),sep='');
		}
		if(fileSkipColumns == 0) {
			cnt = 0;
			for(sl in 1:nSlices()) {
				nr = nrow(getSlice(sl));
				rowNameSlices[[sl]] <<- paste('Row_',cnt + (1:nr),sep='');
				cnt = cnt + nr;
			}
		}		
		rowNameSlices <<- rowNameSlices[1:curSliceId];
		cat(sprintf('Rows read: %d done.\n',nRows()));
	},
	nRows = function() {
		if(nSlices()==0) {
			return(0);
		} else {
			s = 0;
			for(sl in 1:nSlices()) {
				s = s + nrow(getSlice(sl));
			}
			return( s )
		}
	},
	nCols = function() {
		if(nSlices()==0) {
			return(0);
		} else {
			return( ncol(getSlice(1)) )
		}
	},
	nSlices = function() {
		return( nSlices1 )
	},
	Clear = function() {
		# dataSlices1 <<- list();
		for(sl in 1:nSlices()) {
			# setSlice(sl,numeric())
			rm(list = paste(sl), envir = dataEnv)
		}
		nSlices1 <<- 0;
		rowNameSlices <<- list();
		columnNames <<- character(0);
	},
	SetNanRowMean = function() {
		for(sl in 1:nSlices()) {
			slice = getSlice(sl);
			if(any(!is.finite(slice))) {
				rowmean = rowMeans(slice, na.rm = TRUE);
				rowmean[!is.finite(rowmean)] = 0;
				for(j in 1:ncol(slice)) {
					where1 = is.na(slice[,j]);
					slice[where1,j] = rowmean[where1];
				}
				setSlice(sl,slice);
			}
		}
	},
	RowStandardizeCentered = function(){
		for(sl in 1:nSlices()) {
			slice = getSlice(sl);
			div = sqrt(rowSums(slice^2));
			div[div==0] = 1;
			setSlice(sl, slice/div);
		}
	},
	CombineInOneSlice = function(){
		nc = nCols();
		nr = nRows();
		zzz = list();
		for(sl in 1:nSlices()) {
			zzz[[sl]] = t(getSlice(sl));
			setSlice(sl, numeric());
		}
		newData = c(zzz,recursive=TRUE);
		rm(zzz);
		dim(newData) = c(nc, nr);
		newData = t(newData);
		
		nSlices1 <<- 1;
		setSlice(1, newData);
		rm(newData);
		
		newrowNameSlices = c(rowNameSlices,recursive=TRUE);
		rowNameSlices <<- vector('list',1);
		rowNameSlices[[1]] <<- newrowNameSlices;
	},
	Clone = function(){
		clone = SlicedData$new();
		for(sl in 1:nSlices()) {
			clone$setSlice(sl,getSlice(sl));
		}
		clone$rowNameSlices = rowNameSlices;
		clone$columnNames = columnNames;
		clone$fileDelimiter = fileDelimiter;
		clone$fileSkipColumns = fileSkipColumns;
		clone$fileSkipRows = fileSkipRows;
		clone$fileSliceSize = fileSliceSize;
		clone$fileOmitCharacters = fileOmitCharacters;
		return(clone);		
	},
	RowMatrixMultiply = function(multiplier) {
		for(sl in 1:nSlices()) {
			setSlice(sl, getSlice(sl) %*% multiplier);
		}
	},
	ColumnSubsample = function(subset) {
		for(sl in 1:nSlices()) {
			setSlice(sl, getSlice(sl)[ ,subset]);
		}
		columnNames <<- columnNames[subset];
	},
	RowRemoveZeroEps = function(){
		for(sl in 1:nSlices()) {
			slice = getSlice(sl);
			amean = rowMeans(abs(slice));
			remove = (amean < .Machine$double.eps*nCols());
			if(any(remove)) {
				rowNameSlices[[sl]] <<- rowNameSlices[[sl]][!remove];
				setSlice(sl, slice[!remove, ]);
			}
		}
	}
	)
)

Matrix_eQTL_engine = function(snps, gene, cvrt = SlicedData$new(), output_file_name, pvOutputThreshold = 1e-5,useModel = modelLINEAR, errorCovariance = numeric(), verbose=FALSE){
	
	lastTime = 0;
	status <- function(text) {
		gc();
		newTime = proc.time()[3];
		if(lastTime != 0) {
			cat('Task finished in ', newTime-lastTime, ' seconds\n');
		}
		cat(text,'\n');
		lastTime <<- newTime;
		unused = flush.console();
	}
	if(!verbose) {status = function(text){}}

	# Check dimensions
	status('Checking input data dimensions');
	if(snps$nCols() != gene$nCols()) {
		stop('Different number of samples in the genotype and gene expression files');
	}
	
	if(cvrt$nRows()>0) {
		if(snps$nCols() != cvrt$nCols()) {
			stop('Wrong number of samples in the file with covariates');
		}	
		cvrt$SetNanRowMean();
		cvrt$CombineInOneSlice();
	}

	if(length(errorCovariance)>0) {
		status('Processing the errorCovariance matrix');
		errorCovariance = as.matrix(errorCovariance);
		if(nrow(errorCovariance) != ncol(errorCovariance)) {
			stop('The covariance matrix is not square');
		}	
		if(nrow(errorCovariance) != snps$nCols()) {
			stop('The covariance matrix size does not match the data');
		}
		# test for symmetry
		if(!all(errorCovariance==t(errorCovariance))) {
			stop('The covariance matrix is not symmetric');
		}
		eig = eigen(errorCovariance, symmetric = TRUE)
		d = eig$values;
		v = eig$vectors;
		#  errorCovariance == v %*% diag(d) %*% t(v)
		#  errorCovariance^0.5 == v*sqrt(d)*v'
		#  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v'
		if(any(d<=0)) {
			stop('The covariance matrix is not positive definite');
		}
		correctionMatrix = v %*% diag(1./sqrt(d)) %*% t(v);
		rm(eig,v,d)
	}
	
	# Add constant as a covariate
	if(cvrt$nRows()>0) {
		cvrt$setSlice(1, rbind(matrix(1,1,snps$nCols()),cvrt$getSlice(1)));
	} else {
		cvrt$setSlice(1, matrix(1,1,snps$nCols()));
	}
	
	# Correct for the error covariance structure
	if(length(errorCovariance)>0) {
		status('Rotating cvrt based on the errorCovariance matrix');
		cvrt$RowMatrixMultiply(correctionMatrix);
	}
	
	OrthonormalizeRows = function(x){
		for(i in 1:nrow(x)) {
			if(i > 1) {
				for(j in 1:(i-1)) {
					x[i, ] = x[i, ] - crossprod(x[i, ],x[j, ]) * x[j, ];
				}
		}
			div = sqrt(crossprod(x[i, ]));
			if(div < .Machine$double.eps*ncol(x)) {
					stop('Colinear or zero covariates detected');
			}
			x[i, ] = x[i, ] / div;
		}
		return(x)
	}	

	# Orthonormalize covariates
	status('Orthonormalizing covariates');
	cvrt$setSlice(1, OrthonormalizeRows(cvrt$getSlice(1)));
	
	if(useModel == modelLINEAR) {
		status('Imputing missing genotype for linear model (with average)');
		snps$SetNanRowMean();
		snps_list = list();
		snps_list[[1]] = snps$Clone();
		snps$Clear();
	} else if(useModel == modelANOVA ) {
		# split into 2 dummy variables
		status('Splitting genotype variable into dummy variables (for ANOVA)');
		snps_list = list();
		snps_list[[1]] = snps$Clone();
		snps_list[[2]] = snps$Clone();
		# mostrepeated <- function(x) as(names(which.max(table(x))), mode(x))
		for(sl in 1:snps$nSlices()){
			slice = snps$getSlice(sl);
			
			uniq = unique(c(slice));
			uniq = uniq[!is.na(uniq)];
			
			if(length(uniq)>3) {
				stop('More than three genotype categories');
			}
			if(length(uniq) == 2) {
				uniq = c(uniq, min(uniq)-1);
			}
			if(length(uniq) == 1) {
				uniq = c(uniq, min(uniq)-1,min(uniq)-2);
			}
			
			
			freq = matrix(0,nrow(slice),length(uniq));
			for(i in 1:length(uniq)) {
				freq[ ,i] = rowSums(slice==uniq[i],na.rm = TRUE);
			}
			
			md = apply(freq,1,which.max);
			freq[cbind(1:nrow(slice),md)] = -1;
			
			md = apply(freq,1,which.max); # min(freq[cbind(1:nrow(slice),md)] - rowSums(select,na.rm = TRUE ))
			new_slice1 = (slice == uniq[md]);
			new_slice1[is.na(new_slice1)] = 0;
			freq[cbind(1:nrow(slice),md)] = -1;
					
			md = apply(freq,1,which.max);
			new_slice2 = (slice == uniq[md]);
			new_slice2[is.na(new_slice2)] = 0;

			snps_list[[1]]$setSlice(sl, new_slice1);
			snps_list[[2]]$setSlice(sl, new_slice2);
			snps$setSlice(sl, numeric());
		}
	} else stop('Unknown value of useModel');
	
	
	# Correct for the error covariance structure
	if(length(errorCovariance)>0) {
		status('Rotating snps/dummies based on the errorCovariance matrix');
		for(d in 1:length(snps_list)) {
			snps_list[[d]]$RowMatrixMultiply(correctionMatrix);
		}
	}
	
	# Orthogonolize the SNPs w.r.t. covariates
	status('Orthogonolizing snps/dummies w.r.t. covariates');
	for(d in 1:length(snps_list)) {
		for(sl in 1:snps_list[[d]]$nSlices()) {
			slice = snps_list[[d]]$getSlice(sl);
			slice = slice - tcrossprod(slice,cvrt$getSlice(1)) %*% cvrt$getSlice(1);
			snps_list[[d]]$setSlice(sl, slice);
		}
	}
	
	# Orthogonolize the dummies (SNPs) w.r.t. each other
	status('Orthonormalizing snps/dummies');
	snps_list[[1]]$RowStandardizeCentered();
	# crossprod(snps_list[[1]]$dataSlices[[1]][3,])
	if(length(snps_list) > 1) {
		d1 = 2; d2 = 1;
		for(sl in 1:snps_list[[d1]]$nSlices()) {
			slice1 = snps_list[[d1]]$getSlice(sl);
			slice2 = snps_list[[d2]]$getSlice(sl);
			new_slice1 = slice1 - rowSums(slice1*slice2)*slice2;
			snps_list[[d1]]$setSlice(sl, new_slice1);
			# snps_list[[d1]]$dataSlices[[sl]] = snps_list[[d1]]$dataSlices[[sl]] -
			#		rowSums(snps_list[[d1]]$dataSlices[[sl]]*snps_list[[d2]]$dataSlices[[sl]]) * 
			#		snps_list[[d2]]$dataSlices[[sl]];
		}
		snps_list[[2]]$RowStandardizeCentered();
	}
	
	# snps_list[[1]]$dataSlices[[1]]
	
	# Impute gene expression
	status('Imputing missing expression');
	gene$SetNanRowMean();
	# Correct for the error covariance structure
	if(length(errorCovariance)>0) {
		status('Rotating expression based on the errorCovariance matrix');
		gene$RowMatrixMultiply(correctionMatrix);
	}

	gene$RowStandardizeCentered();
	
	# Orthogonolize expression w.r.t. covariates
	status('Orthogonolizing expression w.r.t. covariates');
	for(sl in 1:gene$nSlices()) {
		slice = gene$getSlice(sl);
		slice = slice - tcrossprod(slice,cvrt$getSlice(1)) %*% cvrt$getSlice(1);
		gene$setSlice(sl, slice);
	}
	gene$RowRemoveZeroEps();

	
	status('Standardizing expression');
	gene$RowStandardizeCentered();
	status('Running garbage collector');
	gc();
	
	##
	nSamples = snps_list[[1]]$nCols();
	nGenes = gene$nRows();
	nSnps  = snps_list[[1]]$nRows();
	nCov = cvrt$nRows();
	nVarTested = length(snps_list);
	# dfNull = nSamples - nCov;
	dfFull = nSamples - nCov - nVarTested;

	snameSlices = snps_list[[1]]$rowNameSlices;
	gnameSlices = gene$rowNameSlices;

	if(useModel == modelLINEAR) {
		if(pvOutputThreshold>=1) {
			rThresh = 0;
		} else {
			tThresh = -qt(pvOutputThreshold/2,dfFull);
			rThresh = sqrt(  tThresh^2 / (dfFull + tThresh^2)  );
		}
	} else if(useModel == modelANOVA) {
		if(pvOutputThreshold>=1) {
			r2Thresh = 0;
		} else {
			fThresh = qf(1-pvOutputThreshold, nVarTested,dfFull);
			r2Thresh = fThresh * nVarTested / (dfFull + fThresh*nVarTested);
		}
	}
	
	totalCount = nGenes*nSnps;
	curCount = 0;
	dumpCount = 0;
	
	status('Performing eQTL analysis');
	tic = proc.time()[3];
	fid = file(description = output_file_name, open = 'wt', blocking = FALSE, raw = FALSE)
	for(sc in 1:snps_list[[1]]$nSlices()) {
		for(gc in 1:gene$nSlices()) {
			if(useModel == modelLINEAR) {
				cor = (tcrossprod(snps_list[[1]]$getSlice(sc),gene$getSlice(gc)));
				select = which(abs(cor) >= rThresh, arr.ind=TRUE );
				rsub  = cor[select];
				test = rsub *sqrt( dfFull / (1-rsub ^2));
				pv = pt(-abs(test),dfFull)*2;
				rm(cor);
			} else if(useModel == modelANOVA) {
				r2 = tcrossprod(snps_list[[1]]$getSlice(sc),gene$getSlice(gc))^2;
				for(d in 2:nVarTested) {
					r2 = r2 + tcrossprod(snps_list[[d]]$getSlice(sc),gene$getSlice(gc))^2;
				}
				select = which(r2 >= r2Thresh, arr.ind=TRUE );
				rsub = r2[select];
				test = rsub/(1-rsub) * (dfFull/nVarTested);
				pv = pf(1./test, dfFull, nVarTested);					
				rm(r2);
			}
			snames = snameSlices[[sc]][select[,1]];
			gnames = gnameSlices[[gc]][select[,2]];
			dump2 = data.frame(snames,gnames,pv, test, row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = FALSE)
			write.table(dump2, file = fid, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE);
			
			curCount = curCount + nrow(snps_list[[1]]$getSlice(sc))*nrow(gene$getSlice(gc));
			dumpCount = dumpCount + length(rsub);
			cat( floor(curCount/totalCount*1000)/10, '% done, ',dumpCount,' eQTLs found.\n');
			flush.console();
		}
		# Remove data no longer needed
		for(d in 1:nVarTested) {
			snps_list[[d]]$setSlice(sc,numeric());
		}		
	}
	close(fid);
	toc = proc.time()[3];
	cat('eQTL time: ',toc-tic,' sec\n');
}