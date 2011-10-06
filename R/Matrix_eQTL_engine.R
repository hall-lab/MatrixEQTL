# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

modelLINEAR = 117348;
modelANOVA  = 47074;

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
				rowtag[i] = temp[1];#paste(temp,collapse=' ');
				rowvals[[i]] = scan(file = fid, what = double(), nlines = 1, quiet = TRUE, sep = fileDelimiter, na.strings = fileOmitCharacters);
				if(length(rowvals[[i]])==0) {
					if(i==1) {
						rowtag = matrix(0,0,0);
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
			if(nSlices()>0)
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
		if(nSlices()>0)
			for(sl in 1:nSlices()) {
				# setSlice(sl,numeric())
				rm(list = paste(sl), envir = dataEnv)
			}
		nSlices1 <<- 0;
		rowNameSlices <<- list();
		columnNames <<- character(0);
	},
	SetNanRowMean = function() {
		if(nSlices()==0)
			return;
		if(nCols()==0)
			return;
		for(sl in 1:nSlices()) {
			slice = getSlice(sl);
			if(any(is.na(slice))) {
				rowmean = rowMeans(slice, na.rm = TRUE);
				rowmean[is.na(rowmean)] = 0;
				for(j in which(!complete.cases(slice))) {
					where1 = is.na(slice[j, ]);
					slice[j,where1] = rowmean[j];
				}
				setSlice(sl,slice);
			}
		}
	},
	RowStandardizeCentered = function(){
		if(nSlices()==0)
			return;
		for(sl in 1:nSlices()) {
			slice = getSlice(sl);
			div = sqrt(rowSums(slice^2));
			div[div==0] = 1;
			setSlice(sl, slice/div);
		}
	},
	CombineInOneSlice = function(){
		if(nSlices()==0)
			return;
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
		if(nSlices()>0)
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
		if(nSlices()==0)
			return;
		for(sl in 1:nSlices()) {
			setSlice(sl, getSlice(sl) %*% multiplier);
		}
	},
	ColumnSubsample = function(subset) {
		if(nSlices()==0)
			return;
		for(sl in 1:nSlices()) {
			setSlice(sl, getSlice(sl)[ ,subset]);
		}
		columnNames <<- columnNames[subset];
	},
	RowRemoveZeroEps = function(){
		if(nSlices()==0)
			return;
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

.impute_row_mean = function(x) {
	x[ is.na(x) ] = NaN;
	rowmean = rowMeans(x, na.rm = TRUE);
	rowmean[is.na(rowmean)] = 0;
	for(j in which(!complete.cases(x))) {  # 1:ncol(x)
		where1 = is.na(x[j, ]);
		x[j,where1] = rowmean[j];
	}
	rez = list();
	rez[[1]] = x;
	return(rez);
}

.snps_split_for_ANOVA = function(x) {
		# split into 2 dummy variables
		# status('Splitting genotype variable into dummy variables (for ANOVA)');

		uniq = unique(c(x));
		uniq = uniq[!is.na(uniq)];
		
		if( length(uniq) > 3 ) {
			stop('More than three genotype categories');
		} else if ( length(uniq) < 3 ) {
			uniq = c(uniq, min(uniq)-(1:(3-length(uniq))));
		}
		
		freq = matrix(0,nrow(x),length(uniq));
		for(i in 1:length(uniq)) {
			freq[ ,i] = rowSums(x==uniq[i],na.rm = TRUE);
		}
		
		md = apply(freq,1,which.max);
		freq[cbind(1:nrow(x),md)] = -1;
		
		md = apply(freq,1,which.max); # min(freq[cbind(1:nrow(slice),md)] - rowSums(select,na.rm = TRUE ))
		new_slice1 = (x == uniq[md]);
		new_slice1[is.na(new_slice1)] = 0;
		freq[cbind(1:nrow(x),md)] = -1;
				
		md = apply(freq,1,which.max);
		new_slice2 = (x == uniq[md]);
		new_slice2[is.na(new_slice2)] = 0;
		rez = list();
		rez[[1]]=new_slice1;
		rez[[2]]=new_slice2;
		return(rez);
}

#.SetNanRowMean = function(x) {
	#if(any(is.na(x))) {
		#rowmean = rowMeans(x, na.rm = TRUE);
		#rowmean[is.na(rowmean)] = 0;
		#for(j in which(!complete.cases(x))) {
			#where1 = is.na(x[j, ]);
			#x[j,where1] = rowmean[j];
		#}
	#}
	#return(x);
#}

.RowStandardizeCentered = function(x) {
	div = sqrt(rowSums(x^2));
	div[div==0] = 1;
	return(x/div);
}

.CovariatesProcessing = function(cvrt, nsamples, correctionMatrix) {
	# Add constant as a covariate, impute and combine in one slice
	if( cvrt$nRows()>0 ) {
		cvrt$SetNanRowMean();
		cvrt$CombineInOneSlice();
		cvrt = rbind(matrix(1,1,nsamples),cvrt$getSlice(1));
	} else {
		cvrt = matrix(1,1,nsamples);
	}
	
	# Correct for the error covariance structure
	if( length(correctionMatrix)>0 ) {
		# status('Rotating cvrt based on the errorCovariance matrix');
		cvrt = cvrt * correctionMatrix;
	}
	
	# Orthonormalize covariates
	# status('Orthonormalizing covariates');
	q = qr(t(cvrt));
	if( min(abs(diag(qr.R(q)))) < .Machine$double.eps*nsamples) {
		stop('Colinear or zero covariates detected');
	}
	return( t(qr.Q(q)) );
}

Matrix_eQTL_engine = function(snps, gene, cvrt = SlicedData$new(), output_file_name, pvOutputThreshold = 1e-5, useModel = modelLINEAR, errorCovariance = numeric(), verbose=FALSE){
	
	gene = gene$Clone();
	snps = snps$Clone();
	cvrt = cvrt$Clone();
	
	## the timing routing
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
	if( !verbose ) status = function(text){}

	if( useModel == modelLINEAR ) {
		snps_process = .impute_row_mean;
		statistic_name = 't-stat';
		nVarTested = 1;
	} else if( useModel == modelANOVA ) {
		snps_process = .snps_split_for_ANOVA;
		statistic_name = 'F-test';
		nVarTested = 2;
	} else stop('Unknown value of useModel');

	# Check dimensions
	status('Checking input data dimensions');
	if( snps$nRows()*snps$nCols() == 0 )
		stop('Empty genotype dataset');
	if( gene$nRows()*gene$nCols() == 0 )
		stop('Empty expression dataset');
	if( snps$nCols() != gene$nCols() )
		stop('Different number of samples in the genotype and gene expression files');
	
	if( cvrt$nRows()>0 )
		if( snps$nCols() != cvrt$nCols() )
			stop('Wrong number of samples in the file with covariates');

	################################# error covariance processing #################################
	
	if(length(errorCovariance)>0) {
		status('Processing the error covariance matrix');
		errorCovariance = as.matrix(errorCovariance);
		if( nrow(errorCovariance) != ncol(errorCovariance) ) {
			stop('The covariance matrix is not square');
		}	
		if( nrow(errorCovariance) != snps$nCols() ) {
			stop('The covariance matrix size does not match the data');
		}
		# test for symmetry
		if( !all(errorCovariance == t(errorCovariance)) ) {
			stop('The covariance matrix is not symmetric');
		}
		eig = eigen(errorCovariance, symmetric = TRUE)
		d = eig$values;
		v = eig$vectors;
		#  errorCovariance == v %*% diag(d) %*% t(v)
		#  errorCovariance^0.5 == v*sqrt(d)*v'
		#  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v'
		if( any(d<=0) ) {
			stop('The covariance matrix is not positive definite');
		}
		correctionMatrix = v %*% diag(1./sqrt(d)) %*% t(v);
		rm(eig,v,d,errorCovariance)
	} else {
		rm(errorCovariance);
		correctionMatrix = numeric();
	}
	
	################################# covariates processing #################################
	
	cvrt = .CovariatesProcessing(cvrt, snps$nCols(), correctionMatrix);
		
	################################# gene expression processing #################################
	
	# Impute gene expression
	status('Imputing missing expression');
	gene$SetNanRowMean();
	# Correct for the error covariance structure
	if( length(correctionMatrix)>0 ) {
		status('Rotating expression based on the errorCovariance matrix');
		gene$RowMatrixMultiply(correctionMatrix);
	}
	
	# Orthogonolize expression w.r.t. covariates
	gene$RowStandardizeCentered();
	status('Orthogonolizing expression w.r.t. covariates');
	for(sl in 1:gene$nSlices()) {
		slice = gene$getSlice(sl);
		slice = slice - tcrossprod(slice,cvrt) %*% cvrt;
		gene$setSlice(sl, slice);
	}
	gene$RowRemoveZeroEps();
	
	status('Standardizing expression');
	gene$RowStandardizeCentered();
	
	################################# Prepare for main loop    #################################
	
	nSamples = snps$nCols();
	nGenes = gene$nRows();
	nSnps  = snps$nRows();
	nCov = nrow(cvrt);
	# nVarTested = length(snps_list); # set in case(useModel)
	# dfNull = nSamples - nCov;
	dfFull = nSamples - nCov - nVarTested;

	if(useModel == modelLINEAR) {
		rThreshold = function(pth) {
			if(pth >= 1) return(0);
			tThresh = qt(pth/2, dfFull, lower.tail = FALSE);
			return( sqrt(  tThresh^2 / (dfFull + tThresh^2) ) );
		}
		afun = function(x) {return(abs(x))};
		testfun = function(x) { return( x * sqrt( dfFull / (1 - x^2)));	}
		pvfun = function(x) { return( pt(-abs(x),dfFull)*2 ); }
		
		theThresh = rThreshold(pvOutputThreshold);
	} else if(useModel == modelANOVA) {
		r2Threshold = function(pth) {
			if(pth >= 1) return(0);
			fThresh = qf(pth, nVarTested, dfFull, lower.tail = FALSE);
			return( fThresh * nVarTested / (dfFull + fThresh*nVarTested) );
		}
	
		afun = function(x) { return(x) };
		testfun = function(x) { return( x/(1-x) * (dfFull/nVarTested) ); }
		pvfun = function(x) { return( pf(x, nVarTested, dfFull, lower.tail = FALSE) ); }
	
		theThresh = r2Threshold(pvOutputThreshold);
	}
	
	gene_names= c(gene$rowNameSlices	, recursive=TRUE ); 
	snps_names= c(snps$rowNameSlices	, recursive=TRUE ); 
	FDR_collection = vector("list", gene$nSlices()*snps$nSlices())
	FDR_total_count = 0;
	
	totalCount = nGenes*nSnps;
	dumpCount = 0;

	################################# Main loop #################################

	status('Performing eQTL analysis');
	snps_offset = 0;
	FDR_pos = 1;
	for(sc in 1:snps$nSlices()) {
		gene_offset = 0;
		
		################################# prepare SNPs slice #################################
		# get, impute, split in dummies
		cursnps = snps_process( snps$getSlice(sc) );
		
		for(p in 1:length(cursnps)) {
			if(length(correctionMatrix)>0)
				cursnps[[p]] = cursnps[[p]]*correctionMatrix;
			cursnps[[p]] = cursnps[[p]] - tcrossprod(cursnps[[p]],cvrt) %*% cvrt;
			if( p>1 ) {
				for(w in 1:(p-1))
					cursnps[[p]] = cursnps[[p]] - rowSums(cursnps[[p]]*cursnps[[w]]) * cursnps[[w]]
			}
			cursnps[[p]] = .RowStandardizeCentered(cursnps[[p]]);
		}
		
		nrcs = nrow(cursnps[[1]])

		for(gc in 1:gene$nSlices()) {
			curgene = gene$getSlice(gc);
			nrcg = nrow(curgene)
			if(useModel == modelLINEAR) {
				mat = tcrossprod(cursnps[[1]],curgene);
			} else if(useModel == modelANOVA) {
				mat = tcrossprod(cursnps[[1]],curgene)^2;
				for(d in 2:nVarTested) {
					mat = mat + tcrossprod(cursnps[[d]],curgene)^2;
				}
			}
			select = which( afun(mat) >= theThresh, arr.ind=TRUE );
			rsub  = mat[select];
			rm(mat);
			test = testfun(rsub);
			pv = pvfun(test);
			FDR_collection[[FDR_pos]] = 
				t(cbind(snps_offset+select[,1], gene_offset+select[,2], test, pv));
				
			dumpCount = dumpCount + length(pv);
			
			gene_offset = gene_offset + nrcg;
			FDR_total_count = FDR_total_count + + nrcs*nrcg;
			FDR_pos = FDR_pos + 1;
			
			cat( floor(FDR_total_count/totalCount*1000)/10, '% done, ',dumpCount,' eQTLs found.\n');
			flush.console();
		}
		snps_offset = snps_offset + nrow(cursnps[[1]]);
	}

	fid = file(description = output_file_name, open = 'wt', blocking = FALSE, raw = FALSE)
	writeLines(paste('SNP\tgene\t',statistic_name,'\tp-value\tFDR', sep = ''), fid)
	status('Saving results');
	.Calc_FDR_and_save(FDR_collection, FDR_total_count, snps_names, gene_names, fid, -1)
	close(fid);
	status('');
}


Matrix_eQTL_engine_cis = function(snps, 
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
								cisDist = 1e6 ) {
	
	gene = gene$Clone();
	snps = snps$Clone();
	cvrt = cvrt$Clone();
	
	## the timing routing
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
	if( !verbose ) status = function(text){}

	if( useModel == modelLINEAR ) {
		snps_process = .impute_row_mean;
		statistic_name = 't-stat';
		nVarTested = 1;
	} else if( useModel == modelANOVA ) {
		snps_process = .snps_split_for_ANOVA;
		statistic_name = 'F-test';
		nVarTested = 2;
	} else stop('Unknown value of useModel');

	# Check dimensions
	status('Checking input data dimensions');
	if( snps$nRows()*snps$nCols() == 0 )
		stop('Empty genotype dataset');
	if( gene$nRows()*gene$nCols() == 0 )
		stop('Empty expression dataset');
	if( snps$nCols() != gene$nCols() )
		stop('Different number of samples in the genotype and gene expression files');
	
	if( cvrt$nRows()>0 )
		if( snps$nCols() != cvrt$nCols() )
			stop('Wrong number of samples in the file with covariates');

	################################# error covariance processing #################################
	
	if(length(errorCovariance)>0) {
		status('Processing the error covariance matrix');
		errorCovariance = as.matrix(errorCovariance);
		if( nrow(errorCovariance) != ncol(errorCovariance) ) {
			stop('The covariance matrix is not square');
		}	
		if( nrow(errorCovariance) != snps$nCols() ) {
			stop('The covariance matrix size does not match the data');
		}
		# test for symmetry
		if( !all(errorCovariance == t(errorCovariance)) ) {
			stop('The covariance matrix is not symmetric');
		}
		eig = eigen(errorCovariance, symmetric = TRUE)
		d = eig$values;
		v = eig$vectors;
		#  errorCovariance == v %*% diag(d) %*% t(v)
		#  errorCovariance^0.5 == v*sqrt(d)*v'
		#  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v'
		if( any(d<=0) ) {
			stop('The covariance matrix is not positive definite');
		}
		correctionMatrix = v %*% diag(1./sqrt(d)) %*% t(v);
		rm(eig,v,d,errorCovariance)
	} else {
		rm(errorCovariance);
		correctionMatrix = numeric();
	}
	
	################################# covariates processing #################################
	
	cvrt = .CovariatesProcessing(cvrt, snps$nCols(), correctionMatrix);
	
	################################# gene expression processing #################################
	
	# Impute gene expression
	status('Imputing missing expression');
	gene$SetNanRowMean();
	# Correct for the error covariance structure
	if( length(correctionMatrix)>0 ) {
		status('Rotating expression based on the errorCovariance matrix');
		gene$RowMatrixMultiply(correctionMatrix);
	}
	
	# Orthogonolize expression w.r.t. covariates
	gene$RowStandardizeCentered();
	status('Orthogonolizing expression w.r.t. covariates');
	for(sl in 1:gene$nSlices()) {
		slice = gene$getSlice(sl);
		slice = slice - tcrossprod(slice,cvrt) %*% cvrt;
		gene$setSlice(sl, slice);
	}
	gene$RowRemoveZeroEps();
	
	status('Standardizing expression');
	gene$RowStandardizeCentered();
	
	################################# matching gene and SNPs locations #################################
	
	status('Matching data files and location files')
	
	# names in the input data	
	gene_names = c(gene$rowNameSlices	, recursive=TRUE ); 
	snps_names = c(snps$rowNameSlices	, recursive=TRUE ); 

	# match with the location data
	usedgene = matrix(FALSE, nrow(genepos),1);
	genematch = match( gene_names, genepos[ ,1],  nomatch = 0);
	usedgene[ genematch ] = TRUE;
	if(length(genematch) == 0)
		stop('Gene names do not match those in the gene location file.');
	cat(sum(genematch>0),'of',length(gene_names),' genes matched\n');
	
	usedsnps = matrix(FALSE, nrow(snpspos),1);
	snpsmatch = match( snps_names, snpspos[ ,1],  nomatch = 0);
	usedsnps[ snpsmatch ] = TRUE;
	if(length(genematch) == 0)
		stop('SNP names do not match those in the SNP location file.');
	cat(sum(snpsmatch>0),'of',length(snps_names),' genes matched\n');
	
	# find accessed chr names
	chrNames = unique(c( as.character(unique(snpspos[usedsnps,2])), as.character(unique(genepos[usedgene,2])) ))
	
	# match chr names
	genechr = match(genepos[,2],chrNames);
	snpschr = match(snpspos[,2],chrNames);
	
#	chrMax = matrix(0,length(chrNames),1);
#	for (i in 1:length(chrNames)){
#		chrMax[i] = max(   snpspos[usedsnps & (snpschr == i),3],   genepos[usedgene& (genechr == i),4])
#	}
#	chrOffset = c(0,cumsum(chrMax+sicDist));
	
	# max length of chromosome
	chrMax = max( snpspos[usedsnps, 3], genepos[usedgene, 4])
	
	## set single number location
	genepos2 = genepos;
 	genepos2[ ,3:4] = genepos2[ ,3:4] + (genechr-1)*(chrMax+cisDist);
	
	snpspos2 = snpspos;
 	snpspos2[ ,3  ] = snpspos2[ ,3  ] + (snpschr-1)*(chrMax+cisDist);
	
	snps_pos = matrix(0,length(snps_names),1);
	snps_pos[snpsmatch>0] = snpspos2[snpsmatch,3];
	snps_pos[snps_pos==0] = (length(chrNames)+1)*(chrMax+cisDist);
	
	gene_pos = matrix(0,length(gene_names),2);
	gene_pos[genematch>0,] = as.matrix(genepos2[genematch,3:4]);
	gene_pos[gene_pos==0] = (length(chrNames)+2)*(chrMax+cisDist);
	
	rm(genematch, usedgene, snpsmatch, usedsnps, chrNames, genechr, snpschr, chrMax, genepos2, snpspos2)
	
	# Slice it back.
	geneloc = vector("list", gene$nSlices())
	gene_offset = 0;
	for(gc in 1:gene$nSlices()) {
		nr = length(gene$rowNameSlices[[gc]]);
		geneloc[[gc]] = gene_pos[gene_offset + (1:nr), ];
		gene_offset = gene_offset + nr;	
	}
	rm(gc, gene_offset, gene_pos);
	
	snpsloc = vector("list", snps$nSlices())
	snps_offset = 0;
	for(sc in 1:snps$nSlices()) {
		nr = length(snps$rowNameSlices[[sc]]);
		snpsloc[[sc]] = snps_pos[snps_offset + (1:nr), ];
		snps_offset = snps_offset + nr;	
	}
	rm(sc, snps_offset, snps_pos);
	
	################################# Prepare for main loop    #################################
	
	nSamples = snps$nCols();
	nGenes = gene$nRows();
	nSnps  = snps$nRows();
	nCov = nrow(cvrt);
	# nVarTested = length(snps_list); # set in case(useModel)
	# dfNull = nSamples - nCov;
	dfFull = nSamples - nCov - nVarTested;
	
	if(useModel == modelLINEAR) {
		rThreshold = function(pth) {
			if(pth >= 1) return(0);
			tThresh = qt(pth/2, dfFull, lower.tail = FALSE);
			return( sqrt(  tThresh^2 / (dfFull + tThresh^2) ) );
		}
		afun = function(x) {return(abs(x))};
		testfun = function(x) { return( x * sqrt( dfFull / (1 - x^2)));	}
		pvfun = function(x) { return( pt(-abs(x),dfFull)*2 ); }		
		theThresh_cis = rThreshold(pvOutputThreshold_cis);
		theThresh_tra = rThreshold(pvOutputThreshold_tra);
	} else if(useModel == modelANOVA) {
		r2Threshold = function(pth) {
			if(pth >= 1) return(0);
			fThresh = qf(pth, nVarTested, dfFull, lower.tail = FALSE);
			return( fThresh * nVarTested / (dfFull + fThresh*nVarTested) );
		}
	
		afun = function(x) { return(x) };
		testfun = function(x) { return( x/(1-x) * (dfFull/nVarTested) ); }
		pvfun = function(x) { return( pf(x, nVarTested, dfFull, lower.tail = FALSE) ); }	
		theThresh_cis = r2Threshold(pvOutputThreshold_cis);
		theThresh_tra = r2Threshold(pvOutputThreshold_tra);
	}
	
	FDR_collection_cis = vector("list", gene$nSlices()*snps$nSlices())
	FDR_count_cis = 0;
	FDR_collection_tra = vector("list", gene$nSlices()*snps$nSlices())
	FDR_count_tra = 0;
	FDR_total_count = 0;
	
	totalCount = nGenes*nSnps;
	dumpCount_cis = 0;
	dumpCount_tra = 0;

	################################# Main loop #################################

	status('Performing eQTL analysis');
	snps_offset = 0;
	FDR_pos = 1;
	for(sc in 1:snps$nSlices()) { #snps$nSlices()
		gene_offset = 0;
		
		################################# prepare SNPs slice #################################
		# get, impute, split in dummies
		cursnps = snps_process( snps$getSlice(sc) );
		
		for(p in 1:length(cursnps)) {
			if(length(correctionMatrix)>0) {
				cursnps[[p]] = cursnps[[p]]*correctionMatrix;
			}
			cursnps[[p]] = cursnps[[p]] - tcrossprod(cursnps[[p]],cvrt) %*% cvrt;
			if( p>1 ) {
				for(w in 1:(p-1))
					cursnps[[p]] = cursnps[[p]] - rowSums(cursnps[[p]]*cursnps[[w]]) * cursnps[[w]]
			}
			cursnps[[p]] = .RowStandardizeCentered(cursnps[[p]]);
		}
		
		nrcs = nrow(cursnps[[1]])	
		
		for(gc in 1:gene$nSlices()) {
			curgene = gene$getSlice(gc);
			nrcg = nrow(curgene)
			
			srep = rep( snpsloc[[sc]], each = nrcg);
			iscis = ( srep > (geneloc[[gc]][ ,1] - cisDist)) &
					( srep < (geneloc[[gc]][ ,2] + cisDist));
			rm(srep);
			dim(iscis) = c(nrcg, nrcs);
			
			curcis = sum(iscis);
			FDR_count_cis = FDR_count_cis + curcis;
			FDR_count_tra = FDR_count_tra + (nrcs*nrcg - curcis);
			
			if(useModel == modelLINEAR) {
				mat = tcrossprod(curgene, cursnps[[1]]);
			} else if(useModel == modelANOVA) {
				mat = tcrossprod(curgene, cursnps[[1]])^2;
				for(d in 2:nVarTested) {
					mat = mat + tcrossprod(curgene, cursnps[[d]])^2;
				}
			}
			amat = afun(mat);
			select_cis = which((amat >= theThresh_cis) &  iscis, arr.ind=TRUE );
			select_tra = which((amat >= theThresh_tra) & !iscis, arr.ind=TRUE );
			rm(amat)
			rsub_tra  = mat[select_tra];
			rsub_cis  = mat[select_cis];
			rm(mat);
			test = testfun(rsub_tra);
			pv = pvfun(test);
			dumpCount_tra = dumpCount_tra + length(pv);
			FDR_collection_tra[[FDR_pos]] = 
				t(cbind(snps_offset+select_tra[,2], gene_offset+select_tra[,1], test, pv));
			
			test = testfun(rsub_cis);
			pv = pvfun(test);
			dumpCount_cis = dumpCount_cis + length(pv);
			FDR_collection_cis[[FDR_pos]] = 
				t(cbind(snps_offset+select_cis[,2], gene_offset+select_cis[,1], test, pv));

			
			gene_offset = gene_offset + nrcg;
			FDR_total_count = FDR_total_count + nrcs*nrcg;
			FDR_pos = FDR_pos + 1;
			
			flush.console();
		}
		cat( floor(FDR_total_count/totalCount*1000)/10, '% done, ',dumpCount_cis,'+',dumpCount_tra,' gene-snp pairs significant from ',FDR_count_cis,'+',FDR_count_tra,' total\n');
		snps_offset = snps_offset + nrow(cursnps[[1]]);
	}

	fid = file(description = output_file_name, open = 'wt', blocking = FALSE, raw = FALSE)
	writeLines(paste('SNP\tgene\t',statistic_name,'\tp-value\tFDR\tis_cis-', sep = ''), fid)
	status('Processing cis-eQTLs');
	.Calc_FDR_and_save(FDR_collection_cis, FDR_count_cis, snps_names, gene_names, fid, 1)
	status('Processing trans-eQTLs');
	.Calc_FDR_and_save(FDR_collection_tra, FDR_count_tra, snps_names, gene_names, fid, 0)
	close(fid);
	status('');
}

.Calc_FDR_and_save = function(FDR_collection, FDR_total_count, snps_names, gene_names, fid, cis_indicator) {
	FDR_collection = c(FDR_collection, recursive=TRUE ); 
	if(length(FDR_collection) == 0) {
		return();
	}
	dim(FDR_collection) = c(4, length(FDR_collection)/4);
	
	order1 = sort.list(FDR_collection[4,]);
	FDR_collection = FDR_collection[,order1, drop = FALSE];
	rm(order1)
	
	FDR = FDR_collection[4,]*FDR_total_count/(1:ncol(FDR_collection));
	FDR[length(FDR)] = min(FDR[length(FDR)], 1);
	FDR = rev(cummin(rev(FDR)))
	
	#if(length(FDR)>1)
 	#	for(i in (length(FDR)-1):1)
	#		FDR[i] = min(FDR[i],FDR[i+1]);
	#status('Saving results');
	#fid = file(description = output_file_name, open = 'wt', blocking = FALSE, raw = FALSE)
	#writeLines(paste('SNP\tgene\t',statistic_name,'\tp-value\tFDR', sep = ''), fid)
	step = 100000;
	for(part in 1:ceiling(length(FDR)/step)) {
		ndx = ((part-1)*step + 1) :  min(part*step,length(FDR))
		dump2 = data.frame(snps_names[FDR_collection[1,ndx]],gene_names[FDR_collection[2,ndx]],FDR_collection[3,ndx], FDR_collection[4,ndx], FDR[ndx], row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = FALSE)
		if(cis_indicator>=0)
			dump2[['cis']] = cis_indicator;
		write.table(dump2, file = fid, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE);
	}
}