# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

modelLINEAR = 117348L;
modelANOVA  = 47074L;
modelLINEAR_CROSS = 1113461L;

SlicedData <- setRefClass( "SlicedData",
	fields = list( 
		dataEnv = "environment",
		nSlices1 = "numeric",
		rowNameSlices = "list",
		columnNames = "character",
		fileDelimiter = "character",
		fileSkipColumns = "numeric",
		fileSkipRows = "numeric",
		fileSliceSize = "numeric",
		fileOmitCharacters = "character"
	),
	methods = list(
	initialize = function( mat = NULL ) {
		dataEnv <<- new.env(hash = TRUE, size = 29L);
		nSlices1 <<- 0L;
		if(!is.null(mat)) {
			CreateFromMatrix(mat);
		}
		fileSliceSize <<- 1000;
		fileDelimiter <<- "\t";
		fileSkipColumns <<- 1L;
		fileSkipRows <<- 1L;
		fileOmitCharacters <<- "NA"
		return(invisible(.self));
	},
	CreateFromMatrix = function( mat ) {
		stopifnot( class(mat) == "matrix" );
		setSliceRaw( 1L ,mat );
		rns = rownames( mat );
		if( is.null(rns) ) {
			rns = paste( "Row_",(1:nrow(mat)), sep="" );
		}
		rowNameSlices <<- list(rns);
		cns = colnames( mat );
		if( is.null(cns) ){
			cns = paste( "Col_",(1:ncol(mat)), sep="" );
		}
		columnNames <<- cns;
		return(invisible(.self));
	},
	getSlice = function(sl) {
		value = get(paste(sl), dataEnv);
		if( is.raw(value) ) {
			storage.mode(value) = "double";
			value[value == 255] = NA;
		}
		return( value  )	
	},
	getSliceRaw = function(sl) {
		return( get(paste(sl), dataEnv) )	
	},
	setSliceRaw = function(sl, value) {
		assign( paste(sl), value, dataEnv )
		if( nSlices1 < sl ) {
			nSlices1 <<- sl;
		}
	},
	setSlice = function(sl, value) {
		if( length(value) > 0 ) {
			if( all(as.integer(value) == value, na.rm = TRUE) ) {
				if( (min(value, na.rm = TRUE) >= 0 ) && 
                            (max(value, na.rm = TRUE) < 255) )
				{
					nv = value;
					suppressWarnings({storage.mode(nv) = "raw"});
					nv[ is.na(value)] = as.raw(255);
					value = nv;
				} else {
					storage.mode(value) = "integer";
				}
			}
		}
		setSliceRaw(sl, value);
	},	
	nSlices = function() {
		return( nSlices1 );
	},
	LoadFile = function(filename, skipRows = NULL, skipColumns = NULL, sliceSize = NULL, omitCharacters = NULL, delimiter = NULL) {
		if( !is.null(skipRows) ) {
			fileSkipRows <<- skipRows;
		}
		if( !is.null(skipColumns) ) {
			fileSkipColumns <<- skipColumns;
		}
		if( !is.null(omitCharacters) ) {
			fileOmitCharacters <<- omitCharacters;
		}
		if( !is.null(sliceSize) ) {
			fileSliceSize <<- sliceSize;
		}
		if( !is.null(delimiter) ) {
			fileDelimiter <<- delimiter;
		}
		fid = file(description = filename, open = "rt", blocking = FALSE, raw = FALSE)
		# clean object if file is open
		Clear(); 
		lines = readLines(con = fid, n = max(fileSkipRows,1L), ok = TRUE, warn = TRUE)
		line1 = tail(lines,1);
		splt = strsplit(line1, split = fileDelimiter, fixed = TRUE);
		if( fileSkipRows > 0L ) {
			columnNames <<- splt[[1]]; # [ -(1:fileSkipColumns) ];
		} else {
			seek(fid, 0)
		}		
		
		rm( lines, line1, splt );
		
		rowNameSlices <<- vector("list", 15);

		curSliceId = 0L;
		repeat
		{
			# preallocate names and data
			if(length(rowNameSlices) < curSliceId) {
				rowNameSlices[[2L*curSliceId]] <<- NULL;
			}
			curSliceId = curSliceId + 1L;
			
			# read sliceSize rows
			rowtag = vector("character",fileSliceSize);
			rowvals = vector("list",fileSliceSize);
			for(i in 1:fileSliceSize) {
				temp = "";
				if( fileSkipColumns > 0L ) {
					temp = scan(file = fid, what = character(), n = fileSkipColumns, quiet = TRUE,sep = fileDelimiter);
				}
				rowtag[i] = temp[1];#paste(temp,collapse=" ");
				rowvals[[i]] = scan(file = fid, what = double(), nlines = 1, quiet = TRUE, sep = fileDelimiter, na.strings = fileOmitCharacters);
				if( length(rowvals[[i]]) == 0L ) {
					if(i==1L) {
						rowtag = matrix(0, 0, 0);
						rowvals = character(0);
					} else 	{
						rowtag  = rowtag[  1:(i-1) ];
						rowvals = rowvals[ 1:(i-1) ];
					}
					break;			
				}
			}
			if( length(rowtag) == 0L ) {
				curSliceId = curSliceId - 1L;
				break;
			}
			rowNameSlices[[curSliceId]] <<- rowtag;
			data = c(rowvals, recursive = TRUE);
			dim(data) = c(length(rowvals[[1]]), length(rowvals));
			data = t(data);
			setSlice(curSliceId, data);
			if( length(rowtag) < fileSliceSize ) {
				break;
			}
			numtxt = formatC(curSliceId*fileSliceSize, big.mark=",", format = "f", digits = 0)
			cat( "Rows read: ", numtxt, "\n");
			flush.console()
		}
		close(fid)
		if( fileSkipRows == 0 ) {
			columnNames <<- paste("Col_", (1:nCols()), sep="");
		} else {
			columnNames <<- tail(columnNames, ncol(getSliceRaw(1)));
		}
		if( fileSkipColumns == 0 ) {
			cnt = 0L;
			for( sl in 1:nSlices() ) {
				nr = length(getSliceRaw(sl));
				rowNameSlices[[sl]] <<- paste("Row_",cnt + (1:nr),sep="");
				cnt = cnt + nr;
			}
		}
		rowNameSlices <<- rowNameSlices[1:curSliceId];
		cat("Rows read: ", nRows(), " done.\n");
		return(invisible(.self));
	},
	SaveFile = function(filename) {
		if( nSlices() == 0 ) {
			warning("No data to save");
			return();
		}
		fid = file(filename,"wt");
		for( sl in 1:nSlices() ) {
			z = getSlice(sl);
			rownames(z) = rowNameSlices[[sl]];
			colnames(z) = columnNames;
			write.table(z, file = fid, sep = "\t", 
				col.names = (if(sl == 1){NA}else{FALSE}));
		}
		close(fid);
	},
	nRows = function() {
		if( nSlices() == 0L ) {
			return(0L);
		} else {
			s = 0L;
			for(sl in 1:nSlices()) {
				s = s + nrow(getSliceRaw(sl));
			}
			return( s )
		}
	},
	nCols = function() {
		if( nSlices() == 0L ) {
			return(0L);
		} else {
			return( ncol(getSliceRaw(1L)) )
		}
	},
	Clear = function() {
		if( nSlices() > 0L ) {
			for( sl in 1:nSlices() ) {
				rm(list = paste(sl), envir = dataEnv)
			}
		}
		nSlices1 <<- 0L;
		rowNameSlices <<- list();
		columnNames <<- character();
		return(invisible(.self));
	},
	IsCombined = function() {
		return( nSlices() <= 1L );
	},
	GetAllRowNames = function() {
		return( c(rowNameSlices, recursive=TRUE) );
	},
	SetNanRowMean = function() {
		if( (nSlices() == 0L) || (nCols() == 0L) ) {
			return(invisible(.self));
		}
		for( sl in 1:nSlices() ) {
			slice = getSlice(sl);
			if( any(is.na(slice)) ) {
				rowmean = rowMeans(slice, na.rm = TRUE);
				rowmean[is.na(rowmean)] = 0L;
				for( j in which(!complete.cases(slice)) ) {
					where1 = is.na(slice[j, ]);
					slice[j, where1] = rowmean[j];
				}
				setSlice(sl, slice);
			}
		}
		return(invisible(.self));
	},
	RowStandardizeCentered = function() {
		if( nSlices() == 0L ) {
			return(invisible(.self));
		}
		for(sl in 1:nSlices()) {
			slice = getSlice(sl);
			div = sqrt( rowSums(slice^2) );
			div[ div == 0 ] = 1;
			setSlice(sl, slice/div);
		}
		return(invisible(.self));
	},
	CombineInOneSlice = function() {
		if( nSlices() <= 1L ) {
			return(invisible(.self));			
		}
		nc = nCols();
		nr = nRows();
		datatypes = c("raw","integer","double");
		datafuns = c(as.raw, as.integer, as.double);
		datatype = character(nSlices());
		for(sl in 1:nSlices()) {
			datatype[sl] = typeof(getSliceRaw(sl));
		}
		mch = max(match(datatype,datatypes,nomatch = length(datatypes)));
		datafun = datafuns[[mch]];
		newData = matrix(datafun(0), nrow = nr, ncol = nc);
		offset = 0;
		for(sl in 1:nSlices()) {
			if(mch==1) {
				slice = getSliceRaw(sl);
			} else {
				slice = getSlice(sl);
			}
			newData[ offset + (1:nrow(slice)),] = datafun(slice);
			setSlice(sl, numeric());
			offset = offset + nrow(slice);
		}
		
		nSlices1 <<- 1L;
		setSliceRaw(1L, newData);
		rm(newData);
		
		newrowNameSlices = GetAllRowNames();
		rowNameSlices <<- list(newrowNameSlices)
		return(invisible(.self));
	},
	ResliceCombined = function(sliceSize = -1) {
		if( sliceSize > 0L ) {
			fileSliceSize <<- sliceSize;
		}
		if( fileSliceSize <= 0 ) {
			fileSliceSize <<- 1000;
		}
		if( IsCombined() ) {
			nRows1 = nRows();
			if(nRows1 == 0L) {
				return(invisible(.self));
			}
			newNSlices = floor( (nRows1 + fileSliceSize - 1)/fileSliceSize );
			oldData = getSliceRaw(1L);
			#oldNames = rowNameSlices[[1]];
			newNameslices = vector("list",newNSlices)
			for( sl in 1:newNSlices ) {
				range = (1+(sl-1)*fileSliceSize) : (min(nRows1,sl*fileSliceSize));
				newpart = oldData[range, ,drop = FALSE];
				if( is.raw(oldData) ) {
					setSliceRaw( sl, newpart);
				} else {
					setSlice( sl, newpart);
				}
				newNameslices[[sl]] = rowNameSlices[[1]][range];
			}
			rowNameSlices <<- newNameslices ;
		} else {
			stop("Reslice of a sliced matrix is not supported yet. Use CombineInOneSlice first.");
		}
		return(invisible(.self));
	},
	Clone = function() {
		clone = SlicedData$new();
		if( nSlices() > 0L ) {
			for(sl in 1:nSlices()) {
				clone$setSliceRaw(sl,getSliceRaw(sl));
			}
		}
		clone$rowNameSlices = rowNameSlices;
		clone$columnNames = columnNames;
		clone$fileDelimiter = fileDelimiter;
		clone$fileSkipColumns = fileSkipColumns;
		clone$fileSkipRows = fileSkipRows;
		clone$fileSliceSize = fileSliceSize;
		clone$fileOmitCharacters = fileOmitCharacters;
		return( clone );		
	},
	RowMatrixMultiply = function(multiplier) {
		if( nSlices() == 0L ) {
			return(invisible(.self));
		}
		for(sl in 1:nSlices()) {
			setSlice(sl, getSlice(sl) %*% multiplier);
		}
		return(invisible(.self));
	},
	ColumnSubsample = function(subset) {
		if( nSlices() == 0 ) {
			return(invisible(.self));
		}
		for( sl in 1:nSlices() ) {
			setSliceRaw(sl, getSliceRaw(sl)[ ,subset, drop = FALSE]);
		}
		columnNames <<- columnNames[subset];
		return(invisible(.self));
	},
	RowReorderSimple = function(ordr) {
		# had to use an inefficient and dirty method
		# due to horrible memory management in R
		if( (typeof(ordr) == "logical") && all(ordr) ) {
			return(invisible(.self));
		}
		if( (length(ordr) == nRows()) && all(ordr == (1:length(ordr))) ) {
			return(invisible(.self));
		}
		CombineInOneSlice();
		gc();
		setSliceRaw( 1L, getSliceRaw(1L)[ordr, ] );
		rowNameSlices[[1]] <<- rowNameSlices[[1]][ordr];
		gc();
		ResliceCombined();
		gc();
		return(invisible(.self));
	},
	RowReorder = function(ordr) {
		# transform logical into indices 
		if( typeof(ordr) == "logical" ) {
			if( length(ordr) == nRows() ) {
				ordr = which(ordr);
			} else {
				stop("Parameter \"ordr\" has wrong length")
			}
		}
		## first, check that anything has to be done at all
		if( (length(ordr) == nRows()) && all(ordr == (1:length(ordr))) ) {
			return(invisible(.self));
		}
		## check bounds
		#if( (min(ordr) < 1) || (max(ordr) > nRows()) ) {
		#	stop("Parameter \"ordr\" is out of bounds");
		#}
		## slice the data into individual rows
		all_rows = vector("list", nSlices())
		for( i in 1:nSlices() ) {
			slice = getSliceRaw(i)
			all_rows[[i]] = split(slice, 1:nrow(slice))
			setSliceRaw(i,numeric())
		}
		gc();
		all_rows = unlist(all_rows, recursive=FALSE, use.names = FALSE);
		## Reorder the rows
		all_rows = all_rows[ordr];
		## get row names
		all_names = GetAllRowNames();
		## erase the set
		rowNameSlices <<- list();
		## sort names
		all_names = all_names[ordr];
		##
		## Make slices back
		nrows = length(all_rows);
		nSlices1 <<- as.integer((nrows+fileSliceSize-1)/fileSliceSize);
		##cat(nrows, " ", nSlices1);
		rowNameSlices1 = vector("list", nSlices1);
		for( i in 1:nSlices1 ) {
			fr = 1 + fileSliceSize*(i-1);
			to = min( fileSliceSize*i, nrows);

			subset = all_rows[fr:to];
			types = unlist(lapply(subset,typeof));
			israw = (types == "raw")
			if(!all(israw == israw[1])) {
				# some raw and some are not
				subset = lapply(subset, function(x){if(is.raw(x)){x=as.integer(x);x[x==255] = NA;return(x)}else{return(x)}});
			}
			subset = unlist(subset);
			dim(subset) = c( length(all_rows[[fr]]) , to - fr + 1)
			#subset = matrix(subset, ncol = (to-fr+1));
			if(is.raw(subset)) {
				setSliceRaw(i, t(subset)); 
			} else {
				setSlice(i, t(subset)); 
			}
			rowNameSlices1[[i]] = all_names[fr:to];
			all_rows[fr:to] = 0;
			all_names[fr:to] = 0;
		}
		rowNameSlices <<- rowNameSlices1;
		gc();
		return(invisible(.self));
	},
	RowRemoveZeroEps = function(){
		if( nSlices() == 0 ) {
			return(invisible(.self));
		}
		for(sl in 1:nSlices()) {
			slice = getSlice(sl);
			amean = rowMeans(abs(slice));
			remove = (amean < .Machine$double.eps*nCols());
			if(any(remove)) {
				rowNameSlices[[sl]] <<- rowNameSlices[[sl]][!remove];
				setSlice(sl, slice[!remove, , drop = FALSE]);
			}
		}
		return(invisible(.self));
	},
	FindRow = function(rowname) {
		for( sl in 1:nSlices() ) {
			mch = match(rowname,rowNameSlices[[sl]], nomatch = 0);
			if( mch > 0 )
			{
				row = getSlice(sl)[mch[1], , drop=FALSE];
				rownames(row) = rowname;
				colnames(row) = columnNames;
				return( list(slice = sl, item = mch, row = row) );
			}
		}
		return( NULL );
	}
	)
)

#setGeneric("nrow")
setMethod("nrow", "SlicedData",	function(x) {
		return( x$nRows() );
	})
#setGeneric("NROW")
setMethod("NROW", "SlicedData",	function(x) {
		return( x$nRows() );
	})
#setGeneric("ncol")
setMethod("ncol", "SlicedData",	function(x) {
		return( x$nCols() );
	})
#setGeneric("NCOL")
setMethod("NCOL", "SlicedData",	function(x) {
		return( x$nCols() );
	})
#setGeneric("dim")
setMethod("dim", "SlicedData",	function(x) {
		return( c(x$nRows(),x$nCols()) );
	})
#setGeneric("colnames")
setMethod("colnames", "SlicedData",	function(x) {
		return( x$columnNames );
	})
#setGeneric("rownames")
setMethod("rownames", "SlicedData",	function(x) {
		return( x$GetAllRowNames() );
	})
setMethod("[[", "SlicedData",	function(x,i) {
		return( x$getSlice(i) );
	})
#setGeneric("length")
setMethod("length", "SlicedData",	function(x) {
		return( x$nSlices() );
	})
setMethod("[[<-", "SlicedData",	function(x,i,value) {
		x$setSlice(i, value);
		return(x);
})
#setGeneric("summary")
setMethod("summary", "SlicedData",	function(object) {
		z = c(nCols = object$nCols(), nRows = object$nRows(), nSlices = object$nSlices());
		return(z);
	})
#setGeneric("show", standardGeneric("show"))
setMethod("show", "SlicedData",	function(object) {
		cat("SlicedData object. For more information type: ?SlicedData\n");
		cat("Number of columns:", object$nCols(), "\n");
		cat("Number of rows:", object$nRows(), "\n");
		cat("Data is stored in", object$nSlices(), "slices\n");
		if(object$nSlices()>0) {
			z = object$getSlice(1);
			z = z[1:min(nrow(z),10), 1:min(ncol(z),10), drop = FALSE];
			rownames(z) = object$rowNameSlices[[1]][1:nrow(z)];
			colnames(z) = object$columnNames[1:ncol(z)];
			cat("Top left corner of the first slice (up to 10x10):\n");
			show(z)
		}		
	})
setMethod("as.matrix", "SlicedData", function(x) {
		if(x$nSlices() == 0) {
			return( matrix(0,0,0) );
		}
		if(x$nSlices() > 1) {
			copy = x$Clone();
			copy$CombineInOneSlice();
		} else {
			copy = x;
		}
		mat = copy$getSlice(1L);
		rownames(mat) = rownames(copy);
		colnames(mat) = colnames(copy);
		return( mat );
	})
setMethod("colnames<-", "SlicedData", function(x,value) {
		stopifnot( class(value) == "character" );
		stopifnot( length(value) == x$nCols() );
		x$columnNames = value;
		return(x);
	})
setMethod("rownames<-", "SlicedData", function(x,value) {
		stopifnot( class(value) == "character" );
		stopifnot( length(value) == x$nRows() );
		start = 1;
		newNameSlices = vector("list", x$nSlices());
		for( i in 1:x$nSlices() ) {
			nr = length( x$rowNameSlices[[i]] );
			newNameSlices[[i]] = value[ start:(start+nr-1) ];
			start = start + nr;
		}
		x$rowNameSlices = newNameSlices; 
		return(x);
	})

.OutputSaver_direct <- setRefClass(".OutputSaver_direct",
	fields = list( 
		dataEnv = "environment",
		gene_names = "character",
		snps_names = "character",
		gene1 = "SlicedData",
		snps1 = "SlicedData"
	),
	methods = list(
	initialize = function() {
		dataEnv <<- new.env(hash = FALSE, size = 1L);
		gene_names <<- character(0);
		snps_names <<- character(0);
		gene1 <<- SlicedData$new();
		snps1 <<- SlicedData$new();
		return(.self);
	},
	Start = function(filename, gene, snps, statistic_name) {
		# I hope the program stops if it fails to open the file
		fid = file(description = filename, open = "wt", blocking = FALSE, raw = FALSE);
		writeLines(paste("SNP\tgene\t", statistic_name, "\tp-value", sep = ""), fid);
		assign("fid", fid, envir=dataEnv);
		gene1 <<- gene;
		snps1 <<- snps;
	},
	WriteBlock = function(block) {
		if( length(block) == 0 ) {
			return();
		}
		if( length(gene_names) == 0 ) {
			gene_names <<- gene1$GetAllRowNames();
			gene1 <<- SlicedData$new();
			snps_names <<- snps1$GetAllRowNames();
			snps1 <<- SlicedData$new();
		}
		fid = get("fid", envir=dataEnv);
		dump2 = data.frame(snps_names[block[ ,1]],gene_names[block[ ,2]],block[ ,3], block[ ,4], row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = FALSE);
		write.table(dump2, file = fid, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);
	},
	Commit = function(unused) {
		fid = get("fid", dataEnv);
		close(fid);
		assign("fid", NULL, dataEnv)
		gene_names <<- character(0);
		snps_names <<- character(0);
		return(NULL);
	}
	)
)

.OutputSaver_FRD <- setRefClass(".OutputSaver_FRD",
	fields = list( 
		dataEnv = "environment",
		nBlocks = "numeric",
		gene_names = "character",
		snps_names = "character",
		gene1 = "SlicedData",
		snps1 = "SlicedData"
	),
	methods = list(
	initialize = function () {
		dataEnv <<- new.env(hash = TRUE);
		gene_names <<- character();
		snps_names <<- character();
		gene1 <<- SlicedData$new();
		snps1 <<- SlicedData$new();
		nBlocks <<- 0L;
		return(.self);
	},
	getBlock = function(bl) {
		return( get(paste(bl), dataEnv) )	
	},
	setBlock = function(bl, value) {
		assign(paste(bl), value, dataEnv)
		if( nBlocks < bl ) {
			nBlocks <<- bl;
		}
	},	
	Start = function(filename, gene, snps, statistic_name) {
		# I hope that if file is not open the program halts
		fid = file(description = filename, open = "wt", blocking = FALSE, raw = FALSE)
		writeLines( paste("SNP\tgene\t",statistic_name,"\tp-value\tFDR", sep = ""), fid)
		assign("fid", fid, envir=dataEnv)
		gene1 <<- gene;
		snps1 <<- snps;
	},
	WriteBlock = function(block) {
		if( length(block) == 0 ) {
			return();
		}
		setBlock( nBlocks + 1L, block );
	},
	CombineInOneBlock = function() {
		if(nBlocks <= 1) {
			return();
		}
		FDR_collection = vector("list", nBlocks)
		for(i in 1:nBlocks) {
			FDR_collection[[i]] = t(getBlock(i));
			setBlock(i,NULL);
		}
		FDR_collection = c(FDR_collection, recursive = TRUE );
		dim(FDR_collection) = c(4, length(FDR_collection)/4);
		FDR_collection = t(FDR_collection);
		nBlocks <<- 1L;	
		setBlock(1, FDR_collection);
	},
	RemoveCis = function( cis.class ) {
		# do nothing if there is nothing to do
		if( (nBlocks <= 0) || (cis.class$nBlocks <= 0) ) {
			return();
		}
		CombineInOneBlock();
		cis.class$CombineInOneBlock();
		cis.gene.snps = cis.class$getBlock(1);
		my.block = getBlock(1L);
		max.gene = max( cis.gene.snps[ ,2], my.block[ ,2] );
		
		remove = (	my.block[ ,1] * max.gene + my.block[ ,2] 
						) %in% ( 
					cis.gene.snps[ ,1] * max.gene + cis.gene.snps[ ,2] );
		setBlock( 1, my.block[!remove, , drop = FALSE] );
	},
	Commit = function( FDR_total_count ) {
		fid = get("fid", envir=dataEnv);
		
		gene_names <<- gene1$GetAllRowNames();
		gene1 <<- SlicedData$new();
		snps_names <<- snps1$GetAllRowNames();
		snps1 <<- SlicedData$new();

		CombineInOneBlock();		
		
		FDR = matrix(0,0,1);
		FDR_collection = matrix(0,0,4);
		
		
		if(nBlocks > 0) {
			FDR_collection = getBlock(1);
			
			order1 = sort.list(FDR_collection[ , 4]);
			FDR_collection = FDR_collection[order1, , drop = FALSE];
			rm(order1)
			
			FDR = FDR_collection[ , 4] * FDR_total_count / (1:nrow(FDR_collection));
			FDR[length(FDR)] = min(FDR[length(FDR)], 1);
			FDR = rev(cummin(rev(FDR)))
			
			step = 100000;
			if(length(FDR_collection) > 0) {
				for( part in 1:ceiling(length(FDR)/step) ) {
					ndx = ((part-1)*step + 1) :  min(part*step,length(FDR))
					dump2 = data.frame(	snps_names[FDR_collection[ndx, 1]],
										gene_names[FDR_collection[ndx, 2]],
										FDR_collection[ndx, 3], 
										FDR_collection[ndx, 4], 
										FDR[ndx], 
										row.names = NULL, 
										check.rows = FALSE, 
										check.names = TRUE, 
										stringsAsFactors = FALSE);
					write.table(dump2, file = fid, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);
				}
			}
		} else {
			cat("No significant associations were found.\n", file = fid);
		}
		close(fid);
		assign("fid", NULL, dataEnv)
		FDR_collection = data.frame(
			snps = snps_names[FDR_collection[ , 1]],
			gene = gene_names[FDR_collection[ , 2]],
			statistic = FDR_collection[ , 3],
			pvalue = FDR_collection[ , 4],
			FDR = FDR);
		gene_names <<- character(0);
		snps_names <<- character(0);
		return( FDR_collection );
	}
	)
)


.snps_split_for_ANOVA = function(x) {
		# split into 2 dummy variables

		uniq = unique(c(x));
		uniq = uniq[!is.na(uniq)];
		
		if( length(uniq) > 3 ) {
			stop("More than three genotype categories is not handled by ANOVA");
		} else if ( length(uniq) < 3 ) {
			uniq = c(uniq, min(uniq)-(1:(3-length(uniq))));
		}
		
		freq = matrix(0, nrow(x), length(uniq));
		for(i in 1:length(uniq)) {
			freq[ ,i] = rowSums(x==uniq[i], na.rm = TRUE);
		}
		
		md = apply(freq, 1, which.max);
		freq[ cbind(1:nrow(x),md) ] = -1;
		
		md = apply(freq, 1, which.max); # min(freq[cbind(1:nrow(slice),md)] - rowSums(select,na.rm = TRUE ))
		new_slice1 = (x == uniq[md]);
		new_slice1[is.na(new_slice1)] = 0;
		freq[ cbind(1:nrow(x),md) ] = -1;
				
		md = apply(freq,1,which.max);
		new_slice2 = (x == uniq[md]);
		new_slice2[ is.na(new_slice2) ] = 0;
		rez = vector("list", 2);
		rez[[1]] = new_slice1;
		rez[[2]] = new_slice2;
		return( rez );
}

.my.pmin = function(x, val) {
	# minimum "pmin" function that can handle empty array
	if(length(x) == 0) {
		return(x)
	} else {
		return(pmin.int(x,val));
	}	
}

.my.pmax = function(x, val) {
	# minimum "pmin" function that can handle empty array
	if(length(x) == 0) {
		return(x)
	} else {
		return(pmax.int(x,val));
	}	
}

.pv.nz = function(x){return( .my.pmax(x,.Machine$double.xmin) )}
 
.SetNanRowMean = function(x) {
	if( any(is.na(x)) ) {
		rowmean = rowMeans(x, na.rm = TRUE);
		rowmean[ is.na(rowmean) ] = 0;
		for( j in which(!complete.cases(x)) ) {
			where1 = is.na( x[j, ] );
			x[j,where1] = rowmean[j];
		}
	}
	return(x);
}

.impute_row_mean = function(x) {
	return( list(.SetNanRowMean(x)) );
}

.RowStandardizeCentered = function(x) {
	div = sqrt( rowSums(x^2) );
	div[ div == 0 ] = 1;
	return( x/div );
}

Matrix_eQTL_engine = function(
						snps, 
						gene, 
						cvrt = SlicedData$new(), 
						output_file_name, 
						pvOutputThreshold = 1e-5, 
						useModel = modelLINEAR, 
						errorCovariance = numeric(), 
						verbose = TRUE,
 						pvalue.hist = FALSE ) {
	rez = Matrix_eQTL_main(
				snps = snps, 
				gene = gene, 
				cvrt = cvrt, 
				output_file_name = output_file_name, 
				pvOutputThreshold = pvOutputThreshold,
				useModel = useModel, 
				errorCovariance = errorCovariance, 
				verbose = verbose,
 				pvalue.hist = pvalue.hist );
	return( rez );
}

.Matrix_eQTL_engine_cis = function(
						snps, 
						gene, 
						cvrt = SlicedData$new(), 
						output_file_name,
 						output_file_name.cis = paste(output_file_name,".cis.txt",sep=""),
						pvOutputThreshold.cis = 1e-3,
						pvOutputThreshold.tra = 1e-6,
						useModel = modelLINEAR, 
						errorCovariance = numeric(), 
						verbose = TRUE, 
						snpspos, 
						genepos,
						cisDist = 1e6,
						pvalue.hist = FALSE ) {
	rez = Matrix_eQTL_main(
				snps = snps, 
				gene = gene, 
				cvrt = cvrt, 
				output_file_name = output_file_name, 
				pvOutputThreshold = pvOutputThreshold.tra,
				useModel = useModel, 
				errorCovariance = errorCovariance, 
				verbose = verbose, 
				output_file_name.cis = output_file_name.cis, 
				pvOutputThreshold.cis = pvOutputThreshold.cis,
				snpspos = snpspos, 
				genepos = genepos,
				cisDist = cisDist,
 				pvalue.hist = pvalue.hist);
	return( rez );
}

Matrix_eQTL_main = function(	
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
 						pvalue.hist = FALSE) {
	################################# Basic variable checks #################################
 	{				
		# status("Performing basic checks of the input variables");
		stopifnot( "SlicedData" %in% class(gene) );
		stopifnot( "SlicedData" %in% class(snps) );
		stopifnot( "SlicedData" %in% class(cvrt) );
		
		# Check dimensions
		if( min(snps$nRows(),snps$nCols()) == 0 )
			stop("Empty genotype dataset");
		if( min(gene$nRows(),gene$nCols()) == 0 )
			stop("Empty expression dataset");
		if( snps$nCols() != gene$nCols() )
			stop("Different number of samples in the genotype and gene expression files");
		if( cvrt$nRows()>0 ) {
			if( snps$nCols() != cvrt$nCols() )
				stop("Wrong number of samples in the matrix of covariates");
		}

		stopifnot( class(pvOutputThreshold) == "numeric" );
		stopifnot( length(pvOutputThreshold) == 1 );
		stopifnot( pvOutputThreshold >= 0 );
		stopifnot( pvOutputThreshold <= 1 );
		if( pvOutputThreshold > 0 ) {
			stopifnot( class(output_file_name) == "character" );
			stopifnot( length(output_file_name) == 1 );
		}
		stopifnot( class(pvOutputThreshold.cis) == "numeric" );
		stopifnot( length(pvOutputThreshold.cis) == 1 );
		stopifnot( pvOutputThreshold.cis >= 0 );
		stopifnot( pvOutputThreshold.cis <= 1 );
		stopifnot( !((pvOutputThreshold > 0) & (pvOutputThreshold.cis > 0) & (pvOutputThreshold > pvOutputThreshold.cis)) );
		stopifnot( (pvOutputThreshold > 0) | (pvOutputThreshold.cis > 0) );
		
		stopifnot( class(useModel) == class(modelLINEAR) );
		stopifnot( length(useModel) == 1 );
		stopifnot( useModel %in% c(modelLINEAR, modelANOVA, modelLINEAR_CROSS) );
		if( useModel == modelLINEAR_CROSS ) {
			if( cvrt$nRows() == 0 ) {
				stop( "Model \"modelLINEAR_CROSS\" requires at least one covariate" );
			}
		}
		
		stopifnot( class(verbose) == "logical" );
		stopifnot( length(verbose) == 1 );
		
		if( pvOutputThreshold.cis > 0 ) {
			stopifnot( class(output_file_name.cis) == "character" );
			stopifnot( length(output_file_name.cis) == 1 );
			stopifnot( class(snpspos) == "data.frame" );
			stopifnot( ncol(snpspos) == 3 );
			stopifnot( class(genepos) == "data.frame" );
			stopifnot( ncol(genepos) == 4 );
			stopifnot( nzchar(output_file_name.cis) )
		}
		if( pvOutputThreshold > 0) {
			stopifnot( nzchar(output_file_name) )
		}
		
		stopifnot( class(errorCovariance) %in% c("numeric", "matrix") );
		errorCovariance = as.matrix(errorCovariance);
		if(length(errorCovariance)>0) {
			if( nrow(errorCovariance) != ncol(errorCovariance) ) {
				stop("The covariance matrix is not square");
			}	
			if( nrow(errorCovariance) != snps$nCols() ) {
				stop("The covariance matrix size does not match the number of samples");
			}
			if( !all(errorCovariance == t(errorCovariance)) ) {
				stop("The covariance matrix is not symmetric");
			}
		}
	}

	# preserve the expression and covariate classes from changing
	params = list(
		output_file_name = output_file_name, 
		pvOutputThreshold = pvOutputThreshold,
		useModel = useModel, 
		errorCovariance = errorCovariance , 
		verbose = verbose, 
		output_file_name.cis = output_file_name.cis, 
		pvOutputThreshold.cis = pvOutputThreshold.cis,
		cisDist = cisDist ,
 		pvalue.hist = pvalue.hist);

	gene = gene$Clone();
	#snps = snps$Clone();
	cvrt = cvrt$Clone();
		
	## the timing routing
	if( verbose ) {
		lastTime = 0;
		status <- function(text) {
			# gc();
			newTime = proc.time()[3];
			if(lastTime != 0) {
				cat("Task finished in ", newTime-lastTime, " seconds\n");
			}
			cat(text,"\n");
			lastTime <<- newTime;
			unused = flush.console();
		}
	} else {
		status = function(text){}
	}
	start.time = proc.time()[3];
	## set model-specific variables

	################################# Create the Saver class(es) #################################
	{
		status("Creating output file(s)");
		if( pvOutputThreshold > 0 ) {
			Saver = .OutputSaver_FRD$new();
		}
		if( pvOutputThreshold.cis > 0 ) {
			Saver.cis = .OutputSaver_FRD$new();
		}
		if( pvOutputThreshold > 0 ) {
			if( pvOutputThreshold * gene$nRows() * snps$nRows() > 1000000 ) {
				cat("Warning: pvOutputThreshold is not small enough.\n");
				if( pvOutputThreshold.cis > 0 ) {
					stop("Lower pvOutputThreshold or turn off cis analysis to proceed");
				} else {
					Saver = .OutputSaver_direct$new();
				}
			}
		}
		if( (useModel == modelLINEAR) || (useModel == modelLINEAR_CROSS) ){
			statistic_name = "t-stat";
		} else if( useModel == modelANOVA ) {
			statistic_name = "F-test";
		}
		if( pvOutputThreshold > 0 ) {
			Saver$Start(output_file_name, gene, snps, statistic_name);
		}
		if( pvOutputThreshold.cis > 0 ) {
			Saver.cis$Start(output_file_name.cis, gene, snps, statistic_name);
		}
		rm( statistic_name );
	}
	################################# Error covariance matirx processing #################################
	{
		if( length(errorCovariance) > 0 ) {
			status("Processing the error covariance matrix");
			eig = eigen(errorCovariance, symmetric = TRUE)
			d = eig$values;
			v = eig$vectors;
			#  errorCovariance == v %*% diag(d) %*% t(v)
			#  errorCovariance^0.5 == v*sqrt(d)*v" (no single quotes anymore)
			#  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v"
			if( any(d<=0) ) {
				stop("The covariance matrix is not positive definite");
			}
			correctionMatrix = v %*% diag(1./sqrt(d)) %*% t(v);
			rm( eig, v, d, errorCovariance )
		} else {
			rm( errorCovariance );
			correctionMatrix = numeric();
		}
	}
	################################# Covariates processing #################################
	{	
		status("Processing covariates");
		if( useModel == modelLINEAR_CROSS ) {
			last.covariate = tail( cvrt$getSlice(cvrt$nSlices()), n = 1);
		}		
		if( cvrt$nRows()>0 ) {
			cvrt$SetNanRowMean();
			cvrt$CombineInOneSlice();
			cvrt = rbind(matrix(1,1,snps$nCols()),cvrt$getSlice(1));
		} else {
			cvrt = matrix(1,1,snps$nCols());
		}
		# Correct for the error covariance structure
		if( length(correctionMatrix)>0 ) {
			cvrt = cvrt %*% correctionMatrix;
		}
		# Orthonormalize covariates
		# status("Orthonormalizing covariates");
		q = qr(t(cvrt));
		if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * snps$nCols() ) {
			stop("Colinear or zero covariates detected");
		}
		cvrt = t( qr.Q(q) );
		rm( q );
	}
	################################# Gene expression processing #################################
	{
		status("Processing gene expression data (imputation, residualization, etc.)");
		# Impute gene expression
		gene$SetNanRowMean();
		# Correct for the error covariance structure
		if( length(correctionMatrix)>0 ) {
			gene$RowMatrixMultiply(correctionMatrix);
		}
		# Orthogonolize expression w.r.t. covariates
		gene$RowStandardizeCentered();
		# status("Orthogonolizing expression w.r.t. covariates");
		for( sl in 1:gene$nSlices() ) {
			slice = gene$getSlice(sl);
			slice = slice - tcrossprod(slice,cvrt) %*% cvrt;
			gene$setSlice(sl, slice);
		}
		rm( sl, slice );
		gene$RowRemoveZeroEps();
		# status("Standardizing expression");
		gene$RowStandardizeCentered();
	}
	################################# matching gene and SNPs locations #################################
	if( pvOutputThreshold.cis > 0 ) {
		status("Matching data files and location files")
		# names in the input data	
		gene_names = gene$GetAllRowNames();
		snps_names = snps$GetAllRowNames();

		# match with the location data
		usedgene = matrix(FALSE, nrow(genepos),1); # genes in "genepos" that are matching  "gene_names"
		genematch = match( gene_names, genepos[ ,1],  nomatch = 0L);
		usedgene[ genematch ] = TRUE;
		if( !any(genematch) ) {
			stop("Gene names do not match those in the gene location file.");
		}
		cat( sum(genematch>0), "of", length(gene_names), " genes matched\n");
		
		usedsnps = matrix(FALSE, nrow(snpspos),1);
		snpsmatch = match( snps_names, snpspos[ ,1],  nomatch = 0L);
		usedsnps[ snpsmatch ] = TRUE;
		if( !any(snpsmatch) ) {
			stop("SNP names do not match those in the SNP location file.");
		}
		cat( sum(snpsmatch>0), "of", length(snps_names), " SNPs matched\n");
		
		# list used chr names
		chrNames = unique(c( as.character(unique(snpspos[usedsnps,2])), as.character(unique(genepos[usedgene,2])) ))
		chrNames = chrNames[ sort.list( suppressWarnings(as.integer(chrNames)), method = "radix", na.last = TRUE ) ];
		# match chr names
		genechr = match(genepos[,2],chrNames);
		snpschr = match(snpspos[,2],chrNames);
		
		# max length of a chromosome
		chrMax = max( snpspos[usedsnps, 3], genepos[usedgene, 4], na.rm = TRUE) + cisDist;
		
		# Single number location for all rows in "genepos" and "snpspos"
 		genepos2 = as.matrix(genepos[ ,3:4, drop = FALSE] + (genechr-1)*chrMax);
 		snpspos2 = as.matrix(snpspos[ ,3  , drop = FALSE] + (snpschr-1)*chrMax);
		
		# the final location arrays;
		snps_pos = matrix(0,length(snps_names),1);
		snps_pos[snpsmatch>0, ] = snpspos2[snpsmatch, , drop = FALSE];
		snps_pos[rowSums(is.na(snps_pos))>0, ] = 0;
		snps_pos[snps_pos==0] = (length(chrNames)+1) * (chrMax+cisDist);
		rm(snps_names, snpsmatch, usedsnps, snpschr, snpspos2)
		
		gene_pos = matrix(0,length(gene_names),2);
		gene_pos[genematch>0, ] = genepos2[genematch, , drop = FALSE];
		gene_pos[rowSums(is.na(gene_pos))>0, ] = 0;
		gene_pos[gene_pos==0] = (length(chrNames)+2) * (chrMax+cisDist);
		rm(gene_names, genematch, usedgene, genechr, genepos2)
		rm(chrNames, chrMax);

		if( is.unsorted(snps_pos) ) {
			status("Reordering SNPs\n");
			ordr = sort.list(snps_pos);
			snps$RowReorder(ordr);
			snps_pos = snps_pos[ordr, , drop = FALSE];
			rm(ordr);
		}
		if( is.unsorted(rowSums(gene_pos)) ) {
			status("Reordering genes\n");
			ordr = sort.list(rowSums(gene_pos));
			gene$RowReorder(ordr);
			gene_pos = gene_pos[ordr, , drop = FALSE];
			rm(ordr);
		}
		
		# Slice it back.
		geneloc = vector("list", gene$nSlices())
		gene_offset = 0;
		for(gc in 1:gene$nSlices()) {
			nr = length(gene$rowNameSlices[[gc]]);
			geneloc[[gc]] = gene_pos[gene_offset + (1:nr), , drop = FALSE];
			gene_offset = gene_offset + nr;	
		}
		rm(gc, gene_offset, gene_pos);
		
		snpsloc = vector("list", snps$nSlices())
		snps_offset = 0;
		for(sc in 1:snps$nSlices()) {
			nr = length(snps$rowNameSlices[[sc]]);
			snpsloc[[sc]] = snps_pos[snps_offset + (1:nr), , drop = FALSE];
			snps_offset = snps_offset + nr;	
		}
		rm(nr, sc, snps_offset, snps_pos);
	}
	################################# Prepare for main loop    #################################
	
	nSamples = snps$nCols();
	nGenes = gene$nRows();
	nSnps  = snps$nRows();
	nCov = nrow(cvrt);
	# nVarTested = length(snps_list); # set in case(useModel)
	# dfNull = nSamples - nCov;
	# d.f. of the full model
	
	if( useModel == modelLINEAR ) {
		snps_process = .impute_row_mean;
		nVarTested = 1;
		dfFull = nSamples - nCov - nVarTested;
		statistic.fun = function(mat_list) {
			return( mat_list[[1]] );
		}
		afun = function(x) {return(abs(x))};
		threshfun = function(pv) {
			thr = qt(pv/2, dfFull, lower.tail = FALSE);
			thr = thr^2;
			thr = sqrt(  thr / (dfFull + thr) );
			thr[pv >= 1] = 0;
			thr[pv <= 0] = 1;
			return( thr );
		}
		testfun = function(x) { return( x * sqrt( dfFull / (1 - .my.pmin(x^2,1))));	}
		pvfun = function(x) { return( .pv.nz(pt(-abs(x),dfFull)*2)); }
		thresh.cis = threshfun(pvOutputThreshold.cis);
		thresh = threshfun(pvOutputThreshold);
	} else if( useModel == modelANOVA ) {
		snps_process = .snps_split_for_ANOVA;
		nVarTested = 2;
		dfFull = nSamples - nCov - nVarTested;
		statistic.fun = function(mat_list) {
			return( mat_list[[1]]^2 + mat_list[[2]]^2 );
		}
		afun = identity;
		threshfun = function(pv) {
			thr = qf(pv, nVarTested, dfFull, lower.tail = FALSE);
			thr = thr / (dfFull/nVarTested + thr);
			thr[pv >= 1] = 0;
			thr[pv <= 0] = 1;
			return( thr );
		}
		testfun = function(x) { return( x / (1 - .my.pmin(x,1)) * (dfFull/nVarTested) ); }
		pvfun = function(x) { return( .pv.nz(pf(x, nVarTested, dfFull, lower.tail = FALSE)) ); }
		thresh.cis = threshfun(pvOutputThreshold.cis);
		thresh = threshfun(pvOutputThreshold);
	} else if( useModel == modelLINEAR_CROSS ) {
		last.covariate = as.vector( last.covariate );
		snps_process = function(x) {
			out = vector("list", 2);
			out[[1]] = .SetNanRowMean(x);
			out[[2]] = t( t(out[[1]]) * last.covariate );
			return( out );
		}
		nVarTested = 1;
		dfFull = nSamples - nCov - nVarTested - 1;
		statistic.fun = function(mat_list) {
			return( mat_list[[2]] / sqrt(1 - mat_list[[1]]^2) );
		}
		afun = function(x) {return(abs(x))};
		threshfun = function(pv) {
			thr = qt(pv/2, dfFull, lower.tail = FALSE);
			thr = thr^2;
			thr = sqrt(  thr / (dfFull + thr) );
			thr[pv >= 1] = 0;
			thr[pv <= 0] = 1;
			return( thr );
		}
		testfun = function(x) { return( x * sqrt( dfFull / (1 - .my.pmin(x^2,1))));	}
		pvfun = function(x) { return( .pv.nz(pt(-abs(x),dfFull)*2 )); }		
		thresh.cis = threshfun(pvOutputThreshold.cis);
		thresh = threshfun(pvOutputThreshold);				
	}

	################################# Some useful functions #####################
	
	orthonormalize.snps = function(cursnps) {
		for(p in 1:length(cursnps)) {
			if(length(correctionMatrix)>0) {
				cursnps[[p]] = cursnps[[p]] %*% correctionMatrix;
			}
			cursnps[[p]] = cursnps[[p]] - tcrossprod(cursnps[[p]],cvrt) %*% cvrt;
			if( p>1 ) {
				for(w in 1:(p-1))
					cursnps[[p]] = cursnps[[p]] - rowSums(cursnps[[p]]*cursnps[[w]]) * cursnps[[w]]
			}
			cursnps[[p]] = .RowStandardizeCentered(cursnps[[p]]);
		}
		return(cursnps);
	}
	is.cis.pair = function(gg,ss) {
		return(!( ( snpsloc[[ss]][1, 1] - tail( geneloc[[gg]][ , 2], n = 1L) > cisDist) |
			    ( geneloc[[gg]][1, 1] - tail( snpsloc[[ss]]      , n = 1L) > cisDist) ) );
	}
	
	################################# Prepare counters and histogram bins ##########################
	pvbins = NULL; # bin edges for p-values
	statbins = 0;  # bin edges for the test statistic (|t| or F)
	do.hist = FALSE;
	if( length(pvalue.hist) == 1 ) {
		if(pvalue.hist == "qqplot") {
			lowlim = c(pvOutputThreshold,pvOutputThreshold.cis);
			lowlim = lowlim[lowlim >0];
			lowlim = min(lowlim);
			pvbins = c(0, 10^seq(log10(lowlim)-1, 0, 0.05));
			rm(lowlim);
		} else
		if( is.numeric(pvalue.hist) ) {
			pvbins = seq(from = 0, to = 1, length.out = pvalue.hist+1);
		} else
		if( pvalue.hist == TRUE ) {
			pvbins = seq(from = 0, to = 1, length.out = 100+1);
		}
	} else
	if( is.numeric(pvalue.hist) && (length(pvalue.hist) > 1) ) {
		pvbins = pvalue.hist;
	}
	if( is.null(pvbins) && (pvalue.hist != FALSE) ) {
		stop("Wrong value of pvalue.hist. Must be FALSE, TRUE, \"qqplot\", or numerical");
	}
	do.hist = !is.null(pvbins);
	if( do.hist ) {
		pvbins = sort(pvbins);
		statbins = threshfun(pvbins);
	}
	if( pvOutputThreshold > 0) {
		counter     = .eQTL_summary$new(snps, gene, pvbins, statbins);
	}
	if( pvOutputThreshold.cis > 0) {
		counter.cis = .eQTL_summary$new(snps, gene, pvbins, statbins);
	}
	rm( pvbins, statbins);

	################################# Main loop #################################

	status("Performing eQTL analysis");
	
	snps_offset = 0;
	for(ss in 1:snps$nSlices()) {
	# for(ss in 1:min(15,snps$nSlices())) { #for debug
		cursnps = NULL;
		nrcs = nrow(snps$getSliceRaw(ss));
		
		gene_offset = 0;
		for(gg in 1:gene$nSlices()) {
			curgene = gene$getSlice(gg);
			nrcg = nrow(curgene);
			
			rp = "";
			
			statistic = NULL;
			## do trans/all analysis
			if(pvOutputThreshold>0) {
				if( is.null(cursnps) ) {
					cursnps = orthonormalize.snps( snps_process( snps$getSlice(ss) ) );
				}
				mat = vector("list", length(cursnps));
				for(d in 1:length(cursnps)) {
					mat[[d]] = tcrossprod(curgene, cursnps[[d]]);
				}
				statistic = statistic.fun( mat );
				#rm(mat);
				select = which( afun(statistic) >= thresh, arr.ind = TRUE );
				statistic.select = statistic[ select ];
				test = testfun( statistic.select );
				pv = pvfun( test );
				Saver$WriteBlock( cbind( snps_offset + select[ , 2], gene_offset + select[ , 1], test, pv) );
				counter$Update(gg, ss, select, pv, n.tests = nrcs*nrcg, if(do.hist) afun(statistic) )
				rp = paste(rp, ", ", formatC(counter$neqtls, big.mark=",", format = "f", digits = 0), " eQTLs found", sep = "")
			}
			## do cis analysis
			if( (pvOutputThreshold.cis > 0) && ( is.cis.pair(gg, ss) ) ) {
				if( is.null(cursnps) ) {
					cursnps = orthonormalize.snps( snps_process( snps$getSlice(ss) ) );
				}
				if( is.null( statistic ) ) {
					mat = vector("list", length(cursnps));
					for(d in 1:length(cursnps)) {
						mat[[d]] = tcrossprod(curgene, cursnps[[d]]);
					}
					statistic = statistic.fun( mat );
				}
				
				srep = rep( snpsloc[[ss]], each = nrcg);
				is.cis = ( (geneloc[[gg]][ ,1] - srep) < cisDist) &
						( (srep - geneloc[[gg]][ ,2]) < cisDist);
				#rm(srep);
				dim(is.cis) = c(nrcg, nrcs);
				select = which( (afun(statistic) >= thresh.cis) &  is.cis, arr.ind=TRUE );

				statistic.select  = statistic[ select ];
				#rm( statistic );
				test = testfun( statistic.select );
				pv = pvfun(test);
				Saver.cis$WriteBlock( cbind(snps_offset + select[ , 2], gene_offset + select[ , 1], test, pv) );
				counter.cis$Update(gg, ss, select, pv, n.tests = sum(is.cis), if(do.hist) afun(statistic[is.cis]) )
				rp = paste(rp, ", ", formatC(counter.cis$neqtls, big.mark=",", format = "f", digits = 0), " cis-eQTLs found", sep = "")
			}
			gene_offset = gene_offset + nrcg;
			per = 100*(gg/gene$nSlices() + ss-1) / snps$nSlices();
			if( !is.null(statistic) ) {
				cat( formatC(per, format = "f", width = 5, digits = 2), "% done" , rp, "\n", sep = "");
			}
			flush.console();
		}
		snps_offset = snps_offset + nrow(cursnps[[1]]);
	}
	rez = list(time.in.sec = proc.time()[3] - start.time);
	rez$param = params;
	
	if( pvOutputThreshold > 0 ) {
		rez$all = counter$getResults();
		
		if( (pvOutputThreshold > 0) & (pvOutputThreshold.cis > 0) ) {
			Saver$RemoveCis(Saver.cis);
			rez$trans$eqtls = Saver$Commit(counter$ntests);
		} else {
			rez$all$eqtls = Saver$Commit(counter$ntests);
		}
	}
	if( pvOutputThreshold.cis > 0 ) {
		eqtls = Saver.cis$Commit(counter.cis$ntests);
		rez$cis = counter.cis$getResults();
		if( !is.null(eqtls) ) {
			rez$cis$eqtls = eqtls;
		}
		rm(eqtls);
	}
	if( (pvOutputThreshold > 0) && (pvOutputThreshold.cis > 0) ) {
		rez$trans$hist.counts = rez$all$hist.counts - rez$cis$hist.counts;
		rez$trans$hist.bins = rez$all$hist.bins;
		rez$trans$ntests = rez$all$ntests - rez$cis$ntests;
	}
	class(rez) = "MatrixEQTL";
	status("");
	return(rez);
}

.eQTL_summary <- setRefClass(".eQTL_summary",
	fields = list(
		dataEnv = "environment",
		ntests = "numeric",
		neqtls = "numeric",
		snSlices = "numeric",
		gnSlices = "numeric",
		pvbins1 = "numeric",
		statbins1 = "numeric",
		hist.count = "numeric"
	),
	methods = list(
	initialize = function (snps, gene,  pvbins, statbins) {
		dataEnv <<- new.env(hash = TRUE);
		for( ss in 1:snps$nSlices() ) {
			assign(paste( ss), rep(0L, length(snps$rowNameSlices[[ss]])), dataEnv );
		}
		snSlices <<- snps$nSlices();
		for( gg in 1:gene$nSlices() ) {
			assign(paste(-gg), rep(0L, length(gene$rowNameSlices[[gg]])), dataEnv );
		}
		gnSlices <<- gene$nSlices();
		ntests <<- 0;
		neqtls <<- 0;
		if(length(pvbins)) {
			pvbins1 <<- rev(pvbins);
			statbins1 <<- rev(statbins);
			hist.count <<- double(length(pvbins)-1);
		} else {
			pvbins1 <<- numeric(0);
		}
		return(.self);
	},
	Update = function(gg, ss, select, pv, n.tests, stats.for.hist) {
		neqtls <<- neqtls + length(pv);
		ntests <<- ntests + n.tests;
		if( length(pv) > 0 ) {
			snpseqtl = get( paste( ss), dataEnv );
			snpseqtl = snpseqtl + tabulate(select[ , 2], nbins = length(snpseqtl)); #rowSums(select.mat);
			assign( paste( ss), snpseqtl, dataEnv );
			geneeqtl = get( paste(-gg), dataEnv );
			geneeqtl = geneeqtl + tabulate(select[ , 1], nbins = length(geneeqtl)); #colSums(select.mat);
			assign( paste(-gg), geneeqtl, dataEnv );
		}
		if( !is.null(stats.for.hist) ) {
			 h = counts <- .C("bincount",
					stats.for.hist, as.integer(length(stats.for.hist)),
					statbins1, as.integer(length(statbins1)),
					counts = integer(length(statbins1) - 1), right = as.logical(TRUE),
					include = as.logical(TRUE), naok = FALSE, NAOK = FALSE,
					DUP = FALSE, PACKAGE = "base")$counts;
			hist.count <<- hist.count + h;
		}
	},
	getResults = function() {
		gene.with.eqtls = 0;
		for( gg in 1:gnSlices ) {
			geneeqtl = get( paste(-gg), dataEnv );
			gene.with.eqtls = gene.with.eqtls + sum(geneeqtl>0);
		}
		rm( gg, geneeqtl); 
		snps.with.eqtls = 0;
		for( ss in 1:snSlices ) {
			snpseqtl = get( paste( ss), dataEnv );
			snps.with.eqtls = snps.with.eqtls + sum(snpseqtl>0);
		}
		rm( ss, snpseqtl);
		rez = list( ntests = ntests, neqtls = neqtls, snps.with.eqtls = snps.with.eqtls, gene.with.eqtls = gene.with.eqtls);
		if( length(pvbins1) > 0 ) {
			rez$hist.bins = rev(pvbins1);
			rez$hist.counts = rev(hist.count);
		}
		return(rez);
	}
	)
)

.histme = function(m, name1, name2, ...) {
	cnts = m$hist.counts;
	bins = m$hist.bins;
	ntst = m$ntests;
	centers = 0.5 * (bins[-1L] + bins[-length(bins)]);
	density = 0.5 / (bins[-1L] - centers) * cnts / ntst;
	ntext = paste("P-value distribution for ", name1, formatC(ntst, big.mark=",", format = "f", digits = 0), name2, " gene-SNP pairs ",sep="");
	r = structure(list(breaks = bins, counts = cnts, density = density,
	      mids = centers, equidist = FALSE), class = "histogram");
	plot(r, main = ntext, ylab = "Density", xlab = "P-values", ...)
	abline( h = 1, col = "blue");
	return(invisible());
}

.qqme = function(m, lcol, cex, pch, ...) {
	cnts = m$hist.counts;
	bins = m$hist.bins;
	ntst = m$ntests;
	ypvs = m$eqtls$pvalue;
	xpvs = if(length(ypvs)>0) {1:length(ypvs) / ntst;} else {numeric();}

	if(length(ypvs) > 1000) {
		# need to filter a bit, make the plotting faster
		levels = log(xpvs);
		levels = (levels - levels[1])/(tail(levels,1) - levels[1]);
		levels = floor(levels *2000);
		keep = c(TRUE, levels[-1] != levels[-length(levels)]);
		ypvs = ypvs[keep];
		xpvs = xpvs[keep];
	}
	cusu = cumsum(cnts) / ntst;
	ypos = bins[-1][is.finite(cusu)];
	xpos = cusu[is.finite(cusu)];
	points(xpvs, ypvs, col = lcol, pch = pch, cex = cex, ...);
	lines(xpos, ypos, col = lcol, ...);
}

plot.MatrixEQTL = function(x, cex = 0.5, pch = 19, ...) {
	if( x$param$pvalue.hist == FALSE ) {
		warning("Cannot plot p-value distribution: the information was not recorded.\nUse pvalue.hist!=FALSE.");
		return(invisible());
	}
	if( x$param$pvalue.hist == "qqplot" ) {
		xmin = 1/max(x$cis$ntests, x$all$ntests);
		ymin = min( x$cis$eqtls$pvalue[1], x$cis$hist.bins[c(FALSE,x$cis$hist.counts>0)][1],
				x$all$eqtls$pvalue[1], x$all$hist.bins[c(FALSE,x$all$hist.counts>0)][1],
				x$trans$eqtls$pvalue[1], x$trans$hist.bins[c(FALSE,x$trans$hist.counts>0)][1],
				na.rm = TRUE)/1.5;
		if(ymin == 0) {
			ymin == .Machine$double.xmin
		}
		plot(numeric(),numeric(), xlab = "p-value, theoretical",
			ylab = "p-value, actual",
			log = "xy",
			xlim = c(1, xmin/1.5),
			ylim = c(1, ymin),
			xaxs="i", yaxs="i", ...);
		lines(c(1,xmin), c(1,xmin), col = "gray");
		if((x$param$pvOutputThreshold > 0) && (x$param$pvOutputThreshold.cis > 0)) {
			.qqme( x$cis, "red", cex, pch, ...);
			.qqme( x$trans, "blue", cex, pch, ...);
			title(paste("Q-Q plot for",
				formatC(x$cis$ntests, big.mark=",", format = "f", digits = 0),
				"local and",
				formatC(x$trans$ntests, big.mark=",", format = "f", digits = 0),
				"distant gene-SNP p-values"));
			lset = c(1,2,4);
		} else
		if(x$param$pvOutputThreshold.cis > 0) {
			.qqme(x$cis, "red", cex, pch, ...);
			title(paste("Q-Q plot for",
				formatC(x$cis$ntests, big.mark=",", format = "f", digits = 0),
				"local gene-SNP p-values"));
			lset = c(1,4);
		} else {
			.qqme(x$all, "blue", cex, pch, ...);
			title(paste("Q-Q plot for all",
				formatC(x$all$ntests, big.mark=",", format = "f", digits = 0),
				"gene-SNP p-values"));
			lset = c(3,4);
		}
		legend("topleft",
			c("Local p-values","Distant p-values","All p-values","diagonal")[lset],
			col =      c("red","blue","blue","gray")[lset],
			text.col = c("red","blue","blue","gray")[lset],
			pch = 20, lwd = 1, pt.cex = c(1,1,1,0)[lset])
	} else {
		if((x$param$pvOutputThreshold > 0) && (x$param$pvOutputThreshold.cis > 0)) {
			par(mfrow=c(2,1));
			.histme(x$cis, "", " local", ...);
			tran = list(hist.counts = x$all$hist.counts - x$cis$hist.counts,
					hist.bins = x$all$hist.bins,
					ntests =  x$all$ntests - x$cis$ntests);
			.histme(x$trans,""," distant", ...);
			par(mfrow=c(1,1));
		} else
		if(x$param$pvOutputThreshold.cis > 0) {
			.histme(x$cis, "", " local", ...);
		} else {
			.histme(x$all, "all ", ""  , ...);
		}
	}
	return(invisible());
}

