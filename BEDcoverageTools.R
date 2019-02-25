library("Biostrings");
library("BSgenome");
library("rtracklayer");
library("Rsamtools");
library("ShortRead");
library("chipseq")
library("data.table")

#global constants
gSigDigits = 3;
gNormReads = 1e6

#global plot parameters
gNorm = "rpm"
gheight = 1.5;
gwidth = 1.5
gmai = c(0.3, 0.3, 0.1, 0.1);
gmgp = c(0.6, 0.1, 0);
gtcl = -0.2
gcex = 0.8
gcex.lab = 0.6
gcex.axis = 0.6
glwd = 0.5
glas = 1;
gbty = "n"
gNorm = "rpm"

#global Hist constants
gFragSize = NA;
gHistBin = 10;
gHistOutFile = "paste(bamFiles, '.', gsub('^.*/', '', bedFile), '.range', range, '.bin', bin, '.csv', sep = '')"
glim = 50

#function to make histogram from bam files from a bed file of coordinates
makeBamHist = function(bamFiles, bed, range, outFile = c(gHistOutFile), bin = c(gHistBin), fragSizes = c(gFragSize), normReads = c(gNormReads), sigDigits = c(gSigDigits), sampleNames = c(NA), flag = c(NA), readCounts = c(NA), ...) {

	#Debug
	#bamFiles = bamFile; bed = motifr; range = 100; outFile = NA; bin = 1; fragSizes = c(gFragSize); normReads = gNormReads; sampleNames = NA; sigDigits = gSigDigits; flag = NA; readCounts = NA;
	#normReads = 1e6 / files$frip[i]; fragSizes = 1; outFile = statFile

	#Determine mean width of bed range
	meanWidth = round(mean(width(bed)))

	#binned region to calculate enrichment
	bins = seq(from = -range, to = range + meanWidth, by = bin)

	#Create bin labels, normalize distance within bed region 
	binlabs = bins
	binlabs[binlabs > 0 & binlabs <= meanWidth] = binlabs[binlabs > 0 & binlabs <= meanWidth] / meanWidth
	binlabs[binlabs > meanWidth] = binlabs[binlabs > meanWidth] - meanWidth + 1

	if (!is.null(bed@elementMetadata$name) & length(bed@elementMetadata$name) == length(unique(bed@elementMetadata$name))) {
		histRows = bed@elementMetadata$name
	} else histRows = paste(as.character(bed@seqnames), start(bed), sep = "_");

	if (outFile == gHistOutFile)  outFile = eval(parse(text = outFile))

	#bedCenters = round(apply(cbind(start(bed), end(bed)), MARGIN = 1, mean))

	#Make hist matrix index
	histIdx = matrix(bins, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);
	histIdx = histIdx + start(bed)
	histIdx[histIdx <= 0] = NA

	#make hist matrix
	hist = matrix(0, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);

	hists = list();
	for (i in 1:length(bamFiles)) hists[[i]] = hist;

	if (all(is.na(readCounts))) {
		readCounts = list();
		for (i in 1:length(bamFiles)) readCounts[[i]] = countBam(bamFiles[i], param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isDup = F)))$records
	}
	
	for (chr in levels(bed@seqnames)) {

		print(paste("chr", chr, Sys.time()))

		#Recursively build coverage for each bamFile 
		for (i in 1:length(bamFiles)) {

			si = seqinfo(BamFile(bamFiles[i]))	

			#Check to make sure the bam file has reads on chromosome chr 
			if (chr %in% unique(getBamChrs(as.character(bamFiles[i])))) {

				bedChr = bed[bed@seqnames == chr];
				
				#Read in chromosome from each bam file, extend reads to fragment size
				sbp = ScanBamParam(which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))

				if (is.na(fragSizes[i])) {
					bam = granges(readGAlignments(as.character(bamFiles[i]), param = sbp));
				} else 	bam = suppressWarnings(resize(granges(readGAlignments(as.character(bamFiles[i]), param = sbp)), width = fragSizes[i], fix = "start"));

				#Calculate base level coverage
				cv = coverage(bam)[[chr]];
		
				#Calculate the running mean to put into bins
				means = runmean(cv, k = bin)

				#Get the names of rows for chromosome chr
				chrRows = rownames(hists[[i]])[as.character(seqnames(bed)) == chr];

				#Replace regions that are larger than than the genome with NA
				hists[[i]][chrRows, ][hists[[i]][chrRows, ] > si@seqlengths[si@seqnames == chr]] = NA

				#Replace hist rows with mean
				hists[[i]][chrRows, ] = as.numeric(means)[histIdx[chrRows, ]]

				rm(bam); rm(cv); rm(means); gc();
			}
		}
	}

	#add up hists for each bam file
	histSum = NA;
	for (i in 1:length(bamFiles)) if (all(is.na(histSum))) histSum = hists[[i]] else histSum = histSum + hists[[i]];
		
	#Adjust for strandness
	histSum[as.character(bed@strand) == "-", ] = histSum[as.character(bed@strand) == "-", rev(seq_len(ncol(histSum)))];

	colnames(histSum) = binlabs; #label bins correctly	
	
	#Normalize and round
	histAvg = histSum * normReads / sum(unlist(readCounts));
	histAvg = cbind(row = row.names(histAvg), round(histAvg, sigDigits))
	
	if (!is.na(outFile)) write.table(histAvg, file = outFile, sep = ",", row.names = F, quote = F) 
	invisible(histSum)

	#rm(histSum); rm(histAvg); gc();
}

#Function to return chromosomes included in a Bam file
getBamChrs = function(bamFile) return(names(scanBamHeader(bamFile)[[1]][1]$targets));

#Function to return chromosome lengths included in a Bam file
getBamChrLengths = function(bamFile) {

	#Debug: bamFile = "";

	h = scanBamHeader(bamFile)
	cl = h[[1]][2]$text[names(h[[1]][2]$text) == "@SQ"]
	
	chrs = gsub("SN:", "", unlist(cl)[grepl("SN:", unlist(cl))]);
	len = as.numeric(gsub("LN:", "", unlist(cl)[grepl("LN:", unlist(cl))]));
	names(len) = chrs;

	return(len);
}


#Function to summarize makeBamHist output br bp for histogram plotting
sumBamHistCol = function(bedName, manifest, range, outDir, sumFile) {

	#Debug: bedName = bedName; manifest = files; range = 100, outDir = homeDir;

	print(paste("Summarizing", bedName, Sys.time()))

	#summarize data by row (bp)
	
	stats = list()	
	for (j in 1:dim(manifest)[1]) {
		print(paste(manifest$sample[j], Sys.time()))

		outFile = paste0(outDir, manifest$sample[j], ".", bedName, ".range.", range, "bp.csv")

		if (file.exists(outFile)) {
			stat = fread(outFile, header = T)
			stats[[j]] = colMeans(subset(stat, , colnames(stat)[colnames(stat) != "row"]), na.rm = T)
		}
	}

	#convert to data frame and label columns
	sum = data.frame(stats)
	colnames(sum) = paste0(manifest$sample, ".", bedName)
	sum = cbind(rownames(sum), sum)
	colnames(sum)[1] = "Distance"
	
	#write file
	write.table(sum, sumFile, sep = "\t", quote = F, row.names = F)

}

#Function to plot histogram from sumBamHistCol data
plotBedHist = function(manifest, sumFile = c(NA), lim = c(glim), norm = c(gNorm), colFun = c(NA), ylim = c(NA), height = c(gheight), width = c(gwidth), mai = c(gmai), mgp = c(gmgp), tcl = c(gtcl), cex = c(gcex), cex.lab = c(gcex.lab), cex.axis = c(gcex.axis), lwd = c(glwd), las = c(glas), bty = c(gbty), ...) {

	#Debug: bedName = bedName; manifest = files; lim = 50; outDir = outDir;

	#read in summary file and limit to range set above
	sum = read.table(file = sumFile, header = T, sep = "\t")
	sum.lim = sum[which(sum$Distance < lim+1 & sum$Distance > -lim-1), ]

	#get max y-axis values
	sum.plot = sum.lim[!grepl(c("Distance"), names(sum.lim))]
	ymax = max(unlist(lapply(sum.plot, max)))

	#plot by cellType	
	cairo_pdf(paste0(outDir, bedName,".", lim, "bpLimits.hist.pdf"), height = height, width = width)
	par(mai = mai, mgp = mgp, tcl = tcl, bty = bty)
	plot(NA, xlim = c(0, dim(sum.lim)[1]), ylim = c(0, ymax), cex.lab = cex.lab, cex.axis = cex.axis, ylab = norm, xlab = NA, xaxt = "n")

	#cycle through all sample groups and plot
	for (group in unique(manifest$group)) {

		print(paste(group, Sys.time()));

		#map group to color
		group.col <- unlist(lapply(group, colFun))

		#subset data by group
		filesGrp = manifest[manifest$group == group, ]
		sum.subset = sum.lim[,grepl(paste0(filesGrp$sample, collapse = "|"), names(sum.lim))]

		#find row means
		ave = rowMeans(sum.subset)
		
		#add line to plot
		lines(ave, col = group.col, lwd = lwd)
		}

	#add axis labels	
	axis.at = which(abs(as.numeric(sum.lim$Distance)) == lim | as.numeric(sum.lim$Distance) == 0 | as.numeric(sum.lim$Distance) == 1)
	axis(side = 1, at = axis.at, labels = sum.lim$Distance[axis.at], cex.axis = cex.axis, lwd = lwd)
	dev.off()
}

#Function to summarize makeBamHist output br bp for histogram plotting
sumBamHistRow = function(bedName, manifest, range, lim = c(glim), outDir, sampFile, grpFile, SumGrP = c(FALSE)) {

	#Debug: bedName = bedName; manifest = files; range = 100, lim = 50; outDir = homeDir;

	#combine coverage data for all samples and generate summary file by row for each bed feature
	print(paste("Summarizing", bedName, Sys.time()))

	#summarize data for each sample
	stats = list()
	for (j in 1:dim(manifest)[1]) {
		print(paste(manifest$sample[j], Sys.time()))

		outFile = paste0(outDir, manifest$sample[j], ".", bedName, ".range.", range, "bp.csv")

		if (file.exists(outFile)) {
			stat = fread(outFile, header = T)
		
			#subset stat file to match base pair range set with lim variable above
			sub = as.numeric(colnames(stat)) < lim+1 & as.numeric(colnames(stat)) > -lim-1 & colnames(stat) != "row"
			stat.lim = data.frame(stat)[, sub]
		
			#calculate row sums
			stats[[j]] = rowSums(stat.lim, na.rm = T)
		}
	}
	#convert to data frame and label columns
	stats = data.frame(stats)
	colnames(stats) = paste0(manifest$sample, ".", bedName)

	#write file
	write.table(stats, sampFile, sep = "\t", quote = F, row.names = F)

	if (SumGrP) {
		#summarize data by group
		grps = unique(manifest$group)

		#specify group order..work in progress
		#grpOrder = c("HC_rN", "SLE_rN", "HC_T3", "SLE_T3", "HC_SM", "SLE_SM", "HC_aN", "SLE_aN", "HC_DN", "SLE_DN")
		#set order
		#grps = factor(grpOrder, levels = grpOrder, ordered = T)

		#summarize data for each sample group
		means = list();
		for (i in 1:length(grps)) {

			print(paste(grps[i], Sys.time()))

			filesGrp = manifest[manifest$group == grps[i], ]
			#subset and take row means for each group
			means[[i]] = rowMeans(stats[, grepl(paste0(filesGrp$sample, collapse = "|"), names(stats))], na.rm = T)
		}

		#convert to data frame and label columns
		means = data.frame(means)
		colnames(means) = paste0(grps, ".", bedName)

		#write file
		write.table(means, grpFile, sep = "\t", quote = F, row.names = F)
	}
}

#Function to plot histogram from sumBamHistCol data
plotBedBox = function(manifest, grpFile, norm = c(gNorm), colFun = c(NA), CalcStats = c(TRUE), ylim = c(NA), height = c(gheight), width = c(gwidth), mai = c(0.5, 0.4, 0.1, 0.1), mgp = c(1.2, 0.2, 0), tcl = c(gtcl), cex = c(gcex), cex.lab = c(gcex.lab), cex.axis = c(1), lwd = c(glwd), las = c(3), bty = c(gbty), order=c(NA), ...) {

	#establish group order for boxplot
	if(!is.na(order[1])){
		grps = order
	} else {
		grps = unique(manifest$group)
	}

	#establish color scheme
	cols = unlist(lapply(grps, colFun))

	#read in data
	means = read.table(grpFile, sep = "\t", header = T)

	#plot by boxplot	
	cairo_pdf(gsub(".txt", ".boxPlot.pdf", grpFile), height = height, width = width)
	par(mai = mai, mgp = mgp, tcl = tcl, bty = bty)
	if(!is.na(order[1])) means=means[,pmatch(order, colnames(means))]
	boxplot(means, notch = T, outline = F, col = cols, names = grps, las = las, ylim = ylim, ylab = norm, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, lwd = lwd)
	dev.off()

	#calc stats
	if (CalcStats) {

		#establish groups for comparisons
		grps = unique(manifest$group)
		grpPairs = combn(grps, 2)
	
		#establish data frame to capture stats
		stats = data.frame(Comp = paste0(grpPairs[2,], ".v.", grpPairs[1,]))

		#loop through and perform Ttest for all comparisons
		for (i in 1:dim(grpPairs)[2]) {

			sig = t.test(means[,match(paste0(grpPairs[2,i], ".", bedName), colnames(means))], means[,match(paste0(grpPairs[1,i], ".", bedName), colnames(means))])
			stats$Ttest.Pval[i] = sig$p.value
		}
		
		#write file
		write.table(stats, gsub(".txt", ".Ttest.Stats.txt", grpFile), sep = "\t", quote = F, row.names = F)

	}
}