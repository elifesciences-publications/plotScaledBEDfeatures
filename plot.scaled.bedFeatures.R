#Plot boxPlots or Histogram read depth coverage of a bed file
#Written primarily for Transcription Factor motifs but can be used for any bed data

#set working directory
homeDir = "~/TF.footprinting/";
setwd(homeDir);

#load histogram Function
source("BEDcoverageTools.R")

#bed Name and location of bed file
bedName = "IRF4"

#uses a HOMER generated bed file of locations/motifs such as generated from the annotatePeaks.pl command.  Example below
#annotatePeaks.pl ATAC.peaks.bed mm10 -size given -noann -m irf4.motif -mbed IRF4.motifs.bed
bedFile = "IRF4.motifs.bed"

#read in the sample manifest file
manifestDir = "~/TF.footprinting/" #set directory if not current
manifestFile = "Sample.manifest.txt"
files = read.table(paste0(manfiestDir, manifestFile), header = T, sep = "\t", comment.char = "", quote = "");

#limit samples to use from manifest
files = files[files$include, ]

#bp range around bed coordinates to gather data
range = 100
#base pair range to plot, can be different from range of data gathered above
lim = 50

#remove these chromosomes from the analysis
removeChr = "random|chrY|dm6|Un"

#set output directory and file names
outDir = paste0(homeDir, bedName, "/");
#create new directory if necessary
if (!file.exists(outDir)) dir.create(outDir);
sumFile = paste0(outDir, bedName, ".hist.", range, ".bp.summary.txt")
sampFile = paste0(outDir, bedName, ".", lim, ".bp.rowSums.bySample.txt")
grpFile = paste0(outDir, bedName, ".", lim, ".bp.rowSums.byGroup.txt")

#choose normalization method
norm = "rpm"
#norm = "rppm" #frip score column required in manifest

#color scheme
nB.col = rgb(112, 229, 249, maxColorValue = 255)
PC.col = rgb(151, 247, 151, maxColorValue = 255)

#color mapping function...this can be customized and expanded for any group/sample size
color.map.group = function(name) {if(grepl("nB", name)) nB.col else if (grepl("PC", name)) PC.col}

#####
#Calculate bed coverage from sample BAM files
#####

	#read in bed file
	bed = read.table(bedFile, sep = "\t", skip = 1); 

	#bed file column names
	bedCols = c("chr", "start", "end", "motif", "score", "strand")
	colnames(bed) = bedCols

	#limit chromosomes
	bedLim = subset(bed, !grepl(removeChr, chr))

	#establish Granges object for bed file
	bedr = GRanges(seqnames = bedLim$chr, ranges = IRanges(start = bedLim$start, end = bedLim$end))

	#loop through each sample and calculate coverage
	for (j in 1:dim(files)[1]) {

		print(paste(files$sample[j], Sys.time()))

		#set file paths
		bamFile = paste0(files$dir[j], files$bamFile[j])
		outFile = paste0(outDir, files$sample[j], ".", bedName, ".range.", range, "bp.csv")

		if (norm == "rpm") {
			makeBamHist(bamFile, bedr, range = range, outFile = outFile, bin = 1, fragSizes = 1);
		}
		if (norm == "rppm") {
			makeBamHist(bamFile, bedr, range = range, outFile = outFile, bin = 1, normReads = (1e6 / files$frip[j]), fragSizes = 1);
		}
		#system(paste("pigz -f", statFile))
	}

######
#Plot histogram from data
######

	#summarize coverage data for all samples  at each bp
	sumBamHistCol(bedName = bedName, manifest = files, range = range, outDir = outDir, sumFile = sumFile)

	#plot histogram
	plotBedHist(manifest = files, sumFile = sumFile, lim = lim, norm = norm, colFun = color.map.group) 

######
#Plot boxplot from data
######

	#summarize coverage data at all bp for each sample
	sumBamHistRow(bedName = bedName, manifest = files, range = range, lim = lim, outDir = outDir, sampFile = sampFile, grpFile = grpFile, SumGrP = TRUE)

	#set y-axis limits, this should be modified for each bed file right now
	ylim = c(0,5)

	#plot by boxplot	
	plotBedBox(manifest = files, grpFile = grpFile, norm = norm, colFun = color.map.group, CalcStats = TRUE, ylim = ylim)

