# pipelineSetup.R

# high level wrappers for a variety of pipeline steps

library(DuffyTools)
library(DuffyNGS)
nCores <- multicore.setup(max.cores=1)

annT <- readAnnotationTable("Annotation.txt")
setCurrentSpecies("MT_H37")
bowtie2Par.defaults("Options.txt")
cat("\nBowtie 2 Version: ", getBowtie2Version(), "\n")

# Just recalcuate the transcription of a sample...
`runTranscript` <- function(sampleSet=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,
				speciesID=NULL, ...)
{
	if (is.null(sampleSet)) {
		sampleSet <- readAnnotationTable(annotationFile)$SampleID
	}
	
	nCores <- multicore.currentCoreCount()

	if (nCores == 1) {
		for (sampleID in sampleSet) pipe.Transcriptome(sampleID, annotationFile, optionsFile,
					results.path=results.path, speciesID=speciesID, ...)
	} else {
		cat("\n\nParallel execution over ", nCores, " cores...\n")
		ans <- system.time(multicore.ans <- multicore.lapply(sampleSet, pipe.Transcriptome,
			annotationFile=annotationFile, optionsFile=optionsFile, 
			results.path=results.path, speciesID=speciesID, ...))
		cat("\nDone.\n\nParallel timing:\n")
		print(ans)
	}
}



# Differential Expression between all samples after all transcripts are done...
`runDiffExpression` <- function(sampleSet=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, 
				speciesID=NULL, ...)
{
	if (is.null(sampleSet)) {
		sampleSet <- readAnnotationTable(annotationFile)$SampleID
	}

	nCores <- multicore.currentCoreCount()

	if (nCores == 1) {
		pipe.DiffExpression(sampleSet, annotationFile, optionsFile, results.path=results.path,
				speciesID=speciesID, ...)
	} else {
		cat("\n\nParallel execution over ", nCores, " cores...\n")
		pairedSamples <- multicore.samplePairs(sampleSet)
		ans <- system.time(multicore.ans <- multicore.lapply(pairedSamples,
			pipe.DiffExpression, annotationFile=annotationFile, optionsFile=optionsFile, 
			results.path=results.path, speciesID=speciesID, ...))
		cat("\nDone.\n\nParallel timing:\n")
		print(ans)
	}
}


# Meta Results: Consensus of RoundRobin, Rank Product, and SAM for all samples after DiffExpression is done...
`runMetaResults` <- function(folderName="", sampleSet=NULL, speciesID=getCurrentSpecies(),
			annotationFile="Annotation.txt", optionsFile="Options.txt", 
			results.path=NULL, groupColumn="Group", 
			doNonGenes=FALSE, tailWidth=1000, ...)
{
	if (nchar(folderName) < 1) stop("Meta Results needs an explicit 'folderName' for its results.")


	if (is.null(sampleSet)) {
		annT <- readAnnotationTable(annotationFile)
		sampleSet <- annT$SampleID
	}
	
	pipe.MetaResults(sampleSet, speciesID=speciesID, annotationFile=annotationFile,
				folderName=folderName, optionsFile=optionsFile, 
				results.path=results.path, groupColumn=groupColumn, 
				useMultiHits=TRUE, keepIntergenics=doNonGenes,
				tailWidth=tailWidth, verbose=TRUE, ...)
}


runQuickQC <- function(sampleSet=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt") {
	if (is.null(sampleSet)) sampleSet <- annT$SampleID

	for (i in 1:length(sampleSet)) {
		sampleID <- sampleSet[i]
		pipe.QuickQC(sampleID, annotationFile, optionsFile,
			banner=paste("QuickQC:  ", sampleID),
			chunkSize=1000000, maxReads=5000000, nUSRkeep=500000, 
			mode="all", pause=0)
		if (i < length(sampleSet)) Sys.sleep(60)
	}
	return()
}

