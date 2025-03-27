# rnaseqfin.R

setup_path <- "@SETUP_PATH@"

library(DuffyTools)
library(DuffyNGS)
source(setup_path)

# gather all the main alignment metrics into one summary file, written
# to the current working directory

pipsum = pipe.ExtractPipelineSummaryDetails(annT$SampleID)
write.csv(pipsum,'Pipeline.Summary.Details.csv')

# turn all the transcriptome files into one matrix of gene expression
# that you can then send to cluster tools, PCA, etc.

fileset <- file.path(
	"results/transcript/",
	paste(
		annT$SampleID, # annT defined in source(setup_path)
		"MTB.Transcript.txt",
		sep="."
	)
)

expdata <- expressionFileSetToMatrix(fnames=fileset, fids=annT$SampleID, verbose=T)
write.csv(expdata,'Expression.MTB.GeneData.csv')

plot(expressionCluster(expdata))

# pipe.MetaResults(
# 	annT$SampleID,
# 	folder="Induced.v.Uninduced",
# 	groupColumn="Group.Induced",
# 	colorColumn="Color.Induced",
# 	PLOT.FUN=plotGeneExpression,
# )

