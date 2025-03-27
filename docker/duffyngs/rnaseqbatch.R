# rnaseqbatch.R

setup_path <- "@SETUP_PATH@"
sample_set <- c("@SAMPLE_SET@")

library(DuffyTools)
library(DuffyNGS)
source(setup_path)

for (sample_id in sample_set) pipeline(sample_id)
