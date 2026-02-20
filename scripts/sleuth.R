# sleuth.R
# Using sleuth and dplyr packages, performs differential expression analysis between the 2dpi and 6dpi samples.
# Creates a tab-delimited table of significant transcripts to be written to final pipeline report.

library(sleuth)
library(dplyr)

#this makes sure it will use the input table generated during the pipeline
args = commandArgs(trailingOnly=TRUE)
stab = read.table(args[1], header=TRUE)

#initialize sleuth object using sleuth_prep function from sleuth library
so = sleuth_prep(stab)
#fit a model comparing the two conditions (in this case, 2dpi and 6dpi)
so = sleuth_fit(so, ~condition, 'full')
#fit the reduced model to compare in the likelihood ratio test
so = sleuth_fit(so, ~1, 'reduced')
#likelihood ratio test for differential expression between conditions 
so = sleuth_lrt(so, 'reduced', 'full')

#extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 
#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 


#write to table,tab-delimited, with header row
write.table(dplyr::select(sleuth_significant, target_id, test_stat, pval, qval), file="results/sleuth_results.txt",quote = FALSE,row.names = FALSE,sep="\t")