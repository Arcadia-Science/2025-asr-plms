seqfile = 060525_run/isomaltase_dereplicated_final_recoded_mafft.ph
treefile = 060525_run/isomaltase_dereplicated_final_recoded_mafft_rooted.txt
outfile = codeml_results.txt

noisy = 9
verbose = 2
runmode = 0

seqtype = 2  # Protein sequences
aaRatefile = /Users/isabel/miniconda3/envs/asr/dat/lg.dat
model = 2

fix_alpha = 0  # Estimate alpha for rate heterogeneity
alpha = 0.5  # Initial value for gamma distribution
ncatG = 4

RateAncestor = 2  # Ancestral sequence reconstruction
cleandata = 0 # do not remove gaps
method = 0

fix_blength = 2 # fix tree branch lengths
