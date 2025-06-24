seqfile = /Users/isabel/Documents/ASR/old_run/online_generated_phylip2.txt
treefile = /Users/isabel/Documents/ASR/old_run/paml_input_tree.txt
outfile = codeml_results.txt

noisy = 9
verbose = 2
runmode = 0

seqtype = 2  # Protein sequences
aaRatefile = /Users/isabel/tools/paml-dev/dat/jones.dat  # substitution model
model = 2  

fix_alpha = 0  # Estimate alpha for rate heterogeneity
alpha = 0.5  # Initial value for gamma distribution
ncatG = 8  

RateAncestor = 2  # Ancestral sequence reconstruction
cleandata = 0 # do not remove gaps
method = 0

fix_blength = 2 # fix tree branch lengths