# Setup

	git clone https://github.com/kumarlabgit/ESL ESL
	cd ESL
	bash setup.sh

Additionally, the phylogeny testing pipeline requires python>=3.8 to be installed, as well as the following python packages:

	biopython
	matplotlib

# Phylogeny Testing Pipeline

The pipeline is implemented as a python script which takes a set of gene alignments and a newick tree with at least one named internal node then, for each named internal node, creates a model predicting the chance that a given sequence descends from that node, and finally applies each predictive model to each input sequence and generates a table of predictive values grouped by gene for each sequence.

Usage:
	python ESL_pipeline.py tree_file.nwk alignment_list.txt [--parameter_name parameter_value] ... [--boolean_option] ...

	required arguments:
	  tree_file                       Newick tree file.
	  alignment_list                  List of FASTA alignment files.
	optional parameters:
	  --nodelist <string>             File containing list of named internal nodes to test. If no file is specified, each named internal node of input phylogeny will be tested.
	  --smart_sampling <int>          Mode 1: For each selected positive sample, attempt to select a balanced, phylogenetically informed negative sample.
	  --cladesize_cutoff_lower <int>  Internal nodes with fewer than cladesize_cutoff_lower terminal descendants will not be tested.
	  --cladesize_cutoff_upper <int>  Internal nodes with greater than cladesize_cutoff_upper terminal descendants will not be tested.
	  --method <string>               SGLasso type to use. Options are "logistic", "leastr", "ol_logistic", or "ol_leastr". Defaults to "leastr".
	  --lambda1 <float>               Feature sparsity parameter.
	  --lambda2 <float>               Group sparsity parameter.
	  --ensemble_parts <int>          Build gene-wise ensemble models, splitting the set of genes into ensemble_parts partitions for each run.
	  --ensemble_coverage <int>       Number of ensemble models to build. Each gene will be included in this many individual models.
	  --gene_penalties <string>       File of penalty values (same order as aln_list) to specify penalty score for each gene.
	  --slep_opts <string>            File of tab-separated name-value pairs (one per line) to specify SLEP options.
	  --auto_name_length <int>        Number of characters to take from sequence IDs when automatically generating internal node labels.
	  --gene_display_limit <int>      Limits the number of genes displayed in the generated graph images.
	  --gene_display_cutoff <float>   Limits genes displayed in the generated graph images to those with sum-of-squares greater than gene_display_cutoff.
	  --output <string>               Output directory.
	boolean options:
	  --auto_name_nodes               Assign automatically generated names to unnamed internal nodes, causing them to be tested.
	  --sparsify                      Iteratively increase sparsity until selected set of genes fits in one partition.

sample usage:

	python3 ESL_pipeline.py sample_files/ESL_test.nwk sample_files/angiosperm_100_sample_alns.txt --output sample_output


# C++ Components

## preprocess:

	parameter 1: response matrix file
	parameter 2: file containing list of alignment file paths
	parameter 3: basename of output files
	optional parameters (must be specified after parameters 1-3):
		n: "normalize" feature weights by column
		ct {N}: ignore mutations observed fewer than {N} times (must be an integer)
sample usage:

	cd sample_files
	../bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns.txt angiosperm_input n ct 2
	mv angiosperm_input ..
	cd ..

Notes and Caveats

The paths for the alignment files must be relative to the directory that preprocess is run from.

Any species with a response value of 0 in the predictions file (first argument to preprocess executable) will be left out of the features file entirely.

The preprocessing script is not very well tested. The easiest way to check if it worked correctly is to inspect the feature mapping file - if something went wrong, it will usually have nonsense characters as the alleles for some feature columns.

The list of alignments given to the preprocess program can also generate overlapping groups of input for the overlapping SGLasso algorithm, by specifying multiple comma-separated files/genes on a single line.

For example, for a non-overlapping set of input groups, the file might look like:

	aln_dir/gene1.fas
	aln_dir/gene2.fas
	aln_dir/gene3.fas
	aln_dir/gene4.fas

While input where gene2 shares properties with both gene1 and gene3 might look like:

	aln_dir/gene1.fas,aln_dir/gene2.fas
	aln_dir/gene3.fas,aln_dir/gene2.fas
	aln_dir/gene4.fas


## sg_lasso_leastr

	required inputs:
	  --features_file (-f)          Matrix containing feature set A.
	  --opts_ind_file (-n)          Matrix of indices defining non-overlapping group information.
	  --responses_file (-r)         Vector containing responses y.
	optional inputs:
	  --feature_weights_file (-w)   Output file to write learned feature weights to.
	  --lambda1 (-z)                Feature regularization parameter (z1 >=0). Default value 0.
	  --lambda2 (-y)                Group regularization parameter (z2 >=0). Default value 0.

sample usage:

	bin/sg_lasso_leastr -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights


## overlapping_sg_lasso_leastr

	required inputs:
	  --features_file (-f)          Matrix containing feature set A.
	  --opts_ind_file (-n)          Matrix of indices defining overlapping group information.
	  --field_file (-g)             Vector of feature indices for overlapping groups.
	  --responses_file (-r)         Vector containing responses y.
	optional inputs:
	  --feature_weights_file (-w)   Output file to write learned feature weights to.
	  --lambda1 (-z)                Feature regularization parameter (z1 >=0). Default value 0.
	  --lambda2 (-y)                Group regularization parameter (z2 >=0). Default value 0.

sample usage:

	bin/overlapping_sg_lasso_leastr -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -g angiosperm_input/field_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights

## sg_lasso

	required inputs:
	  --features_file (-f)          Matrix containing feature set A.
	  --opts_ind_file (-n)          Matrix of indices defining non-overlapping group information.
	  --responses_file (-r)         Vector containing responses y.
	optional inputs:
	  --feature_weights_file (-w)   Output file to write learned feature weights to.
	  --lambda1 (-z)                Feature regularization parameter (z1 >=0). Default value 0.
	  --lambda2 (-y)                Group regularization parameter (z2 >=0). Default value 0.

sample usage:

	bin/sg_lasso -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights


## overlapping_sg_lasso_logisticr

	required inputs:
	  --features_file (-f)          Matrix containing feature set A.
	  --opts_ind_file (-n)          Matrix of indices defining overlapping group information.
	  --field_file (-g)             Vector of feature indices for overlapping groups.
	  --responses_file (-r)         Vector containing responses y.
	optional inputs:
	  --feature_weights_file (-w)   Output file to write learned feature weights to.
	  --lambda1 (-z)                Feature regularization parameter (z1 >=0). Default value 0.
	  --lambda2 (-y)                Group regularization parameter (z2 >=0). Default value 0.

sample usage:

	bin/overlapping_sg_lasso_logisticr -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -g angiosperm_input/field_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights


## Parsing raw outputs of C++ components

The weights in the output XML file are in the same order as the lines in the feature mapping file.
Simple commands to remove the XML formatting, merge the two, and remove all features with a weight of zero:

	grep -P "<item>.*</item>" angiosperm_out_feature_weights.xml | sed -re "s/.*<item>(.*)<\/item>.*/\1/" > temp_angiosperm_out_feature_weights.txt
	paste <(sed -e "1d" angiosperm_input/feature_mapping_angiosperm_input.txt) temp_angiosperm_out_feature_weights.txt | grep -v "0.00000000000000000e+00" > angiosperm_out_feature_weights.txt

# Building from source

To install all dependencies (hopefully) for building from source:

	sudo apt -y install g++-8 libopenblas-dev liblapack-dev

To build everything, do the following:

	git clone -b dev https://github.com/kumarlabgit/ESL ESL
	cd ESL
	bash build_script.sh
