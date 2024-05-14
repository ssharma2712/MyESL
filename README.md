# Evolutionary Sparse Learning with MyESL #

## Introduction ##
MyESL is a set of programs written in C++ and Python designed to build an evolutionary sparse learning model for a phylogenetic hypothesis or binary classes of species with or without a phenotype. MyESL uses sparse group lasso regression to optimize weights for each position in the sequence alignments, quantifying the strength of association between a position and the hypothesis/class used in the analysis. MyESL also performs pre-possing of input data and processes ESL model outputs for further biological discoveries. MyESL is currently available as a standalone software `(MyESL.exe)` compatible with Windows and Linux operating systems. Moreover, it offers the option of installation via a Python environment to conduct ESL analysis. 

## Implementation of MyESL ## 
We discuss every required and optional argument for ESL analysis of an example dataset using MyESL software (MyESL.exe). Next, we will discuss installing MyESL in the Python environment.  

### Download ###
You can download `MyESL` from the GitHub repository or use the command 

	git clone -b grid_search https://github.com/kumarlabgit/MyESL MyESL
	cd ESL

### Usage ###

One can perform ESL analysis using `MyESL.exe` after downloading and setting the current working directory to MyESL.  

<br />

```
MyESL.exe alignment_list.txt  --tree phylogenetic_tree.nwk

OR

MyESL.exe alignment_list.txt  --classes classes.txt
```
<br />


### Required arguments:

<br />

```

alignment_list.txt                       : A text file contains a list of paths for all sequence alignments. For example,
                                           angiosperm_alns/7276_C12.fasta
                                           angiosperm_alns/5111_C12.fasta
                                           angiosperm_alns/5507_C12.fasta

--tree <phylogenetic_tree.nwk>           : A phylogenetic tree in newick format with a node ID to construct a hypothesis for the clade of interest.
                                           The hypothesis can also be specified with a separate file using the --classes parameter.
                                           It is highly recommended that the number of species in the clade be equal to or less than those outside of the clade.
                                           It is also recommended to use the smart sampling option (—-class_bal) when the number of species inside the clade is less or greater than the number outside the clade.
OR

--classes <classes.txt>                  : Requires a text file containing a user-defined hypothesis. It has two columns, which are tab-separated. The first column contains species names, and the second column contains the response 
                                          value for the species (+1/-1). A member species in the clade or with a specific phenotype receives +1 and -1 otherwise. This hypothesis is unconstrained by the tree structure. It is highly 
                                          recommended that the number of species with the response +1 equal the number of species with -1. 

```
<br />	

### Optional arguments:

Users can also specify other options in MyESL for processing the input data, building ESL models, and post-processing the ESL models. 

#### Pre-processing input data:

<br />

```
--clade_list <string1, string2,...> : Users can test multiple phylogenetic hypotheses when the input phylogenetic tree contains multiple clade IDs. This option must be used with "--tree" option.

--gen_clade_list <int, int>         : Users can generate multiple hypotheses when the input phylogeny contains no clade ID. The size of the clade is determined by the input integers defining the upper and lower limits
                                      of clade size, respectively.
   
--class_bal <string>                : MyESL performs class balancing, a common practice in supervised machine learning. Class balancing helps balance the number of species inside and outside the focal clade of interest.
                                      Class balancing in MyESL is performed by phylogenetic aware sampling <phylo> when the "--tree" option provides the phylogenetic hypothesis. MyESL also makes a balance between classes 
                                      by using inverse weights using the option <weight>. Upsampling and downsampling between classes can also be performed in MyESL using <up> and <down> options. 

--data_type <string>                : <nucleotide> informs MyESL to treat A, T, C, G, and U as valid characters without case sensitivity, and all other characters will be treated as missing data.
                                      <protein> option treats all unambiguous IUPAC amino acid letters (case insensitive) as valid characters.
                                      <molecular> option provides a way to use both nucleotide and acid letters as valid characters.
                                      <universal> option is used as default, which allows the analysis of presence/absence (0/1) data.

--bit_ct <int>                      : One can drop all bit-columns in which the bit 1 appears fewer than a certain number of times.

```
<br />	

#### ESL model building:

<br />

```
--lambda1 <float>                    : The site sparsity parameter that ranges from 0 to 1. When not specified, the default is 0.1. It is required to build a single-clade model.

--lambda2 <float>                    : The gene sparsity parameter ranges from 0 to 1, and the default is 0.1 when not specified. It is required to build a single-clade model.

--lamda1_grid <min, max, step>       : This option allows users to set the range for the site sparsity parameter. The site sparsity grid is defined by a string of float numbers min, max, and step_size, which range from 0 to 1.
                                       For example, --lamda1_range 0.1, 0.9, 0.1. This option must be used with --lamda2_range.  

--lamda2_grid <min, max, step>       : This option allows users to set the range for the group sparsity parameter. The group sparsity grid is defined by a string of float numbers min, max, and step_size, which range from 0 to 1.
                                       For example, --lamda2_range 0.1, 0.9, 0.1. This option must be used with --lamda1_range. 

--min_groups <int>                   : This option allows users to set the minimum number of genes included in the multi-gene ESL models and helps early stopping in the grid search over the sparsity parameter space.
                                       It takes a value greater than zero (0) and builds models containing more or equal numbers of groups in the model.

--group_wt <filename.txt>            : MyESL uses the square root of the number of bit columns for groups as group weights by default. Users can provide group weights via a text file (<filename.txt>), a tab-separated two-column file.
                                       The first column contains group names, and the second contains corresponding group weights.

--kfold <int>                        : This directive allows users to build multiple ESL models using k-fold cross-validation. For example, if k-fold is set to 5, 80% of rows are used in model training, and 20% of rows are withheld for validation.

--no_group_penalty                   : This directive allows users to perform LASSO regression without group penalty. The LASSO regression without group penalty is also performed when a single FastA file is provided as input. 

```
<br />	

#### Processing ESL-model outputs:

<br />

```
--output <string>                   : The name of the output directory where all results from MyESL analysis will be stored. The program creates an output directory automatically if it is not specified.

--stats_out <BPGHS>                 : MyESL processes the ESL model and outputs, Bit Sparsity Scores (<B>), Position Sparsity Scores (<P>), Gene Sparsity Scores (<G>), Hypothesis Sparsity Scores (<H>), and
                                      Species Sparsity Scores and Prediction Probability (<S>). Users can output multiple files using multiple inputs like <BPS>. Details of these scores can be found in reference #1.  
```
<br />	

### MyESL output files

MyESL produces multiple output files based on different user analysis directives. 

#### Output files from a single ESL model

<br />

```
ESL_model_{clade_ID}.txt   : A tab-separated text file containing allele at each position in a group and non-zero beta values estimated from sparse group lasso analysis. 

BSS_{clade_ID}.txt         : A tab-separated text file containing allele at each position in a group and non-zero bit sparsity score.
PSS_{clade_ID}.txt         : A tab-separated text file containing position ID in a group and non-zero position sparsity scores.
GSS_{clade_ID}.txt         : A tab-separated text file containing group ID and non-zero group sparsity scores.
HSS_{clade_ID}.txt         : A text file containing hypothesis sparsity scores.
SPS_SPP_{clade_ID}.txt     : A tab-separated text file containing species name, regression response, species prediction score, and prediction probability. 

```
<br />

#### Output files from multiple ESL models

When grid options for penalty parameters are used, MyESL produces multiple ESL models and other sparsity scores. In addition to a single model output for each lambda pair, MyESL outputs a summary of all single models. 

<br />

```
ESL_model_{clade_ID}_l1_l2.txt   : A tab-separated text file containing allele at each position in a group and non-zero beta values estimated from sparse group lasso analysis. 

BSS_{clade_ID}_l1_l2.txt         : A tab-separated text file containing allele at each position in a group and non-zero bit sparsity score.
PSS_{clade_ID}_l1_l2.txt         : A tab-separated text file containing position ID in a group and non-zero position sparsity scores.
GSS_{clade_ID}_l1_l2.txt         : A tab-separated text file containing group ID and non-zero group sparsity scores.
HSS_{clade_ID}_l1_l2.txt         : A text file containing hypothesis sparsity scores.
SPS_SPP_{clade_ID}_l1_l2.txt     : A tab-separated text file containing species name, regression response, species prediction score, and prediction probability. 

BSS_{clade_ID}_summary.txt         : A tab-separated text file containing allele at each position in a group and non-zero bit sparsity score.
PSS_{clade_ID}_summary.txt         : A tab-separated text file containing position ID in a group and non-zero position sparsity scores.
GSS_{clade_ID}_summary.txt         : A tab-separated text file containing group ID and non-zero group sparsity scores.
HSS_{clade_ID}_summary.txt         : A text file containing hypothesis sparsity scores.
SPS_SPP_{clade_ID}_summary.txt     : A tab-separated text file containing species name, regression response, species prediction score, and prediction probability. 

Note: l1 and l2 refer to site and group sparsity parameters used for model building in grid search.
```
<br />

### DrPhylo analysis using MyESL

Users can perform DrPhylo analysis using MyESL. Instructions and Outputs for DrPhylo analysis will be found at [DrPhylo](https://github.com/kumarlabgit/MyESL/tree/DrPhylo). 


## Table of Contents ##

1. [Introduction](#introduction)
2. [Implementation of MyESL](#implementation)
3. [Requirements](#requirements)
4. [Setup](#setup)
5. [Usage](#usage)
6. [Output](#output)
7. [C++ Components](#individual-components)




## SGLASSO ##

Applying the sparse group Lasso penalty yields a solution that simultaneously achieves the within- and between-group sparsity. That is, many feature groups are exactly zero (thus not selected) and within the non-zero (thus selected) feature groups, some features are also exactly zero. The simultaneous within- and between- group sparsity makes sparse group Lasso ideal for applications where we are interested in identifying important groups as well as important features within the selected groups.

The algorithm can optimize feature weights under either a least squares or logistic regression model. Additionally, the overlapping variants of those programs can accept input where feature groups overlap.


## Requirements ##

To run ESL, you will need python 3.8 or later installed, as well as the following python packages

- numpy
- biopython
- matplotlib
- pandas

You can install these libraries using pip:

`pip install biopython numpy pandas matplotlib`

## Setup ##
After downloading the 
To setup ESL, perform the following on the command line:

	git clone -b grid_search https://github.com/kumarlabgit/ESL ESL
	cd ESL
	bash setup.sh


## Usage: ##

Once setup has been performed, the main pipeline can be applied to a set of alignment files and a hypothesis to determine the features/feature groups which have the most explanatory power in term of the provided hypothesis.

	python ESL_pipeline.py tree_file.nwk alignment_list.txt [--parameter_name parameter_value] ... [--boolean_option] ...

ESL can also be run in cross-validation or grid-search modes, using the --xval or (--grid_z + --grid_y) parameters, respectively. When in grid-search mode, each regression will be run with various feature- and group- sparsity parameters, and the results for each hypothesis will be aggregated across runs, with options to filter result sets with low accuracy/predictive power.

#### Required arguments:

* `tree_file.nwk`: A newick formatted tree file specifying the phylogenetic groups used to construct the hypothesis. Multiple hypotheses can be contructed from the provided newick tree file, or a single hypothesis may be specified with a separate hypothesis file provided using the --response parameter. To specify a hypothesis using the newick tree, a label should be added to the internal node that will be used as the positive response group. For each labeled internal node, a hypothesis will be generated where all terminal children of that node will be assigned the +1 response value, and all other terminal nodes will be given the -1 response value.

* `alignment_list.txt`: The alignment list file should contain one alignment filename (or comma-separated group of alignment filenames, in the case of overlapping groups) on each line.


#### Optional arguments:

* `--response hypothesis_file.txt`: To specify a single hypothesis to test, unconstrained by tree structure, the --response parameter can be pointed to a tab-separated two-column text file which contains a sequence ID in the first column, and either -1, 1, or 0 in the second column. Sequences with a 0 in the hypothesis file will be ignored completely.

* `--nodelist node_list.txt`: If the user does not wish to test a hypothesis for every labeled internal node, a node_list file may be used to specify a limited set of labeled internal nodes to test. The node_list file should have the name of one internally labeled node on each line.

* `--output, -o output_dir`: Directory to place outputs in.

* `--lambda1, -z <float>`: The feature level sparsity parameter, higher values will result in more feature sparsity, i.e. fewer features with a non-zero weight. Must be a positive, non-zero value less than one. Defaults to 0.1.

* `--lambda2, -y <float>`: The group level sparsity parameter, higher values will result in more group sparsity, i.e. more feature groups in which *every* feature in the group has a value of zero. Must be a positive, non-zero value less than one. Defaults to 0.1.

* `--method [leastr, logistic, ol_leastr, ol_logistic]`: The sparse group LASSO method to use (default: logistic).

* `--xval <int>`: Perform cross validation with N partitions for each hypothesis specified.

* `--grid_z [min,max,step_size]`: Run in grid-search mode with feature sparsity parameter in the specified range using the specified step-size. Must be used with --grid_y option.

* `--grid_y [min,max,step_size]`: Run in grid-search mode with group sparsity parameter in the specified range using the specified step-size. Must be used with --grid_z option.

* `--grid_rmse_cutoff <float>`: Filter models with RMSE above the cutoff when selecting models to aggregate.

* `--grid_acc_cutoff <float>`: Filter models with accuracy below the cutoff when selecting models to aggregate (range:0.0 - 1.0).

* `--grid_gene_threshold <int>`: Skip generating models that are guaranteed to produce less than or equal to this number of non-zero feature groups, i.e., skip generating a new model for any lambda pair where an existing model with L1<=new_L1 and L2<=new_L2 has the threshold number of genes or lower.

* `--auto_name_nodes`: If the user tree contains no labelled internal nodes, the --auto_name_nodes option may be specified to automatically generate names for each internal node, causing all each internal nodes to be used to generate a unique hypothesis for testing. Note: this may result in a very large number of hypotheses to be tested if the user tree is large, in which case it is suggested to also use the --cladesize_cutoff_lower and --cladesize_cutoff_upper parameters to limit the number of hypotheses that will be tested.

* `--cladesize_cutoff_lower <int>`: Internal nodes with fewer than cladesize_cutoff_lower terminal descendants will not be tested.

* `--cladesize_cutoff_upper <int>`: Internal nodes with greater than cladesize_cutoff_upper terminal descendants will not be tested.

* `--smart_sampling <int>`: Instead of assigning a -1 response value to all species not assigned a +1 value, smart sampling will try to select a phylogenetically informed negative sample of equal size to the positive sample. Modes:

* `--gene_penalties penalties.txt`: File of custom penalty values (same order as alignment list) to assign to alignments/groups/genes. Used to balance importance of feature groups of varying size.

#### Usage examples:

Run basic ESL on a tree with 2 internally labeled nodes, producing a result set for each of two hypotheses generated from the input tree:

	python3 ESL_pipeline.py sample_files/ESL_test.nwk sample_files/angiosperm_100_sample_alns.txt --output sample_tree_output

Run basic ESL on a tree with 2 internally labeled nodes, but specify a single hypothesis/response file, ignoring the tree (a tree must still be specified) and producing a single result set for the provided response file:

	python3 ESL_pipeline.py sample_files/ESL_test.nwk sample_files/angiosperm_100_sample_alns.txt --output sample_tree_output --response sample_files/test_pred.txt

Run an ESL grid search with both lambda values in the range 0.05-0.1 and a step size of 0.05, specifying a single hypothesis/response file:

	python3 ESL_pipeline.py sample_files/ESL_test.nwk sample_files/angiosperm_100_sample_alns.txt --grid_z 0.05,0.1,0.05 --grid_y 0.05,0.1,0.05 --output sample_grid2x2_output --response sample_files/test_pred.txt

Run an ESL cross validation with 2 partitions, specifying a single hypothesis/response file:

	python3 ESL_pipeline.py sample_files/ESL_test.nwk sample_files/angiosperm_100_sample_alns.txt --output sample_xval2_output --response sample_files/test_pred.txt --xval 2

Run basic ESL on a tree with 2 internally labeled nodes, producing a result set for each of two hypotheses generated from the input tree, using overlapping feature groups:

	python3 ESL_pipeline.py sample_files/ESL_test.nwk sample_files/angiosperm_100_sample_alns_overlapping.txt --output sample_tree_output_ol_logistic --method ol_logistic


## Output ##

#### Basic Results Files

`{internal_node_label}_{output_dir}_hypothesis_out_feature_weights.txt`: Tab-separated, two-column file of all non-zero feature weights generated by the sparse group LASSO.
	
`{internal_node_label}_{output_dir}_gene_predictions.txt`: Table of predictions produced by model, including the total weight contributed by each feature group to each sample's total score.

`{internal_node_label}_{output_dir}_gene_predictions.png`: Visual representation of *_gene_predictions.txt file.

`{internal_node_label}_{output_dir}_PSS.txt`: Position significance score for each position, calculated as the sum of absolute weights for each one-hot feature referring to that position. Only positions with non-zero significance are displayed. Also includes MIC values.

`{internal_node_label}_{output_dir}_GSS.txt`: Gene/group significance score for each non-zero gene/feature group, calculated as the sum of absolute weights for each feature in that gene.


#### Grid Search Results Files

In addition to containing versions of each of the basic output files with "_{lambda1}_{lambda2}" appended before the file extension (e.g. {internal_node_label}_{output_dir}_gene_predictions_{lambda1}_{lambda2}.txt) for each model generated, the output directory of a grid search process will contain the following results files aggregated across all non-filtered models produced by the grid search:

`{internal_node_label}_GCS_median.txt`: Median of each gene's contribution scores for each species (i.e. the values under the individual gene columns in the *_gene_predictions.txt file) across all non-filtered models.

`{internal_node_label}_GCS_median.png`: Visual representation of *_GCS_median.txt file.

`{internal_node_label}_PSS_median.txt`: Median of position significance scores across all non-filtered models.

`{internal_node_label}_GSS_median.txt`: Median of gene/group significance scores across all non-filtered models.


#### Cross-Validation Results Files

`{internal_node_label}_{output_dir}_gene_predictions_xval.txt`: Table of predictions produced by cross validation (i.e. the weights/predictions displayed for each sample are generated by the model where that sample was in the holdout set), including the total weight contributed by each feature group to each sample's total score.

`{internal_node_label}_{output_dir}_gene_predictions_xval.png`: Visual representation of *_gene_predictions_xval.txt file.


#### Intermediate/auxiliary files

`{internal_node_label}_{output_dir}_hypothesis.txt`: Tab-separated, two-column hypothesis/response file generated by the preprocess program from the provided newick tree.

`feature_{internal_node_label}_{output_dir}_hypothesis.txt`: Matrix of one-hot encoded features generated by the preprocess program from alignments and one hypothesis derived from the input tree. If --response option is used in place of deriving a hypothesis from the tree, {internal_node_label} will reflect the name of the file used for the --response option.

`feature_mapping_{internal_node_label}_{output_dir}_hypothesis.txt`: File mapping each one-hot column in the features file to a feature label specifying the alignment, position, and allele that the one-hot feature column refers to.

`field_{internal_node_label}_{output_dir}_hypothesis.txt`: An intermediate index mapping file needed for overlapping inputs, generated by the preprocess program and only used by overlapping versions of the sparse group LASSO.

`group_indices_{internal_node_label}_{output_dir}_hypothesis.txt`: A 3-row file providing the start coordinate, end coordinate, and weight penalty of each group in the features file. Note, for overlapping inputs, the start and end coordinates refer to ranges in the field file, which contain a set of positions/indexes into the features file. For non-overlapping inputs the start/end coordinates can be indexed directly to the features file.

`{internal_node_label}_{output_dir}_xval_groups.txt`: File indicating cross-validation partition to place each sample in.

`missing_seqs_{output_dir}.txt`: File indicating gene/sequence pairs for which there was no sequence data available in the provided alignments. Used to annotate gene_predictions files to differentiate between feature groups that contributed zero weight vs those that were missing from the input.

`{output_dir}_lambda_list.txt`: File of lambda value pairs to generate models for, used for the grid search mode.

## Individual Components ##

The phylogeny testing pipeline utilizes a set of programs written in C++ to transform inputs and perform individual sparse learning tasks. The individual programs have the following specifications. These programs can be found in the bin directory after running the setup script.

### Preprocess

#### Usage:

 	bin/preprocess hypothesis_file.txt alignment_list.txt output_basename [optional parameters]

#### Required arguments:

* `hypothesis_file.txt`: A tab-separated two-column text file which contains a sequence ID in the first column, and either -1, 1, or 0 in the second column. Sequences with a 0 in the hypothesis file will be ignored completely.

* `alignment_list.txt`: The alignment list file should contain one alignment filename (or comma-separated group of alignment filenames, in the case of overlapping groups) on each line. The paths for the alignment files must be relative to the directory that preprocess is run from. If alignments are in multiple directories, they should each have a unique basename, e.g. having both alns/cohort1/BRCA.fasta and alns/cohort2/BRCA.fasta in your alignment list will cause issues.

* `output_basename`: Basename to use for generated output files.

#### Optional arguments:

* `is`: Ignore singleton mutations.
* `ct <int>`: Ignore mutations observed fewer than N times.
* `useTabs`: Generate tab-separated outputs instead of comma-separated.
* `fuzzIndels`: Fuzz indels.

#### Usage Example
	cd sample_files
	../bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns.txt angiosperm_input ct 2
	mv angiosperm_input ..
	cd ..

### sg_lasso

Sparse group LASSO using logistic regression.

#### Usage:

 	bin/sg_lasso -f features_file.txt -n group_indices.txt -r response_file.txt [optional parameters]

#### Required arguments:

* `features_file.txt`: Matrix containing feature set where rows are samples and columns are features.
* `group_indices.txt`: Matrix of indexes defining non-overlapping group information.
* `response_file.txt`: Single-column response in same sample order as features file.

#### Optional arguments:

* `--output (-w)`: Output file basename to write learned feature weights to.
* `--model_format (-l)`: Specify an output model format of either "xml" or "flat" (defaults to xml).
* `--lambda1 (-z)`: Feature regularization parameter (z1 >=0). Default value 0.
* `--lambda2 (-y)`: Group regularization parameter (z2 >=0). Default value 0.
* `--xval (-x)`: Specify a tab-separated, two-column file of seq_id/index pairs, indicating the cross-validation partition index of each seq_id.
* `--lambda_list (-l)`: Specify a tab-separated, two-column file of lambda1/lambda2 pairs. A model will be generated for each pair in the list.
* `--gene_count_threshold (-c)`: Specify gene selection cutoff for grid-based model building. Only usable in conjunction with lambda_list option.
* `--slep (-s)`: Specify a tab-separated, two-column file of key/value SLEP options.


#### Usage Example
	cd sample_files
	../bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns.txt angiosperm_input ct 2
	mv angiosperm_input ..
	cd ..
 	bin/sg_lasso -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights

### sg_lasso_leastr

Sparse group LASSO using least squares regression.

#### Usage:

 	bin/sg_lasso_leastr -f features_file.txt -n group_indices.txt -r response_file.txt [optional parameters]

#### Required arguments:

* `features_file.txt`: Matrix containing feature set where rows are samples and columns are features.
* `group_indices.txt`: Matrix of indexes defining non-overlapping group information.
* `response_file.txt`: Single-column response in same sample order as features file.

#### Optional arguments:

* `--output (-w)`: Output file basename to write learned feature weights to.
* `--model_format (-l)`: Specify an output model format of either "xml" or "flat" (defaults to xml).
* `--lambda1 (-z)`: Feature regularization parameter (z1 >=0). Default value 0.
* `--lambda2 (-y)`: Group regularization parameter (z2 >=0). Default value 0.
* `--xval (-x)`: Specify a tab-separated, two-column file of seq_id/index pairs, indicating the cross-validation partition index of each seq_id.
* `--lambda_list (-l)`: Specify a tab-separated, two-column file of lambda1/lambda2 pairs. A model will be generated for each pair in the list.
* `--gene_count_threshold (-c)`: Specify gene selection cutoff for grid-based model building. Only usable in conjunction with lambda_list option.
* `--slep (-s)`: Specify a tab-separated, two-column file of key/value SLEP options.


#### Usage Example
	cd sample_files
	../bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns.txt angiosperm_input ct 2
	mv angiosperm_input ..
	cd ..
 	bin/sg_lasso_leastr -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights


### overlapping_sg_lasso_logisticr

Overlapping sparse group LASSO using logistic regression.

#### Usage:

 	bin/overlapping_sg_lasso_logisticr -f features_file.txt -n group_indices.txt -g field.txt -r response_file.txt [optional parameters]

#### Required arguments:

* `features_file.txt`: Matrix containing feature set where rows are samples and columns are features.
* `group_indices.txt`: Matrix of indexes defining non-overlapping group information.
* `field.txt`: Feature indices file for overlapping groups. Index ranges in group indices must be contiguous and non-overlapping so, for overlapping input, group_indices ranges are not direct indexes into the features file, but rather indexes into the field file, which is a list of indexes into the features file. E.g. for a feature file containing 3 features, with 2 overlapping groups each containing 2 features, your group indices would be (1-2, 3-4), and your field file would be (1,2,2,3).
* `response_file.txt`: Single-column response in same sample order as features file.


#### Optional arguments:

* `--output (-w)`: Output file basename to write learned feature weights to.
* `--model_format (-l)`: Specify an output model format of either "xml" or "flat" (defaults to xml).
* `--lambda1 (-z)`: Feature regularization parameter (z1 >=0). Default value 0.
* `--lambda2 (-y)`: Group regularization parameter (z2 >=0). Default value 0.
* `--xval (-x)`: Specify a tab-separated, two-column file of seq_id/index pairs, indicating the cross-validation partition index of each seq_id.
* `--lambda_list (-l)`: Specify a tab-separated, two-column file of lambda1/lambda2 pairs. A model will be generated for each pair in the list.
* `--gene_count_threshold (-c)`: Specify gene selection cutoff for grid-based model building. Only usable in conjunction with lambda_list option.
* `--slep (-s)`: Specify a tab-separated, two-column file of key/value SLEP options.


#### Usage Example
	cd sample_files
	../bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns_overlapping.txt ol_angiosperm_input ct 2
	mv angiosperm_input ..
	cd ..
 	bin/overlapping_sg_lasso_logisticr -f ol_angiosperm_input/feature_ol_angiosperm_input.txt -z 0.1 -y 0.5 -n ol_angiosperm_input/group_indices_ol_angiosperm_input.txt -g ol_angiosperm_input/field_ol_angiosperm_input.txt -r ol_angiosperm_input/response_ol_angiosperm_input.txt -w angiosperm_out_feature_weights


### overlapping_sg_lasso_leastr

Overlapping sparse group LASSO using least squares regression.

#### Usage:

 	bin/overlapping_sg_lasso_leastr -f features_file.txt -n group_indices.txt -g field.txt -r response_file.txt [optional parameters]

#### Required arguments:

* `features_file.txt`: Matrix containing feature set where rows are samples and columns are features.
* `group_indices.txt`: Matrix of indexes defining non-overlapping group information.
* `field.txt`: Feature indices file for overlapping groups. Index ranges in group indices must be contiguous and non-overlapping so, for overlapping input, group_indices ranges are not direct indexes into the features file, but rather indexes into the field file, which is a list of indexes into the features file. E.g. for a feature file containing 3 features, with 2 overlapping groups each containing 2 features, your group indices would be (1-2, 3-4), and your field file would be (1,2,2,3).
* `response_file.txt`: Single-column response in same sample order as features file.


#### Optional arguments:

* `--output (-w)`: Output file basename to write learned feature weights to.
* `--model_format (-l)`: Specify an output model format of either "xml" or "flat" (defaults to xml).
* `--lambda1 (-z)`: Feature regularization parameter (z1 >=0). Default value 0.
* `--lambda2 (-y)`: Group regularization parameter (z2 >=0). Default value 0.
* `--xval (-x)`: Specify a tab-separated, two-column file of seq_id/index pairs, indicating the cross-validation partition index of each seq_id.
* `--lambda_list (-l)`: Specify a tab-separated, two-column file of lambda1/lambda2 pairs. A model will be generated for each pair in the list.
* `--gene_count_threshold (-c)`: Specify gene selection cutoff for grid-based model building. Only usable in conjunction with lambda_list option.
* `--slep (-s)`: Specify a tab-separated, two-column file of key/value SLEP options.


#### Usage Example
	cd sample_files
	../bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns_overlapping.txt ol_angiosperm_input ct 2
	mv angiosperm_input ..
	cd ..
 	bin/overlapping_sg_lasso_leastr -f ol_angiosperm_input/feature_ol_angiosperm_input.txt -z 0.1 -y 0.5 -n ol_angiosperm_input/group_indices_ol_angiosperm_input.txt -g ol_angiosperm_input/field_ol_angiosperm_input.txt -r ol_angiosperm_input/response_ol_angiosperm_input.txt -w angiosperm_out_feature_weights
