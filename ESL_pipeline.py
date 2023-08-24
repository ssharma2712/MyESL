import os
import time
import argparse
import shutil
import math
import statistics
import copy
import pipeline_funcs as pf
import pandas
import gene_contribution_visualizer as gcv
from multiprocessing import Process


def grid_search(args, original_output, input_files):
	hypothesis_file_list = input_files["hypothesis_file_list"]
	slep_opts_file_list = input_files["slep_opts_file_list"]
	features_filename_list = input_files["features_filename_list"]
	groups_filename_list = input_files["groups_filename_list"]
	response_filename_list = input_files["response_filename_list"]
	gene_list = input_files["gene_list"]
	field_filename_list = input_files["field_filename_list"]
	group_list = input_files["group_list"]
	HSS = {}
	missing_seqs = set()
	gcv_files = []
	# features_filename_list, groups_filename_list, response_filename_list, gene_list, field_filename_list, group_list = pf.generate_input_matrices(args.aln_list, hypothesis_file_list, args)
	with open(os.path.join(original_output, "missing_seqs_" + original_output + ".txt"), "r") as file:
		for line in file:
			data = line.strip().split("\t")
			missing_seqs.add((data[1], os.path.splitext(os.path.basename(data[0]))[0]))
	weights_file_list = pf.run_esl(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, args.lambda1, args.lambda2, args.method, slep_opts_file_list, z_ind = args.z_ind, y_ind = args.y_ind)
	pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list)
	os.mkdir(args.output)
	for hypothesis_filename in hypothesis_file_list:
		# shutil.move(hypothesis_filename, args.output)
		shutil.move(hypothesis_filename.replace(".txt","_{}_{}_out_feature_weights.xml".format(args.z_ind, args.y_ind)), args.output)
		shutil.move(hypothesis_filename.replace("hypothesis.txt", "{}_{}_gene_predictions.txt".format(args.z_ind, args.y_ind)), args.output)
		shutil.move(hypothesis_filename.replace("hypothesis.txt", "{}_{}_mapped_feature_weights.txt".format(args.z_ind, args.y_ind)), args.output)
		shutil.move(hypothesis_filename.replace("hypothesis.txt", "{}_{}_PSS.txt".format(args.z_ind, args.y_ind)), args.output)
		shutil.move(hypothesis_filename.replace("hypothesis.txt", "{}_{}_GSS.txt".format(args.z_ind, args.y_ind)), args.output)
	# if args.auto_name_nodes:
	#	shutil.move("auto_named_{}".format(os.path.basename(args.tree)), args.output)
	with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
		file.write("{}\t{}\n".format("Hypothesis", "HSS"))
		for hypothesis_filename in hypothesis_file_list:
			file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), HSS[hypothesis_filename]))
			if args.slep_sample_balance or args.smart_sampling:
				shutil.move(hypothesis_filename.replace("hypothesis.txt", "slep_opts.txt"), args.output)
				shutil.move(hypothesis_filename.replace("hypothesis.txt", "sweights.txt"), args.output)
			if not args.grid_summary_only:
				gcv_files.append(gcv.main(os.path.join(args.output,hypothesis_filename.replace("hypothesis.txt", "{}_{}_gene_predictions.txt".format(args.z_ind, args.y_ind))),gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff))
	for file in gcv_files:
		if os.path.dirname(file)!=os.path.normpath(args.output):
			shutil.move(file, args.output)
	return hypothesis_file_list


def main(args):
	hypothesis_file_list, slep_opts_file_list = pf.generate_hypothesis_set(args)
	HSS = {}
	missing_seqs = set()
	merged_rep_predictions_files = {hypothesis_filename:[] for hypothesis_filename in hypothesis_file_list}
	gcv_files = []
	if args.ensemble_parts is not None and args.ensemble_parts >= 1:
		tempdir_list = []
		merged_parts_prediction_files = {hypothesis_filename:[] for hypothesis_filename in hypothesis_file_list}
		for i in range(0, args.ensemble_coverage):
			partitioned_aln_lists = pf.split_gene_list(args.aln_list, args.ensemble_parts)
			j = 0
			gene_prediction_files = {hypothesis_filename:[] for hypothesis_filename in hypothesis_file_list}
			for part_aln_list in partitioned_aln_lists:
				tempdir = "{}_rep{}_part{}".format(args.output, i+1, j+1)
				tempdir_list.append(tempdir)
				features_filename_list, groups_filename_list, response_filename_list, gene_list, field_filename_list, group_list = pf.generate_input_matrices(part_aln_list, hypothesis_file_list, args)
				with open(os.path.join(args.output, "missing_seqs_" + args.output + ".txt"), "r") as file:
					for line in file:
						data = line.strip().split("\t")
						missing_seqs.add((data[1], os.path.splitext(os.path.basename(data[0]))[0]))
				weights_file_list = pf.run_esl(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, args.lambda1, args.lambda2, args.method, slep_opts_file_list)
				pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list)
				for hypothesis_filename in hypothesis_file_list:
					if i+1 == args.ensemble_coverage and j+1 == args.ensemble_parts and not args.sparsify:
						shutil.move(hypothesis_filename, args.output)
					else:
						shutil.copy(hypothesis_filename, args.output)
					try:
						shutil.move(hypothesis_filename.replace(".txt","_out_feature_weights.xml"), args.output)
					except:
						pass
					shutil.move(hypothesis_filename.replace("hypothesis.txt","gene_predictions.txt"), args.output)
					gene_prediction_files[hypothesis_filename].append(os.path.join(tempdir, hypothesis_filename.replace("hypothesis.txt","gene_predictions.txt")))
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "mapped_feature_weights.txt"), args.output)
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "PSS.txt"), args.output)
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "GSS.txt"), args.output)
				shutil.move(args.output, tempdir)
				j += 1
			if not args.sparsify:
				for hypothesis_filename in hypothesis_file_list:
					merged_parts_prediction_files[hypothesis_filename].append(gcv.merge_predictions(gene_prediction_files[hypothesis_filename],hypothesis_filename.replace("hypothesis.txt","merged_gene_predictions_rep{}.txt".format(i))))
		if not args.sparsify:
			for hypothesis_filename in hypothesis_file_list:
				merged_rep_predictions_files[hypothesis_filename] = gcv.merge_predictions(merged_parts_prediction_files[hypothesis_filename],hypothesis_filename.replace("hypothesis.txt","merged_gene_predictions_final.txt"))
				gcv_files.extend(merged_parts_prediction_files[hypothesis_filename])
				gcv_files.append(merged_rep_predictions_files[hypothesis_filename])
				gcv_files.append(gcv.main(merged_rep_predictions_files[hypothesis_filename], gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff))
		os.mkdir(args.output)
		for tempdir in tempdir_list:
			shutil.move(tempdir, args.output)
		for file in gcv_files:
			shutil.move(file, args.output)
		if args.auto_name_nodes:
			shutil.move("auto_named_{}".format(os.path.basename(args.tree)), args.output)
		with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
			file.write("{}\t{}\n".format("Hypothesis", "HSS"))
			for hypothesis_filename in hypothesis_file_list:
				file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), HSS[hypothesis_filename]))
				if args.slep_sample_balance:
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "slep_opts.txt"), args.output)
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "sweights.txt"), args.output)
		result_files_list = pf.find_result_files(args, hypothesis_file_list)
		weights = pf.parse_result_files(args, result_files_list)
		return pf.analyze_ensemble_weights(args, weights)
	else:
		features_filename_list, groups_filename_list, response_filename_list, gene_list, field_filename_list, group_list = pf.generate_input_matrices(args.aln_list, hypothesis_file_list, args)
		with open(os.path.join(args.output, "missing_seqs_" + args.output + ".txt"), "r") as file:
			for line in file:
				data = line.strip().split("\t")
				missing_seqs.add((data[1], os.path.splitext(os.path.basename(data[0]))[0]))
		weights_file_list = pf.run_esl(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, args.lambda1, args.lambda2, args.method, slep_opts_file_list)
		pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list)
		for hypothesis_filename in hypothesis_file_list:
			shutil.move(hypothesis_filename, args.output)
			try:
				shutil.move(hypothesis_filename.replace(".txt","_out_feature_weights.xml"), args.output)
			except:
				pass
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "gene_predictions.txt"), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "mapped_feature_weights.txt"), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "PSS.txt"), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "GSS.txt"), args.output)
		if args.auto_name_nodes:
			shutil.move("auto_named_{}".format(os.path.basename(args.tree)), args.output)
		with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
			file.write("{}\t{}\n".format("Hypothesis", "HSS"))
			for hypothesis_filename in hypothesis_file_list:
				file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), HSS[hypothesis_filename]))
				if args.slep_sample_balance or args.smart_sampling:
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "slep_opts.txt"), args.output)
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "sweights.txt"), args.output)
				gcv_files.append(gcv.main(os.path.join(args.output,hypothesis_filename.replace("hypothesis.txt", "gene_predictions.txt")),gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff))
		for file in gcv_files:
			if os.path.dirname(file)!=os.path.normpath(args.output):
				shutil.move(file, args.output)
		return hypothesis_file_list


def threaded_main(args, original_output):
	main(args)
	shutil.move(args.output, original_output)


def threaded_grid_search(args, original_output, input_files):
	grid_search(args, original_output, input_files)
	shutil.move(args.output, original_output)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Phylogenetic hypothesis tester.")
	parser.add_argument("tree", help="Input phylogeny to perform testing on.", type=str)
	parser.add_argument("aln_list", help="List of alignment files to extract features from.", type=str)
	parser.add_argument("--nodelist", help="File containing list of named internal nodes to test. If no file is specified, each named internal node of input phylogeny will be tested.",
						type=str, default=None)
	parser.add_argument("--auto_name_nodes", help="Assign automatically generated names to unnamed internal nodes, causing them to be tested.", action='store_true', default=False)
	parser.add_argument("--auto_name_length", help="Number of characters to take from sequence IDs to generate internal node labels.", type=int, default=5)
	parser.add_argument("--cladesize_cutoff_lower", help="Internal nodes with fewer than cladesize_cutoff_lower terminal descendants will not be tested.", type=int, default=0)
	parser.add_argument("--cladesize_cutoff_upper", help="Internal nodes with greater than cladesize_cutoff_upper terminal descendants will not be tested.", type=int,default=None)
	parser.add_argument("--smart_sampling", help="For each selected positive sample, select a balanced, phylogenetically informed negative sample.", type=int, default=None)
	parser.add_argument("--slep_sample_balance", help="Automatically uses the SLEP option to balance sample weights, only available for non-overlapping, logistic sglasso.",
						action='store_true', default=False)
	parser.add_argument("--response", help="File containing list of named node/response value pairs.", type=str, default=None)
	parser.add_argument("-z", "--lambda1", help="Feature sparsity parameter.", type=float, default=0.1)
	parser.add_argument("-y", "--lambda2", help="Group sparsity parameter.", type=float, default=0.1)
	parser.add_argument("--grid_z", help="Grid search sparsity parameter interval specified as 'min,max,steps'", type=str, default=None)
	parser.add_argument("--grid_y", help="Grid search group sparsity parameter interval specified as 'min,max,steps'", type=str, default=None)
	parser.add_argument("--grid_rmse_cutoff", help="RMSE cutoff when selecting models to aggregate.", type=float, default=100.0)
	parser.add_argument("--grid_acc_cutoff", help="Accuracy cutoff when selecting models to aggregate.", type=float, default=0.0)
	parser.add_argument("--grid_threads", help="Number of threads to use when running grid search.", type=int, default=1)
	parser.add_argument("--grid_summary_only", help="Skip generating graphics for individual runs/models.", action='store_true', default=False)
	parser.add_argument("--no_group_penalty", help="Perform mono-level optimization, ignoring group level sparsity penalties.",
						action='store_true', default=False)
	parser.add_argument("-o", "--output", help="Output directory.", type=str, default="output")
	parser.add_argument("--upsample_balance", help="Balance positive and negative response sets by upsampling the underpopulated set.", action='store_true', default=False)
	parser.add_argument("--downsample_balance", help="Balance positive and negative response sets by downsampling the overpopulated set.", action='store_true', default=False)
	parser.add_argument("--fuzz_indels", help="Assign indels in variable sites a one-hot value of 0.5 instead of 0.", action='store_true', default=False)
	parser.add_argument("--ensemble_parts", help="Build gene-wise ensemble models, splitting the set of genes into N partitions for each run.", type=int, default=None)
	parser.add_argument("--ensemble_coverage", help="Number of ensemble models to build. Each gene will be included in this many individual models.", type=int, default=5)
	parser.add_argument("--sparsify", help="Iteratively increase sparsity until selected set of genes fits in one partition.", action='store_true', default=False)
	parser.add_argument("--method", help="SGLasso type to use. Options are \"logistic\", \"leastr\", or \"ol_leastr\". Defaults to \"leastr\".", type=str, default="logistic")
	parser.add_argument("--slep_opts", help="File of tab-separated name-value pairs (one per line) to specify SLEP options.", type=str, default=None)
	parser.add_argument("--gene_penalties", help="File of penalty values (same order as aln_list) to specify penalty score for each gene.", type=str, default=None)
	parser.add_argument("--gene_display_limit", help="Limits the number of genes displayed in the generated graph images.", type=int, default=100)
	parser.add_argument("--gene_display_cutoff", help="Limits genes displayed in the generated graph images to those with sum-of-squares greater than cutoff value.", type=int, default=0.0)
	args = parser.parse_args()
	if args.no_group_penalty:
		args.lambda2 = 0.0000001
	if (args.grid_z is None and args.grid_y is not None) or (args.grid_y is None and args.grid_z is not None):
		raise Exception("Only one grid search parameter specified, --grid_z and --grid_y must be specified together.")
	elif args.grid_z is not None and args.grid_y is not None:
		args_original = copy.deepcopy(args)
		output_folder = args.output
		hypothesis_file_list = []
		# for x in args.grid_z.strip().split(','):
		# 	print(x)
		# 	print(float(x))
		[z_min, z_max, z_steps] = [float(x) for x in args.grid_z.strip().split(',')]
		z_interval = (z_max - z_min) / z_steps
		z_list = [z_min + (x * z_interval) for x in range(0, int(z_steps) + 1)]
		[y_min, y_max, y_steps] = [float(x) for x in args.grid_y.split(',')]
		y_interval = (y_max - y_min) / y_steps
		y_list = [y_min + (x * y_interval) for x in range(0, int(y_steps) + 1)]
		[z_count, y_count] = [0, 0]
		input_files = {}
		input_files["hypothesis_file_list"], input_files["slep_opts_file_list"] = pf.generate_hypothesis_set(args)
		input_files["features_filename_list"], input_files["groups_filename_list"], input_files["response_filename_list"], input_files["gene_list"], \
		input_files["field_filename_list"], input_files["group_list"] = pf.generate_input_matrices(args.aln_list, input_files["hypothesis_file_list"], args)
		# os.mkdir(args_original.output)
		thread_id_list = []
		if args_original.grid_threads > 1:
			for lambda1 in z_list:
				for lambda2 in y_list:
					# Add wait/check loop here for multiprocessing thread limit
					while True:
						completed_threads = []
						for thread_id in thread_id_list:
							if os.path.exists(os.path.join(args_original.output, thread_id)):
								completed_threads.append(thread_id)
						for thread_id in completed_threads:
							thread_id_list.remove(thread_id)
						if len(thread_id_list) < args_original.grid_threads:
							break
						else:
							time.sleep(1)
					args = copy.deepcopy(args_original)
					args.lambda1 = lambda1
					args.lambda2 = lambda2
					args.z_ind = z_count
					args.y_ind = y_count % len(y_list)
					args.output = args_original.output + "_{}_{}".format(z_count, y_count % len(y_list))
					#hypothesis_file_list += main(args)
					# Start new thread for:
					###################
					thread_id_list.append(args_original.output + "_{}_{}".format(z_count, y_count % len(y_list)))
					p = Process(target=threaded_grid_search, args=(args, args_original.output, input_files))
					p.start()
					###################
					y_count += 1
				z_count += 1
			# Add wait/check loop here for multiprocessing step completion
			while True:
				completed_threads = []
				for thread_id in thread_id_list:
					if os.path.exists(os.path.join(args_original.output, thread_id)):
						completed_threads.append(thread_id)
				for thread_id in completed_threads:
					thread_id_list.remove(thread_id)
				if len(thread_id_list) == 0:
					break
				else:
					time.sleep(1)
			for z_ind in range(0, len(z_list)):
				for y_ind in range(0, len(y_list)):
					# Find all hypothesis files by pattern
					excluded = ['response', 'pos_stats', 'group_indices', 'field', 'feature']
					hypothesis_files = [fn for fn in os.listdir(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind))) if fn.endswith("{}_{}_{}_hypothesis.txt".format(args_original.output,z_ind,y_ind)) and not any(fn.startswith(prefix) for prefix in excluded)]
					hypothesis_file_list += hypothesis_files
		else:
			for lambda1 in z_list:
				for lambda2 in y_list:
					# Add wait/check loop here for multiprocessing
					args = copy.deepcopy(args_original)
					args.lambda1 = lambda1
					args.lambda2 = lambda2
					args.z_ind = z_count
					args.y_ind = y_count % len(y_list)
					args.output = args_original.output + "_{}_{}".format(z_count, y_count % len(y_list))
					# hypothesis_file_list += main(args)
					grid_search(args, args_original.output, input_files)
					shutil.move(args.output, args_original.output)
					y_count += 1
				z_count += 1
		pss_vals = {}
		gss_vals = {}
		# gene_predictions = []
		hypothesis_list = list(set(["_".join(x.replace("_" + args_original.output, "").split("_")[0:-1]) for x in input_files["hypothesis_file_list"]]))
		for hypothesis in hypothesis_list:
			gene_predictions = []
			for z_ind in range(0, len(z_list)):
				for y_ind in range(0, len(y_list)):
					# Parse gene_predictions files
					rmse = 0
					acc = 1
					#model = pandas.read_csv(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_gene_predictions.txt".format(hypothesis, args_original.output, z_ind, y_ind)), delim_whitespace=True, index_col=0)
					model = pandas.read_csv(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_gene_predictions.txt".format(hypothesis, args_original.output, z_ind, y_ind)), delim_whitespace=True, index_col=0)
					rmse = ((model.Response - model.Prediction) ** 2).mean() ** 0.5
					acc = sum([1 for x in zip(model.Response, model.Prediction) if (x[0] == -1 and x[1] < 0) or (x[0] == 1 and x[1] > 0)])/float(model.shape[0])
					# Calculate acc and RMSE
					if rmse <= args_original.grid_rmse_cutoff and acc >= args_original.grid_acc_cutoff:
						gene_predictions.append(model)
					# Parse GSS file into GSS vals
					# with open(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_GSS.txt".format(hypothesis, args_original.output, z_ind, y_ind)), 'r') as gss_file:
					with open(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_GSS.txt".format(hypothesis, args_original.output, z_ind, y_ind)), 'r') as gss_file:
						for line in gss_file:
							data = line.strip().split('\t')
							if data[0] == "Gene":
								continue
							if data[0] in gss_vals:
								gss_vals[data[0]].update({"{}_{}".format(z_ind, y_ind): data[1]})
							else:
								gss_vals[data[0]] = {"{}_{}".format(z_ind, y_ind): data[1]}
					# Parse PSS file into PSS vals
					# with open(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_PSS.txt".format(hypothesis, args_original.output, z_ind, y_ind)), 'r') as pss_file:
					with open(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_PSS.txt".format(hypothesis, args_original.output, z_ind, y_ind)), 'r') as pss_file:
						for line in pss_file:
							data = line.strip().split('\t')
							if data[0] == "Position Name" or float(data[1]) == 0.0:
								continue
							if data[0] in pss_vals:
								pss_vals[data[0]].update({"{}_{}".format(z_ind, y_ind): data[1]})
							else:
								pss_vals[data[0]] = {"{}_{}".format(z_ind, y_ind): data[1]}
			# Write aggregated GSS file
			with open(os.path.join(args_original.output, "{}_GSS_median.txt".format(hypothesis)), 'w') as file:
				for gene in gss_vals.keys():
					file.write("{}\t{}\n".format(gene, statistics.median([float(x) for x in gss_vals[gene].values() if x != 0])))
			# Write aggregated PSS file
			with open(os.path.join(args_original.output, "{}_PSS_median.txt".format(hypothesis)), 'w') as file:
				for pos in pss_vals.keys():
					file.write("{}\t{}\n".format(pos, statistics.median([float(x) for x in pss_vals[pos].values() if x != 0])))
			# Write aggregated gene_predictions file
			with open(os.path.join(args_original.output, "{}_GCS_median.txt".format(hypothesis)), 'w') as file:
				file.write("SeqID\t{}\n".format("\t".join([gene for gene in gene_predictions[0].columns if gene not in ["SeqID", "Prediction", "Intercept"]])))
				for species in gene_predictions[0].index:
					file.write("{}\t".format(species))
					for gene in gene_predictions[0].columns:
						if gene in ["SeqID", "Prediction", "Intercept"]:
							continue
						GCS_list = [x[gene][species] for x in gene_predictions if x[gene][species] != 0]
						if len(GCS_list) == 0:
							GCS_list = [0]
						# if gene == "Response":
						# 	print(GCS_list)
						# 	print(statistics.median(GCS_list))
						file.write("{}\t".format(statistics.median(GCS_list)))
					file.write("\n")
			gcv.gcv_median(os.path.join(args_original.output, "{}_GCS_median.txt".format(hypothesis)), gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff)
		for file in input_files["hypothesis_file_list"]:
			shutil.move(file, args_original.output)
	else:
		score_tables = main(args)
	gene_target = None
	if args.sparsify and score_tables is not None:
		args_original = copy.deepcopy(args)
		output_folder = args.output
		aln_list_dir = os.path.dirname(args.aln_list)
		aln_list_basename = os.path.splitext(os.path.basename(args.aln_list))[0]
		aln_file_list = {}
		with open(args.aln_list, 'r') as file:
			for line in file:
				for aln_filename in line.strip().split(","):
					basename = os.path.splitext(os.path.basename(aln_filename.strip()))[0]
					if basename in aln_file_list.keys():
						if aln_file_list[basename] != aln_filename.strip():
							raise Exception("Found multiple alignment files with identical basename {}.".format(basename))
					else:
						aln_file_list[basename] = aln_filename.strip()
		for hypothesis in score_tables.keys():
			args = copy.deepcopy(args_original)
			if gene_target is None:
				gene_target = math.ceil(len(score_tables[hypothesis])/args.ensemble_parts)
			elif gene_target != math.ceil(len(score_tables[hypothesis])/args.ensemble_parts):
				raise Exception("Mismatched gene counts between hypotheses.")
			selected_count = sum([1 for val in score_tables[hypothesis].values() if val[0] > 0])
			counter = 1
			final_round = False
			if selected_count <= gene_target and not final_round:
				final_round = True
			shutil.move("{}_hypothesis.txt".format(hypothesis), "{}.txt".format(hypothesis))
			while selected_count > gene_target or final_round:
				# Generate new aln_list file and point args.aln_list at it
				aln_list_filename = os.path.join(aln_list_dir, "{}_{}_{}.txt".format(aln_list_basename, hypothesis, counter))
				with open(aln_list_filename, 'w') as file:
					with open(args.aln_list, 'r') as original_file:
						stashed_gene_groups = []
						nonzero_genes = [os.path.basename(key) for key in score_tables[hypothesis].keys() if score_tables[hypothesis][key][0] > 0]
						#print(nonzero_genes)
						new_gene_groups = []
						for line in original_file:
							gene_group = line.strip().split(",")
							#print(len(gene_group))
							new_gene_group = sorted([os.path.splitext(os.path.basename(val))[0] for val in gene_group if os.path.splitext(os.path.basename(val))[0] in nonzero_genes])
							print("{}\t::\t{}".format(",".join(gene_group), ",".join(new_gene_group)))
							if tuple(new_gene_group) not in stashed_gene_groups:
								stashed_gene_groups.append(tuple(new_gene_group))
								new_gene_groups.append(new_gene_group)
						#print(new_gene_groups)
						for new_gene_group in new_gene_groups:
							if len(new_gene_group) > 0:
								file.write("{}\n".format(",".join([aln_file_list[val] for val in new_gene_group])))
				args.aln_list = aln_list_filename
				# Also point --response at this hypothesis' response file
				args.response = "{}.txt".format(hypothesis)
				# args.output = os.path.join(output_folder, "{}_sparsify_round{}".format(hypothesis, counter))
				args.output = "{}_sparsify_round{}".format(hypothesis, counter)
				# temp_score_tables = main(args)
				if final_round:
					args.output = "{}_sparsify_final".format(hypothesis)
					args.lambda1 = 10**-6
					args.lambda2 = 10**-6
					args.ensemble_parts = 1
					args.ensemble_coverage = 1
					args.sparsify = False
				score_tables.update(main(args))
				args.sparsify = True
				shutil.move(args.output, output_folder)
				new_selected_count = sum([1 for val in score_tables[hypothesis].values() if val[0] > 0])
				if not final_round and new_selected_count == selected_count:
					# This round didn't increase sparsity, so increase group sparsity parameter for the next round
					print("The number of selected genes was unchanged, increasing group sparsity parameter from {} to {}.".format(args.lambda2, args.lambda2 * 1.2))
					args.lambda2 = args.lambda2 * 1.2
				else:
					selected_count = new_selected_count
				if selected_count <= gene_target and not final_round:
					final_round = True
				else:
					final_round = False
				counter += 1
			if os.path.exists("{}.txt".format(hypothesis)):
				os.remove("{}.txt".format(hypothesis))
			try:
				shutil.move(os.path.join(args_original.output,"{}_sparsify_final".format(hypothesis),"{}_merged_gene_predictions_final.txt".format(hypothesis)), args_original.output)
				shutil.move(os.path.join(args_original.output,"{}_sparsify_final".format(hypothesis),"{}_merged_gene_predictions_final.png".format(hypothesis)), args_original.output)
			except:
				print("Couldn't find summary graphic file at {}.".format(os.path.join(args_original.output,"{}_sparsify_final".format(hypothesis),"{}_merged_gene_predictions_final.png".format(hypothesis))))

	# generate_gene_prediction_table(weights_filename, responses_filename, groups_filename, features_filename, output_filename)
	# if False:
	# 	pf.generate_gene_prediction_table("angiosperm_out_feature_weights.xml", "angiosperm_20spec_pred.txt", "angiosperm_input/group_indices_angiosperm_input.txt", "angiosperm_input/feature_angiosperm_input.txt", "test_output.txt")
	# if False:
	# 	hypothesis_file_list = pf.generate_hypothesis_set("testtree2.nwk", "nodelist.txt")
	# 	features_filename, groups_filename, response_filename_list = pf.generate_input_matrices("sample_files/angiosperm_100_sample_alns.txt", hypothesis_file_list)
	# 	weights_file_list = pf.run_esl(features_filename, groups_filename, response_filename_list)
	# 	pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename, features_filename)


