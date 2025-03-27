import argparse
from datetime import datetime
import subprocess
import os
import re
import shutil
import subprocess
import sys
import warnings

def generate_pipeline_setup(
	pls_path=None, annotation_path=None, options_path=None, n_cores=None,
	work_path='/work', n_cores_default=4, dummy_x11=False,
):
	if not pls_path:
		pls_path = work_path + '/pipelineSetup.R'

	if os.path.isfile(pls_path):
		warnings.warn(f'Pipeline setup path {pls_path} already exists. Backing up existing file.')
		pls_path_bak = pls_path
		while os.path.isfile(pls_path_bak):
			pls_path_bak_match = re.fullmatch(r'.*[.]bak(\d+)', pls_path_bak)
			if not pls_path_bak_match:
				pls_path_bak = pls_path_bak + '.bak0'
			else:
				iter_old = pls_path_bak_match.group(1)
				iter_new = str(int(iter_old) + 1)
				pls_path_bak = pls_path_bak + '.bak' + iter_new
		shutil.copyfile(pls_path, pls_path_bak)

	shutil.copyfile('/home/hippo/pipelineSetup.R', pls_path)

	_set_pipeline_setup(
		pls_path,
		annotation_path=annotation_path if annotation_path else work_path + '/Annotation.txt',
		options_path=options_path if options_path else work_path + '/Options.txt',
		n_cores=n_cores if n_cores else n_cores_default,
		warn=False, dummy_x11=dummy_x11,
	)

	return pls_path

def quick_qc(
	sample_ids, work_path='/work', annotation_path=None, fastq_path=None, options_path=None,
	pipeline_setup_path=None, n_cores=4,
):
	if not sample_ids:
		raise ValueError('Please provide specific sample IDs to run "quickQC".')

	tokens = {
		'@SAMPLE_SET@': '", "'.join(sample_ids)
	}

	os.chdir(work_path)
	_rscript(
		'quickqc', tokens, work_path, annotation_path, fastq_path, options_path,
		pipeline_setup_path, dummy_x11=True, n_cores=n_cores,
	)

def rnaseq_batch(
	sample_ids, work_path='/work', annotation_path=None, fastq_path=None, options_path=None,
	pipeline_setup_path=None, n_cores=4,
):
	if not sample_ids:
		raise ValueError('Please provide specific sample IDs to analyze.')

	tokens = {
		'@SAMPLE_SET@': '", "'.join(sample_ids)
	}

	os.chdir(work_path)
	_rscript(
		'rnaseqbatch', tokens, work_path, annotation_path, fastq_path, options_path,
		pipeline_setup_path, dummy_x11=True, n_cores=n_cores,
	)

def rnaseq_finalize(
	sample_ids, work_path='/work', annotation_path=None, fastq_path=None, options_path=None,
	pipeline_setup_path=None, n_cores=4, pairwise_comparisons=[],
):
	tokens = {
		'@SAMPLE_SET@': '", "'.join(sample_ids)
	}

	# TODO: need to handle differential expression here!
	# slice annotation_file based on the pairwise_comparisons provided
	# also add that info to the tokens somehow

	os.chdir(work_path)
	_rscript(
		'rnaseqfin', tokens, work_path, annotation_path, fastq_path, options_path,
		pipeline_setup_path, dummy_x11=True, n_cores=n_cores,
	)

def _check_annotation(annotation_path, work_path='/work'):
	if not annotation_path:
		annotation_path = work_path + '/Annotation.txt'

	if not os.path.isfile(annotation_path):
		raise ValueError(f'Annotation file {annotation_path} does not exist or is not a file')
	# maybe validate the contents, too??
	return annotation_path

def _generate_options(fastq_path=None, options_path='/work/Options.txt', n_cores=None):
	shutil.copyfile('/home/hippo/Options.txt', options_path)
	if n_cores is not None:
		_set_option(options_path, 'nCores', n_cores)
	if fastq_path is not None:
		_set_option(options_path, 'fastqData.path', fastq_path)
	return options_path, _read_options(options_path)

def _get_options(options_path, fastq_path=None, n_cores=None, n_cores_default=4):
	if not os.path.isfile(options_path):
		n_cores = n_cores if n_cores is not None else n_cores_default
		options_path, options = _generate_options(
			fastq_path=fastq_path, options_path=options_path, n_cores=n_cores,
		)
	else:
		options = _read_options(options_path)
		reset = False

		if n_cores is not None and options['nCores'] != n_cores:
			warnings.warn(f'Changing Options.txt file to match the supplied number of cores')
			_set_option(options_path, 'nCores', n_cores)
			reset = True

		if fastq_path is not None and options['fastqData.path'] != fastq_path:
			warnings.warn(f'Changing Options.txt file to match the supplied FASTQ path')
			_set_option(options_path, 'fastqData.path', fastq_path)
			reset = True

		if reset:
			options = _read_options(options_path) # read again to ensure we report the true options
	return options

def _prep_pipeline_setup(
	pl_setup_path, annotation_path, work_path='/work', options_path=None, n_cores=None,
	dummy_x11=False,
):
	if not os.path.isfile(pl_setup_path):
		generate_pipeline_setup(
			pl_setup_path, annotation_path=annotation_path, options_path=options_path,
			work_path=work_path, n_cores=n_cores, dummy_x11=dummy_x11,
		)
	else:
		_set_pipeline_setup(
			pl_setup_path, annotation_path=annotation_path, options_path=options_path,
			n_cores=n_cores, dummy_x11=dummy_x11,
		)
	return pl_setup_path

def _read_options(options_path):
	with open(options_path, 'r') as options_file:
		return dict(line.rstrip('\n').split('\t') for line in options_file.readlines())

def _rscript(
	name, tokens, work_path, annotation_path, fastq_path, options_path, pipeline_setup_path,
	n_cores=None, dummy_x11=False,
):
	# validate/prep configuration files

	annotation_path = _check_annotation(annotation_path, work_path=work_path)
	options = _get_options(options_path, fastq_path=fastq_path, n_cores=n_cores)
	pipeline_setup_path = _prep_pipeline_setup(
		pipeline_setup_path, annotation_path=annotation_path, work_path=work_path,
		options_path=options_path, n_cores=n_cores, dummy_x11=dummy_x11,
	)

	tokens.update({
		'@ANNOT_PATH@': annotation_path,
		'@OPTS_PATH@': options_path,
		'@SETUP_PATH@': pipeline_setup_path,
	})

	# prep R script

	with open(f'/home/hippo/{name}.R', 'r') as templ_file:
		template = templ_file.read()

	for key, value in tokens.items():
		template = re.sub(key, value, template)

	script_path = f'{work_path}/{name}_{datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f")}.R'
	with open(script_path, 'w') as script_file:
		script_file.write(template)

	# run R script

	completed_process = subprocess.run([
		'R', '--no-save', '--no-restore', '--no-site-file', '--no-init-file', '-f', script_path,
	], check=True)

def _set_option(options_path, key, value):
	subprocess.run(['sed', '-i', f's/{key}\t.*/{key}\t{value}/', options_path])

def _set_pipeline_setup(
	pl_setup_path, annotation_path=None, options_path=None, n_cores=None, warn=True,
	dummy_x11=False,
):
	with open(pl_setup_path, 'r') as pl_setup_file:
		pl_setup_content = pl_setup_file.read()

	if annotation_path is not None:
		annotation_path_match = re.search(r'readAnnotationTable\("(.+)"\)', pl_setup_content)
		if not annotation_path_match:
			raise ValueError('Error in pipeline setup script: no annotation table')
		annotation_path_old = annotation_path_match.group(1)
		if annotation_path_old != annotation_path:
			if warn:
				warnings.warn(
					f'Changing pipeline setup: annotation path from {annotation_path_old} ' + \
						f'to {annotation_path}')
			subprocess.run([
				'sed', '-i',
				f's;"{annotation_path_old}";' + \
					f'"{annotation_path}";',
				pl_setup_path
			])

	if options_path is not None:
		options_path_match = re.search(r'bowtie2Par.defaults\("(.+)"\)', pl_setup_content)
		if not options_path_match:
			raise ValueError('Error in pipeline setup script: no options table')
		options_path_old = options_path_match.group(1)
		if options_path_old != options_path:
			if warn:
				warnings.warn(
					f'Changing pipeline setup: options path from {options_path_old} ' + \
						f'to {options_path}')
			subprocess.run([
				'sed', '-i',
				f's;"{options_path_old}";' + \
					f'"{options_path}";',
				pl_setup_path
			])

	if n_cores is not None:
		n_cores_match = re.search(r'multicore[.]setup\(max[.]cores=([\d+])', pl_setup_content)
		if not n_cores_match:
			raise ValueError('Error in pipeline setup script: no multicore setup')
		n_cores_old = n_cores_match.group(1)
		if int(n_cores_old) != n_cores:
			if warn:
				warnings.warn(f'Changing pipeline setup: n cores from {n_cores_old} to {n_cores}')
			subprocess.run([
				'sed', '-i',
				f's/multicore.setup.max.cores={n_cores_old}/' + \
					f'multicore.setup(max.cores={n_cores}/',
				pl_setup_path
			])

	if dummy_x11:
		with open(pl_setup_path, 'a') as pl_setup_file:
			pl_setup_file.write('''

checkX11 <- function(bg="white", type="dbcairo", width=10, height=7, ...) {
	dev <- png(filename="Rplot%03d.png", type="cairo-png")
	.Last <- function() dev.off()
}
'''
)

#
# main
#

def main(
	command, sample_ids=None, work_path='/work', annotation_path=None, fastq_path=None,
	options_path=None, pipeline_setup_path=None, n_cores=None,
	pairwise_comparisons=[],
):
	if command == 'quickqc':
		quick_qc(
			sample_ids, work_path=work_path, annotation_path=annotation_path, fastq_path=fastq_path,
			options_path=options_path, pipeline_setup_path=pipeline_setup_path, n_cores=n_cores,
		)
	elif command == 'rnaseqbatch':
		rnaseq_batch(
			sample_ids, work_path=work_path, annotation_path=annotation_path, fastq_path=fastq_path,
			options_path=options_path, pipeline_setup_path=pipeline_setup_path, n_cores=n_cores,
		)
	elif command == 'rnaseqfinalize':
		rnaseq_finalize(
			sample_ids, work_path=work_path, annotation_path=annotation_path, fastq_path=fastq_path,
			options_path=options_path, pipeline_setup_path=pipeline_setup_path, n_cores=n_cores,
			pairwise_comparisons=pairwise_comparisons,
		)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description=(
			'Utility for executing NGS alignment and analysis using DuffyNGS.'
		)
	)

	parser.add_argument('command',
		choices=['quickqc', 'rnaseqbatch', 'rnaseqfinalize'],
		help='The DuffyNGS subcommand to execute. Choices are %(choices)s.',
	)

	parser.add_argument('-s', '--sample-ids',
		default=[],
		nargs='*',
		help=(
			'The list of samples to be included in the selected alignment/analysis. Required in '
			'the case of quickQC and batched RNA-seq analysis executions.'
		)
	)

	parser.add_argument('-a', '--annotation-path',
		default='/work/Annotation.txt',
		help=(
			'The file system location (within the docker image) where the main annotation file can '
			'be found, which lists the relevant samples (and their associated metadata). Default '
			'is %(default)s. A valid annotation file is required for the pipeline to proceed.'
		)
	)
	parser.add_argument('-f', '--fastq-path',
		default='/work/fastqs',
		help=(
			'The file system location (within the docker image) where the sequencing FASTQ files '
			'can be found (i.e., the directory containing said files). Default is %(default)s. '
		)
	)
	parser.add_argument('-o', '--options-path',
		default='/work/Options.txt',
		help=(
			'The file system location (within the docker image) where the options config file can '
			'be found, which sets parameters for the alignment software. Default is %(default)s. '
			'If a file is not found at this location, a default one is generated.'
		)
	)
	parser.add_argument('-p', '--pipeline-setup-path',
		default='/work/pipelineSetup.R',
		help=(
			'The file system location (within the docker image) where the pipeline setup script '
			'can be found, which sets some helpful parameters for the R scripts. Default is '
			'%(default)s. If a file is not found at this location, a default one is generated.'
		)
	)
	parser.add_argument('-w', '--work-path',
		default='/work',
		help=(
			'The file system location (within the docker image) where files should be read and '
			'written. Default is %(default)s.'
		)
	)

	parser.add_argument('-n', '--n-cores',
		type=int,
		help=(
			'The number of CPU cores (per job, if applicable) to be used.'
		)
	)

	# pairwise_comparisons not implemented for the moment

	args = parser.parse_args(sys.argv[1:])
	args_dict = vars(args)
	main(**args_dict)
