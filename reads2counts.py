### Boas Pucker ###
### v2.1 ###
### bpucker@cebitec.uni-bielefeld.de ###

import os, re, sys, glob, time, datetime
from shutil import copyfile
import matplotlib.pyplot as plt

# --- end of imports --- #

__usage__ = """python reads2counts2.py\n
				required:
				--fastq_file_dir <FULL_PATH_TO_DIRECTORY>
				--tmp_cluster_dir <FULL_PATH_TO_TEMPORARY_DIRECTORY_ON_CLUSTER_VOLUME>
				--result_dir <FULL_PATH_TO_RESULT_DIRECTORY>\n
				--ref_gff_file <FULL_PATH_TO_GFF_FILE_MATCHING_THE_PROVIDED_GENOME>
				--ref_genome_file <FULL_PATH_TO_GENOME_FILE_MATCHING_THE_PROVIDED_GFF_FILE>
				
				optional:
				--dissimilarity <FLOAT, 1-identity, value between 0.0 and 1.0>[0.05]
				--length_fraction <FLOAT, value between 0.0 and 1.0>[0.9]
				--para_jobs <INTEGER, number of jobs to be processed at the compute cluster at the same time>[50]
				
				NOTE: all subdirectories of given input directory will be searched as well
				
				bug report and feature request: bpucker@cebitec.uni-bielefeld.de
			"""


def generate_STAR_commands( 	STAR_path,  
							genome_dir,
							fw_rna_seq_read_files,
							rv_rna_seq_read_files,
							output_directories,
							max_missmatches_per_read_length="0.05",
							min_matches_per_read_length="0.9",
							para_jobs=10
						):
	"""! @brief submit STAR mapping jobs to cluster
	
		@param STAR_path (string) is full path to STAR binary file
		@param genome_dir (string) full path to STAR reference genomeDir (without / at the end)
		@param rna_seq_read_files (list) contains full paths to all RNA-seq files in fastq.gz format
		@param output_directories (list) contains full paths to all output directories (/ at the end)
		@param max_missmatches_per_read_length (string) percentage of allowed missmatches in read mapping
		@param min_matches_per_read_length (string) percentage of read length that needs to be mapped
		@param para_jobs (integer) number of mapping jobs that can be lunched in parallel
	"""
	
	print "submitting STAR jobs to cluster ... "
	time.sleep( 10 )
	
	IDs_to_check = []
	files_to_delete = []
	batch_ID = str( datetime.datetime.now() )[-4:]
	for idx, filename in enumerate( sorted( fw_rna_seq_read_files ) ):	#iterate over all provided data files (two)
		time.sleep( 1 )
		
		ID = "S_" + batch_ID + '_' + str( idx ).zfill(3)
		IDs_to_check.append( ID )
		
		# --- setting file names --- #
		err_file = output_directories[idx] + ID + ".err"
		out_file = output_directories[idx] + ID + ".out"
		
		sh_file = output_directories[idx] + ID + ".sh"
		
		files_to_delete.append( err_file )
		STAR_cmd = [ 	STAR_path,
						" --genomeDir ",
						genome_dir,
						" --readFilesIn ",
						filename + "," + rv_rna_seq_read_files[ idx ],
						" --readFilesCommand zcat --runThreadN 3 --outFileNamePrefix ",
						output_directories[idx],
						" --limitBAMsortRAM 80000000000 --outBAMsortingThreadN 1",
						" --outSAMtype BAM SortedByCoordinate --twopassMode Basic",
						" --outFilterMismatchNoverLmax ",
						max_missmatches_per_read_length,
						" --outFilterMatchNminOverLread ",
						min_matches_per_read_length
					]
		STAR_cmd = "".join( STAR_cmd )
		
		with open( sh_file, "w" ) as out:
			out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
															STAR_cmd + '"',
															"| qsub -cwd",
															"-N",
															ID,
															"-l vf=10G",
															"-l arch=lx-amd64",
															"-pe multislot 3",
															#"-l idle=1",
															"-P fair_share",
															"-o",
															out_file,
															"-e",
															err_file
														] ) + '\n'
								 )
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		print "job  " + str( idx ) + " submitted to cluster."
		time.sleep( 1 )
		
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "S_" + batch_ID + "_\d{3}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )	
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "S_" + batch_ID + "_\d{3}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 1 )
	
	for each in files_to_delete:
		try:
			os.remove( each )
		except:
			"file cannot be deleted, because it is missing"
	
	print "STAR read mapping done..."


def submit_featureCount_jobs( 	featureCounts_path,
														annotation_file,
														output_directories,
														count_table_names,
														bam_file_names,
														feature_level="gene",
														feature_ID="ID",
														para_jobs=50
												):
	"""! @brief submit featureCounts jobs to cluster
	
	"""
	
	print "submitting featureCount jobs to cluster ... "
	time.sleep( 10 )
	
	IDs_to_check = []
	files_to_delete = []
	batch_ID = str( datetime.datetime.now() )[-4:]
	for idx, filename in enumerate( sorted( bam_file_names ) ):	#iterate over all provided data files (two)
		time.sleep( 1 )
		
		ID = "F_" + batch_ID + '_' + str( idx ).zfill(3)
		IDs_to_check.append( ID )
		
		# --- setting file names --- #
		err_file = output_directories[idx] + ID + ".report"
		out_file = output_directories[idx] + ID + ".out"
		
		sh_file = output_directories[idx] + ID + ".sh"
		
		files_to_delete.append( out_file )
		
		featureCounts_cmd = [ 	featureCounts_path,
								" -t " + feature_level,
								" -g " + feature_ID,
								" -a " + annotation_file,
								" -O -o " + count_table_names[ idx ],
								" " + filename						
							]
		featureCounts_cmd = "".join( featureCounts_cmd )
		
		with open( sh_file, "w" ) as out:
			out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
															featureCounts_cmd + '"',
															"| qsub -cwd",
															"-N",
															ID,
															"-l vf=1G",
															"-l arch=lx-amd64",
															"-P fair_share",
															"-o",
															out_file,
															"-e",
															err_file
														] ) + '\n'
								 )
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		print "job " + str( idx ) + " submitted to cluster."
		time.sleep( 1 )
		
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "F_" + batch_ID + "_\d{3}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )	
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "F_" + batch_ID + "_\d{3}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 1 )
	
	for each in files_to_delete:
		try:
			os.remove( each )
		except:
			"file cannot be deleted, because it is missing"
	
	print "feature counting done..."

def main( parameters ):
	"""! @brief run all functions """
	
	# --- setting path names and file names --- #
	fastq_file_dir = parameters[ parameters.index( '--fastq_file_dir' )+1 ]
	result_file_dir = parameters[ parameters.index( '--tmp_cluster_dir' )+1 ]
	final_result_dir = parameters[ parameters.index( '--result_dir' )+1 ]
	
	ref_gff_file = parameters[ parameters.index( '--ref_gff_file' )+1 ]
	ref_genome_file = parameters[ parameters.index( '--ref_genome_file' )+1 ]
	
	# --- setting options --- #
	if '--dissimilarity' in parameters:
		max_missmatches_per_read_length = parameters[ parameters.index( '--dissimilarity' )+1 ]
	else:
		max_missmatches_per_read_length="0.05"
	
	if '--length_fraction' in parameters:
		min_matches_per_read_length = parameters[ parameters.index( '--length_fraction' )+1 ]
	else:
		min_matches_per_read_length="0.9"
		
	
	feature_level="gene"
	feature_ID="ID"
	
	STAR_path = "STARlong"
	featureCounts_path = "featureCounts"
	
	
	# --- construct references --- #
	ref_genome_dir = result_file_dir + "reference_genome"
	
	if not os.path.exists( ref_genome_dir ):
		os.makedirs( ref_genome_dir )
	
	#genome dir construction:
	print "constructing STAR references ... "
	cmd1 = "".join( [ 	STAR_path,
						" --runMode genomeGenerate",
						" --genomeDir ",
						ref_genome_dir,
						" --genomeFastaFiles ",
						ref_genome_file,
						" --runThreadN 7 --limitGenomeGenerateRAM 70000000000",
						" --genomeSAindexNbases 4"
					])
	os.popen( cmd1 )
	
	# --- get all RNA-seq read files --- #
	fw_compressed_fastq_files = sorted( glob.glob( fastq_file_dir + '*_1.fastq.gz' ) + glob.glob( fastq_file_dir + '*_1.fq.gz' ) )
	rv_compressed_fastq_files = []
	for fw_file in fw_compressed_fastq_files:
		try:
			rv_compressed_fastq_files.append( fw_file.replace( '_1.fastq.gz', '_2.fastq.gz' ) )
		except:
			rv_compressed_fastq_files.append( fw_file.replace( '_1.fq.gz', '_2.fq.gz' ) )
	
	# --- construct output directories for STAR mappings and predicted names of final bam files --- #
	bam_files = []
	report_files = []
	mapping_result_dirs = []
	fw_rna_seq_read_files = []
	rv_rna_seq_read_files = []
	IDs = []
	
	for idx, filename in enumerate( fw_compressed_fastq_files ):	
		fw_rna_seq_read_files.append( filename )
		rv_rna_seq_read_files.append( rv_compressed_fastq_files[ idx ] )
		ID = filename.split('/')[-1].split('.trimmed')[0].split('_trimmed')[0]
		IDs.append( ID )
		mapping_result_dir =  result_file_dir + ID + "/"
		if not os.path.exists( mapping_result_dir ):
			os.makedirs( mapping_result_dir )
		else:
			pass	#sys.exit( "ERROR: mapping results already exist!" )
		mapping_result_dirs.append( mapping_result_dir )
		bam_files.append( mapping_result_dir + "Aligned.sortedByCoord.out.bam" )
		report_files.append( mapping_result_dir + "Log.final.out" )
	
	# --- run STAR mappings on cluster --- #
	generate_STAR_commands( 	STAR_path,  
												ref_genome_dir,
												fw_rna_seq_read_files,
												rv_rna_seq_read_files,
												mapping_result_dirs,
												max_missmatches_per_read_length,
												min_matches_per_read_length
											)
	
	# --- run featureCounts analysis on cluster --- #
	count_table_names = []
	cont_table_summary_names = []
	for bam_file in bam_files:
		count_table_names.append( bam_file + ".count_table" )
		cont_table_summary_names.append( bam_file + ".count_table.summary" )
	
	submit_featureCount_jobs( 	featureCounts_path,
														ref_gff_file,
														mapping_result_dirs,
														count_table_names,
														bam_files,
														feature_level=feature_level,
														feature_ID=feature_ID,
														para_jobs=50
												)
	
	# --- save all result files --- #
	if not os.path.exists( final_result_dir ):
		os.makedirs( final_result_dir )
	
	for idx, ID in enumerate( IDs ):
		try:
			copyfile( report_files[idx], final_result_dir + ID + "_STAR.logfile" )
			copyfile( count_table_names[idx], final_result_dir + ID + ".count_table" )
			copyfile( cont_table_summary_names[idx], final_result_dir + ID + ".count_table.summary" )
		except:
			print "ERROR while transferring results of " + ID


if __name__ == '__main__':
	
	if '--fastq_file_dir'  in sys.argv and '--tmp_cluster_dir' in sys.argv and '--result_dir' in sys.argv and '--ref_gff_file' in sys.argv and '--ref_genome_file' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	print "all done!"
