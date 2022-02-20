### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python construct_DESeq2_input.py
	--counts <FULL_PATH_TO_OUTPUT_COUNT_TABLE>
	--info <FULL_PATH_TO_SAMPLE_INFO_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import glob, sys, re, os
from datetime import date
from operator import itemgetter

# --- end of imports --- #

def load_expression_values( filename, ID ):
	"""! @brief load all expression values from featureCounts result file and return raw counts """
	
	expression_data = {}
	
	with open( filename, "r" ) as f:
		header = f.readline().strip().split('\t')
		idx_of_interest = header.index( ID )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression_data.update( { parts[0]: float( parts[ idx_of_interest ] ) } )
			line = f.readline()
		
	return expression_data


def load_sample_information( sample_information_file, combined_count_table ):
	"""! @brief load information about all samples into list of dictionaries """
	
	sample_information = []
	
	with open( sample_information_file, "r" ) as f:
		f.readline()	#header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			sample_information.append( { 	'name': parts[0],
																'expression':load_expression_values( combined_count_table, parts[0] ),
																'replicate': parts[2],
																'genotype': parts[1],
																
															} )
			line = f.readline()
	print "number of loaded sample information data sets: " + str( len( sample_information ) )
	return sample_information


def construct_sample_table( sample_information, sample_table ):
	"""! @brief construct sample table for DeSeq2 """
	
	# --- construct sample table --- #
	ordered_valid_samples = sorted( sample_information, key=itemgetter( 'genotype', 'replicate' ) )
	with open( sample_table, "w" ) as out:
		out.write( '\t'.join( [ "",  'genotype',  'replicate' ] ) + '\n' )
		for idx, sample in enumerate( ordered_valid_samples ):
			new_line = "\t".join( [ sample['name'],			#sample name
												sample['genotype'],
												sample['replicate']		#replicate
												] ) + '\n'
			out.write( new_line )
	
	return ordered_valid_samples


def construct_data_matrix( ordered_valid_samples, data_matrix_file ):
	"""! @brief construct the data matrix for DESeq2 based on featureCount result files """
	
	with open( data_matrix_file, "w" ) as out:
		
		# --- construct header line (all sample IDs in order) --- #
		ordered_valid_ids = []
		for sample in ordered_valid_samples:
			ordered_valid_ids.append( sample['name'] )
		out.write( "\t".join( [ "" ] + ordered_valid_ids ) + '\n' )
		
		# --- construct data body (expression data per sample in one column) --- #
		gene_IDs = sorted( ordered_valid_samples[0]['expression'].keys() )
		for gene in gene_IDs:			
			new_line = [ gene ]	#add geneID
			for sample in ordered_valid_samples:
				new_line.append( int( sample[ 'expression' ][ gene ] ) )
			out.write( '\t'.join( map( str, new_line ) ) + '\n' )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	combined_count_table = arguments[ arguments.index('--counts')+1 ]	#path to data table
	sample_information_file = arguments[ arguments.index('--info')+1 ]	#path to table with meta information
	
	prefix = arguments[ arguments.index('--out')+1 ]	#prefix for output
	if prefix[ -1 ] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	# --- nothing needs to be changed beyond this point --- #
	if prefix[-1] != '/':
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	sample_table = prefix + "clean_sample_table.txt"
	data_matrix_file = prefix + "clean_data_matrix.txt"
	
	
	sample_information = load_sample_information( sample_information_file, combined_count_table )
	
	ordered_valid_samples = construct_sample_table( sample_information, sample_table )
	
	construct_data_matrix( ordered_valid_samples, data_matrix_file )


if '--counts' in sys.argv and '--info' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
