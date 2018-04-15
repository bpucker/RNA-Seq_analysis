### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python find_housekeeping_genes.py
	--in <FULL_PATH_TO_EXPRESSION_FILE>
	--out <FULL_PATH__TO_OUTPUT_FILE>
	optional:
	--anno <FULL_PATH__TO_ANNOTATION_FILE>
	--cutoff <MINIMAL_EXPRESSION_PER_SAMPLE(INTEGER)>[50]
					"""

from operator import itemgetter
import numpy as np
import sys

# --- end of imports --- #

def load_expression_values( filename, expression_cutoff ):
	"""! @brief load all expression values """
	
	expression_data = []
	with open( filename, "r" ) as f:
		tissues = f.readline().strip().split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression = []
			for idx, each in enumerate( parts[1:] ):
				expression.append( float( parts[ idx+1 ] ) )
			if sum( expression ) > len( expression )*expression_cutoff:
				exp_var = False
				exp_var = np.std( expression )	/ np.mean( expression )
				expression_data.append( { 'id': parts[0], 'expvar': exp_var, 'mean': np.mean( expression ),  'std': np.std( expression ) } )
			line = f.readline()
	return expression_data


def load_annotation( annotation_file ):
	"""! @brief load annotation mapping table """
	
	annotation_mapping_table = {}
	
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			annotation_mapping_table.update( { parts[0]: parts[1] } )
			line = f.readline()
	return annotation_mapping_table


def main( arguments ):
	"""! @brief run all parts of this script """
	
	expression_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	if '--anno' in arguments:
		annotation_file = arguments[ arguments.index( '--anno' )+1 ]
		annotation_mapping_table = load_annotation( annotation_file )
	else:
		annotation_mapping_table = {}
	
	if '--cutoff' in arguments:
		expression_cutoff = int( arguments[ arguments.index( '--cutoff' )+1 ] )
	else:
		expression_cutoff = 50
	
	gene_expression = load_expression_values( expression_file, expression_cutoff )
	
	print "number of considered genes: " + str( len( gene_expression ) )
	sorted_genes = sorted( gene_expression, key=itemgetter('expvar') )
	with open( output_file, "w" ) as out:
		out.write( 'ID\tExpVar\tMean\tStandardDeviation\tAnnotation\n' )
		for gene in sorted_genes:
			try:
				out.write( gene['id'] + '\t' + str( gene['expvar'] ) + '\t' + str( gene['mean'] ) + '\t' + str( gene['std'] ) + '\t' + annotation_mapping_table[ gene['id'] ] +  '\n' )
			except KeyError:
				out.write( gene['id'] + '\t' + str( gene['expvar'] ) + '\t' + str( gene['mean'] ) + '\t' + str( gene['std'] ) + '\tN/A\n' )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
