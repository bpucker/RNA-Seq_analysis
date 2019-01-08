### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """
					python qPCR_ref_gene_finder.py
					--exp <FULL_PATH_TO_EXPRESSION_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""
import numpy as np
from operator import itemgetter
import sys

# --- end of imports --- #

def load_expression_values( data_file ):
	"""! @brief load expression values """
	
	exp_data = {}
	with open( data_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			exp_data.update( { parts[0]: map( float, parts[1:] ) } )
			line = f.readline()
	return exp_data


def main( arguments ):
	"""! @brief runs everything """
	
	data_file = arguments[ arguments.index( '--exp' ) +1 ]
	report_file = arguments[ arguments.index( '--out' ) +1 ]
	
	expression = load_expression_values( data_file )
	
	exp_cutoff = 10
	gene_characteristics = []
	for gene in expression.keys():
		avg = np.median( expression[ gene ] )
		if avg > 0:
			dev = np.std( expression[ gene ] ) / avg
			if avg > exp_cutoff and min(  expression[ gene ] ) > exp_cutoff:
				gene_characteristics.append( { 'id': gene, 'avg': avg, 'dev': dev } )
	gene_characteristics.sort( key=itemgetter( 'avg' ), reverse=True )
	gene_characteristics.sort( key=itemgetter( 'dev' ) )
	
	
	with open( report_file, "w" ) as out:
		out.write( "GeneID\tNormalizedStandardDeviation\tAverageExpression\n" )
		for candidate in gene_characteristics:
			out.write( "\t".join( map( str, [ candidate['id'],  candidate['dev'], candidate['avg'] ] ) ) + '\n' )


if __name__ == '__main__':
	
	if '--exp' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
