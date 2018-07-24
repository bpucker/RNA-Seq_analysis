### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
	python get_mapping_stats.py
	--in <FULL_PATH_TO_INPUT_DIRECTORY>
	--out <FULL_PATH_TO_OUTPUT_FILE>
					"""

import glob, re, sys

# --- end of imports --- #

def get_stats_from_file( filename ):
	"""! @brie get al relevant stats from given mapping report file """
	
	stats = {}
	with open( filename, "r" ) as f:
		line = f.readline()
		while line:
			if 'Number of input reads' in line:
				value = int( re.findall( "\d+", line )[0] )
				stats.update( { 'total': value } )
			elif 'Uniquely mapped reads number' in line:
				value = int( re.findall( "\d+", line )[0] )
				stats.update( { 'unique': value } )
			elif 'Number of reads mapped to multiple loci' in line:
				value = int( re.findall( "\d+", line )[0] )
				stats.update( { 'multi': value } )
			line = f.readline()
	return stats

def main( arguments ):
	"""! @brief runs everything """
	
	input_dir =arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]

	all_files = glob.glob( input_dir + "*.logfile" ) + glob.glob( input_dir + "/*.logfile" ) 
	all_stats = {}
	for filename in all_files:
		ID = filename.split('/')[-1]
		all_stats.update( { ID: get_stats_from_file( filename ) } )

	with open( output_file, "w" ) as out:
		out.write( 'Sample\tTotalReads\tUniquelyMappedReads\tUniquelyMappedReads%\tMultiMappedReads\tMultiMappedReads%\n' )
		for ID in sorted( all_stats.keys() ):
			out.write( "\t".join( map( str, [ 	ID,
																all_stats[ ID ]['total'],
																all_stats[ ID ]['unique'],
																str( ( 100.0*all_stats[ ID ]['unique'] ) / all_stats[ ID ]['total'] )[:5]+"%",
																all_stats[ ID ]['multi'],
																str( ( 100.0*all_stats[ ID ]['multi'] ) / all_stats[ ID ]['total'] )[:5]+"%" ] )
											 ) + '\n' )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
