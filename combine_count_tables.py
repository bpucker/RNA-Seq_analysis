### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
	python combine_count_tables.py
	--in <FULL_PATH_INPUT_FILE_DIRECTORIES(S)>
	(multiple directories can be provided comma-separated)
	--gff <FULL_PATH_TO_GFF3_FILE>
	--out <FULL_PATH_TO_OUTPUT_FILE>
					"""

import glob, re, sys, os
from operator import itemgetter

# --- end of imports --- #

def load_raw_counts( count_table, raw_counts ):
	"""! @brief load all data from given count table """
	
	with open( count_table, "r" ) as f:
		f.readline()
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				values = raw_counts[ parts[0] ]
				del raw_counts[ parts[0] ]
				values.append( int( parts[-1] ) )
				raw_counts.update( { parts[0]: values } )
			except KeyError:
				raw_counts.update( { parts[0]: [ int( parts[-1] ) ] } )
			line = f.readline()
	return raw_counts


def load_TPMs( count_table, TPMs ):
	"""! @brief load all data from given count table """
	
	lib_size_counter = 0
	with open( count_table, "r" ) as f:
		f.readline()
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			lib_size_counter += int( parts[ -1 ] )
			line = f.readline()
	if lib_size_counter > 0:
		factor = 1000000.0  / lib_size_counter
		with open( count_table, "r" ) as f:
			f.readline()
			f.readline()
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				try:
					values = TPMs[ parts[0] ]
					del TPMs[ parts[0] ]
					values.append( int( parts[-1] )*factor )
					TPMs.update( { parts[0]: values } )
				except KeyError:
					TPMs.update( { parts[0]: [ int( parts[-1] )*factor ] } )
				line = f.readline()
	else:
		print count_table
	return TPMs


def get_collapsed_exon_length( exons ):
	"""! @brief collapse overlapping exon fragments to get effective exon length """
	
	pos_to_sort = []
	for exon in exons:
		pos_to_sort.append( { 'start': exon[0], 'end': exon[1] } )
	sorted_positions = sorted( pos_to_sort, key=itemgetter('start') )
	final_exon_positions = []
	for element in sorted_positions:
		if len( final_exon_positions ) == 0:
			final_exon_positions.append( element )
		elif element['start'] < final_exon_positions[ -1 ]['end'] and element['end'] > final_exon_positions[ -1 ]['end']:
			final_exon_positions[-1] = { 'start': final_exon_positions[ -1 ]['start'], 'end': element['end'] }
		elif element['start'] > final_exon_positions[ -1 ]['end']:
			final_exon_positions.append( element )
	total_length = 0
	for exon in final_exon_positions:
		total_length += exon['end']-exon['start']
	return total_length


def get_gene_lengths( gff3_file ):
	"""! @brief get gene lengths """
	
	rna_to_gene = {}
	exons_per_gene = {}
	with open( gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "mRNA", "transcript", "lnc_RNA", "tRNA", "rRNA" ]:
					ID = re.findall( "rna\d+", parts[-1] )[0]
					parent = re.findall( "gene\d+", parts[-1] )[0]
					rna_to_gene.update( { ID: parent } )
				elif parts[2] in  [ "exon", "CDS", "five_prime_UTR", "three_prime_UTR"  ]:
					start, end = int( parts[3] ), int( parts[4] )
					try:
						parent = re.findall( "rna\d+", parts[-1] )[0]
						try:
							exons_per_gene[ rna_to_gene[ parent ] ].append( ( start, end ) )
						except KeyError:
							exons_per_gene.update( { rna_to_gene[ parent ]: [ (start, end) ] } )
					except IndexError:
						parent = re.findall( "gene\d+", parts[-1] )[0]
						try:
							exons_per_gene[ parent ].append( ( start, end ) )
						except KeyError:
							exons_per_gene.update( { parent: [ (start, end) ] } )
			line = f.readline()
	
	# --- calculate exon length per gene by collapsing overlaps --- #
	length_per_gene = {}
	for key in exons_per_gene.keys():
		exons = exons_per_gene[ key ]
		length_per_gene.update( { key: get_collapsed_exon_length( exons ) } )
	return length_per_gene


def main( arguments ):
	"""! @brief run all parts of this script """
	
	input_directories = []
	input_files = arguments[ arguments.index('--in')+1 ]
	if not ',' in input_files:
		input_directories.append( input_files )
	else:
		for each in input_files.split(','):
			if len( each ) > 3:
				input_directories.append( each )
	
	gff3_file = arguments[ arguments.index('--gff')+1 ]
	prefix = arguments[ arguments.index('--out')+1 ]
	
	if prefix[-1] != '/':
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	input_files = []
	for each in input_directories:
		input_files += glob.glob( each + "*.count_table" ) + glob.glob( each + "*/*.count_table" )
		input_files += glob.glob( each + "*.countTable.txt" ) + glob.glob( each + "*/*.countTable.txt" )
	
	print "number of identified count tables: " + str( len( input_files ) )
	
	try:
		gene_lengths = get_gene_lengths( gff3_file )
	except:
		gene_lengths = {}
	
	raw_counts = {}
	TPMs = {}
	samples = []
	for each in input_files:
		samples.append( each.split('/')[-1].split('.')[0] )
		raw_counts = load_raw_counts( each, raw_counts )
		TPMs = load_TPMs( each, TPMs )
	
	sample_names_file = prefix + "SAMPLE_NAMES.txt"
	with open( sample_names_file, "w" ) as out:
		out.write( "\n".join( samples ) )
	
	output_raw_counts_file = prefix + "raw_counts.txt"
	output_TPM_file = prefix + "TPMs.txt"
	output_FPKM_file = prefix + "FPKMs.txt"
	
	with open( output_raw_counts_file, "w" ) as out:
		out.write( "\t".join( [ "gene" ] + samples ) + '\n' )
		for gene in sorted( raw_counts.keys() ):
			out.write( "\t".join( map( str, [ gene ] + raw_counts[ gene ] ) ) + '\n' )
	
	with open( output_TPM_file, "w" ) as out:
		out.write( "\t".join( [ "gene" ] + samples ) + '\n' )
		for gene in sorted( TPMs.keys() ):
			out.write( "\t".join( map( str, [ gene ] + TPMs[ gene ] ) ) + '\n' )
	
	if len( gene_lengths.keys() ) > 0:
		with open( output_FPKM_file, "w" ) as out:
			out.write( "\t".join( [ "gene" ] + samples ) + '\n' )
			for gene in sorted( TPMs.keys() ):
				values = TPMs[ gene ]
				FPKMs = []
				try:
					gene_len = gene_lengths[ gene ]
					for value in values:
						FPKMs.append( ( 1000.0*value ) / gene_len )
					out.write( "\t".join( map( str, [ gene ] + FPKMs ) ) + '\n' )
				except KeyError:
					print gene


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
