### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python construct_gene_expression_plots.py
	--candidates <FULL_PATH_TO_CANDIDATE_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	at least one expression data file is required:
	--counts <FULL_PATH_TO_RAW_COUNT_TABLE>
	--tpms <FULL_PATH_TO_TPM_FILE>
	--fpkms <FULL_PATH_TO_FPKM_FILE>
		
	optional:
	--samples <FULL_PATH_TO_SAMPLE_ORDER_FILE>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import matplotlib.pyplot as plt
import os, sys
import numpy as np

# --- end of imports --- #


def load_expression( count_table ):
	"""! @brief load all data from given count table """
	
	expression = {}
	with open( count_table, "r" ) as f:
		samples = f.readline().strip().split('\t')[1:]
		if len( samples ) != len( list( set( samples ) ) ):
			critical_elements = []
			for element in list( set( samples ) ):
				if samples.count( element ) > 1:
					critical_elements.append( element )
			sys.exit( "ERROR: duplicate sample names detected! => " + str( sorted( list( set( critical_elements ) ) ) ) )
			
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				expression[ parts[0] ]
				sys.exit( "ERROR: duplicate gene entries detected!! => " + parts[0] )
			except KeyError:
				gene_exp = {}
				for idx, value in enumerate( parts[1:] ):
					gene_exp.update( { samples[ idx ]: float( value ) } )
				expression.update( { parts[0]: gene_exp } )
			line = f.readline()
	return expression


def construct_expression_plot( fig_file, expression, exp_type, candidate, sample_order, group_info ):
	"""! @brief construct gene expression plot """
	
	fig, ax = plt.subplots( figsize=(10,5) )
	
	# --- add data points --- #
	val_to_plot = expression[ candidate ]
	all_y_values = []
	for idx, sample in enumerate( sample_order ):
		ax.plot( [ idx, idx ], [ 0, expression[candidate][ sample ] ], linewidth=5, color="blue" )
		all_y_values.append( expression[candidate][ sample ] )
	
	# --- add group infos --- #
	group_off_set = 0
	for group in group_info:
		x = group_off_set+group['len']+0.5
		ax.plot( [ x, x ], [ 0, max( all_y_values ) ], linewidth=3, linestyle=":", color="black" )
		ax.text( x-( group['len']*0.5 ), max( all_y_values ), group['id'], color="black", ha="center", va="center" )
		group_off_set += group['len']
	
	ax.set_title( exp_type + " of " + candidate )
	
	ax.set_ylabel( "gene expresion ("+exp_type+")" )
	
	ax.set_xlim( -1, len( sample_order ) )
	start, stop = ax.get_xlim()
	xticks = np.arange( start, stop, 1 )
	ax.set_xticks( xticks )
	ax.set_xticklabels( [""] + sample_order, rotation=90 )
	
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	
	plt.subplots_adjust( left=0.08, right=0.99, top=0.92, bottom=0.15 )
	fig.savefig( fig_file, dpi=600 )
	plt.close("all")


def main( arguments ):
	"""! @brief runs everything """
	
	candidate_file = arguments[ arguments.index( '--candidates' )+1 ]	#one gene ID per line
	outdir = arguments[ arguments.index( '--out' )+1 ]
	
	if outdir[-1] != '/':
		outdir += "/"
	if not os.path.exists( outdir ):
		os.makedirs( outdir )
	
	# --- load candidate genes --- #
	with open( candidate_file, "r" ) as f:
		candidates = []
		for each in f.read().strip().split('\n'):
			candidates.append( each.strip() )
	
	# ---- load expression values --- #
	data_files = []
	if '--counts' in arguments:
		data_files.append( { 'id': 'rawCounts', 'file': arguments[ arguments.index( '--counts' )+1 ] } )
	if '--tpms' in arguments:
		data_files.append( { 'id': 'TPMs', 'file': arguments[ arguments.index( '--tpms' )+1 ] } )
	if '--fpkms' in arguments:
		data_files.append( { 'id': 'FPKMs', 'file': arguments[ arguments.index( '--fpkms' )+1 ] } )
	
	for idx, entry in enumerate( data_files ):
		data_files[idx].update( { 'expression': load_expression( entry['file'] ) } )
	
	for entry in data_files:
		if len( entry['expression'] ) < 1:
			sys.exit( "ERROR: no expression loaded! => " + entry['id'] )
	
	# --- load sample order --- #
	sample_order = []
	group_order = []
	group_length = {}
	if '--samples' in arguments:
		sample_order_file = arguments[ arguments.index( '--samples' )+1 ]	#all sample names to include in desired order
		with open( sample_order_file, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				try:
					entry['expression'].values()[0][ parts[0] ]
					sample_order.append( parts[0] )
				except KeyError:
					sys.exit( "ERROR: sample name given in sample order file is missing! => " + parts[0] )
				if len( parts ) > 1:
					if parts[1] not in group_order:
						group_order.append( parts[1] )
						group_length.update( { parts[1]: 1 } )
					else:
						group_length[ parts[1] ] += 1
				line = f.readline()
		group_info = []
		for group in group_order:
			group_info.append( { 'id': group, 'len': group_length[ group ] } )
	else:
		sample_order = sorted( entry['expression'].values()[0].keys() )
		
	# --- construct expression plots --- #
	for candidate in candidates:
		for entry in data_files:
			try:
				entry['expression'][ candidate ]
			except KeyError:
				sys.exit( "ERROR: candidate gene not present in expression data set! => " + candidate )
			fig_file = outdir + candidate + '_' + entry['id'] + ".png"
			construct_expression_plot( fig_file, entry['expression'], entry['id'], candidate, sample_order, group_info )


if __name__ == '__main__':
	
	if '--candidates' in sys.argv and '--counts' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	elif '--candidates' in sys.argv and '--tpms' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	elif '--candidates' in sys.argv and '--fpkms' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
