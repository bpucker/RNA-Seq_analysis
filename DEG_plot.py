### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python DEG_plot.py
	--in <FULL_PATH_TO_INPUT_DIRECTORY>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, glob
from operator import itemgetter
import matplotlib.pyplot as plt

# --- end of imports --- #

def load_data( filename ):
	"""! @brief load expression data from DESeq2 output file """
	
	# --- load data from file --- #
	data = []
	with open( filename, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data.append( { 'id': parts[0], 'baseMean': float( parts[1] ), 'logFC': float( parts[2] ), 'padj': float( parts[6] ), 'name': parts[7] } )
			except IndexError:
				data.append( { 'id': parts[0], 'baseMean': float( parts[1] ), 'logFC': float( parts[2] ), 'padj': float( parts[6] ), 'name': parts[0] } )
			line = f.readline()
	
	# --- split into up and down --- #
	up = []
	down = []
	for entry in data:
		if entry['logFC'] > 0:
			up.append( entry )
		elif entry['logFC'] < 0:
			down.append( entry )
	
	up = sorted( up, key=itemgetter('padj') )
	down = sorted( down, key=itemgetter('padj') )
	
	return up, down


def construct_plot( up, down, fig_file, ID ):
	"""! @brief construct plot with selected genes """
	
	fig, ax = plt.subplots( )
	up_values = [ 0 ]
	for idx, gene in enumerate( up ):
		ax.plot( [ 0, gene['logFC'] ], [ idx, idx ], color="green", linewidth=3 )
		up_values.append( gene['logFC'] )
		ax.text( gene['logFC']+0.1, idx, gene['name'], ha="left" )
	down_values = [ 0 ]
	for idx, gene in enumerate( down ):
		ax.plot( [ 0, gene['logFC'] ], [ idx, idx ], color="red", linewidth=3 )
		down_values.append( gene['logFC'] )
		ax.text( gene['logFC']-0.1, idx, gene['name'], ha="right" )
	
	ax.set_title( ID )
	ax.set_xlabel( "log2FC" )
	
	ax.set_xlim( min( down_values )-2, max( up_values )+2 )
	ax.set_ylim( -0.5, max( [ len( down_values ), len( up_values ) ] )-0.5 )
	
	ax.plot( [ 0,0 ], [ -0.5, max( [ len(up), len(down) ] )-0.5 ], color="black", linewidth=3 )	#new Y axis
	
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.get_yaxis().set_ticks([])
	
	plt.subplots_adjust( left=0.01, right=0.99, top=0.9, bottom=0.1 )	
	fig.savefig( fig_file, dpi=300 )
	plt.close("all")


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	input_dir = arguments[ arguments.index( '--in' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	
	if '--gene_num' in arguments:
		gene_number = int( arguments[ arguments.index( '--gene_num' )+1 ] )
	else:
		gene_number = 10
	
	input_files = glob.glob( input_dir + "*.txt" )
	if len( input_files ) == 0:
		sys.exit( "ERROR: no data files detected!" )
	
	for filename in input_files:
		ID = filename.split('/')[-1].split('.')[0]
		up, down = load_data( filename )
		fig_file = output_dir + ID + '.png'
		print "number of upregulated genes: " + str( len( up ) )
		print "number of downregulated genes: " + str( len( down ) )
		construct_plot( up[ :min( [ len( up ), gene_number ] ) ], down[ :min( [ len( down ), gene_number ] ) ], fig_file, ID )


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
		
