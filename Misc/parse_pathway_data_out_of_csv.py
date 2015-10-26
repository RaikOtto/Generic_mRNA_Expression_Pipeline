"""
Take the pathway file, read out all genes in all pathways, optionally limited for bigger amount of genes per pathway.
Read the called genes files
write out how many mutated genes were found in the pathway and write that information into a format, to be read by the r
script BH
"""
epilog="""Raik Otto <raikotto@bayer-ag.de> 26.03.2013"""

import argparse, sys

def get_genes( input_file ):
	"""
	Reads the 'integrated data file format, which is proprietary'
	"""
	g_d = {}
	pos_i = {}

	with open( input_file, "r" ) as i_h:
		for line in i_h:

			line 	= line.strip().split("\t")
			if line[0].startswith("#"): continue
			elif line[0].startswith(">"):
				mut_types = line
				for entry in mut_types:
					pos_i[entry] = line.index(entry)

			else:

				gene = line[pos_i['>gene']].upper()

				g_d[gene] = {}

				for entry in mut_types:

					if entry == ">gene": continue
					g_d[gene][entry] = int(line[pos_i[entry]])

	return g_d

def get_pathways( pathway_file, cut_off_pathway ):

	p_d	=	{}
	m_d	=	{} # m parameter dictionary to count what genes occur

	with open( pathway_file, "r" ) as p_h:
		for line in p_h:

			line = line.strip().split("\t")
			if line[0].startswith("pathway\tsource"): continue

			name 	= line[0]
			source 	= line[1]
			genes	=	line[-1].split(",")

			if len(genes) <= cut_off_pathway: continue

			for gene in genes:	m_d[gene] = True
			#	if not m_d.has_key(gene): m_d[gene] = True

			p_d[name] = {}
			p_d[name]['not_seen'] = line[-1].split(",")
			p_d[name]['seen'] = []
			p_d[name]['source'] = source

	return (p_d, len(m_d.keys()))

def analyze_pathways( gene_file, formated_res_file, pathway_file, cut_off_muts, cut_off_pathway):

	p_d,N = get_pathways( pathway_file, cut_off_pathway )
	g_d = get_genes( gene_file )

	gene_count = 0
	for gene in g_d.keys(): 
		if (g_d[gene]['novel_missense'] >= cut_off_muts):
			for pathway in p_d.keys():
				if (gene in p_d[pathway]['not_seen']): gene_count += 1; break

	header = "q\tm\tn\tk\tpathway\tgenes\tgene_source\r\n"
	with open( formated_res_file, "w+" ) as o_h: o_h.write( header )

	for pathway in p_d.keys():

		for gene in list(p_d[pathway]['not_seen']):

			gene = gene.upper()

			if g_d.has_key(gene):

				if (int(g_d[gene]['novel_missense']) < cut_off_muts): continue

				if gene in p_d[pathway]['not_seen']:

					p_d[pathway]['seen'].append(gene)

					if len(p_d[pathway]['not_seen']) > 1:	p_d[pathway]['not_seen'].remove(gene)
					else:									p_d[pathway]['not_seen'] = []

		q = str(len(p_d[pathway]['seen']))
		m = len(p_d[pathway]['seen']) + len(p_d[pathway]['not_seen'])
		n = str(N - m); m = str(m)
		k = str(gene_count+1)
		genes = ",".join(p_d[pathway]['seen'])

		with open( formated_res_file, "a" ) as o_h:
			if (len(p_d[pathway]['seen'])>0):
				o_h.write("\t".join([q,m,n,k,pathway,genes,p_d[pathway]['source']])+"\r\n")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=__doc__, epilog=epilog)
	parser.add_argument('-i',  '--input_file' , nargs=1, 	 	help='input files', required=True)
	parser.add_argument('-o',  '--output_file' , nargs=1, 	 	help='file to contain the counts per pathway', required=True)
	parser.add_argument("-p", "--pathway_file", nargs=1, help="file that contains the pathway information", required=True)
	parser.add_argument("-cm", "--cut_off_muts", nargs=1, help="just consider genes with more than or equal to x mutations", required=False, default = 1, type=int)
	parser.add_argument("-cp", "--cut_off_pathway", nargs=1, help="just consider pahtways with more than x members", required=False, default = 0, type=int)

	parser = parser.parse_args()
	analyze_pathways( parser.input_file[0], parser.output_file[0], parser.pathway_file[0], parser.cut_off_muts, parser.cut_off_pathway)
