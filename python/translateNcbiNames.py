import sys
import gzip
import csv

input_gxf=sys.argv[1]
output_gxf=sys.argv[2]
ncbi_to_ens_names=sys.argv[3]


#The NCBI GTF uses names based on RefSeq, so translate:
ncbiIdMap = dict()
with open(ncbi_to_ens_names, 'r') as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		ncbiIdMap[row[0]] = row[1]
		
		
noAlias = set()

def openPossiblyGz(fname): 
	return gzip.open(fname, mode = 'rt') if fname.endswith('.gz') else open(fname, 'r')

with openPossiblyGz(input_gxf) as input:
	with open(output_gxf, 'w') as output:
		for row in input:
			if not row.startswith('#'):
				#Translate NCBI ID into chromosome #, since NCBI's GTF is named based on accession:
				row = row.split('\t')
				if row[0] in ncbiIdMap.keys():
					row[0] = ncbiIdMap[row[0]]
				else:
					noAlias.add(row[0])
					
				row = '\t'.join(row)
				
			output.write(row)
			

if len(noAlias) > 0:
	print('NCBI contigs without accession -> name mapping (this might be ok): ' + ','.join(noAlias))
else:
	print('All NCBI contig names were aliased to Ensembl names')