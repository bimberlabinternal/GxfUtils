import gffutils
import sys
import os
from natsort import natsort_keygen

input_gxf=sys.argv[1]
output_gxf=sys.argv[2]

transcriptKey='transcript_id'
geneKey='gene_id'
sortFeatureType='gene'
dbOut='./tempDb'

if os.path.exists(dbOut):
	os.remove(dbOut)

print('Reading Input: ' + input_gxf)
db = gffutils.create_db(input_gxf, dbOut, 
	gtf_transcript_key=transcriptKey, 
	gtf_gene_key=geneKey, 
	verbose=False, 
	disable_infer_transcripts=True,
	disable_infer_genes=True
)

nsKey = natsort_keygen()

def sortFeatures(feat):
	sortPriority = 3
	if feat.featuretype.lower() == 'gene':
		sortPriority = 1
	elif feat.featuretype.lower() == 'transcript' or feat.featuretype.lower() == 'mrna':
		sortPriority = 2
		
	return ( nsKey(feat.seqid), feat.start, sortPriority, feat.strand )

print('Writing Output')
i = 0
with open(output_gxf, 'w') as output:
	for feat in sorted(db.all_features(), key = sortFeatures):
		output.write(str(feat) + '\n')
		i += 1
		if i % 10000 == 0:
			print(str(i) + ' lines written')

if os.path.exists(dbOut):
	os.remove(dbOut)
	
print('Done')