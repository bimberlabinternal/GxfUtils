import gffutils
import sys

input_gxf=sys.argv[0]
output_gxf=sys.argv[1]

transcriptKey='transcript_id'
geneKey='gene_id'
sortFeatureType='gene'


db = gffutils.create_db(input_gxf, ":memory:", 
	gtf_transcript_key=transcriptKey, 
	gtf_gene_key=geneKey, 
	verbose=False, 
	disable_infer_transcripts=True,
	disable_infer_genes=True
)

def sortFeatureGroup(feats):
	return ( feats[0].seqid, feats[0].start, feats[0].strand )

with open(output_gxf, 'w') as output:
	for feats in sorted(db.iter_by_parent_childs(featuretype=sortFeatureType, order_by = ( 'seqid', 'start', 'strand')), key = sortFeatureGroup):
		parent = feats[0]
		children = feats[1:]
		
		for c in children:
			output.write(str(c) + '\n')
