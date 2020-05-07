import gffutils
import sys
import os
import csv

ncbi_gtf=sys.argv[1]
ensembl_gtf=sys.argv[2]
gene_intersect_bedtools=sys.argv[3]
ncbi_to_ens_mapping=sys.argv[4]
outputBase=sys.argv[5]
outDir=sys.argv[6] or './'

if outDir != './' and not os.path.exists(outDir):
	os.makedirs(outDir)

mergedGtfOut = outDir + outputBase + '.gtf'
mergedGffOut = outDir + outputBase + '.gff'

#build map of NCBI <-> Ensembl by transcript
ncbiGeneIdToEnsemblIdMap = dict()
ensemblGeneToNcbiIdMap = dict()
ncbiTranscriptIdtoEnsemblMap = dict()
ensemblTranscriptIdtoNcbiMap = dict()
with open(ncbi_to_ens_mapping, 'r') as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		ncbiGeneIdToEnsemblIdMap[row[1]] = row[2]
		ensemblGeneToNcbiIdMap[row[2]] = row[1]

		#Drop version (i.e. .1) from name:
		ensT = row[4].split('.')[0]

		ncbiTranscriptIdtoEnsemblMap[row[3]] = ensT
		ensemblTranscriptIdtoNcbiMap[ensT] = row[3]

print('Total NCBI->Ensembl mappings from NCBI: ' + str(len(ncbiGeneIdToEnsemblIdMap)))

#Now based on reciprocal overlap (bedtools)
ncbiToEnsOverlapByGeneMap = dict()
ensToNcbiOverlapByGeneMapToGene = dict()
ensToNcbiOverlapByGeneMapToId = dict()
with open(gene_intersect_bedtools, 'r') as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		n = gffutils.feature.feature_from_line('\t'.join(row[0:9]) + '\n')
		e = gffutils.feature.feature_from_line('\t'.join(row[9:18])+ '\n')
		
		#NCBI Gene Name to EnsemblId
		ncbiToEnsOverlapByGeneMap[n.attributes['gene_id'][0]] = e.attributes['gene_id'][0]
		ensToNcbiOverlapByGeneMapToGene[e.attributes['gene_id'][0]] = n.attributes['gene_id'][0]
		
		#NCBI Numeric Gene Id to EnsemblId
		if 'db_xref' in n.attributes.keys():
			ncbiId = n.attributes['db_xref'][0].split(':')[1]
			
			ncbiToEnsOverlapByGeneMap[ncbiId] = e.attributes['gene_id'][0]
			ensToNcbiOverlapByGeneMapToId[e.attributes['gene_id'][0]] = ncbiId
		else:
			print('No db_xref')
			print(n)

results = {
	'TotalNCBIGenes': set(),
	'TotalNCBIMappedToEns': set(),
	'TotalNCBINotMappedToEns': set(),
	'TotalNCBINotMappedToEnsWithOverlapping': set(),
	'TotalEnsGenes': set(),
	'TotalEnsMappedToNCBI': set(),
	'TotalEnsNotMappedToNCBI': set(),
	'TotalEnsNotMappedToNCBIWithOverlapping': set()
}

rowCounter = {
	'ens': 0,
	'ncbi': 0
}

ncbiGeneIdToGeneName = {}

def transformNCBI(d):
	# Unclear why this happens, but the parser picks up an emtpy attribute some of the time:
	if '' in d.attributes.keys():
		del d.attributes['']

	#Store this for sorting later.
	rowCounter['ncbi'] = rowCounter['ncbi'] + 1
	d.attributes['sort_order'] = str(rowCounter['ncbi'])

	#More explicitly store this ID, since it is used heavily within the GTF
	ncbiId = "-1"
	if 'db_xref' in d.attributes.keys():
		ncbiId = d.attributes['db_xref'][0].split(':')[1]
		d['ncbi_gene_dbid'] = ncbiId
		ncbiGeneIdToGeneName[ncbiId] = d['gene_id'][0]
	
	d.attributes['source_gtf'] = 'ncbi'
	d.source = 'ncbi'
	
	#Normalize attribute name with Ensembl:
	if 'gene' in d.attributes.keys():
		d.attributes['gene_name'] = d.attributes['gene']
		del d.attributes['gene']

	#Only both cross-annotating the parent records:	
	if d.featuretype == 'gene' or d.featuretype == 'transcript':
		#gene_id is generally the gene name in NCBI
		if d['gene_id'][0] in ncbiToEnsOverlapByGeneMap.keys():
			d.attributes['overlapping_gene_id'] = ncbiToEnsOverlapByGeneMap[d['gene_id'][0]]
		
		if ncbiId != "-1" and ncbiId in ncbiToEnsOverlapByGeneMap.keys():
			d.attributes['overlapping_gene_dbid'] = ncbiToEnsOverlapByGeneMap[ncbiId]
		
		if ncbiId != "-1" and ncbiId in ncbiGeneIdToEnsemblIdMap.keys():
			d.attributes['ensembl_geneid'] = ncbiGeneIdToEnsemblIdMap[ncbiId]

		results['TotalNCBIGenes'].add(d['gene_id'][0])
		if 'ensembl_geneid' in d.attributes.keys():
			results['TotalNCBIMappedToEns'].add(d['gene_id'][0])
		elif 'overlapping_gene_id' in d.attributes.keys():
			results['TotalNCBINotMappedToEnsWithOverlapping'].add(d['gene_id'][0])
		else:
			results['TotalNCBINotMappedToEns'].add(d['gene_id'][0])
	else:
		if 'transcript_id' in d.attributes.keys():
			d['ncbi_transcript_id'] = d['transcript_id']
			if d['transcript_id'][0] in ncbiTranscriptIdtoEnsemblMap.keys():
				d.attributes['ensembl_transcript_id'] = ncbiTranscriptIdtoEnsemblMap[d['transcript_id'][0]]
		elif d.featuretype != 'gene':
			print('No transcript ID')

	return d

def transformEns(d):
	d.attributes['source_gtf'] = 'ensembl'
	d.source = 'ensembl'
	
	#Store this for sorting later.
	rowCounter['ens'] = rowCounter['ens'] + 1
	d.attributes['sort_order'] = str(rowCounter['ens'])

	#Only both cross-annotating the parent records:
	if d.featuretype == 'gene' or d.featuretype == 'transcript':
		if d['gene_id'][0] in ensemblGeneToNcbiIdMap.keys():
			d.attributes['ncbi_gene_dbid'] = ensemblGeneToNcbiIdMap[d['gene_id'][0]]
			d.attributes['ncbi_geneid'] = ncbiGeneIdToGeneName[ensemblGeneToNcbiIdMap[d['gene_id'][0]]]

		if d['gene_id'][0] in ensToNcbiOverlapByGeneMapToGene.keys():
			d.attributes['overlapping_gene_id'] = ensToNcbiOverlapByGeneMapToGene[d['gene_id'][0]]

		if d['gene_id'][0] in ensToNcbiOverlapByGeneMapToId.keys():
			d.attributes['overlapping_gene_dbid'] = ensToNcbiOverlapByGeneMapToId[d['gene_id'][0]]
			
		results['TotalEnsGenes'].add(d['gene_id'][0])
		if 'ncbi_gene_dbid' in d.attributes.keys():
			results['TotalEnsMappedToNCBI'].add(d['gene_id'][0])
		elif 'overlapping_gene_id' in d.attributes.keys():
			results['TotalEnsNotMappedToNCBIWithOverlapping'].add(d['gene_id'][0])
		else:
			results['TotalEnsNotMappedToNCBI'].add(d['gene_id'][0])
	else:
		if 'transcript_id' in d.attributes.keys():
			d.attributes['ensembl_transcript_id'] = d['transcript_id']
			if d['transcript_id'][0] in ensemblTranscriptIdtoNcbiMap.keys():
				d.attributes['ncbi_transcript_id'] = ensemblTranscriptIdtoNcbiMap[d['transcript_id'][0]]
		elif d.featuretype != 'gene':
			print('No transcript ID')

	return d

#Note: NCBI must be imported first, to build ID->Name map
print('Parsing NCBI GTF: ' + ncbi_gtf)
ncbi = gffutils.create_db(ncbi_gtf, ":memory:", 
	id_spec={'transcript': 'transcript_id', 'gene': 'gene_id', 'Name': 'gene_name'}, 
	gtf_transcript_key='transcript_id', 
	gtf_gene_key='gene_id', 
	verbose=False, 
	transform=transformNCBI, 
	disable_infer_transcripts=True,
	disable_infer_genes=True
)
print('Done')

print('Parsing Ensembl GTF: ' + ensembl_gtf)
ensdb = gffutils.create_db(ensembl_gtf, ":memory:",
	id_spec={'transcript': 'transcript_id', 'gene': 'gene_id', 'Name': 'gene_name'},
	gtf_transcript_key='transcript_id',
	gtf_gene_key='gene_id',
	verbose=False,
	transform=transformEns,
	disable_infer_transcripts=True,
	disable_infer_genes=True
)
print('Done')

#with open('ncbi.annotated.gtf', 'w') as fout:
#	for f in ncbi.all_features():
#		fout.write(str(f) + '\n')

#with open('ensembl.annotated.gtf', 'w') as fout:
#	for f in ensdb.all_features():
#		fout.write(str(f) + '\n')

with open(outDir + 'GenesWithOverlapNotJoined.txt', 'w') as out: 
	out.write('NCBI_GeneId\tEns_GeneId\tEns_GeneName\n')
	for name in results['TotalNCBINotMappedToEnsWithOverlapping']:
		geneName = ''
		f = ensdb[ncbiToEnsOverlapByGeneMap[name]]
		if 'gene_name' in f.attributes.keys():
			geneName = f.attributes['gene_name'][0]

		out.write(name + '\t' + ncbiToEnsOverlapByGeneMap[name] + '\t' + geneName + '\n')

	for name in results['TotalEnsNotMappedToNCBIWithOverlapping']:
		geneName = ''
		f = ensdb[name]
		if 'gene_name' in f.attributes.keys():
			geneName = f.attributes['gene_name'][0]
		
		out.write(ensToNcbiOverlapByGeneMapToGene[name] + '\t' + name + '\t' + geneName + '\n')

print('Summary:')
for key in results.keys():
	print(key + ': ' + str(len(results[key])))

def sortByOrder(record):
	return int(record.attributes['sort_order'][0])
	
#Iterate existing Ensembl genes.  If a gene is mapped to NCBI, iterate transcripts on both sides and merge.
transcriptsToAdd = []
transcriptNamesToAdd = set()
transcriptsMappingUsingNcbi = 0
transcriptWithIdenticalExons = 0
with open(outDir + 'EnsemblGenesMergedWithNcbi.txt', 'w') as mergeOut:
	mergeOut.write('\t'.join(['EnsGeneId', 'EnsGeneName', 'NCBI_GeneId', 'TotalEnsTranscripts', 'TotalNCBITranscript', 'TotalMerged', 'MergedNames']) + '\n')
	
	for records in ensdb.iter_by_parent_childs(featuretype='gene'):
		parent = records[0]
		children = records[1:]
		
		#Build map of transcripts and exon patterns:
		transcripts = {}
		parentTranscriptIds = set()
		for record in children:
			if 'transcript_id' in record.attributes.keys():
				parentTranscriptIds.add(record.attributes['transcript_id'][0])
				
			if record.featuretype == 'exon':
				if 'transcript_id' not in record.attributes.keys():
					raise 'Record lacks transcript_id: ' + str(record)
				
				tName = record.attributes['transcript_id'][0]
				if tName not in transcripts.keys():
					transcripts[tName] = []
					
				transcripts[tName].append(str(record.start) + '-' + str(record.end))

		parentExonPatterns = {}
		for key in transcripts.keys():
			parentExonPatterns[';'.join(transcripts[key])] = key

		#If we have a cognate NCBI record, retrieve it and compare exons/CDS:
		if 'ncbi_geneid' in parent.attributes.keys():
			name = parent.attributes['ncbi_geneid'][0]
			cognate = ncbi[name]
			cognateChildren = ncbi.children(name)
			cognateChildren = sorted(cognateChildren, key = sortByOrder)
			
			#Iterate the transcripts, append any with 
			childTranscripts = {}
			childTranscriptFeatures = {}
			for record in cognateChildren:
				if record.featuretype == 'exon':
					#indicates this matches
					if 'ensembl_transcript_id' in record.attributes.keys() and record.attributes['ensembl_transcript_id'][0] in parentTranscriptIds:
						transcriptsMappingUsingNcbi += 1
						continue
				
					if 'transcript_id' not in record.attributes.keys():
						raise 'Record lacks transcript_id: ' + str(record)
						
					tName = record.attributes['transcript_id'][0]
					if tName not in childTranscripts.keys():
						childTranscripts[tName] = []
						childTranscriptFeatures[tName] = []

					childTranscripts[tName].append(str(record.start) + '-' + str(record.end))
					childTranscriptFeatures[tName].append(record)
			
			addedForGene = set()
			for tName in childTranscripts.keys():
				joined = ';'.join(childTranscripts[tName])
				if joined not in parentExonPatterns.keys():
					transcriptNamesToAdd.add(tName)
					addedForGene.add(tName)
					
					#Actually add:
					feats = childTranscriptFeatures[tName]
					for f in feats:
						f.attributes['gene_id'] = parent.attributes['gene_id']
						transcriptsToAdd.append(f)
					
				else:
					transcriptWithIdenticalExons += 1
			
			parentGeneName = ''
			if 'gene_name' in parent.attributes.keys():
				parentGeneName = parent.attributes['gene_name'][0]
			else:
				parentGeneName = 'MISSING'
				
			mergeOut.write('\t'.join([parent.attributes['gene_id'][0], parentGeneName, cognate.attributes['gene_id'][0], str(len(transcripts)), str(len(childTranscripts)), str(len(addedForGene)), ','.join(addedForGene)]) + '\n')

print('Transcripts ignored because of NCBI->Ens mapping: ' + str(transcriptsMappingUsingNcbi))
print('Transcripts ignored because of itentical exons: ' + str(transcriptWithIdenticalExons))

ensdb.update(transcriptsToAdd, 
	id_spec={'transcript': 'transcript_id', 'gene': 'gene_id', 'Name': 'gene_name'}, 
	gtf_transcript_key='transcript_id', 
	gtf_gene_key='gene_id', 
	disable_infer_transcripts = True, 
	disable_infer_genes = True, 
	verbose=False
)

print('Total NCBI transcripts added to an existing Ensembl gene: '  + str(len(transcriptNamesToAdd)))

# These genes from NCBI have no mate in Ensembl, and no overlapping gene.  Add them outright:
print('Adding ' + str(len(results['TotalNCBINotMappedToEns'])) + ' new genes that lack mapped or overlapping feature')
toAdd = []
for geneId in results['TotalNCBINotMappedToEns']:
	children = ncbi.children(ncbi[geneId])
	children = sorted(children, key = sortByOrder)
	for c in children:
		toAdd.append(c)

with open(outDir + 'FeaturesToAdd.gtf', 'w') as fout:
	for f in toAdd:
		fout.write(str(f) + '\n')
	for f in transcriptsToAdd:
		fout.write(str(f) + '\n')
		
ensdb.update(toAdd, 
	id_spec={'transcript': 'transcript_id', 'gene': 'gene_id', 'Name': 'gene_name'}, 
	gtf_transcript_key='transcript_id', 
	gtf_gene_key='gene_id', 
	disable_infer_transcripts = True, 
	disable_infer_genes = True, 
	verbose=False
)

def sortFeatureGroup(feats):
	return ( feats[0].seqid, feats[0].start, feats[0].strand )

geneMismatch = []
with open(mergedGtfOut, 'w') as gtfOut, open(mergedGffOut, 'w') as gffOut:
	gtfOut.write('##gtf-version 2\n')
	gtfOut.write('#Ensembl Input: ' + ensembl_gtf + '\n')
	gtfOut.write('#NCBI Input: ' + ncbi_gtf + '\n')
	
	gffOut.write('##gff-version 3\n')
	gffOut.write('#Ensembl Input: ' + ensembl_gtf + '\n')
	gffOut.write('#NCBI Input: ' + ncbi_gtf + '\n')

	for feats in sorted(ensdb.iter_by_parent_childs(featuretype='gene', order_by = ( 'seqid', 'start', 'strand')), key = sortFeatureGroup):
		parent = feats[0]
		children = feats[1:]

		children = sorted(children, key = sortByOrder)

		for c in children:
			del c.attributes['sort_order']
			if 'overlapping_gene_id' in c.attributes.keys():
				del c.attributes['overlapping_gene_id']
				
			if 'overlapping_gene_dbid' in c.attributes.keys():
				del c.attributes['overlapping_gene_dbid']
			
			if 'db_xref' in c.attributes.keys():
				del c.attributes['db_xref']

			if 'ncbi_genename' in c.attributes.keys() and 'gene_name' in c.attributes.keys() and c.attributes['ncbi_genename'][0] != c.attributes['gene_name'][0]:
				geneMismatch.append([c.attributes['gene_id'][0], c.attributes['gene_name'][0], c.attributes['ncbi_genename'][0]])
			
			if 'ncbi_genename' in c.attributes.keys() and 'gene_name' not in c.attributes.keys():
				c.attributes['gene_name'] = c.attributes['ncbi_geneid']
				geneMismatch.append([c.attributes['gene_id'][0], 'MISSING', c.attributes['ncbi_genename'][0]])
				
			gtfOut.write(str(c) + '\n')
			
			c.dialect = {
				'field separator': ';',
				'fmt': 'gff3',
				'keyval separator': '=',
				'leading semicolon': False,
				'multival separator': ',',
				'quoted GFF2 values': False,
				'repeated keys': False,
				'trailing semicolon': False
			}
			
			gffOut.write(str(c) + '\n')

with open(outDir + 'GeneNameMismatch.txt', 'w') as fout:
	fout.write('Ensembl_GeneName' + '\t' + 'NCBI_GeneName' + '\n')
	for l in geneMismatch:
		fout.write('\t'.join(l) + '\n')
		
print('Done')