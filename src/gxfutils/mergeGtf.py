import gffutils
from gffutils import Feature
import sys
import os
import csv
import copy

ncbi_gtf=sys.argv[1]
ensembl_gtf=sys.argv[2]
gene_intersect_bedtools=sys.argv[3]
ncbi_to_ens_mapping=sys.argv[4]
outputBase=sys.argv[5]
outDir=sys.argv[6] or './'

gtfInMemory = False

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

def getSortPriority(d):
	if d.featuretype == 'gene':
		return "1"
	elif d.featuretype == 'transcript' or d.featuretype == 'mRNA':
		return "2"
	else:
		return "3"
		
def transformNCBI(d):
	# Unclear why this happens, but the parser picks up an emtpy attribute some of the time:
	if '' in d.attributes.keys():
		del d.attributes['']

	#Store this for sorting later.
	rowCounter['ncbi'] = rowCounter['ncbi'] + 1
	d.attributes['sort_order'] = str(rowCounter['ncbi'])
	d.attributes['sort_priority'] = getSortPriority(d)

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
	d.attributes['sort_priority'] = getSortPriority(d)

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
ncbiDbFile = outDir + 'ncbiDB'
if not gtfInMemory and os.path.exists(ncbiDbFile):
	os.remove(ncbiDbFile)
	
ncbi = gffutils.create_db(ncbi_gtf, (":memory:" if gtfInMemory else ncbiDbFile), 
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
ensDbFile = outDir + 'ensDB'
if not gtfInMemory and os.path.exists(ensDbFile):
	os.remove(outDir + 'ensDB')
	
ensdb = gffutils.create_db(ensembl_gtf, (":memory:" if gtfInMemory else ensDbFile),
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

def sortFeatureGroup(record):
	return ( int(record.attributes['sort_order'][0]), int(record.attributes['sort_priority'][0]) )
	
#Iterate existing Ensembl genes.  If a gene is mapped to NCBI, iterate transcripts on both sides and merge.
print('Comparing/merging transcripts')
transcriptsToAdd = []
transcriptNamesToAdd = set()
transcriptsMappingUsingNcbi = 0
transcriptWithIdenticalExons = 0

i = 0
with open(outDir + 'NCBITranscriptsMergedToEnsemblGenes.txt', 'w') as mergeOut:
	mergeOut.write('\t'.join(['EnsGeneId', 'EnsGeneName', 'NCBI_GeneId', 'TotalEnsTranscripts', 'TotalNCBITranscript', 'TotalMerged', 'MergedNames']) + '\n')
		
	for records in ensdb.iter_by_parent_childs(featuretype='gene'):
		i += 1
		if i % 2500 == 0:
			print('processed: ' + str(i) + ' records')

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
			cognateChildren = sorted(cognateChildren, key = sortFeatureGroup)
			
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
	verbose=False,
	make_backup=False
)

transcriptFeaturesAdded = []
uniqueTranscriptIds = set()
transcriptNameUpdates = set()

def ensureUniqueTranscriptName(transcriptName, uniqueTranscriptIds):
	if transcriptName not in uniqueTranscriptIds:
		return transcriptName
	
	transcriptNameOrig = transcriptName
	i = 0
	while transcriptName in uniqueTranscriptIds:
		i += 1
		transcriptName = transcriptNameOrig + '.' + str(i)
		
	return transcriptName

def possiblyAddTranscriptRecord(parent, children):
	transcriptFeaturesFound = {}
	for f in children:
		if 'transcript_id' in f.attributes.keys() and f.attributes['transcript_id'][0] not in transcriptFeaturesFound.keys():
			transcriptFeaturesFound[f.attributes['transcript_id'][0]] = False
		if f.featuretype == 'transcript':
			transcriptFeaturesFound[f.attributes['transcript_id'][0]] = True
			uniqueTranscriptIds.add(f.attributes['transcript_id'][0])
	
	transcriptUpdatesNeeded = {}
	for transcriptName in transcriptFeaturesFound.keys():
		if transcriptFeaturesFound[transcriptName] == False:
			minStart = 0
			maxEnd = 0
			for f in children:
				if 'transcript_id' in f.attributes.keys() and f.attributes['transcript_id'][0] == transcriptName:
					if minStart == 0 or minStart > f.start:
						minStart = f.start
						
					if maxEnd < f.end:
						maxEnd = f.end
			
			props = {}
			props['seqid'] = parent.seqid
			props['start'] = minStart
			props['end'] = maxEnd
			props['strand'] = parent.strand
			props['source'] = parent.source
			props['featuretype'] = 'transcript'
			props['attributes'] = copy.deepcopy(parent.attributes)
			feat = Feature(dialect = parent.dialect, **props)
			
			feat.attributes['sort_priority'] = "2"
			transcriptNameToUse = transcriptName
			if 'gene_id' in feat.attributes.keys() and transcriptName == feat.attributes['gene_id'][0]:
				transcriptNameToUse = 'transcript-' + transcriptName

			transcriptNameToUse = ensureUniqueTranscriptName(transcriptNameToUse, uniqueTranscriptIds)
			if transcriptNameToUse != transcriptName:
				transcriptNameUpdates.add( (feat.attributes['gene_id'][0], transcriptName, transcriptNameToUse) )
				transcriptUpdatesNeeded[transcriptName] = transcriptNameToUse
				
			uniqueTranscriptIds.add(transcriptNameToUse)
			
			feat.attributes['transcript_id'] = transcriptNameToUse
			feat.attributes['new_transcript'] = "Y"	
			children.append(feat)
			transcriptFeaturesAdded.append(feat)
	
	if len(transcriptUpdatesNeeded) > 0:
		for c in children:
			if 'transcript_id' in c.attributes.keys() and c.attributes['transcript_id'][0] in transcriptUpdatesNeeded.keys():			
				c.attributes['transcript_id'] = transcriptUpdatesNeeded[c.attributes['transcript_id'][0]]
	
	children = sorted(children, key = sortFeatureGroup)	
	
	return children

	
print('Total NCBI transcripts added to an existing Ensembl gene: '  + str(len(transcriptNamesToAdd)))
results['NCBITranscriptsAddedToEnsemblGene'] = len(transcriptNamesToAdd)

# These genes from NCBI have no mate in Ensembl, and no overlapping gene.  Add them outright:
print('Adding ' + str(len(results['TotalNCBINotMappedToEns'])) + ' new genes that lack mapped or overlapping feature')
toAdd = []
for geneId in results['TotalNCBINotMappedToEns']:
	children = ncbi.children(ncbi[geneId])
	children = sorted(children, key = sortFeatureGroup)

	parent = None
	for c in children:
		if c.featuretype == 'gene':
			parent = c
			break
			
	if parent == None:
		raise 'No gene record found for group: ' + geneId

	# Create a transcript/mRNA record if necessary:
	children = possiblyAddTranscriptRecord(parent, children)
	
	for c in children:
		toAdd.append(c)

with open(outDir + 'FeaturesAdded.gtf', 'w') as fout:
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
	verbose=False,
	make_backup=False
)

def sortFeatureGroupByPosition(feats):	
	return ( feats[0].seqid, feats[0].start, feats[0].strand )
	
print('Writing final GTF/GFFs')
geneNameMismatch = []
transcriptsLackingId = []
uniqueIds = set()
duplicateIds = []
with open(mergedGtfOut, 'w') as gtfOut, open(mergedGffOut, 'w') as gffOut:
	gtfOut.write('##gtf-version 2\n')
	gtfOut.write('#Ensembl Input: ' + ensembl_gtf + '\n')
	gtfOut.write('#NCBI Input: ' + ncbi_gtf + '\n')
	
	gffOut.write('##gff-version 3\n')
	gffOut.write('#Ensembl Input: ' + ensembl_gtf + '\n')
	gffOut.write('#NCBI Input: ' + ncbi_gtf + '\n')

	i = 0
	for feats in sorted(ensdb.iter_by_parent_childs(featuretype='gene', order_by = ( 'seqid', 'start', 'strand')), key = sortFeatureGroupByPosition):
		parent = feats[0]
		children = feats[1:]

		# A GFF file will require a parent 'transcript' feature for each. Prepare this here:
		children = possiblyAddTranscriptRecord(parent, children)
		children = sorted(children, key = sortFeatureGroup)
				
		for c in children:
			del c.attributes['sort_order']
			del c.attributes['sort_priority']
			
			if 'overlapping_gene_id' in c.attributes.keys():
				del c.attributes['overlapping_gene_id']
				
			if 'overlapping_gene_dbid' in c.attributes.keys():
				del c.attributes['overlapping_gene_dbid']
			
			if 'db_xref' in c.attributes.keys():
				del c.attributes['db_xref']

			if 'ncbi_genename' in c.attributes.keys() and 'gene_name' in c.attributes.keys() and c.attributes['ncbi_genename'][0] != c.attributes['gene_name'][0]:
				geneNameMismatch.append([c.attributes['gene_id'][0], c.attributes['gene_name'][0], c.attributes['ncbi_genename'][0]])
			
			if 'ncbi_genename' in c.attributes.keys() and 'gene_name' not in c.attributes.keys():
				c.attributes['gene_name'] = c.attributes['ncbi_geneid']
				geneNameMismatch.append([c.attributes['gene_id'][0], 'MISSING', c.attributes['ncbi_genename'][0]])
			
			gtfOut.write(str(c) + '\n')
			i += 1
			if i % 10000 == 0:
				print(str(i) + ' lines written')
			
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
			
			#Specifically for GFF, create ID and Parent:
			if 'transcript' == c.featuretype:
				# GFFs also seem to expect mRNA as feature type:
				c.featuretype = 'mRNA'
				if 'transcript_id' not in c.attributes.keys() or c.attributes['transcript_id'][0] == '':
					transcriptsLackingId.append(c)
				else:
					featId = c.attributes['transcript_id'][0]
					if featId in uniqueIds:
						duplicateIds.append(c)
						
					uniqueIds.add(featId)
					
					c.attributes['ID'] = 'transcript:' + featId
					c.attributes['Parent'] = 'gene:' + c.attributes['gene_id'][0]
			elif 'gene' == c.featuretype:
				c.attributes['ID'] = 'gene:' + c.attributes['gene_id'][0]
			elif 'gene' != c.featuretype:
				if 'transcript_id' not in c.attributes.keys() or c.attributes['transcript_id'] == '':					
					transcriptsLackingId.append(c)
				else:
					c.attributes['Parent'] = 'transcript:' + c.attributes['transcript_id'][0]

			gffOut.write(str(c) + '\n')

print('Transcript features added: ' + str(len(transcriptFeaturesAdded)))
print('Transcript IDs update to make unique: ' + str(len(transcriptNameUpdates)))
results['TranscriptFeaturesAdded'] = len(transcriptFeaturesAdded)

with open(outDir + 'TranscriptsAdded.gtf', 'w') as fout:	
	for f in transcriptFeaturesAdded:
		fout.write(str(f) + '\n')

results['TranscriptNamesUpdates'] = len(transcriptNameUpdates)
if len(transcriptNameUpdates) > 0:
	with open(outDir + 'TranscriptNamesUpdates.gtf', 'w') as fout:
		fout.write('GeneId\tOrigTranscriptName\tUniquifiedTranscriptName\n')
		for c in transcriptNameUpdates:
			fout.write(c[0] + '\t' + c[1] + '\t' + c[2] + '\n')
			
results['TranscriptsLackingId'] = len(transcriptsLackingId)
if len(transcriptsLackingId) > 0:
	with open(outDir + 'TranscriptsLackingTranscriptId.gtf', 'w') as transcriptOut:
		for c in transcriptsLackingId:
			transcriptOut.write(str(c) + '\n')

if len(duplicateIds) > 0:
	with open(outDir + 'DuplicateIds.gtf', 'w') as fout:
		for c in duplicateIds:
			fout.write(str(c) + '\n')
			
if len(geneNameMismatch) > 0:
	with open(outDir + 'GeneNameMismatch.txt', 'w') as fout:
		fout.write('Ensembl_GeneName' + '\t' + 'NCBI_GeneName' + '\n')
		for l in geneNameMismatch:
			fout.write('\t'.join(l) + '\n')

with open(outDir + 'Summary.txt', 'w') as output:
	for key in results.keys():
		if hasattr(results[key], '__len__'):
			output.write(key + '\t' + str(len(results[key])) + '\n')
		else:
			output.write(key + '\t' + str(results[key]) + '\n')

del ncbi
if not gtfInMemory and os.path.exists(ncbiDbFile):
	os.remove(ncbiDbFile)

del ensdb
if not gtfInMemory and os.path.exists(ensDbFile):
	os.remove(ensDbFile)
	
print('Done')