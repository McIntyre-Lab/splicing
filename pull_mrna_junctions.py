#!/usr/bin/env python

#	DESCRIPTION: This program uses the transcript FASTA file and extracts the junctions made by each pair of exons. 
#	The output is a bed file with all of these junctions. 
#
#	AUTHOR: Chelsea Tymms


import gffutils
import itertools


def start_pos(x):
    """ Little function to let me sort gffutils objects by start position """
    return x.start

def createExonArray(exonList,mergedExons):
    """ Group overlapping exons together as an list within a list """
    exonArray = []

    for mergedExon in mergedExons:
        exonSubset = []

        for exon in exonList:
            if exon.start >= mergedExon.start and exon.stop <= mergedExon.stop:
                exonSubset.append(exon)

        exonArray.append(exonSubset)

    return exonArray

def createJunctionArray(exonArray):
    """ Using itertools create all possible combination between groups of exons """
    iterProduct = []
    s=len(exonArray)
    for i in range(0,s-1):
#        for j in range (i+1,len(exonArray)):
        j=i+1
        for k in itertools.product(exonArray[i],exonArray[j]):
            iterProduct.append(k)

    return iterProduct

def defineRegion(index,junctionArray):
    """ Instead of keeping the entire exon, only keep up to 75bp on either side of the junction """
    if len(junctionArray[index][0]) >= 38:
        exon1Start = junctionArray[index][0].stop - 37
        exon1Stop = junctionArray[index][0].stop
    else:
        exon1Start = junctionArray[index][0].start
        exon1Stop = junctionArray[index][0].stop

    if len(junctionArray[index][1]) >= 38:
        exon2Start = junctionArray[index][1].start 
        exon2Stop = junctionArray[index][1].start + 37
    else:
        exon2Start = junctionArray[index][1].start 
        exon2Stop = junctionArray[index][1].stop

    return (exon1Start,exon1Stop,exon2Start,exon2Stop)


def main():
    # Get the database
    #fname='../gff/AceView.GRCh37.hg19.gtf'
    #db_fn='AceView.GRCh37.hg19.gtf.db'
    #db_fn='AceView.spliced.GRCh37.hg19.gtf.db'
    db_fn='../../useful_dmel_data/flybase557/gff/dmel-all-no-analysis-r5.57.gff.db'
    # gffutils.create_db(fname, db_fn, merge_strategy='merge', id_spec=["Name", "ID"])
    db=gffutils.FeatureDB(db_fn)
    rnas=db.features_of_type('mRNA')

    # Open the output BED file
    with open('../generated_files/mRNA_annotated_junctions.tsv','wb') as outputFile:
        for rna in rnas:

            # Get all exons
            exonList = list(db.children(rna,featuretype='exon'))
            
            # Skip genes that only have 1 exon, or genes that don't have an exon e.g. miRNA
            if len(exonList) == 0 or len(exonList) == 1:
                continue

            # Convert to BED 0-based coordinate system. Only have to change start.
            for exon in exonList:
                exon.start = exon.start - 1 

            # Create list of overlapping exons
            # Note: mod(mdg4) has exons on both strands, so ignore_strand needs
            # to be set for it
            #mergedExons = list(db.merge_features(exonList,ignore_strand=True))
            mergedExons = list(db.merge(exonList,ignore_strand=True))

            # Sort merged exons by start postiion
            mergedExons.sort(key=start_pos)

            # Using the merged exons to identify overlapping exons create an array where
            # overlapping exons are grouped together
            exonArray = createExonArray(exonList,mergedExons)

            # Now create all possible junctions 
            junctionArray = createJunctionArray(exonArray)

            # Create BED file
            for i in range(0,len(junctionArray)):

                # Only keep up to 75bp on each side of the junction
                (exon1Start, exon1Stop, exon2Start, exon2Stop) = defineRegion(i,junctionArray)

                # Construct various parts of the junction BED file
                totalStart=exon1Start # Take exon1 start
                totalStop=exon2Stop  # Take exon2 stop
                chrom=junctionArray[i][0].chrom
                strand=junctionArray[i][0].strand
                name="".join(junctionArray[i][0].attributes['Name']) + '|' + str(junctionArray[i][0].chrom) + ':' + str(junctionArray[i][0].start) + '-' + str(junctionArray[i][0].stop) + ':' + junctionArray[i][0].strand + '|' + "".join(junctionArray[i][1].attributes['Name']) + str(junctionArray[i][1].chrom) + ':' + str(junctionArray[i][1].start) + '-' + str(junctionArray[i][1].stop) + ':' + junctionArray[i][1].strand
                cat_gene="|".join([gene.id for gene in db.parents(junctionArray[i][0], featuretype='gene', level=2)])
                cat_xscript="|".join(junctionArray[i][0].attributes['Parent'])
                xscript_gene="|".join([gene.id for gene in db.parents(rna, featuretype='gene')])
                xscripter="|".join(rna.attributes['ID'])
                exact_jnct=str(junctionArray[i][0].chrom) + ':' + str(junctionArray[i][0].stop) + ':' + str(junctionArray[i][1].start) + ':' + str(junctionArray[i][0].strand)
                bedArray=[chrom,totalStart,totalStop,name,strand,totalStart,totalStop,cat_xscript,cat_gene,exact_jnct,xscripter,xscript_gene]
                # Output to BED file
                outputFile.write("\t".join(str(i) for i in bedArray)+"\n")
                        
if __name__ == '__main__':
    main()

