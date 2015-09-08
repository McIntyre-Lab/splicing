#!/usr/bin/env python

#	DESCRIPTION: This program uses the transcript FASTA file and extracts the junctions made by each pair of exons. 
#	The output is a bed file with all of these junctions. 
#
#	AUTHOR: Chelsea Tymms

# Built-in packages
import argparse
import os

# Add-on packages
import gffutils
import itertools


def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--gff", dest="gffInput", action='store', required=True, help="Input GFF or GTF file")
    parser.add_argument("--output", dest="outCSV", action='store', required=True, help="Output CSV file name")
    
    args = parser.parse_args()
    return args

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

def main():
    # Get the database
    db_fn=args.gffInput + '.db'
    db=gffutils.FeatureDB(db_fn)
    
    rnas=db.features_of_type(['mRNA','ncRNA'])
    
    csvOut = open(args.outCSV, 'w')
    
    # Write column headers for TSV
    csvHeader = ','.join(['chrom','start','stop','strand','exon_id','transcript_id','gene_id']) + "\n"
    csvOut.write(csvHeader)


    # Open the output BED file
    with open('aceview_hg19_transcript_junctions.tsv','wb') as outputFile:
        for rna in rnas:

            # Get all exons
            exonList = list(db.children(rna,featuretype='exon',order_by=('start','end')))
            
            # Skip genes that only have 1 exon, or genes that don't have an exon e.g. miRNA
            if len(exonList) == 0 or len(exonList) == 1:
                continue

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

                # Construct various parts of the junction BED file
                junctionID=''.join(junctionArray[i][0].attributes['Name']) + '|' + ''.join(junctionArray[i][1].attributes['Name']) 
                xscriptGenes=list(db.parents(rna, featuretype='gene'))
                xscriptGene=[]
                for gene in xscriptGenes:
                    xscriptGene.append(gene.id)
                transcriptID=''.join(rna.id)
                geneID=''.join(xscriptGene)
                juncCoord=str(junctionArray[i][0].chrom) + ':' + str(junctionArray[i][0].stop) + ':' + str(junctionArray[i][1].start) + ':' + str(junctionArray[i][0].strand)
                csvArray=[junctionID, juncCoord, transcriptID, geneID]
                # Output to BED file
                csvOut.write(",".join(str(i) for i in csvArray)+"\n")
    
    csvOut.close()
                        
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()


