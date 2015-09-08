#!/usr/bin/env python

# Built-in packages
import argparse
import os

# Add-on packages
import gffutils


def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--gff", dest="gffInput", action='store', required=True, help="Input GFF or GTF file")
    parser.add_argument("--obed", dest="outBED", action='store', required=True, help="Output BED file name")
    parser.add_argument("--otable", dest="outCSV", action='store', required=True, help="Output CSV file name")
    
    args = parser.parse_args()
    return args

def main():
    # Get the database
    db_fn=args.gffInput + '.db'
    db=gffutils.FeatureDB(db_fn)
    
    # Open output files for writing
    bedOut = open(args.outBED, 'w')
    csvOut = open(args.outCSV, 'w')
    
    # Write column headers for TSV
    csvHeader = ','.join(['chrom','start','stop','strand','exon_id','transcript_id','gene_id']) + "\n"
    csvOut.write(csvHeader)

    exons = db.features_of_type('exon')
    for exon in exons:
        chrom=exon.chrom
        start=exon.start-1
        stop=exon.stop
        strand=exon.strand
        name=''.join(exon.attributes['Name'])
        #aceview doesn't assign a unique id for each exon, so we can make one based on it's coordinates and gene_id
        #exon_coord=''.join(exon.attributes['gene_id']) + '|' + str(chrom) + ':' + str(start) + '-' + str(stop) + ':' + exon.strand
        xscripts=list(db.parents(exon, featuretype='ncRNA'))
        xscript_cat=[]
        for xscript in xscripts:
            xscript_cat.append(xscript.id)
        genes=list(db.parents(exon, featuretype='gene'))
        gene_cat=[]
        for gene in genes:
            gene_cat.append(gene.id)
        score='.'
        xscriptList='|'.join(xscript_cat)
        geneList='|'.join(gene_cat)
        csvEntry=[chrom, start, stop, strand, name, xscriptList, geneList]
        bedEntry=[chrom, start, stop, name, score, strand]
        csvOut.write(",".join(str(x) for x in csvEntry) + "\n")
        bedOut.write("\t".join(str(x) for x in bedEntry) + "\n")
    
    bedOut.close()
    csvOut.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()


