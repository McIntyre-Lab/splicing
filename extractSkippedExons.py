#!/usr/bin/env python

### Testing the concatentation of skipped exons
### Might make this later output each skipped exon to a separate list, but for now let's try this


# Built-in packages
import argparse

# Add-on packages
import gffutils

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--gff", dest="gffInput", action='store', required=True, help="Input GFF or GTF file")
    parser.add_argument("--otable", dest="outCSV", action='store', required=True, help="Output CSV file name")
    parser.add_argument("--olist", dest="outList", action='store', required=True, help="Output CSV file name")
    
    args = parser.parse_args()
    return args

def main():
    # Get the database
    db_fn=args.gffInput + '.db'
    db=gffutils.FeatureDB(db_fn)
    
    csvOut = open(args.outCSV, 'w')
    listOut = open(args.outList, 'w')
    
    # Write column headers for CSV
    csvHeader = ','.join(['junction_id', 'flag_exonskip', 'num_skipped_exons', 'cat_skipped_exons']) + "\n"
    csvOut.write(csvHeader)

    # Write column headers for List
    listHeader = ','.join(['junction_id', 'skipped_exon_id', 'flag_exonskip']) + "\n"
    listOut.write(listHeader)

    db_fn = args.gffInput + '.db'
    db = gffutils.FeatureDB(db_fn)

    # Get genes - we only want to concern ourselves with junctions within a gene
    genes = db.features_of_type('gene', order_by=('seqid','start','end'))

    # Not worrying about strand, since we only want the list of exons within a gene
    for gene in genes:
        exonList = list(db.children(gene, featuretype='exon', order_by=('start','end')))
        numExons = len(exonList)
        if numExons > 1:                            # Only do for genes with two or more exons
            for i in range(0,numExons-1):           # All possible donor exons
                for j in range(1,numExons):         # All possible acceptor exons
                    
                    # Logic check - if exon pair is not possible then start on next exon pair
                    donorStop = exonList[i].end
                    acceptorStart = exonList[j].start-1 # BED output is 0-based, want to remain consistent!
                    donorName = ''.join(exonList[i].attributes['Name'])
                    acceptorName = ''.join(exonList[j].attributes['Name'])
                    junctionID = donorName + "|" + acceptorName
                    if acceptorStart <= donorStop:      # if exon2 overlaps exon1
                        continue                        # move to next exon pair
                    else:                               # if exon2 does not overlap exon1
                        skippedExons = []
                        for k in range(0,numExons):     # Iterate through all exons for gene
                            exonStart = exonList[k].start
                            exonStop = exonList[k].end
                            if (exonStart > donorStop and exonStop < acceptorStart): # If exon falls between donor and acceptor
                                skippedID = ''.join(exonList[k].attributes['Name'])
                                skippedExons.append(skippedID)
                                isSkipped = '1'
                                listEntry=[junctionID, skippedID, isSkipped]
                                listOut.write(",".join(str(i) for i in listEntry)+"\n")
                    if len(skippedExons) == 0:          # If no skipped exons, then output this
                        isSkipped = '0'                 # I think we still need this info!
                        listEntry=[junctionID, '', isSkipped]
                        listOut.write(",".join(str(i) for i in listEntry)+"\n")
                    
                    numSkipped = len(skippedExons)
                    if numSkipped > 0: # Flag if exon skipping junction
                        doesSkip = '1'
                    else:
                        doesSkip = '0'
                    skippedExonList='|'.join(skippedExons)
                    csvEntry = [junctionID, doesSkip, str(numSkipped), skippedExonList]
                    csvOut.write(",".join(str(i) for i in csvEntry)+"\n")
    
    csvOut.close()
    listOut.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
