#!/usr/bin/env python

##        DESCRIPTION: This program takes a sorted GTF and reformats it to a GFF3 file with a similar structure
##        to FlyBase GFF3 files. It outputs a GFF3 file containing entries for chromosomes, genes, mRNAs,
##        exons, CDSs and introns. A known limitation is that the chromosome end position is 1000 bp
##        after the last gene on a chromosome and so does not likely reflect the true size of a chromsome.
##        As GTF files typically don't include chromosome entries this is acceptable workaround. The output
##        GFF should work with current McLab scripts.
##        NOTE: this script will take a long time for a complete GTF file. A recommended workaround is to split
##        your input GTF file by chromosome and convert each individually, then merge the final GFF3 outputs.
##        
##        Author: Jeremy R. B. Newman

### FLYBASE FORMAT
## FIELDS
# chrom, source, feature_type, start, end, score, strand, frame, attribute
## IMPORTANT ATTRIBUTE FIELD:
# Genes:	ID, Name
# Transcripts:	ID, Name, Parent
# Exons:	Name (Gene:Num), Parent, parent_type=mRNA (multiple parent mRNAs)
# Introns:	Name (Gene-in), Parent, parent_type_mRNA (one parent mRNA)
# CDS:		Name (Gene-cds), Parent, parent_type (one parent mRNA)

# Built-in packages
import argparse
from operator import itemgetter

# Add-on packages
import gffutils


def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--input", dest="gffInput", action='store', required=True, help="Input GTF file, sorted by chromosome, start and stop positions")
    parser.add_argument("--output", dest="outputFile", action='store', required=True, help="Output GFF3 file name")

    args = parser.parse_args()
    return args

def main():

    ### Write GFF output
    
    with open(args.outputFile, 'wb') as outputFile:
        outputFile.write("##gff-version 3" + "\n")

        # Import database
        db_fn= args.gffInput + '.db'
        db = gffutils.FeatureDB(db_fn)
        
        # Get Number of genes - doing this further down doesn't seem to let the entries print so I've put it up here
        genes = db.features_of_type('gene')
        geneNum=0
        geneLastIndex=len(list(genes))
        
        # Set GFF output array
        
        gffList=[]
        geneEntry=[]
        xscriptEntry=[]
        cdsEntry=[]
        exonEntry=[]
        intronEntry=[]
         
        # Iterate through each gene
        # For each gene do this
        genes = db.features_of_type('gene', order_by=('seqid','start','end'))
        for gene in genes:
            # For getes we need ID, Name, start, end, chrom, strand
            # Can use the chrom and strand information for all child objects
            geneID=gene.id
            geneChrom=gene.chrom
            geneStart=gene.start
            geneStop=gene.end
            geneStrand=gene.strand
            geneSource=str(gene.source) + '_gff3_converted'
            
            # Format GFF entry and append to GFF list
            geneEntry=[geneChrom, geneSource, 'gene', str(geneStart), str(geneStop), '.', geneStrand, '.', "ID=" + str(gene.id) + ";Name=" + str(gene.id)]
            outputFile.write("\t".join(geneEntry) + "\n")
            
            # Format transcripts and CDS
            xscriptList = list(db.children(gene, featuretype='transcript', order_by=('start','end')))
            # For transcripts we need ID, Name, start, end
            for xscript in xscriptList:
                xscriptID=xscript.id
                xscriptStart=xscript.start
                xscriptStop=xscript.end
                
                # Format GFF entry and append to GFF list
                xscriptEntry=[geneChrom, geneSource, 'mRNA', str(xscriptStart), str(xscriptStop), '.', geneStrand, '.', "ID=" + xscriptID + ";Name=" + xscriptID + ";Parent=" + geneID]
                outputFile.write("\t".join(xscriptEntry) + "\n")
                
                # For CDS we need start, end, frame
                cdsList = list(db.children(xscript, featuretype='CDS', order_by='start'))
                for cds in cdsList:
                    cdsStart=cds.start
                    cdsStop=cds.end
                    cdsFrame=cds.frame
                    cdsID=str(gene.id) + "-cds"
                    
                    #Format GFF entry and append to GFF list
                    cdsEntry=[geneChrom, geneSource , 'CDS' , str(cdsStart), str(cdsStop), '.', geneStrand, str(cdsFrame), "Name=" + cdsID + ";Parent=" + xscriptID + ";parent_type=mRNA"]
                    outputFile.write("\t".join(cdsEntry) + "\n")
                
                # For Introns we need start, end, frame
                # Doing this by transcript rather than gene. Could change this later if necessary
                intronList = list(db.children(xscript, featuretype='intron', order_by=('start','end')))
                intronIndex=1
                for intron in intronList:
                    intronStart=intron.start
                    intronStop=intron.end
                    intronID=str(gene.id) + "-in" + str(intronIndex)
                    intronIndex=intronIndex+1 # Append an intron number for each intron in a transcript
                    
                    # Format GFF entry and append to GFF list
                    intronEntry=[geneChrom, geneSource, 'intron', str(intronStart), str(intronStop), '.', geneStrand, '.', "Name=" + intronID + ";Parent=" + xscriptID + ";parent_type=mRNA"]
                    outputFile.write("\t".join(intronEntry) + "\n")
            
            # Merge duplicate exons and format exons
            exonList = list(db.children(gene, featuretype='exon', order_by=('start','end')))
            
            # Adds in the last exon again so that it actually outputs
            # Could probably change this to handle similar to chromosome but this works right now
            exonList.append(exonList[len(exonList)-1])
            index=1
            exonIndex=0
            lastExon=len(exonList)-1 ## Check if last exon!
            for exon in exonList:
                
                if index == 1: # If first exon in a gene get all the required pieces
                    
                    #Get ID start stop parents, etc.
                    exonID=str(gene.id) + ":" + str(index)
                    exonStart=exon.start
                    exonStop=exon.end
                    parentList = list(db.parents(exon, featuretype='transcript', order_by='start'))
                    exonParents = []
                    for parent in parentList:
                        parentID=str(parent.id)
                        exonParents.append(parentID)
                    currExonStart=exon.start # These entries are used as a reference point for checking
                    currExonStop=exon.end    # if subsequent exons are duplicates or unique 
                    currExonID=str(gene.id) + ":" + str(index)
                    currExonParents=exonParents
                    index=index+1
                    exonIndex=exonIndex+1
                
                else:
                    
                    #Check to see if exon matched exon currently stored
                    exonStart=exon.start
                    exonStop=exon.end
                    
                    if (exonStart == currExonStart and exonStop == currExonStop):
                        # if the same coordinates as the current exon then only need to append parent transcript to list
                        parentList = list(db.parents(exon, featuretype='transcript', order_by='start'))
                        for parent in parentList:
                            parentID=str(parent.id)
                            exonParents.append(parentID)
                        exonIndex=exonIndex+1
                    
                    else:
                        # if new exon format GFF entry and append to GFF list
                        # Then proceed with building the new exon entry                
                        currExonParentsString=','.join(currExonParents)
                        exonEntry=[geneChrom, geneSource, 'exon', str(currExonStart), str(currExonStop), '.', geneStrand, '.', "Name=" + currExonID + ";Parent=" + currExonParentsString + ";parent_type=mRNA"]
                        outputFile.write("\t".join(exonEntry) + "\n")
                        
                        exonID=str(gene.id) + ":" + str(index)
                        parentList = list(db.parents(exon, featuretype='transcript', order_by='start'))
                        exonParents = []
                        for parent in parentList:
                            parentID=str(parent.id)
                            exonParents.append(parentID)
                        currExonStart=exon.start
                        currExonStop=exon.end
                        currExonID=str(gene.id) + ":" + str(index)
                        currExonParents=exonParents
                        index=index+1
                        exonIndex=exonIndex+1
                
                # If last exon in gene format GFF entry and append
                if exonIndex == lastExon:
                    currExonParentsString=','.join(currExonParents)
                    exonEntry=[geneChrom, geneSource, 'exon', str(currExonStart), str(currExonStop), '.', geneStrand, '.', "Name=" + currExonID + ";Parent=" + currExonParentsString + ";parent_type=mRNA"]
                    outputFile.write("\t".join(exonEntry) + "\n")
            
            # Adding in chromosome entries
            # Check that if gene chromosome is the same as the one stored in memory then print the stored chromosome entry
            # Otherwise build new extry
            # Because most GTF/GFF files don't seem to store chromosome as an entry we will have to build this manually
            # For now chromosome ranges are 1 to 1000bp 3' of last gene on chromosome. This should be sufficient for most things
            
            # If first gene set chromosome number and stop position
            if geneNum == 0:
                chromName=gene.chrom
                chromStop=gene.end
                currChromName=gene.chrom
                currChromStop=gene.end
                geneNum=geneNum+1
            
            # If not first gene then check with stored chromosome information
            else:
                chromName=gene.chrom
                chromStop=gene.end
                
                if chromName == currChromName: #if same chromosome
                    if chromStop > currChromStop: #if stop position is later than stored stop
                        currChromStop=chromStop #set new stop as chromosome stop
                    geneNum=geneNum+1
                
                else: # if new chromosome
                    chromEnd=currChromStop+1000
                    chromEntry=[currChromName, geneSource, 'chromosome', '1', str(chromEnd), '.' ,'.' ,'.', "ID=" + str(currChromName) + ";Name=" + str(currChromName)]
                    outputFile.write("\t".join(chromEntry) + "\n")
                    geneNum=geneNum+1
                    currChromName=chromName
                    currChromStop=chromStop
            
            if geneNum == geneLastIndex:
                chromEnd=currChromStop+1000
                chromEntry=[currChromName, geneSource, 'chromosome', '1', str(chromEnd), '.' ,'.' ,'.', "ID=" + str(currChromName) + ";Name=" + str(currChromName)]
                outputFile.write("\t".join(chromEntry) + "\n")
        
   

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()




