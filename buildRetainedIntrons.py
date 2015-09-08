#!/usr/bin/env python

##        DESCRIPTION: This program takes a GTF or GFF file and creates putative intron retention events
##        from exons. It outputs two files: a 12-column BED file for intron retention events, and a TSV
##        describing exon(s) associated with an event. Intron retention events are only made for introns
##        that are unambiguously intronic sequence. Events are not made for intron sequences that overlap
##        with exonic regions, For plus-strand genes, retained introns are made from the end of exons. For
##        minus-strand genes, retained introns are made from the start of exons. Events are not made for the
##        last exons for a gene (genomic 3' last for plus-strand genes, genomic 5' first for minus-strand genes).
##        Information for last exons are output to only the TSV file with the event ID exonID|no_event" for 
##        later reference if necessary. The intron_position column in the TSV is the start positon (plus strand)
##        or stop position (minus strand) of the intron used in making the intron retention event.
##        
##        Author: Jeremy R. B. Newman

# Built-in packages
import argparse

# Add-on packages
import gffutils


def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--gff", dest="gffInput", action='store', required=True, help="Input GFF or GTF file")
    parser.add_argument("--obed", dest="outBED", action='store', required=True, help="Output BED file name")
    parser.add_argument("--otable", dest="outCSV", action='store', required=True, help="Output CSV file name")
    parser.add_argument("--size", dest="eventSize", action='store', type=int, required=True, help="Size of intron retention event parts (event size will be twice this)")

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
    csvHeader = ','.join(['gene_id','event_id','chr','strand','intron_position','exon_id','exon_cat','flag_lastexon']) + "\n"
    csvOut.write(csvHeader)
    
    # Make events
    genes=db.features_of_type('gene')
    
    for gene in genes:
        if gene.strand == '+':
            # If plus strand, order exons by ascending start position
            exonList = list(db.children(gene, featuretype='exon', order_by='start')) 
        elif gene.strand == '-':
            # If plus strand, order exons by descending stop position
            exonList = list(db.children(gene, featuretype='exon', order_by='end', reverse=True))
        exonList.append(exonList[len(exonList)-1])
        index=1
        exonIndex=0
        lastExon=len(exonList)-1 # Check for last exon
        exonNameList=[]
        for exon in exonList:
            if index == 1:
                #Get ID start stop parents, etc.
                exonStart=exon.start
                exonStop=exon.end
                exonID=''.join(exon.attributes['Name'])
                currExonStart=exon.start
                currExonStop=exon.end
                currExonID=''.join(exon.attributes['Name'])
                exonNameList.append(exonID)
                index=index+1
                exonIndex=exonIndex+1
            else:
                #Check with previous
                exonStart=exon.start
                exonStop=exon.end
                exonID=''.join(exon.attributes['Name'])
                
                # Overlapping exons check
                if gene.strand == '+':
                    if (exonStart <= currExonStop):     # If exon overlaps with current exon
                        # if exon overlaps check to see what has the bigger coordinates
                        if (exonStop >= currExonStop):  # If exon ends later than current exon
                            currExonStart=exonStart     # Make exon current exon
                            currExonStop=exonStop
                            currExonID=exonID
                        exonNameList.append(exonID)     # Append exonID to list of exonIDs
                        exonIndex=exonIndex+1
                    else:                               # Otherwise if new exon write outputs for current exon
                        exonLength=currExonStop-currExonStart+1
                        if exonLength > args.eventSize:             # Need to make this number customizable
                            donorStart=currExonStop - args.eventSize
                        else:
                            donorStart=currExonStop - exonLength
                        donorStop=currExonStop
                        acceptorStart=currExonStop+1
                        acceptorStop=currExonStop+args.eventSize
                        intronPos=acceptorStart
                        
                        # Write output for BED file
                        eventChrom=gene.chrom
                        totalStart=donorStart
                        totalStop=acceptorStop
                        eventName=currExonID + "|intron"
                        eventScore="."
                        eventStrand=gene.strand
                        color="255,0,0"
                        blocks="1"
                        blockLength=str(totalStop-totalStart)
                        blockStarts="0"
                        exonListOut='|'.join(exonNameList)
                        exonEntry='\t'.join([str(eventChrom), str(totalStart), str(totalStop), eventName, eventScore, eventStrand, color, blocks, blockLength, blockStarts]) + "\n"
                        bedOut.write(exonEntry)

                        # Write output for TSV file
                        eventChrom=gene.chrom
                        eventGene=gene.id
                        eventName=currExonID + "|intron"
                        isLastExon='0'
                        eventStrand=gene.strand
                        exonListOut='|'.join(exonNameList)
                        csvEntry=','.join([eventGene, eventName, str(eventChrom), eventStrand, str(intronPos), currExonID, exonListOut, isLastExon]) + "\n"
                        csvOut.write(csvEntry)
                        
                        # Get pieces for next exon
                        exonNameList=[]
                        exonID=''.join(exon.attributes['Name'])
                        currExonStart=exon.start
                        currExonStop=exon.end
                        currExonID=''.join(exon.attributes['Name'])
                        exonNameList.append(exonID)
                        index=index+1
                        exonIndex=exonIndex+1
                
                if gene.strand == '-':
                    if (exonStop >= currExonStart):     # If exon overlaps with current exon
                        # if exon overlaps check to see what has the bigger coordinates
                        if (exonStart <= currExonStart): # If exon start earlier than current exon
                            currExonStart=exonStart      # Make exon current exon
                            currExonStop=exonStop
                            currExonID=exonID
                        exonNameList.append(exonID)
                        exonIndex=exonIndex+1
                    else:                               # Otherwise if new exon write outputs for current exon
                        exonLength=currExonStop-currExonStart+1
                        if exonLength > args.eventSize:             # Need to make this number customizable
                            acceptorStop=currExonStart + args.eventSize
                        else:
                            acceptorStop=currExonStart + exonLength
                        acceptorStart=currExonStart
                        donorStart=currExonStart-args.eventSize
                        donorStop=currExonStart-1
                        intronPos=donorStop
                        
                        # Write output for BED file
                        eventChrom=gene.chrom
                        totalStart=donorStart
                        totalStop=acceptorStop
                        eventName=currExonID + "|intron"
                        eventScore="."
                        eventStrand=gene.strand
                        color="255,0,0"
                        blocks="1"
                        blockLength=str(totalStop-totalStart)
                        blockStarts="0"
                        exonListOut='|'.join(exonNameList)
                        exonEntry='\t'.join([str(eventChrom), str(totalStart), str(totalStop), eventName, eventScore, eventStrand, color, blocks, blockLength, blockStarts]) + "\n"
                        bedOut.write(exonEntry)

                        # Write output for TSV file
                        eventChrom=gene.chrom
                        eventGene=gene.id
                        eventName=currExonID + "|intron"
                        isLastExon='0'
                        eventStrand=gene.strand
                        exonListOut='|'.join(exonNameList)
                        csvEntry=','.join([eventGene, eventName, str(eventChrom), eventStrand, str(intronPos), currExonID, exonListOut, isLastExon]) + "\n"
                        csvOut.write(csvEntry)
                        
                        # Get pieces for next exon
                        exonNameList=[]
                        exonID=''.join(exon.attributes['Name'])
                        currExonStart=exon.start
                        currExonStop=exon.end
                        currExonID=''.join(exon.attributes['Name'])
                        exonNameList.append(exonID)
                        index=index+1
                        exonIndex=exonIndex+1
            if exonIndex == lastExon:
                #Get ID start stop parents, etc.
                exonStart=exon.start
                exonStop=exon.end
                eventName=currExonID + "|no_event"
                if gene.strand == '+':
                    intronPos=exonStop+1
                elif gene.strand == '-':
                    intronPos=exonStart-1
                exonID=''.join(exon.attributes['Name'])
                currExonStart=exon.start
                currExonStop=exon.end
                currExonID=''.join(exon.attributes['Name'])
                #exonNameList.append(exonID)
                intronPos=-9999999
                
                # Write output for TSV file
                eventChrom=gene.chrom
                eventGene=gene.id
                eventName=currExonID + "|no_event"
                isLastExon='1'
                eventStrand=gene.strand
                exonListOut='|'.join(exonNameList)
                csvEntry=','.join([eventGene, eventName, str(eventChrom), eventStrand, str(intronPos), currExonID, exonListOut, isLastExon]) + "\n"
                csvOut.write(csvEntry)
                
                index=index+1
                exonIndex=exonIndex+1
    
    bedOut.close()
    csvOut.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()


