#!/usr/bin/env python

import os
import gffutils



def MAIN():
    db_fn='../gff/AceView.sorted.GRCh37.hg19.gtf.db'
    out='../created_files/hg19_aceview_exons.bed'
    db=gffutils.FeatureDB(db_fn)
    exons=db.features_of_type('exon')
    with open(out, 'w') as OUT:
        for exon in exons:
            mylist=[exon.chrom, exon.start-1, exon.stop, exon.id]
            OUT.write("\t".join(str(x) for x in mylist) + "\n")

if __name__=='__MAIN__':
    MAIN()
