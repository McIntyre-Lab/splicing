#!/usr/bin/env python

import gffutils

gff_fn='aceview.gff'
db_fn='aceview.gff.db'

gffutils.create_db(gff_fn, db_fn,merge_strategy='create_unique')
#gffutils.create_db(gff_fn, db_fn)

