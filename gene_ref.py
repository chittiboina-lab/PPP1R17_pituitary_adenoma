#!/usr/bin/env python

import pyensembl as ensembl

ensembldb = ensembl.EnsemblRelease(release = 86)

genes = ensembldb.genes()
with open('genes.bed', 'w') as fout:
    for gene in genes: 
        if gene.biotype == 'protein_coding' and gene.contig != 'MT':
            fout.write('chr%s\t%s\t%s\t%s\n' % (gene.contig, gene.start, gene.end, gene.gene_name))

