#!/usr/bin/env python
names = ['P1_T', 'P1_N', 'P2_T', 'P2_N', 'P3_T', 'P4_T']
for name in names:
    with open('../results/%s/%s_cort_barcodes.csv' % (name, name), 'r') as fin:
        with open('../results/%s/%s_cort_barcodes_formatted.csv' % (name, name), 'w') as fout:
            for line in fin:
                line = line.replace('"', '').strip('\n').split(',')
                if line[1] == 'x':
                    pass
                else:
                    fout.write('%s\n' % (line[1]))