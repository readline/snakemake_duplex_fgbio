#!/usr/bin/env python
import sys
import pandas as pd

def samplesheet(sspath):
    df = pd.read_csv(sspath, sep='\t', index_col=1)
    sample = {}
    run = {}
    run2sample = {}
    for i in df.index:
        ss = df.loc[i,'sample']
        if ss not in sample:
            sample[ss] = []
        if i not in run:
            sample[ss].append(i)
            run2sample[i] = ss
            run[i] = {'r1':df.loc[i,'read1'], 'r2':df.loc[i,'read2']}
    return sample, run, run2sample

if __name__ == '__main__':
    print(samplesheet(sys.argv[1]))