#!/usr/bin/env python3
"""Write tiny DAG inputs only. The placeholder index/image must NEVER be executed."""
import argparse
import gzip
import json
from pathlib import Path


def create(directory):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    index = directory/'reference/index'
    index.mkdir(parents=True, exist_ok=True)
    for suffix in ('ctab', 'refinfo', 'ssi'):
        (index/f'piscem_idx.{suffix}').touch()
    (index/'simpleaf_index.json').write_text('{"index_type":"piscem"}\n')
    (index/'t2g_3col.tsv').write_text('tx1\tENSG1\tS\ntx2\tENSG1\tU\n')
    (index/'gene_id_to_name.tsv').write_text('ENSG1\tGENE1\n')
    samples=[]
    for sample in ('WTA1', 'WTA2'):
        for mate in (1,2):
            with gzip.open(directory/f'{sample}_R{mate}.fq.gz','wt') as handle:
                handle.write('@fixture\nACGT\n+\nIIII\n')
        samples.append({'name':sample,'uses':{
            'cb_umi_fq':[str(directory/f'{sample}_R1.fq.gz')],
            'cdna_fq':[str(directory/f'{sample}_R2.fq.gz')],
            'allowedlist':384, 'diversity_insets':'yes'}})
    # Exercise optional upstream concatenation/downsampling in the dry-run DAG.
    sample2=samples[1]['uses']
    sample2['cb_umi_fq'] *= 2
    sample2['cdna_fq'] *= 2
    sample2['downsample']=50
    config = dict(working_dir=str(directory/'results'),nthreads=2,bead_chemistry='enhanced_v3',
                  export={'python':'python3'},
                  simpleaf={'index':str(index.parent),'resolution':'usa'},samples=samples)
    path=directory/'config.json'
    path.write_text(json.dumps(config,indent=2)+'\n')
    return path


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    args=parser.parse_args()
    print(create(args.directory))
