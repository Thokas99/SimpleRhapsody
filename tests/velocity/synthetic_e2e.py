#!/usr/bin/env python3
"""Execute the real workflow on tiny synthetic reads and a temporary reference."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'workflow/scripts'))
from bd_barcode_conversion import encode, load_keys
from validate_velocity import validate_outputs
from simpleaf_quant import inspect_reference
from velocity_common import equal_sparse

with tempfile.TemporaryDirectory(prefix='rhapsodist-e2e-') as temporary:
    directory = Path(temporary)
    home = directory/'af_home'
    home.mkdir()
    (home/'chemistries.json').write_text('{}')
    env = dict(os.environ, ALEVIN_FRY_HOME=str(home), XDG_CACHE_HOME=str(directory/'cache'))
    def run(*command):
        try:
            subprocess.run(list(map(str, command)), env=env, check=True, cwd=directory)
        except subprocess.CalledProcessError:
            for log in (directory/'results/logs').glob('*.log'):
                print(log.name, log.read_text(), file=sys.stderr)
            for path in (directory/'results').rglob('quants_mat_cols.txt'):
                print(path.read_text(), file=sys.stderr)
            raise
    run('simpleaf', 'set-paths', '--piscem', shutil.which('piscem'),
        '--alevin-fry', shutil.which('alevin-fry'))
    rng = random.Random(42)
    sequences = [''.join(rng.choices('ACGT', k=1000)) for _ in range(2)]
    fasta = directory/'reference.fa'
    fasta.write_text(''.join(f'>tx{i}\n{seq}\n' for i, seq in enumerate(sequences)))
    reference = directory/'reference'
    run('simpleaf', 'index', '--ref-seq', fasta, '--output', reference,
        '--threads', '2', '--dict', 'sshash', '--work-dir', directory/'index_work')
    index = reference/'index'
    (index/'t2g_3col.tsv').write_text('tx0\tENSG000001\tS\ntx1\tENSG000001\tU\n')
    checked=inspect_reference(index)
    manifest={key:checked[key] for key in ('reference_t2g_sha256','reference_index_meta_sha256')}
    manifest.update(version=1,index_type='spliceu',genome='synthetic',GENCODE_release='not applicable: synthetic reference')
    (index/'reference_manifest.json').write_text(json.dumps(manifest))
    barcode = encode(901, load_keys(ROOT/'workflow/data/whitelist_384x3'))
    reads = [[], []]
    for i in range(40):
        umi = ''.join(rng.choices('ACGT', k=8))
        r1 = 'A'+barcode[:9]+'GTGA'+barcode[9:18]+'GACA'+barcode[18:]+umi+'T'*25
        r2 = sequences[i % 2][100:200]
        for mate, seq in enumerate((r1, r2)):
            reads[mate].append(f'@r{i}\n{seq}\n+\n'+ 'I'*len(seq)+'\n')
    for mate in range(2):
        with gzip.open(directory/f'R{mate+1}.fq.gz', 'wt') as handle:
            handle.write(''.join(reads[mate]))
    config = dict(working_dir=str(directory/'results'), nthreads=2,
                  export={'python':sys.executable}, simpleaf={'index':str(index),'resolution':'usa'},
                  samples=[{'name':'synthetic','uses':{
                      'cb_umi_fq':[str(directory/'R1.fq.gz')],
                      'cdna_fq':[str(directory/'R2.fq.gz')],
                      'allowedlist':384,'diversity_insets':'yes'}}])
    # Two lanes preserve the exact original read order; downsampling stays sample-wide.
    for lane in range(2):
        for mate in range(2):
            with gzip.open(directory/f'lane{lane}_R{mate+1}.fq.gz','wt') as handle:
                handle.write(''.join(reads[mate][lane*20:(lane+1)*20]))
    for name,fraction in (('split',100),('sampled',50)):
        config['samples'].append({'name':name,'uses':{
            'cb_umi_fq':[str(directory/f'lane{lane}_R1.fq.gz') for lane in range(2)],
            'cdna_fq':[str(directory/f'lane{lane}_R2.fq.gz') for lane in range(2)],
            'allowedlist':384,'diversity_insets':'yes','downsample':fraction}})
    config_path = directory/'config.json'
    config_path.write_text(json.dumps(config))
    run('snakemake', '--snakefile', ROOT/'workflow/Snakefile', '--configfile', config_path,
        '--cores', '2', '--printshellcmds')
    sample = directory/'results/synthetic'
    completed = [p for p in (directory/'results').rglob('*') if p.is_file()]
    before_resume = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in completed}
    run('snakemake', '--snakefile', ROOT/'workflow/Snakefile', '--configfile', config_path,
        '--cores', '2')
    assert before_resume == {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in completed}, 'Completed outputs recomputed'
    metrics = validate_outputs(sample/'synthetic_velocity.h5ad', sample/'simpleaf/af_quant',
                               ROOT/'workflow/data/whitelist_384x3')
    assert not (sample/'raw').exists() and not (sample/'usa_counts').exists()
    assert metrics['cross_format'] == 'PASS', metrics
    import anndata
    adata = anndata.read_h5ad(sample/'synthetic_velocity.h5ad')
    assert adata.obs_names.tolist() == ['901'], adata.obs_names
    assert adata.n_vars == 1 and adata.layers['spliced'].sum() > 0 and adata.layers['unspliced'].sum() > 0
    assert adata.uns['pipeline']['genome']=='synthetic'
    provenance=json.loads((sample/'simpleaf/pipeline_provenance.json').read_text())
    assert provenance['workflow_lock_sha256']==hashlib.sha256((ROOT/'pixi.lock').read_bytes()).hexdigest()
    assert provenance['preparation_qc_sha256']==hashlib.sha256((sample/'qc/preparation_qc.json').read_bytes()).hexdigest()
    for field in ('workflow_lock_sha256','workflow_config_sha256','preparation_qc_sha256'):
        assert adata.uns['pipeline'][field]==provenance[field]
        assert len(provenance[field])==64
    assert list(adata.uns['pipeline']['workflow_source_paths'])==provenance['workflow_source_paths']
    assert list(adata.uns['pipeline']['workflow_source_sha256'])==provenance['workflow_source_sha256']
    for relative,digest in zip(provenance['workflow_source_paths'],provenance['workflow_source_sha256']):
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==digest
    split=anndata.read_h5ad(directory/'results/split/split_velocity.h5ad')
    assert split.obs_names.tolist()==adata.obs_names.tolist()
    assert split.var_names.tolist()==adata.var_names.tolist()
    for layer in ('spliced','unspliced'):
        assert equal_sparse(split.layers[layer],adata.layers[layer]), layer
    assert not (directory/'results/data/fastq/combined/split_R1.fastq.gz').exists()
    assert (directory/'results/data/fastq/combined/sampled_R1.fastq.gz').exists()
    sampled=json.loads((directory/'results/sampled/qc/preparation_qc.json').read_text())
    assert sampled['input_read_pairs']==40
    assert 0 < sampled['prepared_read_pairs'] < 40
    assert not list((directory/'results/data/fastq/units').rglob('*.gz'))
    for name in ('split','sampled'):
        assert validate_outputs(directory/f'results/{name}/{name}_velocity.h5ad',
            directory/f'results/{name}/simpleaf/af_quant',ROOT/'workflow/data/whitelist_384x3')['cross_format']=='PASS'
    print('PASS: single, split-lane and downsampled synthetic FASTQs → SimpleAF USA → H5AD; manifest and resume checks')
