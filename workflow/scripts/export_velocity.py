#!/usr/bin/env python3
"""Export raw alevin-fry USA blocks without gene/cell filtering or densification."""
import argparse
import csv
import json
import os
from pathlib import Path
import warnings
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from velocity_common import read_usa, validate, write_metrics
from bd_barcode_conversion import convert_barcodes


def export(quant, output, sample, qc, keys, provenance, gene_map=None, reference=None):
    qc = Path(qc)
    qc.mkdir(parents=True, exist_ok=True)
    barcodes, genes, matrices = read_usa(quant)
    indexes, names, report = convert_barcodes(
        barcodes, keys, qc/'barcode_conversion.tsv', reference,
        qc/'bd_reference_differences.tsv')
    S, U, A = matrices
    obs = pd.DataFrame({'barcode_27nt': barcodes,
                        'bd_cell_index': np.array(indexes, dtype=np.int64),
                        'sample': sample},
                       index=pd.Index(names, name='cell_id'))
    var = pd.DataFrame({'gene_id': genes}, index=pd.Index(genes, name='gene_id'))
    if gene_map:
        mapping = {}
        with Path(gene_map).open() as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if len(row) != 2 or not row[0] or not row[1]:
                    raise ValueError('Gene map must contain two nonempty TSV columns, no header')
                if row[0] in mapping and mapping[row[0]] != row[1]:
                    raise ValueError(f'Conflicting gene symbol: {row[0]}')
                mapping[row[0]] = row[1]
        if any(g not in mapping for g in genes):
            raise ValueError('Reference gene-symbol mapping does not cover quantified IDs')
        # Reference names that are Ensembl IDs are placeholders, not gene symbols.
        var['gene_symbol'] = ['' if mapping[g].startswith('ENSG') else mapping[g] for g in genes]
    else:
        var['gene_symbol'] = ''
    combined = (S + A).tocsr()
    combined.eliminate_zeros()
    adata = ad.AnnData(X=combined, obs=obs, var=var,
                      layers={'spliced': combined.copy(), 'unspliced': U.copy()})
    adata.uns['pipeline'] = json.loads(Path(provenance).read_text())
    adata.uns['pipeline'].update(barcode_conversion='ordered BD B384 1-based mixed radix',
                                 bd_naming=report['naming'],
                                 bd_reference_compatibility=report['reference_compatibility'],
                                 gene_symbols='reference mapping; Ensembl placeholders empty' if gene_map else 'unavailable; empty strings')
    metrics = validate(adata)
    metrics.update(total_S=float(S.sum()), total_U=float(U.sum()), total_A=float(A.sum()))
    adata.uns['pipeline']['fractional_counts'] = metrics['fractional_counts']
    if metrics['fractional_counts']:
        warnings.warn('Fractional raw count estimates preserved without rounding')
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + '.tmp.h5ad')
    try:
        adata.write_h5ad(temporary, compression='gzip')
        # Reopen the serialized artifact before publishing it.
        validate(ad.read_h5ad(temporary))
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    metrics['detected_quantified_cells'] = len(barcodes)
    metrics['corrected_distinct_barcodes'] = adata.uns['pipeline'].get('corrected_distinct_barcodes', 'unavailable')
    if (Path(quant)/'generate_permit_list.json').is_file():
        metrics['permit_list_metadata'] = 'available in af_quant/generate_permit_list.json'
    write_metrics(qc/'quantification_summary.tsv', metrics)
    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('quant', 'output', 'sample', 'qc', 'keys', 'provenance'):
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--gene-map')
    parser.add_argument('--bd-reference')
    args = parser.parse_args()
    export(args.quant, args.output, args.sample, args.qc, args.keys, args.provenance,
           args.gene_map, args.bd_reference)
