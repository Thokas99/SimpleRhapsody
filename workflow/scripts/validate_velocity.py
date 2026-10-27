#!/usr/bin/env python3
"""Independently reopen and validate a serialized velocity-count companion."""
import argparse
import anndata as ad
from pathlib import Path
from bd_barcode_conversion import load_keys, encode, decode
from velocity_common import read_usa, equal_sparse, validate, write_metrics


def _validate_barcodes(adata, keys):
    panels = load_keys(Path(keys))
    inverse = [{s: i for i, s in enumerate(p)} for p in panels]
    for barcode, index in zip(adata.obs['barcode_27nt'], adata.obs['bd_cell_index']):
        if decode(barcode, inverse) != str(index) or encode(index, panels) != barcode:
            raise ValueError('Serialized barcode/index mapping disagrees')
    if adata.uns['pipeline']['bd_naming'] == 'decimal Cell_index':
        if list(adata.obs_names) != adata.obs['bd_cell_index'].astype(str).tolist():
            raise ValueError('Serialized BD naming mismatch')


def validate_outputs(h5ad_path, quant, keys):
    adata = ad.read_h5ad(h5ad_path)
    metrics = validate(adata)
    _validate_barcodes(adata, keys)
    barcodes, genes, matrices = read_usa(quant)
    if barcodes != adata.obs['barcode_27nt'].tolist() or genes != list(adata.var_names):
        raise ValueError('Native SimpleAF cell/gene identifier order differs from H5AD')
    S, U, A = matrices
    for state, matrix in (('spliced', S + A), ('unspliced', U)):
        if not equal_sparse(matrix, adata.layers[state]):
            raise ValueError(f'Native SimpleAF {state} counts differ from H5AD')
    metrics.update(native_matrix_orientation='cells x USA gene states', cross_format='PASS')
    return metrics


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--quant', required=True)
    parser.add_argument('--keys', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    write_metrics(args.output, validate_outputs(args.input, args.quant, args.keys))
