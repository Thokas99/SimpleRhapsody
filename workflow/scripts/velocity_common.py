"""Small shared I/O and sparse count checks for the velocity companion."""
import csv
import json
import os
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.io import mmread


def write_metrics(path, metrics):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['metric', 'value'])
        writer.writerows(metrics.items())
    os.replace(temporary, path)


def check_counts(matrix):
    if not sparse.issparse(matrix) or matrix.format != 'csr':
        raise ValueError('Count matrix must be sparse CSR')
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError('Counts must be finite and nonnegative')
    return bool(np.equal(matrix.data, np.floor(matrix.data)).all())


def read_usa(directory):
    directory = Path(directory)
    meta = json.loads((directory/'quant.json').read_text())
    if meta.get('usa_mode') is not True:
        raise ValueError('quant.json must confirm usa_mode=true')
    matrix = mmread(directory/'alevin/quants_mat.mtx', spmatrix=True)
    if not sparse.issparse(matrix):
        raise ValueError('Expected sparse Matrix Market coordinate data')
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError('Counts must be finite and nonnegative')
    matrix = matrix.tocsr()
    check_counts(matrix)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    genes = (directory/'alevin/quants_mat_cols.txt').read_text().splitlines()
    cells = (directory/'alevin/quants_mat_rows.txt').read_text().splitlines()
    n_states = meta.get('num_genes')
    if not isinstance(n_states, int) or n_states <= 0 or n_states % 3:
        raise ValueError('USA num_genes must describe three complete state blocks')
    n_genes = n_states // 3
    if matrix.shape != (len(cells), n_states):
        raise ValueError('USA matrix dimensions disagree with metadata/cell order')
    if len(genes) == n_states:
        base = genes[:n_genes]
        repeated = base * 3
        labelled = base + [g+'-U' for g in base] + [g+'-A' for g in base]
        if genes != repeated and genes != labelled:
            raise ValueError('S/U/A gene blocks have different gene order')
        genes = genes[:n_genes]
    if len(genes) != n_genes or len(set(genes)) != n_genes or any(not gene for gene in genes):
        raise ValueError('Missing/duplicate gene IDs or incorrect gene count')
    if not all(gene.startswith('ENSG') for gene in genes):
        raise ValueError('Expected stable GENCODE Ensembl gene IDs, not symbols')
    return cells, genes, [matrix[:, i*n_genes:(i+1)*n_genes].tocsr() for i in range(3)]


def equal_sparse(left, right):
    return left.shape == right.shape and (left != right).nnz == 0


def validate(adata):
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError('Empty count matrix')
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError('Duplicate cell or gene identifiers')
    integer = check_counts(adata.X)
    if (set(adata.layers) - {None}) != {'spliced', 'unspliced'}:
        raise ValueError('Expected only spliced and unspliced layers')
    for layer in ('spliced', 'unspliced'):
        matrix = adata.layers[layer]
        if matrix.shape != adata.shape:
            raise ValueError(f'Wrong {layer} dimensions')
        integer = check_counts(matrix) and integer
    if not equal_sparse(adata.layers['spliced'], adata.X):
        raise ValueError('spliced differs from S+A')
    for field in ('barcode_27nt', 'bd_cell_index'):
        if field not in adata.obs or not adata.obs[field].is_unique:
            raise ValueError(f'Missing or duplicated {field}')
    if list(adata.var['gene_id']) != list(adata.var_names):
        raise ValueError('Gene identity/order mismatch')
    return dict(status='PASS', cells=adata.n_obs, genes=adata.n_vars,
                integer_valued=integer, fractional_counts=not integer,
                total_spliced=float(adata.layers['spliced'].sum()),
                total_U=float(adata.layers['unspliced'].sum()),
                X_sparsity=1-adata.X.nnz/(adata.n_obs*adata.n_vars))
