"""Fast development checks only. No mapping, index building, or workflow execution."""
import ast
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'workflow/scripts'))
sys.path.insert(0, str(ROOT/'workflow/src'))
import bd_barcode_conversion as bd
import simpleaf_quant as quant
from export_velocity import export
from validate_velocity import validate_outputs
from velocity_common import equal_sparse, read_usa
import required_helpers
import prepare_reads as preparation
import anndata as ad
import numpy as np
from scipy import sparse
from scipy.io import mmwrite

PANELS = ROOT/'workflow/data/whitelist_384x3'


def rule_block(text, name):
    return re.search(r'^rule '+name+r':\n.*?(?=^rule |^if config|\Z)', text, re.M|re.S).group()


def original():
    return subprocess.check_output(['git', 'show', '7991f870fdfe38e12648673d330cbfd4839e755d:workflow/Snakefile'], cwd=ROOT, text=True)


def fixture(directory, fractional=False):
    directory = Path(directory)
    (directory/'alevin').mkdir(parents=True)
    S = sparse.csr_matrix([[1, 0], [2, 3]], dtype=float)
    U = sparse.csr_matrix([[0, 4], [5, 0]], dtype=float)
    A = sparse.csr_matrix([[6, 0], [0, 7]], dtype=float)
    if fractional:
        S.data[0] = 0.25
    mmwrite(directory/'alevin/quants_mat.mtx', sparse.hstack([S, U, A], format='csr'))
    keys = bd.load_keys(PANELS)
    cells = [bd.encode(n, keys) for n in (1, 384**3)]
    (directory/'alevin/quants_mat_rows.txt').write_text('\n'.join(cells)+'\n')
    (directory/'alevin/quants_mat_cols.txt').write_text('ENSG000001.1\nENSG000002.9\n')
    (directory/'quant.json').write_text(json.dumps({'usa_mode': True, 'num_genes': 6}))
    return S, U, A


class Conversion(unittest.TestCase):
    def setUp(self):
        self.keys = bd.load_keys(PANELS)
        self.inverse = [{s: i for i, s in enumerate(k)} for k in self.keys]

    def test_official_panel_order_and_known_pairs(self):
        source = ROOT/'workflow/data/bd_conversion/rhapsody_cell_label.py'
        tree = ast.parse(source.read_text())
        namespace = {}
        for statement in tree.body:
            if isinstance(statement, ast.Assign):
                try:
                    value = ast.literal_eval(statement.value)
                except (ValueError, TypeError):
                    continue
                for target in statement.targets:
                    if isinstance(target, ast.Name):
                        namespace[target.id] = value
        self.assertEqual(self.keys, [list(namespace[f'B384_cell_key{i}']) for i in (1,2,3)])
        for statement in tree.body:
            if isinstance(statement, ast.FunctionDef) and statement.name in ('label_sections_to_index', 'index_to_sequence'):
                exec(compile(ast.Module(body=[statement], type_ignores=[]), '<official BD functions>', 'exec'), namespace)
        for index in (1, 901, 4748181, 52923177, 384**3):
            sequence = namespace['index_to_sequence'](index, 'EnhV2')
            diversity = namespace['EnhV2_diversity_inserts'][self.keys[0][(index-1)//384**2]]
            offset = len(diversity)
            canonical = sequence[offset:offset+9]+sequence[offset+13:offset+22]+sequence[offset+26:offset+35]
            self.assertEqual(bd.encode(index, self.keys), canonical)
            self.assertEqual(bd.decode(canonical, self.inverse), str(index))
        self.assertEqual(bd.encode(1, self.keys), 'TGTGTTCGCTTGTGTTGTTTGTGGCTG')

    def test_roundtrip_every_panel_position_and_boundaries(self):
        indexes = {1, 384**3, 384, 385, 384**2, 384**2+1}
        for p in range(3):
            indexes.update(i*384**p+1 for i in range(384))
        for index in indexes:
            self.assertEqual(bd.decode(bd.encode(index, self.keys), self.inverse), str(index))
        for index in (0, -1, 384**3+1, '1.0', 'x', '1_suffix'):
            with self.assertRaises(ValueError):
                bd.encode(index, self.keys)
        for sequence in ('A'*26, 'N'*27, 'A'*27, 'T'*28):
            with self.assertRaises(ValueError):
                bd.decode(sequence, self.inverse)

    def test_real_bd4_identifier_representation_examples(self):
        # These strings were read from the real WTA1 H5AD observation index only.
        names=['91','543','1123','1325','2001']
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            reference=d/'bd_ids.txt'
            reference.write_text('\n'.join(names)+'\n')
            _, converted, report=bd.convert_barcodes([bd.encode(n,self.keys) for n in names],
                PANELS,d/'report.tsv',reference)
            self.assertEqual(converted,names)
            self.assertEqual(report['bd_reference_matches'],5)
            self.assertEqual(report['roundtrip_failures'],0)

    def test_duplicate_rejection_and_exact_reference_differences(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            barcodes = [bd.encode(n, self.keys) for n in (1,384**3)]
            reference = d/'reference.txt'
            reference.write_text('1\n901\n')
            indexes, names, report = bd.convert_barcodes(barcodes, PANELS, d/'qc.tsv', reference, d/'diff.tsv')
            self.assertEqual(names, ['1', str(384**3)])
            self.assertEqual(report['bd_reference_matches'], 1)
            self.assertEqual(report['bd_reference_unmatched'], 1)
            self.assertEqual(report['reference_only'], 1)
            self.assertIn('reference_only\t901', (d/'diff.tsv').read_text())
            with self.assertRaisesRegex(ValueError, 'Duplicate'):
                bd.convert_barcodes(barcodes*2, PANELS, d/'duplicate.tsv')
            self.assertIn('\t2\t', (d/'duplicate.tsv').read_text())
            with self.assertRaisesRegex(ValueError, 'Noncanonical'):
                bd.convert_barcodes([barcodes[0].lower()], PANELS, d/'bad.tsv')
            _, missing_names, missing_report = bd.convert_barcodes(barcodes, PANELS, d/'missing.tsv', d/'absent.txt')
            self.assertEqual(missing_names, ['1', str(384**3)])
            self.assertEqual(missing_report['reference_compatibility'], 'unverified: missing reference file')
            reference.write_text('bd_cell_index\tbd_cell_id\n1\t1_runA\n56623104\t56623104_runA\n')
            _, names, report = bd.convert_barcodes(barcodes, PANELS, d/'full.tsv', reference)
            self.assertEqual(names, ['1_runA', '56623104_runA'])


class Export(unittest.TestCase):
    def test_missing_gene_symbols_preserve_gene_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            fixture(directory/'usa')
            provenance = directory/'provenance.json'
            provenance.write_text('{}')
            output = directory/'WTX_velocity.h5ad'
            export(directory/'usa', output, 'WTX', directory/'qc', PANELS, provenance)
            result = ad.read_h5ad(output)
            self.assertEqual(result.var_names.tolist(), ['ENSG000001.1','ENSG000002.9'])
            self.assertEqual(result.var['gene_symbol'].tolist(), ['', ''])
            self.assertEqual(result.uns['pipeline']['gene_symbols'], 'unavailable; empty strings')

    def test_sparse_export_reopen_and_fractional_preservation(self):
        for fractional in (False, True):
            with self.subTest(fractional=fractional), tempfile.TemporaryDirectory() as d:
                d = Path(d)
                S, U, A = fixture(d/'usa', fractional)
                provenance = d/'provenance.json'
                provenance.write_text(json.dumps({'rhapsodist_commit': 'fixture', 'quantification_resolution':'USA'}))
                symbols = d/'genes.tsv'
                symbols.write_text('ENSG000001.1\tA\nENSG000002.9\tB\n')
                output = d/'WTX_velocity.h5ad'
                export_ad = SimpleNamespace(AnnData=ad.AnnData, read_h5ad=ad.read_h5ad)
                with patch('export_velocity.ad', export_ad), patch.object(
                        export_ad, 'AnnData', wraps=ad.AnnData) as constructor:
                    if fractional:
                        with self.assertWarnsRegex(UserWarning, 'Fractional raw count estimates'):
                            export(d/'usa', output, 'WTX', d/'qc', PANELS, provenance, gene_map=symbols)
                    else:
                        export(d/'usa', output, 'WTX', d/'qc', PANELS, provenance, gene_map=symbols)
                constructor.assert_called_once()
                self.assertEqual(set(constructor.call_args.kwargs['layers']) - {None},
                                 {'spliced', 'unspliced'})
                a = ad.read_h5ad(output)
                self.assertTrue(a.obs_names.is_unique and a.var_names.is_unique)
                self.assertEqual(list(a.obs_names), ['1', str(384**3)])
                self.assertEqual(list(a.var_names), ['ENSG000001.1', 'ENSG000002.9'])
                self.assertEqual(a.obs['bd_cell_index'].tolist(), [1, 384**3])
                self.assertEqual((set(a.layers) - {None}), {'spliced', 'unspliced'})
                for name, expected in [('spliced', S+A),('unspliced', U)]:
                    self.assertTrue(sparse.isspmatrix_csr(a.layers[name]))
                    self.assertTrue(equal_sparse(a.layers[name], expected))
                self.assertTrue(equal_sparse(a.X, S+A))
                result = validate_outputs(output, d/'usa', PANELS)
                self.assertEqual(result['fractional_counts'], fractional)
                self.assertEqual(result['cross_format'], 'PASS')
                self.assertTrue((d/'qc/barcode_conversion.tsv').is_file())

    def test_native_labels_and_counts_preserved_without_duplicate_exports(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            native = directory/'usa'
            fixture(native)
            genes = ['ENSG000001.1', 'ENSG000002.9', 'ENSG000003.2']
            S = sparse.csr_matrix([[1, 2, 0], [0, 0, 0]])
            U = sparse.csr_matrix((2, 3))
            A = sparse.csr_matrix([[0, 3, 0], [0, 0, 0]])
            mmwrite(native/'alevin/quants_mat.mtx', sparse.hstack([S, U, A]))
            (native/'alevin/quants_mat_cols.txt').write_text('\n'.join(genes)+'\n')
            (native/'quant.json').write_text(json.dumps({'usa_mode':True, 'num_genes':9}))
            before = {p.relative_to(native): p.read_bytes() for p in native.rglob('*') if p.is_file()}
            provenance = directory/'provenance.json'
            provenance.write_text('{}')
            output = directory/'WTX_velocity.h5ad'
            export(native, output, 'WTX', directory/'qc', PANELS, provenance)
            adata = ad.read_h5ad(output)
            self.assertEqual(adata.shape, (2, 3))
            self.assertEqual(adata.obs_names.tolist(), ['1', str(384**3)])
            self.assertEqual(adata.var_names.tolist(), genes)
            self.assertEqual((set(adata.layers) - {None}), {'spliced', 'unspliced'})
            for state, expected in (('spliced', S+A), ('unspliced', U)):
                self.assertTrue(equal_sparse(adata.layers[state], expected))
            self.assertEqual(validate_outputs(output, native, PANELS)['cross_format'], 'PASS')
            self.assertEqual(before, {p.relative_to(native):p.read_bytes() for p in native.rglob('*') if p.is_file()})
            self.assertFalse((directory/'raw').exists())
            self.assertFalse((directory/'usa_counts').exists())
            labels = native/'alevin/quants_mat_rows.txt'
            labels.write_text('\n'.join(reversed(labels.read_text().splitlines()))+'\n')
            with self.assertRaisesRegex(ValueError, 'identifier order'):
                validate_outputs(output, native, PANELS)

    def test_alevin_fry_labelled_usa_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            expected = fixture(directory/'usa')
            genes = ['ENSG000001.1', 'ENSG000002.9']
            columns = genes + [g+'-U' for g in genes] + [g+'-A' for g in genes]
            (directory/'usa/alevin/quants_mat_cols.txt').write_text('\n'.join(columns)+'\n')
            _, actual_genes, matrices = read_usa(directory/'usa')
            self.assertEqual(actual_genes, genes)
            for actual, target in zip(matrices, expected):
                self.assertTrue(equal_sparse(actual, target))

    def test_bad_usa_metadata_counts_and_gene_order(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            fixture(d/'usa')
            metadata_path = d/'usa/quant.json'
            metadata = json.loads(metadata_path.read_text())
            metadata.pop('num_genes')
            metadata_path.write_text(json.dumps(metadata))
            with self.assertRaisesRegex(ValueError, 'USA num_genes'):
                read_usa(d/'usa')
            metadata['num_genes'] = 6
            metadata_path.write_text(json.dumps(metadata))
            cols = d/'usa/alevin/quants_mat_cols.txt'
            cols.write_text('ENSG000001.1\nENSG000001.1\n')
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                read_usa(d/'usa')
            cols.write_text('ENSG1\nENSG2\nENSG2\nENSG1\nENSG1\nENSG2\n')
            with self.assertRaisesRegex(ValueError, 'gene order'):
                read_usa(d/'usa')
            cols.write_text('ENSG1\nENSG2\n')
            bad = sparse.csr_matrix([[1, -1, 0, 0, 0, 0],[0,0,0,0,0,0]], dtype=float)
            for value in (-1, np.nan):
                bad.data[-1] = value
                mmwrite(d/'usa/alevin/quants_mat.mtx', bad)
                with self.assertRaisesRegex(ValueError, 'finite and nonnegative'):
                    read_usa(d/'usa')


class Interface(unittest.TestCase):
    def test_copied_index_dry_run_without_optional_gene_symbols(self):
        binary = os.environ.get('VELOCITY_SNAKEMAKE')
        if not binary:
            self.skipTest('Set VELOCITY_SNAKEMAKE for a copied-index DAG regression check')
        from make_dryrun_fixture import create
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            config_path = create(directory)
            config = json.loads(config_path.read_text())
            copied = directory/'reference/copied_spliceu'
            (directory/'reference/index').rename(copied)
            (copied/'gene_id_to_name.tsv').unlink()
            config['simpleaf']['index'] = str(copied)
            config_path.write_text(json.dumps(config))
            self.assertEqual(quant.index_directory(copied), copied)
            result = subprocess.run([binary, '-s', str(ROOT/'workflow/Snakefile'),
                                     '--configfile', str(config_path), '--cores', '2',
                                     '--dry-run', '--printshellcmds'], cwd=ROOT,
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertNotIn('--gene-map ', result.stdout)
            self.assertIn(str(copied/'t2g_3col.tsv'), result.stdout)

    def test_cli_uses_velocity_config_and_existing_containers(self):
        from rhapsodist import cli
        with patch.object(sys,'argv',['rhapsodist','--cores','2','-n']), patch.object(cli.subprocess,'run',return_value=subprocess.CompletedProcess([],0)) as call:
            with self.assertRaises(SystemExit) as result:
                cli.main()
            self.assertEqual(result.exception.code,0)
            command=call.call_args.args[0]
        self.assertNotIn('--use-apptainer',command)
        self.assertNotIn('--use-conda',command)
        self.assertIn('configs/velocity.yaml',command)
        self.assertIn('--dry-run',command)

    def test_installed_cli_geometry_and_flags(self):
        binary = os.environ.get('VELOCITY_SIMPLEAF')
        if not binary:
            self.skipTest('Set VELOCITY_SIMPLEAF to inspect the installed executable')
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d/'chemistries.json').write_text('{}')
            env = dict(os.environ, ALEVIN_FRY_HOME=str(d))
            help_text = subprocess.check_output([binary,'quant','--help'],env=env,text=True)
            for flag in ('--unfiltered-pl','--min-reads','--chemistry','--t2g-map'):
                self.assertIn(flag,help_text)
            subprocess.run([binary,'chemistry','add','--name','bd_standardized','--geometry',
                            quant.GEOMETRY,'--expected-ori','both'],env=env,check=True,
                           stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            registered=json.loads((d/'chemistries.json').read_text())
            self.assertEqual(registered['bd_standardized']['geometry'],quant.GEOMETRY)

    def test_command_and_reference_layout(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            index = d/'index'
            index.mkdir()
            (index/'simpleaf_index.json').write_text('{"index_type":"piscem"}')
            (index/'t2g_3col.tsv').write_text('tx1\tENSG1\tS\ntx2\tENSG1\tU\n')
            for suffix in ('ctab', 'refinfo', 'ssi'):
                (index/f'piscem_idx.{suffix}').touch()
            resolved, t2g = quant.reference_files(d)
            command = quant.quant_command(resolved, t2g, 'r1', 'r2', 'wl', 'out', 16)
            self.assertEqual(command[command.index('--chemistry')+1], '1{b[27]u[8]}2{r:}')
            self.assertEqual(command[command.index('--resolution')+1], 'cr-like')
            self.assertEqual(command[command.index('--min-reads')+1], '1')
            self.assertEqual(command[command.index('--t2g-map')+1], str(t2g))
            self.assertNotIn('usa', command)
            self.assertIn('--unfiltered-pl', command)
            with self.assertRaises(ValueError):
                quant.quant_command(index,t2g,'r1','r2','wl','out',1)

    def test_quant_driver_atomic_publication_with_mocked_external_tools(self):
        # Mock all external processes: never invoke quantification or a container.
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            index=d/'reference/index'; index.mkdir(parents=True)
            (index/'simpleaf_index.json').write_text('{"index_type":"piscem"}')
            (index/'t2g_3col.tsv').write_text('tx1\tENSG1\tS\ntx2\tENSG1\tU\n')
            for suffix in ('ctab','ssi','refinfo'):
                (index/f'piscem_idx.{suffix}').touch()
            for name in ('image.sif','R1.gz','R2.gz','wl.txt'):
                (d/name).touch()
            provenance=d/'source.json'
            provenance.write_text(json.dumps({'upstream_commit':'fixture','sha256':{}}))
            def capture(command, **kwargs):
                if command[-1]=='--help':
                    return '--unfiltered-pl --min-reads --chemistry --t2g-map'
                return command[-2]+' '+quant.EXPECTED[command[-2]]
            def run(command, **kwargs):
                if 'quant' in command:
                    output=Path(command[command.index('--output')+1])
                    fixture(output/'af_quant')
                return subprocess.CompletedProcess(command,0)
            with patch.object(quant.subprocess,'check_output',side_effect=capture), patch.object(quant.subprocess,'run',side_effect=run) as mocked:
                quant.run_quant(index,d/'R1.gz',d/'R2.gz',d/'wl.txt',
                                d/'sample/simpleaf',2,d/'quant.log',provenance)
                self.assertEqual(mocked.call_count,2)
            self.assertTrue((d/'sample/simpleaf/af_quant/quant.json').is_file())
            metadata=json.loads((d/'sample/simpleaf/pipeline_provenance.json').read_text())
            self.assertEqual(metadata['quantification_resolution'],'USA')
            self.assertEqual(metadata['genome'],'unknown')
            self.assertEqual(metadata['GENCODE_release'],'unknown')
            self.assertFalse(list((d/'sample').glob('.simpleaf-*')))
            with self.assertRaises(FileExistsError):
                quant.run_quant(index,d/'R1.gz',d/'R2.gz',d/'wl.txt',
                                d/'sample/simpleaf',2,d/'quant.log',provenance)
            with patch.object(quant.subprocess,'check_output',side_effect=capture), patch.object(quant.subprocess,'run',side_effect=subprocess.CalledProcessError(1,'mock')):
                with self.assertRaises(subprocess.CalledProcessError):
                    quant.run_quant(index,d/'R1.gz',d/'R2.gz',d/'wl.txt',
                                    d/'failed/simpleaf',2,d/'failed.log',provenance)
            self.assertFalse((d/'failed/simpleaf').exists())
            self.assertFalse(list((d/'failed').glob('.simpleaf-*')))
            with patch.object(quant.subprocess,'check_output',side_effect=capture), patch.object(quant.subprocess,'run',side_effect=run), patch.object(
                    quant,'execution_provenance',side_effect=[{'snapshot':'before'},{'snapshot':'after'}]):
                with self.assertRaisesRegex(ValueError,'changed during quantification'):
                    quant.run_quant(index,d/'R1.gz',d/'R2.gz',d/'wl.txt',
                                    d/'changed/simpleaf',2,d/'changed.log',provenance)
            self.assertFalse((d/'changed/simpleaf').exists())
            self.assertFalse(list((d/'changed').glob('.simpleaf-*')))

    def test_pair_and_geometry_guards(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            for mate in ('r1','r2'):
                (d/mate).write_text('@read\n'+'A'*35+'\n+\n'+'I'*35+'\n')
            self.assertEqual(quant.paired_count(d/'r1',d/'r2',True),1)
            (d/'r2').write_text('@other\nACGT\n+\nIIII\n')
            with self.assertRaisesRegex(ValueError,'out of order'):
                quant.paired_count(d/'r1',d/'r2')
            (d/'r2').write_text('')
            with self.assertRaisesRegex(ValueError,'mate counts'):
                quant.paired_count(d/'r1',d/'r2')
            (d/'r1').write_text('@r\nACGT\n+\nIIII\n')
            shutil.copy(d/'r1',d/'r2')
            with self.assertRaisesRegex(ValueError,'exactly 27'):
                quant.paired_count(d/'r1',d/'r2',True)


class Preparation(unittest.TestCase):
    def test_retained_rules_and_dependencies_are_verbatim(self):
        upstream = original()
        retained = (ROOT/'workflow/src/prepare_fastq.snmk').read_text()
        for name in ('concat_fastqs','downsample_fastqs'):
            actual = re.sub(r'    resources:\n        mem_mb=.*\n', '', rule_block(retained, name))
            actual = re.sub(r',\n        reference=REFERENCE_CHECK\n', '\n', actual)
            self.assertEqual(re.sub(r'    conda:\n        [^\n]+\n', '', rule_block(upstream,name)), actual)
        for file in [f'workflow/data/whitelist_384x3/BD_CLS{i}.txt' for i in (1,2,3)]:
            self.assertEqual((ROOT/file).read_bytes(),subprocess.check_output(['git','show',f'7991f870fdfe38e12648673d330cbfd4839e755d:{file}'],cwd=ROOT))

    def test_required_helpers_match_upstream_for_active_fastq_config(self):
        upstream = {}
        helper_source = subprocess.check_output(
            ['git', 'show', '7991f870fdfe38e12648673d330cbfd4839e755d:workflow/src/workflow_functions.py'],
            cwd=ROOT, text=True)
        exec(compile(helper_source, '<upstream workflow helpers>', 'exec'), upstream)
        config = {'working_dir': '/tmp/work', 'downsample': 100,
                  'samples': [{'name': 'pilot', 'uses': {
                      'cb_umi_fq': ['lane1_R1.gz', 'lane2_R1.gz'],
                      'cdna_fq': ['lane1_R2.gz', 'lane2_R2.gz'],
                      'allowedlist': 384, 'diversity_insets': 'yes', 'downsample': 50}}]}
        upstream.update(config=config, op=os.path)
        required_helpers.config = config
        required_helpers.op = os.path
        for name, args in (
            ('get_sample_names', ()), ('get_cbumi_inputs', ('pilot',)),
            ('get_cdna_inputs', ('pilot',)), ('_raw_cbumi_path', ('pilot',)),
            ('_raw_cdna_path', ('pilot',)), ('get_downsample_fraction', ('pilot',)),
            ('get_cbumi_by_name', ('pilot',)), ('get_cdna_by_name', ('pilot',)),
            ('get_declared_allowedlist', ('pilot',)),
            ('get_declared_diversity_insets', ('pilot',)),
            ('get_bead_type_by_name', ('pilot',)),
            ('get_barcode_whitelist_by_name', ('pilot',))):
            self.assertEqual(upstream[name](*args), getattr(required_helpers, name)(*args), name)
        upstream['validate_fastq_lists']()
        required_helpers.validate_fastq_lists()
        with tempfile.TemporaryDirectory() as directory:
            reads = Path(directory)/'linkers.fastq.gz'
            with gzip.open(reads, 'wt') as handle:
                for offset in range(4):
                    sequence = 'A'*(9+offset) + 'GTGA' + 'C'*9 + 'GACA' + 'T'*20
                    handle.write(f'@r{offset}\n{sequence}\n+\n' + 'I'*len(sequence) + '\n')
            self.assertEqual(upstream['scan_r1_linkers'](reads),
                             required_helpers.scan_r1_linkers(reads))
            second=Path(directory)/'second.gz'
            second.write_bytes(reads.read_bytes())
            combined=Path(directory)/'combined.gz'
            combined.write_bytes(reads.read_bytes()+second.read_bytes())
            for limit in (2,5,10000):
                self.assertEqual(upstream['scan_r1_linkers'](combined,limit),
                                 required_helpers.scan_r1_linkers([reads,second],limit))


    def test_cutadapt_fixture_matches_original_and_expected_streams(self):
        python = os.environ.get('VELOCITY_CUTADAPT_PYTHON')
        binary_dir = os.environ.get('VELOCITY_PREPROCESS_BIN')
        if not python or not binary_dir:
            self.skipTest('Set VELOCITY_CUTADAPT_PYTHON and VELOCITY_PREPROCESS_BIN for Cutadapt 4.9 fixture')
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            tools = d/'bin'; tools.mkdir()
            launcher = tools/'cutadapt'
            launcher.write_text('#!/bin/sh\nexec '+python+' -m cutadapt "$@"\n')
            launcher.chmod(0o755)
            env = dict(os.environ, PATH=str(tools)+':'+binary_dir+':'+os.environ['PATH'])
            keys = bd.load_keys(PANELS)
            expected=[]; r1=[]; r2=[]
            for i, diversity in enumerate(('', 'A', 'GT', 'TCA')):
                bc = ''.join(k[i] for k in keys)
                umi = 'ACGTACGT'
                seq = diversity+bc[:9]+'GTGA'+bc[9:18]+'GACA'+bc[18:]+umi+'T'*25
                name = f'@r{i}' if i < 3 else '@A00579:123:HT6LFDSX7:1:1101:1000:1000'
                r1.append(f'{name}\n{seq}\n+\n'+ 'I'*len(seq)+'\n')
                r2.append(f'{name}\nACGTTGCA\n+\nIIIIIIII\n')
                expected.append((name[:9]+name[13:22]+name[26:],bc+umi))
            seq='A'*43
            r1.append('@bad\n'+seq+'\n+\n'+'I'*43+'\n')
            r2.append('@bad\nACGTTGCA\n+\nIIIIIIII\n')
            for name, reads in [('r1',r1),('r2',r2)]:
                with gzip.open(d/(name+'.fq.gz'),'wt') as handle:
                    handle.write(''.join(reads))
            streams=[]
            for label, source in [('original',original()),('retained',(ROOT/'workflow/src/prepare_fastq.snmk').read_text())]:
                out=d/label; out.mkdir()
                if label == 'original':
                    body=re.search('    shell:\n        """(.*?)        """',rule_block(source,'standardize_cb_umis_cutadapt'),re.S).group(1)
                    outputs=SimpleNamespace(temp_cb_umi=str(out/'temp.gz'),cdna=str(out/'r2.gz'),standardized_cb_umi=str(out/'r1.gz'))
                    command=body.format(params=SimpleNamespace(path=out,bead_type='enhanced_v2',max_errors=0),threads=2,
                                        input=SimpleNamespace(cb_umi=str(d/'r1.fq.gz'),cdna=str(d/'r2.fq.gz')),
                                        output=outputs,log=str(out/'cutadapt.log'))
                    subprocess.run(command,shell=True,executable='/bin/bash',env=env,check=True)
                else:
                    with patch.dict(os.environ,env):
                        preparation.standardize(d/'r1.fq.gz',d/'r2.fq.gz',out/'unit1.gz',out/'unit2.gz',
                                                out/'unit.json',out/'cutadapt.log','fixture','0000',2)
                    preparation.combine_prepared([out/'unit1.gz'],[out/'unit2.gz'],[out/'unit.json'],
                        [d/'r1.fq.gz'],[d/'r2.fq.gz'],out/'r1.gz',out/'r2.gz',out/'cutadapt.json','fixture',['0000'])
                    self.assertFalse(list(out.glob('.prepare-*')))
                    self.assertFalse((out/'temp.gz').exists())
                self.assertEqual(list(quant.records(out/'r1.gz')),expected)
                self.assertEqual(quant.paired_count(out/'r1.gz',out/'r2.gz',True),4)
                streams.append([gzip.open(out/m,'rt').read() for m in ('r1.gz','r2.gz')])
            self.assertEqual(streams[0],streams[1])
            # Execute original whitelist body, including its exact validation/filtering.
            block=rule_block(original(),'derive_observed_whitelist')
            body=block.split('    run:\n',1)[1]
            import textwrap
            wl=d/'wl.txt'
            exec(textwrap.dedent(body),dict(input=SimpleNamespace(cb1=str(PANELS/'BD_CLS1.txt'),
                 cb2=str(PANELS/'BD_CLS2.txt'),cb3=str(PANELS/'BD_CLS3.txt'),
                 cb_umi_fq=str(d/'retained/r1.gz')),output=SimpleNamespace(observed_whitelist=str(wl)),
                 log=[str(d/'whitelist.log')]))
            self.assertEqual(wl.read_text().splitlines(), sorted(seq[:27] for _,seq in expected))
            quant.fastq_qc([d/'r1.fq.gz'],[d/'r2.fq.gz'],d/'retained/r1.gz',d/'retained/r2.gz',wl,d/'fastq_qc.tsv')
            self.assertIn('5\t4\t80.0\t4', (d/'fastq_qc.tsv').read_text())
            optimized_wl=d/'optimized_whitelist.txt'
            receipt=d/'receipt.json'
            panels=[PANELS/f'BD_CLS{i}.txt' for i in (1,2,3)]
            quant.prepare_qc([d/'r1.fq.gz'],[d/'r2.fq.gz'],d/'retained/r1.gz',d/'retained/r2.gz',
                             panels,optimized_wl,d/'retained/cutadapt.json',receipt,'fixture')
            self.assertEqual(wl.read_bytes(),optimized_wl.read_bytes())
            kwargs=dict(receipt=receipt,sample='fixture',panels=panels,cutadapt_json=d/'retained/cutadapt.json')
            args=([d/'r1.fq.gz'],[d/'r2.fq.gz'],d/'retained/r1.gz',d/'retained/r2.gz',optimized_wl,d/'optimized.tsv')
            with patch.object(quant,'records',side_effect=AssertionError('QC must not read FASTQs')):
                quant.fastq_qc(*args,**kwargs)
            self.assertEqual((d/'fastq_qc.tsv').read_bytes(),(d/'optimized.tsv').read_bytes())
            quant.fastq_qc(*args,**kwargs,exhaustive=True)
            with self.assertRaisesRegex(ValueError,'wrong-sample'):
                quant.fastq_qc(*args,**dict(kwargs,sample='other'))
            # The actual paired Cutadapt reader rejects malformed/mismatched raw input.
            for bad in ('@other\nACGT\n+\nIIII\n', '@r0\nACGT\n+\n'):
                (d/'bad_r2.fq').write_text(bad)
                result=subprocess.run([str(launcher), '-o', str(d/'bad1.fq'), '-p', str(d/'bad2.fq'),
                                       str(d/'r1.fq.gz'), str(d/'bad_r2.fq')],
                                      env=env,capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
            (d/'retained/r1.gz').write_bytes((d/'retained/r1.gz').read_bytes()[:-8])
            with self.assertRaises((EOFError,OSError)):
                quant.prepare_qc([d/'r1.fq.gz'],[d/'r2.fq.gz'],d/'retained/r1.gz',d/'retained/r2.gz',
                                 panels,optimized_wl,d/'retained/cutadapt.json',receipt,'fixture')
            optimized_wl.write_text('A'*27+'\n')
            with self.assertRaisesRegex(ValueError,'Stale'):
                quant.fastq_qc(*args,**kwargs)



class FastqQC(unittest.TestCase):
    def test_record_failures(self):
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary)
            good='@r/1\n'+'A'*35+'\n+\n'+'I'*35+'\n'
            (d/'r1').write_text(good)
            for bad, message in [('@r/2\nAAAA\n+\nIII\n','Malformed'),
                                 ('@r/2\nAAAA\n+\n','Malformed'),
                                 ('@other/2\nAAAA\n+\nIIII\n','out of order'),
                                 ('','mate counts')]:
                (d/'r2').write_text(bad)
                with self.assertRaisesRegex(ValueError,message):
                    quant.paired_count(d/'r1',d/'r2')
            with gzip.open(d/'truncated.gz','wt') as handle:
                handle.write(good)
            data=(d/'truncated.gz').read_bytes()
            (d/'truncated.gz').write_bytes(data[:-8])
            with self.assertRaises((EOFError,OSError)):
                list(quant.records(d/'truncated.gz'))
            (d/'r1').write_text('@r\n'+'Z'*35+'\n+\n'+'I'*35+'\n')
            (d/'r2').write_text('@r\nAAAA\n+\nIIII\n')
            with self.assertRaisesRegex(ValueError,'non-DNA'):
                quant.paired_count(d/'r1',d/'r2',True)

    def test_whitelist_filter_counts_and_stale_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary)
            panels=[PANELS/f'BD_CLS{i}.txt' for i in (1,2,3)]
            barcode=bd.encode(1,bd.load_keys(PANELS))
            # Ambiguous UMIs are legal; invalid barcode panels remain excluded, not filtered reads.
            sequences=[barcode+'NNNNNNNN', 'N'*27+'ACGTACGT', barcode+'ACGTACGT']
            for name in ('raw1','raw2','r1','r2'):
                (d/name).write_text(''.join('@r'+str(i)+'\n'+seq+'\n+\n'+'I'*len(seq)+'\n'
                                          for i,seq in enumerate(sequences)))
            stats=d/'cutadapt.json'
            stats.write_text(json.dumps(dict(cutadapt_version='4.9',input=dict(paired=True,
                path1=str(d/'raw1'),path2=str(d/'raw2')),command_line_arguments=['-p',str(d/'r2')],read_counts=dict(input=3,output=3))))
            args=([d/'raw1'],[d/'raw2'],d/'r1',d/'r2',panels,d/'wl',stats,d/'receipt','sample')
            with patch.object(quant,'records',wraps=quant.records) as scanned:
                quant.prepare_qc(*args)
                self.assertEqual([call.args[0] for call in scanned.call_args_list], [d/'r1'])
            self.assertEqual((d/'wl').read_text(),barcode+'\n')
            kwargs=dict(receipt=d/'receipt',sample='sample',panels=panels,cutadapt_json=stats)
            qcargs=([d/'raw1'],[d/'raw2'],d/'r1',d/'r2',d/'wl',d/'qc.tsv')
            quant.fastq_qc(*qcargs,**kwargs)
            quant.fastq_qc(*qcargs,**kwargs,exhaustive=True)
            self.assertIn('3\t3\t100.0\t1\tunavailable', (d/'qc.tsv').read_text())
            # An unchanged count cannot hide an ordinary same-size edit.
            (d/'r2').write_text((d/'r2').read_text().replace('@r0','@x0'))
            with self.assertRaisesRegex(ValueError,'Stale'):
                quant.fastq_qc(*qcargs,**kwargs)
            with self.assertRaisesRegex(ValueError,'out of order'):
                quant.fastq_qc(*qcargs)
            (d/'r2').write_text((d/'raw2').read_text())
            stats.write_text(stats.read_text().replace('"output": 3','"output": 2'))
            with self.assertRaisesRegex(ValueError,'disagree'):
                quant.prepare_qc(*args)
            stats.write_text(stats.read_text().replace('"output": 2','"output": 3'))
            # Multi-lane and downsampled configs keep raw checks and exact original denominator.
            for downsampled in (False,True):
                with patch.object(quant,'paired_count',wraps=quant.paired_count) as count:
                    quant.prepare_qc(*args,downsampled=downsampled)
                    self.assertEqual(count.call_count, 1 if downsampled else 0)
            quant.prepare_qc([d/'raw1']*2,[d/'raw2']*2,d/'r1',d/'r2',panels,d/'wl',stats,d/'receipt','sample',True)
            self.assertEqual(json.loads((d/'receipt').read_text())['input_read_pairs'],6)
            for bad in ('', barcode+'\n'+barcode+'\n', 'A'*26+'\n'):
                (d/'wl').write_text(bad)
                with self.assertRaisesRegex(ValueError,'whitelist'):
                    quant.fastq_qc(*qcargs)

class ExecutionProvenance(unittest.TestCase):
    def test_hashes_are_stable_sensitive_and_require_complete_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary)
            source=d/'source.py';source.write_text('print(1)\n')
            lock=d/'pixi.lock';lock.write_text('fixture lock\n')
            receipt=d/'receipt.json';receipt.write_text('{"sample":"fixture"}')
            with patch.object(quant,'WORKFLOW_ROOT',d),patch.object(quant,'WORKFLOW_SOURCES',[source]):
                before=quant.execution_provenance('a'*64,receipt)
                self.assertEqual(before,quant.execution_provenance('a'*64,receipt))
                for path,field in ((source,'workflow_source_sha256'),(lock,'workflow_lock_sha256'),
                                   (receipt,'preparation_qc_sha256')):
                    original=path.read_bytes()
                    path.write_bytes(original+b' ')
                    self.assertNotEqual(before[field],quant.execution_provenance('a'*64,receipt)[field])
                    path.write_bytes(original)
                self.assertNotEqual(before['workflow_config_sha256'],
                                    quant.execution_provenance('b'*64,receipt)['workflow_config_sha256'])
                self.assertEqual(quant.execution_provenance()['workflow_config_sha256'],'unavailable')
                with self.assertRaisesRegex(ValueError,'both'):
                    quant.execution_provenance('a'*64)
                with self.assertRaisesRegex(ValueError,'Invalid'):
                    quant.execution_provenance('not-a-sha',receipt)


class ReferenceIdentity(unittest.TestCase):
    def test_manifest_binding_unknown_identity_and_stale_preflight(self):
        from make_dryrun_fixture import create
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary)
            create(d)
            index=d/'reference/index'
            checked=quant.inspect_reference(index)
            self.assertEqual((checked['genome'],checked['GENCODE_release']),('unknown','unknown'))
            manifest={key:checked[key] for key in ('reference_t2g_sha256','reference_index_meta_sha256')}
            manifest.update(version=1,index_type='spliceu',genome='GRCh37',GENCODE_release='v19')
            path=index/'reference_manifest.json'
            path.write_text(json.dumps(manifest))
            checked=quant.inspect_reference(index)
            self.assertEqual((checked['genome'],checked['GENCODE_release']),('GRCh37','v19'))
            self.assertIn('hashes matched',checked['reference_identity_status'])
            preflight=d/'checked.json'
            preflight.write_text(json.dumps(checked))
            self.assertEqual(quant.load_reference_validation(index,preflight),checked)
            (index/'piscem_idx.ctab').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'Stale'):
                quant.load_reference_validation(index,preflight)
            for key in ('reference_t2g_sha256','reference_index_meta_sha256'):
                bad=dict(manifest,**{key:'0'*64})
                path.write_text(json.dumps(bad))
                with self.assertRaisesRegex(ValueError,'hash mismatch'):
                    quant.inspect_reference(index)
            path.write_text('{}')
            with self.assertRaisesRegex(ValueError,'Invalid reference'):
                quant.inspect_reference(index)
            path.unlink()
            (index/'t2g_3col.tsv').write_text('tx\tENSG1\tS\n')
            with self.assertRaisesRegex(ValueError,'both S and U'):
                quant.inspect_reference(index)


class PreparationUpgrades(unittest.TestCase):
    def setUp(self):
        self.binary_dir=os.environ.get('VELOCITY_PREPROCESS_BIN')
        if not self.binary_dir:
            self.skipTest('Set VELOCITY_PREPROCESS_BIN for actual streaming preparation tests')
        self.env=dict(os.environ,PATH=self.binary_dir+':'+os.environ['PATH'],
                      XDG_CACHE_HOME=str(ROOT/'.cache'))
        self.barcode=bd.encode(1,bd.load_keys(PANELS))

    def write_pair(self,d,prefix,mismatch=False):
        paths=[d/(prefix+f'_R{mate}.fq.gz') for mate in (1,2)]
        sequence=self.barcode[:9]+'GTGA'+self.barcode[9:18]+'GACA'+self.barcode[18:]+'ACGTACGT'+'T'*25
        for mate,path in enumerate(paths):
            name='@other' if mismatch and mate==1 else '@shared_read_name'
            seq=sequence if mate==0 else 'ACGTTGCA'
            with gzip.open(path,'wt') as handle:
                handle.write(name+'\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
        return paths

    def test_lane_counts_order_whitelist_and_pipeline_failures(self):
        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ,self.env):
            d=Path(temporary)
            raw=[self.write_pair(d,f'lane{i}') for i in range(2)]
            outputs=[]
            for i,(left,right) in enumerate(raw):
                out=[d/f'unit{i}_R1.gz',d/f'unit{i}_R2.gz',d/f'unit{i}.json']
                preparation.standardize(left,right,*out,d/f'unit{i}.log','fixture',f'{i:04d}',2)
                outputs.append(out)
            final=[d/'final_R1.gz',d/'final_R2.gz',d/'combined.json']
            args=([o[0] for o in outputs],[o[1] for o in outputs],[o[2] for o in outputs],
                  [r[0] for r in raw],[r[1] for r in raw],*final,'fixture',['0000','0001'])
            preparation.combine_prepared(*args)
            self.assertEqual(quant.paired_count(final[0],final[1],True),2)
            for mate in (0,1):
                expected=''.join(gzip.open(o[mate],'rt').read() for o in outputs)
                self.assertEqual(gzip.open(final[mate],'rt').read(),expected)
            panels=[PANELS/f'BD_CLS{i}.txt' for i in (1,2,3)]
            with patch.object(quant,'paired_count',side_effect=AssertionError('No repeated raw lane scan')):
                quant.prepare_qc(args[3],args[4],final[0],final[1],panels,d/'wl',final[2],d/'receipt','fixture')
            quant.fastq_qc(args[3],args[4],final[0],final[1],d/'wl',d/'qc',receipt=d/'receipt',
                           sample='fixture',panels=panels,cutadapt_json=final[2],exhaustive=True)
            self.assertIn('2\t2\t100.0\t1', (d/'qc').read_text())
            with self.assertRaisesRegex(ValueError,'reordered'):
                preparation.combine_prepared(*args[:-1],['0001','0000'])
            bad=self.write_pair(d,'bad',True)
            failed=[d/'failed1.gz',d/'failed2.gz',d/'failed.json']
            with self.assertRaises(subprocess.CalledProcessError):
                preparation.standardize(*bad,*failed,d/'failure.log','fixture','bad',2)
            self.assertTrue(all(not p.exists() for p in failed))
            # Cutadapt alone accepts these mismatches; preserve legacy strict checks in workers.
            for left_name,right_name in (('read1','read2'),('read/2','read/1'),('read.1','read.2')):
                for path,name in zip(raw[0],(left_name,right_name)):
                    text=gzip.open(path,'rt').read()
                    with gzip.open(path,'wt') as handle:
                        handle.write('@'+name+'\n'+text.split('\n',1)[1])
                with self.assertRaises(subprocess.CalledProcessError):
                    preparation.standardize(*raw[0],*failed,d/'failure.log','fixture','bad',4)
                self.assertIn('out of order',(d/'failure.log').read_text())
            self.write_pair(d,'lane0')
            # Fail downstream pipeline processes separately; pipefail must reject both.
            tools=d/'tools';tools.mkdir()
            for executable in ('cut','pigz'):
                fake=tools/executable
                fake.write_text('#!/bin/sh\nexit 7\n');fake.chmod(0o755)
                with patch.dict(os.environ,PATH=str(tools)+':'+self.env['PATH']):
                    with self.assertRaises(subprocess.CalledProcessError):
                        preparation.standardize(*raw[0],*failed,d/'failure.log','fixture','bad',2)
                fake.unlink()
                self.assertTrue(all(not p.exists() for p in failed))
            self.assertFalse(list(d.glob('.prepare-*')))

    def test_early_gates_and_failed_lane_resume(self):
        binary=os.environ.get('VELOCITY_SNAKEMAKE')
        if not binary:
            self.skipTest('Set VELOCITY_SNAKEMAKE for early-gate and failed-lane recovery checks')
        from make_dryrun_fixture import create
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary)
            config_path=create(d)
            config=json.loads(config_path.read_text())
            raw=[self.write_pair(d,'lane0'),self.write_pair(d,'lane1',True)]
            config['samples']=[dict(name='fixture',uses=dict(cb_umi_fq=[str(r[0]) for r in raw],
                              cdna_fq=[str(r[1]) for r in raw],allowedlist=384,diversity_insets='yes'))]
            config_path.write_text(json.dumps(config))
            work=Path(config['working_dir'])
            def run(target,success=True,rerun=False):
                command=[binary,'-s',str(ROOT/'workflow/Snakefile'),'--configfile',str(config_path),
                         '--cores','2',str(target)]
                if rerun:command.append('--rerun-incomplete')
                result=subprocess.run(command,cwd=d,env=self.env,text=True,capture_output=True)
                if success:self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                else:self.assertNotEqual(result.returncode,0,result.stdout+result.stderr)
                return result
            target=work/'data/fastq/fixture_standardized_cb_umi.fq.gz'
            # Reference content fails before any FASTQ processing despite existing index files.
            t2g=d/'reference/index/t2g_3col.tsv';valid=t2g.read_text()
            t2g.write_text('tx\tENSG1\tS\n')
            run(target,False)
            self.assertFalse(list((work/'logs').glob('*unit*.log')))
            t2g.write_text(valid)
            # Non-enhanced/unknown input must fail before Cutadapt, with diagnostics retained in log.
            good=raw[0][0].read_bytes()
            for sequence in ('A'*70,'A'*9+'ACTGGCCTGCGA'+'C'*9+'GGTAGCGGTGACA'+'T'*30):
                with gzip.open(raw[0][0],'wt') as handle:
                    handle.write('@r\n'+sequence+'\n+\n'+'I'*len(sequence)+'\n')
                config['samples'][0]['uses']['cb_umi_fq']=[str(raw[0][0])]
                config['samples'][0]['uses']['cdna_fq']=[str(raw[0][1])]
                config_path.write_text(json.dumps(config))
                run(target,False,True)
                self.assertFalse(list((work/'logs').glob('*unit*.log')))
                self.assertIn('detected bead linker class', (work/'logs/fixture_bead_version_validation.log').read_text())
            raw[0][0].write_bytes(good)
            config['samples'][0]['uses']['cb_umi_fq']=[str(r[0]) for r in raw]
            config['samples'][0]['uses']['cdna_fq']=[str(r[1]) for r in raw]
            config_path.write_text(json.dumps(config))
            good_unit=work/'data/fastq/units/fixture/0000_R1.fq.gz'
            first=run(good_unit,True,True)
            self.assertFalse(target.exists(),first.stdout+first.stderr)
            completed=work/'logs/fixture_unit_0000.json'
            before=completed.stat().st_mtime_ns
            run(target,False)
            self.assertEqual(completed.stat().st_mtime_ns,before)
            raw[1][1].write_bytes(raw[0][1].read_bytes())
            run(target,True,True)
            self.assertEqual(completed.stat().st_mtime_ns,before)
            self.assertEqual(quant.paired_count(target,work/'data/fastq/fixture_standardized_cdna.fq.gz',True),2)
            run(target)
            self.assertEqual(completed.stat().st_mtime_ns,before)
            self.assertFalse((work/'data/fastq/combined').exists())
            self.assertFalse(list((work/'data/fastq/units').rglob('*.gz')))


if __name__ == '__main__':
    unittest.main(verbosity=2)
