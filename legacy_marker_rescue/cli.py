"""Command-line interface for marker-image and genome analyses."""
from __future__ import annotations
import argparse
import datetime as dt
import os
from pathlib import Path
import sys
import traceback
from .io import read_json, write_json, write_manifest, resolve_path, digest

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description='Legacy molecular-marker analysis')
    parser.add_argument('stage', choices=['all', 'genome', 'context', 'detection', 'external', 'selftest'])
    parser.add_argument('--inputs', type=Path, help='Local input JSON; relative paths are resolved against this file')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--config-dir', type=Path, default=ROOT / 'config')
    args = parser.parse_args(argv)
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ.setdefault(name, '1')
    if args.stage == 'selftest':
        import unittest
        suite = unittest.defaultTestLoader.discover(str(ROOT / 'tests'), top_level_dir=str(ROOT))
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        return 0 if result.wasSuccessful() and result.testsRun > 0 and not result.skipped else 1
    if args.output_dir is None: parser.error('--output-dir is required for analysis')
    out = args.output_dir.resolve()
    if out == ROOT or ROOT in out.parents:
        parser.error('Place analysis outputs outside the source repository')
    if out.exists() and any(out.iterdir()): parser.error('Use an empty output directory')
    out.mkdir(parents=True, exist_ok=True)
    conf = args.config_dir.resolve()
    inputs = read_json(args.inputs) if args.inputs else {}
    base = args.inputs.resolve().parent if args.inputs else ROOT
    def path(key, default=None):
        value = inputs.get(key, default)
        if value is None: raise ValueError('Missing local input: ' + key)
        return resolve_path(str(value), base)
    status = {'status': 'RUNNING', 'stages': {}, 'started_utc': dt.datetime.now(dt.timezone.utc).isoformat()}
    write_json(out / 'STATUS.json', status)
    code = 0
    try:
        import numpy, scipy, PIL
        write_json(out / 'software_versions.json', {'python': sys.version, 'numpy': numpy.__version__, 'scipy': scipy.__version__, 'pillow': PIL.__version__})
        config_files = {p.name: digest(p) for p in sorted(conf.glob('*.json'))}
        source_files = {p.relative_to(ROOT).as_posix(): digest(p) for p in sorted((ROOT / 'legacy_marker_rescue').rglob('*.py'))}
        write_json(out / 'analysis_identity.json', {'source_sha256': source_files, 'configuration_sha256': config_files})
        if args.stage in {'all', 'genome'}:
            from .genome.pipeline import run
            status['stages']['genome'] = run(path('genome_manifest'), path('genome_data', ROOT / 'data/genome'),
                      read_json(conf / 'genome.json'), out / 'genome', path('scan_cache', out.parent / 'scan_cache'),
                      path('mash_table') if inputs.get('mash_table') else None)
        if args.stage in {'all', 'context'}:
            from .context.pipeline import run
            source = out / 'genome' if args.stage == 'all' else path('genome_result')
            status['stages']['context'] = run(source, path('genome_manifest'), read_json(conf / 'context.json'), out / 'context')
        if args.stage in {'all', 'detection'}:
            from .gel.pipeline import run
            status['stages']['detection'] = run(path('denoyes_data', ROOT / 'data/denoyes'), conf, out / 'detection')
        if args.stage in {'all', 'external'}:
            from .benchmark.analysis import run
            status['stages']['external'] = run(read_json(conf / 'benchmark.json'), read_json(conf / 'detection.json'),
                                               path('external_cache', out.parent / 'external_cache'), out / 'external')
        status['status'] = 'COMPLETE'
    except BaseException as error:
        if isinstance(error, KeyboardInterrupt): status['status'] = 'INTERRUPTED'; code = 130
        else: status['status'] = 'FAILED'; code = 1
        status['error'] = type(error).__name__ + ': ' + str(error)
        (out / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        print(status['error'], file=sys.stderr, flush=True)
    finally:
        status['ended_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
        write_json(out / 'STATUS.json', status)
        write_manifest(out)
        print('PIPELINE_STATUS: ' + status['status'], flush=True)
    return code
