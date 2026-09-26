#!/usr/bin/env python3
"""Compute statistical analyses and render their outputs from authorized cohort inputs."""
import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def run(*command):
    subprocess.run([str(x) for x in command], check=True)

def validate_outputs(output):
    with (ROOT / 'config/output_contract.csv').open(newline='') as f:
        expected = list(csv.DictReader(f))
    missing = [row['path'] for row in expected if not (output/row['path']).is_file()]
    if missing:
        raise RuntimeError('Incomplete analysis output: ' + ', '.join(missing))
    def rows(name):
        with (output/'tables'/name).open() as f:
            return list(csv.DictReader(f))
    for record in rows('resampling_validity_registry.csv'):
        try:
            attempted, valid, invalid = (int(record[k]) for k in ('attempted','valid','invalid'))
        except (ValueError,KeyError) as e:
            raise RuntimeError('Invalid resampling ledger schema/counts') from e
        if attempted != valid + invalid:
            raise RuntimeError('Resampling ledger does not reconcile')
    if len(rows('table_S5_component_landscape_associations.csv')) != 288:
        raise RuntimeError('Incomplete component comparison family')
    if len(rows('global_1d_landscape_tests.csv')) != 6 or len(rows('global_2d_surface_tests.csv')) != 6:
        raise RuntimeError('Incomplete primary comparison family')


def render(output, anatomy_assets=None):
    if not (output / 'tables/component_selected_statistic_key.csv').is_file():
        raise FileNotFoundError('Statistical stage must finish before rendering: ' + str(output))
    run('Rscript', ROOT / 'R/render_figures.R', '--input-dir', output,
        '--output-dir', output / 'figures')
    if anatomy_assets is not None:
        run(sys.executable, ROOT / 'python/yeo_overlap.py', '--source-dir', anatomy_assets / 'yeo',
            '--output-dir', output / 'figure_data/yeo')
        run(sys.executable, ROOT / 'python/aal3_overlap.py', '--source-dir', anatomy_assets / 'aal3',
            '--output-dir', output)
        run(sys.executable, ROOT / 'python/figure3.py', '--panel-dir', anatomy_assets / 'source_panels',
            '--yeo-csv', output / 'figure_data/yeo/yeo_radar_display.csv', '--output-dir', output / 'figures')
    else:
        print('Anatomical figures not requested: supply --anatomy-assets with authorized external atlas/map assets.')


def external_output_path(path):
    output = Path(path).resolve()
    if output == ROOT or ROOT in output.parents:
        raise ValueError('Output must be outside the code checkout')
    return output


def file_hashes(folder, paths):
    return {str(p.relative_to(folder)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths) if p.is_file()}


def write_manifest(output, manifest):
    # Atomic metadata replacement; a failed render never claims success.
    path = output / 'run_manifest.json'
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(manifest, indent=2) + '\n')
    temporary.replace(path)


def render_with_provenance(output, anatomy_assets=None):
    output = external_output_path(output)
    manifest = json.loads((output / 'run_manifest.json').read_text())
    if manifest.get('mode') not in ('synthetic', 'full'):
        raise ValueError('Unrecognized computed-run manifest')
    hashes = source_hashes()
    assets = file_hashes(anatomy_assets, anatomy_assets.rglob('*')) if anatomy_assets else None
    manifest['render_status'] = 'RUNNING'
    write_manifest(output, manifest)
    try:
        # Retain an earlier successful anatomy render when only statistical
        # figures are requested again. Verify its files before retaining provenance.
        previous_anatomy = manifest.get('anatomy_figures') == 'GENERATED_FROM_EXTERNAL_ASSETS'
        if anatomy_assets is None and previous_anatomy:
            previous = manifest.get('anatomy_output_hashes', {})
            if not previous or file_hashes(output, [output / p for p in previous]) != previous:
                raise RuntimeError('Cannot verify previous anatomy outputs; supply --anatomy-assets to regenerate')
            provenance = manifest.get('anatomy_provenance')
            if provenance is not None:
                if (provenance.get('output_hashes') != previous or
                        provenance.get('asset_hashes') != manifest.get('external_anatomy_asset_hashes', {})):
                    raise RuntimeError('Anatomy provenance does not match retained assets/outputs; regenerate with --anatomy-assets')
            else:
                # Older manifests have only the most recent render record.
                # Recover its generator only if it actually rendered these
                # anatomy files; otherwise explicitly retain the unknown source.
                last = manifest.get('render_provenance', {})
                figure_outputs = {p: h for p, h in previous.items() if p.startswith('figures/')}
                recoverable = (last.get('anatomy_requested_this_render') is True and
                               bool(last.get('source_hashes')) and
                               bool(figure_outputs) and
                               all(last.get('figure_hashes', {}).get(p) == h
                                   for p, h in figure_outputs.items()))
                manifest['anatomy_provenance'] = {
                    'status': 'RECORDED' if recoverable else 'NOT_RECORDED',
                    'source_hashes': last.get('source_hashes') if recoverable else None,
                    'completed_utc': last.get('completed_utc') if recoverable else None,
                    'python': last.get('python') if recoverable else None,
                    'asset_hashes': manifest.get('external_anatomy_asset_hashes', {}),
                    'output_hashes': previous,
                }
        render(output, anatomy_assets)
        if hashes != source_hashes():
            raise RuntimeError('Source changed during rendering')
        if anatomy_assets:
            if assets != file_hashes(anatomy_assets, anatomy_assets.rglob('*')):
                raise RuntimeError('Anatomy assets changed during rendering')
            manifest['anatomy_figures'] = 'GENERATED_FROM_EXTERNAL_ASSETS'
            manifest['external_anatomy_asset_hashes'] = assets
            paths = list((output/'figure_data/yeo').rglob('*')) + list((output/'figure_data/aal3').rglob('*'))
            paths += [p for p in (output/'figures').glob('*')
                      if p.name.startswith(('Figure_3_components_anatomy', 'supplementary_figure3_aal3.'))]
            manifest['anatomy_output_hashes'] = file_hashes(output, paths)
            manifest['anatomy_provenance'] = {
                'status': 'RECORDED',
                'completed_utc': datetime.now(timezone.utc).isoformat(),
                'source_hashes': hashes, 'python': sys.version,
                'asset_hashes': assets,
                'output_hashes': manifest['anatomy_output_hashes'],
            }
        elif not previous_anatomy:
            manifest['anatomy_figures'] = 'NOT_REQUESTED'
        manifest['render_status'] = 'COMPLETE'
        manifest['render_provenance'] = {
            'completed_utc': datetime.now(timezone.utc).isoformat(),
            'source_hashes': hashes, 'python': sys.version,
            'anatomy_requested_this_render': anatomy_assets is not None,
            'figure_hashes': file_hashes(output, (output/'figures').rglob('*')),
        }
        # Existing statistical source/input provenance is deliberately retained.
        write_manifest(output, manifest)
    except Exception:
        manifest['render_status'] = 'FAILED'
        write_manifest(output, manifest)
        raise


def source_hashes():
    files = [x for folder in ('R','python','config') for x in (ROOT/folder).rglob('*')
             if x.is_file() and '__pycache__' not in x.parts]
    return {str(x.relative_to(ROOT)): hashlib.sha256(x.read_bytes()).hexdigest() for x in files}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='mode', required=True)
    for name in ('check', 'run', 'synthetic'):
        parser = sub.add_parser(name)
        parser.add_argument('--config', type=Path, default=ROOT/'config/study.yml')
        if name != 'synthetic':
            parser.add_argument('--input', type=Path, required=True)
            parser.add_argument('--correction', type=Path)
        if name != 'check':
            parser.add_argument('--output', type=Path, required=True)
            parser.add_argument('--cores', type=int, default=1)
            parser.add_argument('--component-map-dir', type=Path)
            parser.add_argument('--anatomy-assets', type=Path, help='Optional authorized external source_panels, yeo and aal3 directories')
        if name == 'run':
            parser.add_argument('--authorize-private-analysis', action='store_true')
    again = sub.add_parser('render', help='Render computed run outputs without inference')
    again.add_argument('--output', type=Path, required=True)
    again.add_argument('--anatomy-assets', type=Path)
    a = p.parse_args()
    anatomy_assets = getattr(a, 'anatomy_assets', None)
    if anatomy_assets is not None:
        anatomy_assets = anatomy_assets.resolve()
        for subdir in ('source_panels','yeo','aal3'):
            if not (anatomy_assets/subdir).is_dir():
                p.error('Missing external anatomy directory: ' + str(anatomy_assets/subdir))
    if a.mode == 'render':
        if not (a.output/'run_manifest.json').is_file():
            p.error('Only this workflow\'s computed run may be rendered')
        render_with_provenance(a.output, anatomy_assets)
        return
    command = ['Rscript', ROOT/'R/run_analysis.R', '--config', a.config.resolve(), '--mode',
               'synthetic' if a.mode == 'synthetic' else 'full']
    if a.mode != 'synthetic':
        command += ['--input', a.input.resolve()]
        if a.correction:
            command += ['--correction', a.correction.resolve()]
    if a.mode == 'check':
        run(*command, '--check-input')
        return
    if a.mode == 'run':
        if not a.authorize_private_analysis:
            p.error('Protected-data inference requires --authorize-private-analysis')
        command += ['--authorize-private-analysis']
    if a.cores < 1:
        p.error('--cores must be positive')
    try:
        output = external_output_path(a.output)
    except ValueError as e:
        p.error(str(e))
    if output.exists():
        p.error('Output must be NEW; previous results are never replaced')
    command += ['--output', output, '--cores', a.cores]
    if a.component_map_dir:
        command += ['--component-map-dir', a.component_map_dir.resolve()]
    starting_hashes = source_hashes()
    try:
        run(*command)
        validate_outputs(output)
        render_with_provenance(output, anatomy_assets)
        hashes = source_hashes()
        if hashes != starting_hashes:
            raise RuntimeError('Source files changed during execution; rerun with unchanged source')
        (output/'code_and_asset_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
        manifest_path = output/'run_manifest.json'
        manifest = json.loads(manifest_path.read_text())
        manifest.update(driver_command=sys.argv, python=sys.version, code_and_asset_hashes=hashes)

        with (output/'tables/method_review_required.csv').open() as f:
            issue_count = len(list(csv.DictReader(f)))
        manifest['method_review_rows'] = issue_count
        manifest['status'] = 'COMPUTED_WITH_METHOD_REVIEW' if issue_count else 'COMPUTED'
        manifest['empirical_validation'] = 'NOT_EXECUTED' if a.mode == 'synthetic' else 'NOT_ASSESSED'
        manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
        (output/'RUN_STATUS.txt').write_text(f"status: {manifest['status']}\nmode: {a.mode}\nmethod_review_rows: {issue_count}\n")
    except Exception:
        if output.is_dir():
            (output/'RUN_STATUS.txt').write_text('status: FAILED\nPartial outputs retained for diagnosis.\n')
        raise

if __name__ == '__main__':
    main()
