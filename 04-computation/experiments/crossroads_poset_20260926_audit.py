"""Replay five poset/Collatz lanes with assertions both enabled and disabled.

--record refreshes raw stdout intentionally. All arithmetic checks in the
lanes use explicit exceptions. The manifest pins LF-normalized content.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
PREFIX = 'crossroads_poset_20260926_'
LANES = ('bridge', 'integer', 'barrier', 'width', 'moran_audit')


def normalized(data):
    return data.replace(b'\r\n', b'\n')


def replay(lane, record):
    script = ROOT/'04-computation'/'experiments'/(PREFIX+lane+'.py')
    output = ROOT/'05-knowledge'/'results'/(PREFIX+lane+'.out')
    normal = subprocess.run([sys.executable,'-B','-X','utf8',str(script)],
                            cwd=ROOT,capture_output=True,check=True).stdout
    optimized = subprocess.run([sys.executable,'-B','-O','-X','utf8',str(script)],
                               cwd=ROOT,capture_output=True,check=True).stdout
    if normalized(normal) != normalized(optimized):
        raise RuntimeError(lane+': optimized mismatch')
    if record:
        output.write_bytes(normalized(normal))
    if normalized(normal) != normalized(output.read_bytes()):
        raise RuntimeError(lane+': retained-output mismatch')
    return lane+': normal = optimized = retained stdout'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record',action='store_true')
    parser.add_argument('--manifest',action='store_true')
    parser.add_argument('--verify-manifest',action='store_true')
    args = parser.parse_args()
    manifest_path = ROOT/'05-knowledge'/'results'/(PREFIX+'manifest.json')
    if args.verify_manifest:
        manifest = json.loads(manifest_path.read_text(encoding='utf8'))
        for item in manifest['files']:
            actual = hashlib.sha256(normalized((ROOT/item['path']).read_bytes())).hexdigest()
            if actual != item['sha256']:
                raise RuntimeError('manifest mismatch: '+item['path'])
        print('PASS: '+str(len(manifest['files']))+' pinned LF-normalized files')
        return
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(replay,lane,args.record) for lane in LANES]
        for future in futures:
            print(future.result(),flush=True)
    if args.manifest:
        files = []
        for folder in ('04-computation/experiments','05-knowledge/results'):
            files.extend(p for p in (ROOT/folder).glob(PREFIX+'*')
                         if p.is_file() and p.suffix in ('.py','.md','.out'))
        files.extend(ROOT/p for p in (
            '01-canon/theorems/THM-4502-collatz-41-tail-growth-family.md',
            '01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md',
            '01-canon/theorems/THM-4504-families-of-27-moran-function-of-the-inverse-tree.md',
            '05-knowledge/results/procgen_family27_20260926_long_orbit_families.md',
            '05-knowledge/hypotheses/HYP-9161-landing-multiplicity-polylog.md'))
        manifest = {'hash_basis':'UTF-8 file bytes with CRLF normalized to LF',
                    'files':[{'path':p.relative_to(ROOT).as_posix(),
                              'sha256':hashlib.sha256(normalized(p.read_bytes())).hexdigest()}
                             for p in sorted(files)]}
        manifest_path.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')
    print('PASS: five exact lanes; submitted global claims remain under review.')


if __name__ == '__main__':
    main()
